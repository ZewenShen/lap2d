classdef LapHssTree < handle
    %LAPHSSTREE A handle class for storing and manipulating the HSS tree.
    %   Detailed explanation goes here
    
    properties(GetAccess='public',SetAccess='public')
        
        %... HSS decomposition factors
        
        interp;
        diag;
        proxy_center_info;
        skeleton_idx;
        
        
        %... inverse factors
        
        Dhat;
        E;
        F;
        G;
        invfact_tree;
        
        
        %... Discretization tree
        
        discret_tree;
        perm;
        
        
        %... input
        
        nnode;
        nlevel;
        npanel;
        nlege;
    end
    
    properties(Constant=true)
        
        %... discret_tree constants
        
        LEFT_CHILD=2;
        RIGHT_CHILD=3;
        PARENT=4;
        START=5;
        END=6;
        DIAG_SIZE=7;
        %height and width of interpolation matrix U 
        %(we enforce the right and left interp mat to be the same)
        U_H=8; 
        U_W=9; 
        SKELETON_SIZE=10;
        
        
        %... proxy_center_info constants
        
        PROXY_X=2;
        PROXY_Y=3;
        PROXY_R=4;
        
        
        %... invfact_tree constants
      
        D_H=2;
        D_W=3;
        E_H=4;
        E_W=5;
        F_H=6;
        F_W=7;
        G_H=8;
        G_W=9;
    end
    
    methods
        function obj = LapHssTree(discret_tree,perm,npanel,nlege,lw)
            %LAPHSSTREE Construct an instance of this class
            
            nnode=size(discret_tree,1);
            obj.nnode=nnode;
            obj.nlevel=log2(nnode+1);
            obj.npanel=npanel;
            obj.nlege=nlege;
            
            obj.discret_tree=discret_tree;
            obj.perm=perm;
            
            obj.proxy_center_info=-1*ones(nnode,6);
            obj.proxy_center_info(:,1)=1:nnode;
            
            obj.skeleton_idx=-1*ones(nnode,lw);
            
            %use zeros since off diagonal blocks should be 0.
            obj.diag=zeros(lw,lw,nnode); 
            
            obj.interp=-ones(lw,lw,nnode);
            
            obj.Dhat=-ones(lw,lw,nnode);
            obj.E=-ones(lw,lw,nnode);
            obj.F=-ones(lw,lw,nnode);
            obj.G=-ones(lw,lw,nnode);
            
            obj.invfact_tree=-ones(nnode,10);
            obj.invfact_tree(:,1)=1:nnode;
        end
        
        
        %... Functions for hssfact
        
        function updateLeafDiag(obj,inode,c_coord,c_norml,c_curv,c_sqrt_whts,alpha)
            %UPDATELEAFDIAG Update the D matrix of the input leaf node.
            %
            %   What it computes is the diagonal block, which does not have
            %   low rank structure and cannot be compressed.
            %
            %       Input parameters:
            %   inode - The index of the leaf node in discret_tree.
            %   c_coord - The coordinates of the Gaussian nodes on the panel.
            %   c_curv - 1 by nlege curvature array
            %   c_sqrt_whts - nlege by 1 sqrt Gaussian weight array.

            obj.diag(1:obj.nlege,1:obj.nlege,inode)=...
                lapintmat_dipole(c_coord,c_coord,c_norml,c_sqrt_whts,...
                c_sqrt_whts,alpha,c_curv);

            obj.discret_tree(inode,LapHssTree.DIAG_SIZE)=obj.nlege;
            
        end
        
        function [proxy_coord,proxy_norml,proxy_sqrt_whts,x0,r]=...
                updateLeafProxy(obj,inode,c_coord,nproxy,lpanel)
            %UPDATELEAFPROXY Update the proxy information of the leaf node.
            
            x0=(c_coord(1,:)+c_coord(end,:))/2;
            [proxy_coord,proxy_norml,r]=lapcircle_full(x0,nproxy,lpanel);
            
            proxy_whts=2*pi*r/nproxy*ones(nproxy,1);
            proxy_sqrt_whts=sqrt(proxy_whts);
            
            %... update proxy center info
            obj.proxy_center_info(inode,LapHssTree.PROXY_X:LapHssTree.PROXY_Y)=x0;
            obj.proxy_center_info(inode,LapHssTree.PROXY_R)=r;
            
        end
        
        function [U,irows] = computeLeafSkel(obj,c_coord,c_norml,...
                c_sqrt_whts,compl_bd_idx,coord,norml,sqrt_whts,ideps,...
                proxy_coord,proxy_norml,proxy_sqrt_whts,x0,r)
            %COMPUTELEAFSKEL Compute the skeleton of leaf panel by ID.
            %
            %       Input parameters:
            %   compl_bd_idx - The index of nodes on neighboring panels.
            %   x0 - Coordinate of proxy circle center.
            
            
            %... compute near field info
            npoint=obj.nlege*obj.npanel;
            
            nmaxnearpoint=obj.nlege*2;
            near_coord=zeros(nmaxnearpoint,2);
            near_norml=zeros(nmaxnearpoint,2);
            near_sqrt_whts=zeros(nmaxnearpoint,1);

            istack=1;
            
            compl_bd_coord=LapHssTree.getper_ran(coord,compl_bd_idx,npoint);
            compl_bd_norml=LapHssTree.getper_ran(norml,compl_bd_idx,npoint);
            compl_bd_sqrt_whts=LapHssTree.getper_ran(sqrt_whts,compl_bd_idx,npoint);
            
            for k=1:size(compl_bd_coord,1)
                if LapHssTree.is_in_circle(compl_bd_coord(k,:),x0,r)
                    near_coord(istack,:)=compl_bd_coord(k,:);
                    near_norml(istack,:)=compl_bd_norml(k,:);
                    near_sqrt_whts(istack)=compl_bd_sqrt_whts(k);
                    istack=istack+1;
                end
            end
            
            nnear=istack-1;
            
            near_coord=near_coord(1:nnear,:);
            near_norml=near_norml(1:nnear,:);
            near_sqrt_whts=near_sqrt_whts(1:nnear);
            
            
            %... compute int mats
            intmat=LapHssTree.compute_intmat(c_coord,c_norml,near_coord,...
                near_norml,proxy_coord,proxy_norml,c_sqrt_whts,...
                near_sqrt_whts,proxy_sqrt_whts);

            %... ID
            [irows,~,~,~,U,~,~,~]=rskeleton(intmat,ideps,'silent'); %assert U*concat_blocks(irows,:)=concat_blocks
            
        end
        
        function updateLevelProxy(obj,inode_list,nproxy,lpanel)
            %UPDATELEVELPROXY Create the proxy information for the whole parent level.
            %
            %       Input parameters:
            %   inode_list - A list of node indices.
            %   lpanel - arclength of the each node on the corresponding
            %            parent level.
            
            
            for inode=inode_list
        
                lchild=obj.discret_tree(inode,LapHssTree.LEFT_CHILD);
                rchild=obj.discret_tree(inode,LapHssTree.RIGHT_CHILD);


                %... new proxy setup
                [lproxy_x,lproxy_y,~]=...
                    LapHssTree.get_proxy_info(lchild,obj.proxy_center_info);
                [rproxy_x,rproxy_y,~]=...
                    LapHssTree.get_proxy_info(rchild,obj.proxy_center_info);
                
                
                x0=[(lproxy_x+rproxy_x)/2 (lproxy_y+rproxy_y)/2];
                obj.proxy_center_info(inode,LapHssTree.PROXY_X:LapHssTree.PROXY_Y)=x0;

                
                cur_npanel=obj.discret_tree(inode,LapHssTree.END)-...
                           obj.discret_tree(inode,LapHssTree.START)+1;
                cur_lpanel=cur_npanel*lpanel;

                
                [~,~,r]=lapcircle_full(x0,nproxy,cur_lpanel);
                obj.proxy_center_info(inode,LapHssTree.PROXY_R)=r;

            end
            
        end
        
        function proxyidx2incircle = computeNearPointTable(obj,inode_list,coord)
            %COMPUTENEARPOINTTABLE Compute points in the near field for each parent node on the given level.
            
            proxy_centers=...
                obj.proxy_center_info(inode_list,LapHssTree.PROXY_X:LapHssTree.PROXY_Y);
            proxy_radius=...
                obj.proxy_center_info(inode_list,LapHssTree.PROXY_R);
            
            proxyidx2incircle=false(length(inode_list),size(coord,1));
            
            for j=1:size(coord,1)
            
                cur_pt=coord(j,:);
                cur_in=LapHssTree.is_in_circle(cur_pt,proxy_centers,proxy_radius);
                proxyidx2incircle(cur_in,j)=true;

            end
            
        end
        
        function [iskel,c_coord,c_norml,c_sqrt_whts,lskel_coord,lskel_norml,...
                lskel_sqrt_whts,rskel_coord,rskel_norml,rskel_sqrt_whts,...
                proxy_coord,proxy_norml,proxy_sqrt_whts] = ...
                getParentInfo(obj,inode,coord,norml,sqrt_whts,lpanel,nproxy)
            
            lchild=obj.discret_tree(inode,LapHssTree.LEFT_CHILD);
            rchild=obj.discret_tree(inode,LapHssTree.RIGHT_CHILD);

            %... Init of skeleton points info
            nlskel=obj.discret_tree(lchild,LapHssTree.SKELETON_SIZE);
            nrskel=obj.discret_tree(rchild,LapHssTree.SKELETON_SIZE);
            ilskel=obj.skeleton_idx(lchild,1:nlskel);
            irskel=obj.skeleton_idx(rchild,1:nrskel);
            iskel=[ilskel irskel];


            c_coord=coord(iskel,:); %locations of all points in the current panel
            c_norml=norml(iskel,:);
            c_sqrt_whts=sqrt_whts(iskel);

            lskel_coord=coord(ilskel,:);
            lskel_norml=norml(ilskel,:);
            lskel_sqrt_whts=sqrt_whts(ilskel);

            rskel_coord=coord(irskel,:);
            rskel_norml=norml(irskel,:);
            rskel_sqrt_whts=sqrt_whts(irskel);
            
            %... when proxy info is not needed (first level)
            if nargin~=5
                [x,y,r]=obj.get_proxy_info(inode,obj.proxy_center_info);
                x0=[x y];
                proxy_whts=2*pi*r/nproxy*ones(nproxy,1);
                proxy_sqrt_whts=sqrt(proxy_whts);

                cur_npanel=obj.discret_tree(inode,LapHssTree.END)-...
                           obj.discret_tree(inode,LapHssTree.START)+1;
                cur_lpanel=cur_npanel*lpanel;

                [proxy_coord,proxy_norml,~]=lapcircle_full(x0,nproxy,cur_lpanel);
            end

        end
        
        function updateParentDiag(obj,inode,lskel_coord,lskel_norml,...
                lskel_sqrt_whts,rskel_coord,rskel_norml,rskel_sqrt_whts)
            %UPDATEPARENTDIAG Update the B matrix of the input parent node.
            %
            %   What it computes is the diagonal block, which does not have
            %   low rank structure and cannot be compressed.
            
            nlskel=size(lskel_coord,1);
            nrskel=size(rskel_coord,1);
            
            obj.diag(1:nlskel,nlskel+1:nlskel+nrskel,inode)=...
                lapintmat_dipole(lskel_coord,rskel_coord,rskel_norml,...
                lskel_sqrt_whts,rskel_sqrt_whts);
            
            obj.diag(nlskel+1:nlskel+nrskel,1:nlskel,inode)=...
                lapintmat_dipole(rskel_coord,lskel_coord,lskel_norml,...
                rskel_sqrt_whts,lskel_sqrt_whts);

            obj.discret_tree(inode,LapHssTree.DIAG_SIZE)=nlskel+nrskel;

        end
        
        function [U,irows] = computeParentSkel(obj,c_coord,c_norml,...
                c_sqrt_whts,coord,norml,sqrt_whts,ideps,proxy_coord,...
                proxy_norml,proxy_sqrt_whts,j,proxyidx2incircle,inode)
            %COMPUTELEAFSKEL Compute the skeleton of parent panel by ID.
        
            %... compute near field points
            compl_bd_idx=true(1,size(coord,1));
            compl_bd_idx(LapHssTree.get_panelpts_idx(...
                obj.discret_tree(inode,LapHssTree.START),...
                obj.discret_tree(inode,LapHssTree.END),...
                obj.nlege)...
                )=false;
            
            inear=compl_bd_idx&proxyidx2incircle(j,:);

            near_coord=coord(inear,:);
            near_norml=norml(inear,:);
            near_sqrt_whts=sqrt_whts(inear);
            
            
            %... compute int mats
            intmat=LapHssTree.compute_intmat(c_coord,c_norml,near_coord,...
                near_norml,proxy_coord,proxy_norml,c_sqrt_whts,...
                near_sqrt_whts,proxy_sqrt_whts);



            %... ID
            [irows,~,~,~,U,~,~,errout]=rskeleton(intmat,ideps,'silent'); %assert U*concat_blocks(irows,:)=concat_blocks

        end
        
        function updateInterpSkel(obj,inode,U,irows,bd_idx)
            %UPDATEINTERPSKEL Update the interpolative skeleton for the given node.
            
            obj.interp(1:size(U,1),1:size(U,2),inode)=U;
            obj.discret_tree(inode,LapHssTree.U_H)=size(U,1);
            obj.discret_tree(inode,LapHssTree.U_W)=size(U,2);

            obj.skeleton_idx(inode,1:length(irows))=bd_idx(irows);
            obj.discret_tree(inode,LapHssTree.SKELETON_SIZE)=length(irows);
            
        end
        
        
        
        %... Functions for hssinv
        
        function Dtilde = getDtilde(obj,inode)
            %GETDTILDE Get the Dtilde matrix given an inode.
            
            if inode>obj.nnode-obj.npanel
                %... Is a leaf.
                Dtilde=...
                    LapHssTree.extract_mat(inode,obj.diag,obj.discret_tree,...
                    [LapHssTree.DIAG_SIZE LapHssTree.DIAG_SIZE]);

            else
                %... Is a parent.
                lchild=obj.discret_tree(inode,LapHssTree.LEFT_CHILD);
                rchild=obj.discret_tree(inode,LapHssTree.RIGHT_CHILD);

                [lDhat,lDhat_h,lDhat_w]=...
                    LapHssTree.extract_mat(lchild,obj.Dhat,obj.invfact_tree,...
                    [LapHssTree.D_H LapHssTree.D_W]);
                [rDhat,rDhat_h,rDhat_w]=...
                    LapHssTree.extract_mat(rchild,obj.Dhat,obj.invfact_tree,...
                    [LapHssTree.D_H LapHssTree.D_W]);


                Dtilde=...
                    LapHssTree.extract_mat(inode,obj.diag,obj.discret_tree,...
                    [LapHssTree.DIAG_SIZE LapHssTree.DIAG_SIZE]);
                Dtilde(1:lDhat_h,1:lDhat_w)=lDhat;
                Dtilde(end-rDhat_h+1:end,end-rDhat_w+1:end)=rDhat;
            end
        end
        
        function updateInvFactors(obj,inode,Dtilde)
            %UPDATEINVFACTORS Update the inverse factors given an inode and its corresponding Dtilde.
            %
            %   Reference: Martinsson paper's algo3. Note that his assignment for G is wrong.

            U=LapHssTree.extract_mat(inode,obj.interp,obj.discret_tree,...
                [LapHssTree.U_H,LapHssTree.U_W]);
            invDtilde=inv(Dtilde);

            tmpDhat=inv(U'*invDtilde*U);
            tmpE=invDtilde*U*tmpDhat;
            tmpF=tmpDhat*U'*invDtilde;
            tmpG=invDtilde-invDtilde*U*tmpDhat*U'*invDtilde;


            %... Store the inverse factors
            
            obj.invfact_tree(inode,[LapHssTree.D_H LapHssTree.D_W])=size(tmpDhat);
            obj.Dhat(1:size(tmpDhat,1),1:size(tmpDhat,2),inode)=tmpDhat;
            
            obj.invfact_tree(inode,[LapHssTree.E_H LapHssTree.E_W])=size(tmpE);
            obj.E(1:size(tmpE,1),1:size(tmpE,2),inode)=tmpE;

            obj.invfact_tree(inode,[LapHssTree.F_H LapHssTree.F_W])=size(tmpF);
            obj.F(1:size(tmpF,1),1:size(tmpF,2),inode)=tmpF;

            obj.invfact_tree(inode,[LapHssTree.G_H LapHssTree.G_W])=size(tmpG);
            obj.G(1:size(tmpG,1),1:size(tmpG,2),inode)=tmpG;
            
        end
    end
    
    
    methods(Static=true,Access=private)
        function [x,y,r] = get_proxy_info(inode,proxy_center_info)
        x=proxy_center_info(inode,LapHssTree.PROXY_X);
        y=proxy_center_info(inode,LapHssTree.PROXY_Y);
        r=proxy_center_info(inode,LapHssTree.PROXY_R);
        end
        
        function idx = get_panelpts_idx(istart,iend,nlege)
        idx=1+nlege*(istart-1):nlege*iend;
        end
        
        function pts = getper_ran(car,range,npoint)
        % for discontinuous range (unlike getper_se)
        
        idx1=range<1;
        range(idx1)=range(idx1)+npoint;
        idx2=range>npoint;
        range(idx2)=range(idx2)-npoint;
        pts=car(range,:);
        
        end
        
        function in = is_in_circle(pt,car_center,radius)
            
        n=size(car_center,1);
        assert(n==size(radius,1),'size(car_center,1)~=size(radius,1)');

        d=sqrt(sum((pt-car_center).^2,2));
        in=d<radius-(1e-15);
        
        end
        
        function intmat = compute_intmat(c_coord,c_norml,near_coord,...
                near_norml,proxy_coord,proxy_norml,c_sqrt_whts,...
                near_sqrt_whts,proxy_sqrt_whts)
            
        intmat_n2t=...
            lapintmat_dipole(c_coord,near_coord,near_norml,c_sqrt_whts,near_sqrt_whts);
        intmat_p2t=...
            lapintmat_dipole(c_coord,proxy_coord,proxy_norml,c_sqrt_whts,proxy_sqrt_whts);
        intmat_t2n=...
            lapintmat_dipole(near_coord,c_coord,c_norml,near_sqrt_whts,c_sqrt_whts);
        intmat_t2p=...
            lapintmat_dipole(proxy_coord,c_coord,c_norml,proxy_sqrt_whts,c_sqrt_whts);
        
        intmat=[intmat_p2t intmat_n2t intmat_t2p' intmat_t2n'];
        
        end
        
        function [mat,h,w] = extract_mat(inode,mat_table,tree,field)
        h=tree(inode,field(1));
        w=tree(inode,field(2));
        mat=mat_table(1:h,1:w,inode);
        end
    end
end

