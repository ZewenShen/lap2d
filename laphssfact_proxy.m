function laphssfact_proxy(coord,norml,curv,whts,alpha,hsstree,...
    npanel,nlege,ideps,lpanel)
%LAPHSSFACT_PROXY Compute the HSS factorization of the structured matrix with the use of proxy.
% 
%       Input parameters:
%
%   coord - A (npanel*nlege) by 2 matrix that represents the coordinates of 
%           Gaussian nodes on the boundary.
%   norml - A (npanel*nlege) by 2 matrix that represents the normal of 
%           Gaussian nodes to the boundary.
%   curv - A 1 by (npanel*nlege) array that represents the curvature of
%          Gaussian nodes on the boundary.
%   whts - A (npanel*lege) by 1 array that represents the Gaussian weights.
%   alpha - A constant coefficient of the identity term. Default = 0.5. 
%   hsstree - A LapHssTree handler that stores the inverse factors and
%             other useful information about the problem.
%   npanel - The number of leaf panels.
%   nlege - the number of Legendre nodes on each leaf panel.
%   ideps - The absolute error tolerance used by rskeleton.
%   lpanel - The arclength of each leaf panel.


%... Init of helper variables
nnode=hsstree.nnode;
nlevel=hsstree.nlevel;
nproxy=50;
sqrt_whts=sqrt(whts);


%... process leaves
for j=1:npanel

    %... Read basic information related to the current leaf.
    bd_idx=get_panelpts_idx(j,j,nlege); %pick the boundary indices that the panel include
    compl_bd_idx=[get_panelpts_idx(j-1,j-1,nlege) get_panelpts_idx(j+1,j+1,nlege)];
    c_sqrt_whts=sqrt_whts(bd_idx);
    c_coord=coord(bd_idx,:); %locations of all points in the current panel
    c_norml=norml(bd_idx,:);
    c_curv=curv(bd_idx);
    inode=nnode-npanel+j; %leave bd_idx in the tree. From nnode-npanel+1 to nnode.


    %... Construction of B matrices
    hsstree.updateLeafDiag(inode,c_coord,c_norml,c_curv,c_sqrt_whts,alpha);

    [proxy_coord,proxy_norml,proxy_sqrt_whts,x0,r]=...
        hsstree.updateLeafProxy(inode,c_coord,nproxy,lpanel);

    %... compute skeleton points of each leave
    [U,irows]=hsstree.computeLeafSkel(c_coord,c_norml,c_sqrt_whts,...
                           compl_bd_idx,coord,norml,sqrt_whts,ideps,...
                           proxy_coord,proxy_norml,proxy_sqrt_whts,x0,r);
                       
   %... update interp and skeleton_idx, plus corresponding sizes in discret_tree
    hsstree.updateInterpSkel(inode,U,irows,bd_idx);

end




%... Process parent nodes

for l=nlevel-1:-1:2 %from second last level to the second level
    
    inode_list=2^(l-1):2^l-1;
    
    hsstree.updateLevelProxy(inode_list,nproxy,lpanel);

    proxyidx2incircle=hsstree.computeNearPointTable(inode_list,coord);
    
    for j=1:length(inode_list)
        
        inode=inode_list(j);
        
        
        [iskel,c_coord,c_norml,c_sqrt_whts,lskel_coord,lskel_norml,...
         lskel_sqrt_whts,rskel_coord,rskel_norml,rskel_sqrt_whts,...
         proxy_coord,proxy_norml,proxy_sqrt_whts]=...
                hsstree.getParentInfo(inode,coord,norml,sqrt_whts,lpanel,nproxy);
            
            
        %... Construction of B matrices
        
        hsstree.updateParentDiag(inode,lskel_coord,lskel_norml,...
                lskel_sqrt_whts,rskel_coord,rskel_norml,rskel_sqrt_whts);
            
        %... ID
        [U,irows]=hsstree.computeParentSkel(c_coord,c_norml,c_sqrt_whts,...
                coord,norml,sqrt_whts,ideps,proxy_coord,...
                proxy_norml,proxy_sqrt_whts,j,proxyidx2incircle,inode);

        
        %... update interp and skeleton_idx, plus corresponding sizes in discret_tree
        hsstree.updateInterpSkel(inode,U,irows,iskel);

    end

end

%... handle the first level: Init


[~,~,~,~,lskel_coord,lskel_norml,lskel_sqrt_whts,...
    rskel_coord,rskel_norml,rskel_sqrt_whts]=hsstree.getParentInfo(1,coord,norml,sqrt_whts);


%... Construction of B matrices
hsstree.updateParentDiag(1,lskel_coord,lskel_norml,lskel_sqrt_whts,...
    rskel_coord,rskel_norml,rskel_sqrt_whts);


end

function idx = get_panelpts_idx(istart,iend,nlege)
idx=1+nlege*(istart-1):nlege*iend;
end
