function q = laphssinv_apply(u,E,F,G,invfact_tree,discret_tree)
%LAPHSSINV_APPLY Compute the weighted density given weighted boundary data and the HSS tree.
%
%   Won't be directly called by the user. Used by laphssdensity.
%
%       Input parameters:
%
%   u - A (npanel*lege) by 1 array that represents weighted boundary 
%       condition data.
%   The rest of inputs are defined in LapHssTree.
%
%       Output parameters:
%
%   q - A (npanel*lege) by 1 array that represents the unweighted
%       density on the boundary.



%... Init helper variables
nnode=size(discret_tree,1);
nlevel=log2(nnode+1);
npanel=2^(nlevel-1);
nlege=discret_tree(end,LapHssTree.DIAG_SIZE);



%... Init helper variables for storing new constants
assert(length(u)==nlege*npanel,"length(u)~=legen*npanel");
uhat=zeros(length(u),nnode);
qhat=zeros(length(u),nnode);
q=zeros(length(u),1);

size_tree_=-ones(nnode,10);
size_tree_(:,1)=1:nnode;
UHAT=2;
QHAT=3;



%... part 1 and 2 of the apply algo. Bottom -> second level.
for imax=nnode:-1:2
    
    %pick the boundary indices that the panel include
    bd_idx=panelpts_idx_new_(imax,discret_tree,nlege); 

    cur_F=extract_(imax,F,invfact_tree,[LapHssTree.F_H LapHssTree.F_W]);
    
    if imax>nnode-npanel
        %... Is a leaf.
        

        uhat(1:size(cur_F,1),imax)=cur_F*u(bd_idx);

    else
        %... Is a parent.
        
        lchild=discret_tree(imax,LapHssTree.LEFT_CHILD);
        rchild=discret_tree(imax,LapHssTree.RIGHT_CHILD);
        
        uhat(1:size(cur_F,1),imax)=...
            cur_F*[uhat(1:size_tree_(lchild,UHAT),lchild); uhat(1:size_tree_(rchild,UHAT),rchild)];
        
    end
    
    size_tree_(imax,UHAT)=size(cur_F,1);
    
end

%... part 3 of the apply algo. Process first level.

G1=extract_(1,G,invfact_tree,[LapHssTree.G_H LapHssTree.G_W]);


tmp_q=G1*[uhat(1:size_tree_(2,UHAT),2); uhat(1:size_tree_(3,UHAT),3)];


qhat(1:size_tree_(2,UHAT),2)=tmp_q(1:size_tree_(2,UHAT));
qhat(1:size_tree_(3,UHAT),3)=tmp_q(size_tree_(2,UHAT)+1:end);
assert(length(tmp_q)-size_tree_(2,UHAT)==size_tree_(3,UHAT));

size_tree_([2 3],QHAT)=[size_tree_(2,UHAT);size_tree_(3,UHAT)];


%... part 4 and 5 of the apply algo. Second -> bottom.

for imax=2:nnode
    
    bd_idx=panelpts_idx_new_(imax,discret_tree,nlege);
    
    cur_E=extract_(imax,E,invfact_tree,[LapHssTree.E_H LapHssTree.E_W]);
    cur_G=extract_(imax,G,invfact_tree,[LapHssTree.G_H LapHssTree.G_W]);
    
    if imax<=nnode-npanel
        %... Is a parent.
        
        lchild=discret_tree(imax,LapHssTree.LEFT_CHILD);
        rchild=discret_tree(imax,LapHssTree.RIGHT_CHILD);
        
        tmp_q=cur_E*qhat(1:size_tree_(imax,QHAT),imax)+...
            cur_G*[uhat(1:size_tree_(lchild,UHAT),lchild);...
                   uhat(1:size_tree_(rchild,UHAT),rchild)];
        
        
        qhat(1:size_tree_(lchild,UHAT),lchild)=tmp_q(1:size_tree_(lchild,UHAT));
        qhat(1:size_tree_(rchild,UHAT),rchild)=tmp_q(size_tree_(lchild,UHAT)+1:end);
        
        size_tree_([lchild rchild],QHAT)=[size_tree_(lchild,UHAT);size_tree_(rchild,UHAT)];
        
    else
        %... Is a leaf.
        
        q(bd_idx)=cur_E*qhat(1:size_tree_(imax,QHAT),imax)+cur_G*u(bd_idx);
        
    end
    
end


end

function [mat,h,w] = extract_(inode,mat_table,tree,field)
h=tree(inode,field(1));
w=tree(inode,field(2));
mat=mat_table(1:h,1:w,inode);
end

function idx = panelpts_idx_private(istart,iend,legen)
idx=1+legen*(istart-1):legen*iend;
end

function idx = panelpts_idx_new_(imax,tree_,legen)
START=5;
END=6;
idx=panelpts_idx_private(tree_(imax,START),tree_(imax,END),legen);
end