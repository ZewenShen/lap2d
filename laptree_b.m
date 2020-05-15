function [discret_tree,perm] = laptree_b(npanel)
%LAPTREE_B Create a hierarchical binary tree of the panels.
%
%       Input parameters:
%   npanel - The number of leaf panels. Has to be a power of 2.
%
%       Output parameters:
%   discret_tree - The discret_tree in LapHssTree.
%   perm - A permutation array. Not used in our project.

nnode=npanel*2-1;
discret_tree=-ones(nnode,10);
perm=zeros(1,npanel);

nlevel=log2(nnode+1);
assert(abs(round(nlevel)-nlevel)<10*eps,'npanel must be in the form 2^j');


istack=0;


for l=1:nlevel
    nnode_per_level=2^(l-1);
    node_size=npanel/nnode_per_level;
    
    for j=1:nnode_per_level % number of nodes on level l
        istack=istack+1;
        discret_tree(istack,1)=istack;

        discret_tree(istack,LapHssTree.START)=(j-1)*node_size+1;
        discret_tree(istack,LapHssTree.END)=j*node_size;

        if l~=nlevel % not leave
            % children info setup
            left=istack*2;
            right=istack*2+1;

            discret_tree(istack,LapHssTree.LEFT_CHILD)=left;
            discret_tree(istack,LapHssTree.RIGHT_CHILD)=right;

            discret_tree(left,LapHssTree.PARENT)=istack;
            discret_tree(right,LapHssTree.PARENT)=istack;
        else % is leave

        end
    end
end

end