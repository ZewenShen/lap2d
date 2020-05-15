function laphssinv(hsstree)
%LAPHSSINV Compute the inverse factors given a HSS decomposition.
% 
%       Input parameters:
% hsstree - A LapHssTree handler.


%... process leaves
for inode=hsstree.nnode:-1:2

    Dtilde=hsstree.getDtilde(inode);
    
    hsstree.updateInvFactors(inode,Dtilde);
    
end

%... Process the first node and compute G1

Dtilde=hsstree.getDtilde(1);

hsstree.G(1:size(Dtilde,1),1:size(Dtilde,2),1)=inv(Dtilde);

hsstree.invfact_tree(1,[LapHssTree.G_H LapHssTree.G_W])=size(Dtilde); 

end