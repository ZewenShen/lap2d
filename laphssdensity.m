function whted_rho = laphssdensity(f,coord,whts,E,F,G,invfact_tree,discret_tree)
%LAPHSSDENSITY Compute the weighted density given boundary condition and the HSS tree.
%
%   A wrapper for laphssinv_apply.m
%
%       Input parameters:
%
%   f - Boundary condition function.
%   coord - A (npanel*nlege) by 2 matrix that represents the coordinates of 
%           Gaussian nodes on the boundary.
%   whts - A (npanel*lege) by 1 array that represents the Gaussian weights.
%   The rest of inputs are defined in LapHssTree.
%
%       Output parameters:
%
%   whted_rho - A (npanel*lege) by 1 array that represents the weighted
%               density on the boundary.

sqrt_whts=sqrt(whts);
b=sqrt_whts.*f(coord(:,1),coord(:,2));
rho=laphssinv_apply(b,E,F,G,invfact_tree,discret_tree);
whted_rho=sqrt_whts.*rho;
end

