function [coord,norml,whts,hsstree,lpanel,nana_w] = lapdiri2d_g_hss_proxy(...
    funcurv,par,rl,nanaSampN,npanel,nlege,naeps,ideps,lw,alpha)
%LAPDIRI2D_G_HSS_PROXY Fast direct solver for 2d Dirichlet Laplace's equation.
%
%   Use HSS scheme proposed by Martinsson plus Gaussian quadrature.
%
%       Input parameters:
%
%   funcurv - A function that parametrizes the boundary of the domain.
%   par - Parameters used by funcurv.
%   rl - The domain of funcurv is [0, rl].
%   nanaSampN - The number of points sampled by nanafast. The larger the
%               boundary is, the more points are required.
%   npanel - The number of leaf panels.
%   nlege - The number of Legendre nodes on each leaf panel.
%   naeps - The absolute error tolerance used by nanafast.
%   ideps - The absolute error tolerance used by rskeleton.
%   lw - The size of working array for LapHssTree. npanel/8 is usually a
%        good choice.
%   alpha - A constant coefficient of the identity term. Default = 0.5. 
%
%       Output parameters:
%
%   coord - A (npanel*nlege) by 2 matrix that represents the coordinates of 
%           Gaussian nodes on the boundary.
%   norml - A (npanel*nlege) by 2 matrix that represents the normal of 
%           Gaussian nodes to the boundary.
%   whts - A (npanel*lege) by 1 array that represents the Gaussian weights.
%   hsstree - A LapHssTree handler that stores the inverse factors and
%             other useful information about the problem.
%   lpanel - The arclength of each leaf panel.
%   nana_w - a large array containing, among other things, nested Legendre
%            expansions of the arc length of the curve as a function of the
%            curve parameter. Can be used as an input to the function
%            nana_fparam (see).

default_alpha=-pi; 


[~,cl,nana_w]=nanafast(funcurv,par,rl,nanaSampN,12,naeps,1000000,'silent');
lpanel=cl/npanel; % the discretized interval length.


if nargin==9
    alpha=default_alpha;
end



[legepts,legewhts]=legeexps(nlege);

coord=[];
norml=[];
curv=[];
whts=[];
for i=1:npanel
    a=(i-1)*lpanel; b=i*lpanel;
    pts=(b-a)*legepts/2+(a+b)/2;
    whts=[whts; (b-a)*legewhts'/2];
    tout=nana_fparam(nana_w,naeps,pts,'silent');
    [x,y,dxdt,dydt,d2xdt2,d2ydt2]=funcurv(tout,par);
    coord=[coord; x' y'];
    norml=[norml; -dydt' dxdt'];
    curv=[curv (dxdt.*d2ydt2-dydt.*d2xdt2)./(dxdt.^2+dydt.^2).^1.5];
end



[tree_,perm]=laptree_b(npanel);

hsstree=LapHssTree(tree_,perm,npanel,nlege,lw);

laphssfact_proxy(coord,norml,curv,whts,alpha,hsstree,npanel,...
    nlege,ideps,lpanel);

laphssinv(hsstree);

end