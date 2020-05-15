function legecoeff = lapdensity_legeexps(rho,nlege,npanel)
%LAPDENSITY_LEGEEXPS Compute the Legendre expansion of the boundary density.
%   
%       Input parameters:
%   rho - an (npanel*nlege) by 1 array that represents the unweighted
%         density.
%   nlege - the number of Legendre nodes on each leaf panel.
%   npanel - the number of leaf panels.
%
%       Output parameters:
%   legecoeff - an (npanel*nlege) by 1 array that represents the Legendre
%               expansion coefficients of the density rho.

assert(length(rho)==nlege*npanel);

[~,u,~,~]=legeexps(nlege);

legecoeff=zeros(length(rho),1);

for j=1:npanel
    
    range=(j-1)*nlege+1:j*nlege;
    
    cur_rho=rho(range);
    
    legecoeff(range)=u*cur_rho;
    
end


end

