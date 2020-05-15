function dG = lapdgreen2d(x,x0,h,curvature)
%LAPDGREEN2D Compute the dipole potential.
%
%   Given a concentrated source at x0, return the first deriv of
%   Green's function along the direction vector that normal to the boundary
%   at x. Can also be regarded as the potential of a unit dipole.
%
%   dG = lapdgreen2d([1 1], [0 0; 3 3], [1 1; -1 -1], [0.5, 1])
%   dG = lapdgreen2d([1 1], [0 0], [-1 -1])
%
%       Input parameters:
%
%   x - A 1 by 2 array that represents the point to evaluate.
%   x0 - An n by 2 matrix that represents the dipole source.
%   h - An n by 2 matrix that represents the dipoles' normal to the boundary 
%       at x0. Not ecessarily ||h|| = 1. Pointing outwards.
%   curvature - a scalar that represents the curvature at x. Only 
%               required when x is very close to x0.
%
%       Output parameters:
%
%   dG - An n by 1 array that represents dipole potential at x. 

h=h./vecnorm(h,2,2);
dG=dot(h,x0-x,2)./vecnorm(x0-x,2,2).^2;

nanidx=sum(abs(x0-x),2)<1e-15;
if any(nanidx,'all')
    dG(nanidx)=-0.5*curvature;
end

end