function intmat = lapintmat_dipole(x,x0,h,sqrt_t_whts,sqrt_s_whts,alpha,curvature)
%LAPINTMAT_DIPOLE Compute the dipole interaction matrix.
% 
%       Input parameters:
%
%   x - An n by 2 matrix that represents the coordinates of target Gaussian
%       nodes on the boundary.
%   x0 - A m by 2 matrix that represents the coordinates of source Gaussian
%        nodes on the boundary.
%   h - A m by 2 matrix that represents the normal of source Gaussian
%        nodes to the boundary.
%   sqrt_t_whts - An n by 1 array that represents the target sqrt
%                 Gaussian weights.
%   sqrt_t_whts - A m by 1 array that represents the source sqrt
%                 Gaussian weights.
%   alpha - A constant coefficient of the identity term. Default = 0.5. 
%   curvature - An n by 1 array that represents the curvature at x. Only 
%               required when x is very close to x0.


if nargin==7
    
    assert(size(x,1)==size(x0,1));
    intmat=alpha*eye(size(x,1));

else
    
    curvature=zeros(1,size(x,1));
    intmat=zeros(size(x,1),size(x0,1));

end
for j=1:size(x,1)
    intmat(j,:)=intmat(j,:)+sqrt_t_whts(j)*(lapdgreen2d(x(j,:),x0,h,curvature(j)).*sqrt_s_whts)';
end
end

