function [coord,norml,proxy_r] = lapcircle_full(x0,nproxy,lpanel,r_ratio)
%LAPCIRCLE_FULL Create a proxy circle around x0 based on the panel length.
%
%       Input parameters:
%
%   x0 - a 1 by 2 array that represents the center of circle.
%   nproxy - the number of proxy points that are uniformly distributed on
%            the circle.
%   lpanel - the length of the panel.
%   r_ratio - the radius of the proxy circle would be lpanel/2*r_ratio.
%   
%       Output parameters:
%   
%   coord - an nproxy by 2 matrix that represents the cartesian coordinate of
%           the proxy points on the circle.
%   norml - an nproxy by 2 matrix that represents the normal of the proxy 
%           points on the circle.
%   proxy_r - the radius of the proxy circle.

if nargin==3
    r_ratio=1.5;
end

r=lpanel/2;

angles=linspace(0,2*pi,nproxy+1);
angles=angles(1:end-1);

car_angles=[cos(angles)' sin(angles)'];
coord=x0+r*r_ratio*car_angles;

norml=-car_angles;

proxy_r=r_ratio*r;