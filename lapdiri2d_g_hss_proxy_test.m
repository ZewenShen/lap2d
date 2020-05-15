%Tests for solving 2D Dirichlet Laplace's equation directly by HSS scheme with proxy.

delete diary;
diary on;



f=@(x,y) 2*x.^2-2*y.^2;
% f = @(x, y) x + y;


%
%... Simple tests
%
point=[0.1 -0.3];
npanel=64;
nanaSampN=1000;
nlege=12;
lw=npanel/8;
fprintf("\nSimple test starts");
prinf('npanel=',npanel);
par{1}=5;
par{2}=4;
rl=2*pi;
naeps=1e-12;
ideps=1e-14;

tic;
[coord,norml,whts,hsstree,lpanel,nana_w]=...
     lapdiri2d_g_hss_proxy(@ellipse_curve,par,rl,nanaSampN,npanel,nlege,naeps,ideps,lw);

whted_rho=laphssdensity(f,coord,whts,hsstree.E,hsstree.F,hsstree.G,...
    hsstree.invfact_tree,hsstree.discret_tree);

u_val=lapfparam(@lapdgreen2d,point,coord,norml,whted_rho);
f_val=f(point(1),point(2));
err_adap=abs(u_val-f_val);
prin2_long('u=',u_val);

prin2_long('sol=',f_val);
prin2('err=',err_adap);

toc





%
%... near boundary evaluations
%
fprintf("\nNear boundary eval test starts\n");
point2=[4.9 0];
adapeps=1e-14;
rho=whted_rho./whts;
legecoeff=lapdensity_legeexps(rho,nlege,npanel);

center_coord=hsstree.proxy_center_info(end-npanel+1:end,...
                              [LapHssTree.PROXY_X LapHssTree.PROXY_Y]);
                          
u_val_adap=lapfparam_adap(@lapdgreen2d,point2,coord,norml,whted_rho,legecoeff,lpanel,...
    nlege,center_coord,@ellipse_curve,par,nana_w,naeps,adapeps);

u_val_naive=lapfparam(@lapdgreen2d,point2,coord,norml,whted_rho);

f_val=f(point2(1),point2(2));
err_adap=abs(u_val_adap-f_val);
err_naive=abs(u_val_naive-f_val);
prin2_long('u_val_adap=',u_val_adap);
prin2_long('u_val_naive=',u_val_naive);
prin2_long('sol=',f_val);
prin2('adap err=',err_adap);
prin2('naive err=',err_naive);





%
%... Wobble
%

test_wobble=true;

if test_wobble
    
    fprintf("\nWobble experiment starts");
    
    par{1}=0.1;
    par{2}=5;
    par{3}=1.3;

    ns=[]; 
    acc=[];
    T=[];
    n=3:6;
    for j=1:length(n)
        tic
        [coord,norml,whts,hsstree]=...
            lapdiri2d_g_hss_proxy(@wobble,par,rl,nanaSampN,2^n(j),nlege,naeps,ideps,lw);
        
        whted_rho=laphssdensity(f,coord,whts,hsstree.E,hsstree.F,hsstree.G,...
            hsstree.invfact_tree,hsstree.discret_tree);
        
        u_val=lapfparam(@lapdgreen2d,point,coord,norml,whted_rho);
        f_val=f(point(1),point(2));
        
        
        acc=[acc abs(u_val-f_val)]; 
        ns=[ns 2^n(j)*nlege];
        T(j)=toc;
    end
    prinf('ns=',ns);
    prin2('acc=',acc);
    prin2('err order=',log2(acc(1:end-1)./acc(2:end)));
    prin2('time=',T);
    prin2('time order=',log2(T(2:end)./T(1:end-1)));
end


diary off;







function [x,y,dxdt,dydt,d2xdt2,d2ydt2]=ellipse_curve(t,par)
%
%   describes an ellipse contour in \R^2
%
%                   input parameters:
%
%  t - curve parameter. between 0 and 2*pi.
%  par - cell array where 
%       a=par{1}
%       b=par{2}
%  And the ellipse is expressed as (x/a)^2 + (y/b)^2 = 1
%
%                   output parameters:
%
%  x,y,dxdt,dydt,d2xdt2,d2ydt2 - information about the curve at t
%
a=par{1};
b=par{2};

x=a*cos(t);
y=b*sin(t);

dxdt=-a*sin(t);
dydt=b*cos(t);

d2xdt2=-a*cos(t);
d2ydt2=-b*sin(t);
end

function [x,y,dxdt,dydt,d2xdt2,d2ydt2]=wobble(t,par)
%
%   describes a wobbly closed contour in \R^2
%
%                   input parameters:
%
%  t - curve parameter. between 0 and 2*pi.
%  par - cell array where 
%       mu=par{1} is the intensity of the wobbles
%       k=par{2} is an integer; the freqences of wobbles
%       beta=par{3} is the eccentricity
%
%                   output parameters:
%
%  x,y,dxdt,dydt - information about the curve at t
%
mu=par{1};
k=par{2};
beta=par{3};

x=beta*cos(t).*(1+mu*sin(2*k*t));
y=sin(t).*(1+mu*sin(2*k*t));

dxdt=-beta*sin(t).*(1+mu*sin(2*k*t)) + beta*mu*2*k*cos(t).*cos(2*k*t);
dydt=cos(t).*(1+mu*sin(2*k*t)) + mu*2*k*sin(t).*cos(2*k*t);

d2xdt2=-beta*cos(t)-beta*mu*sin(2*k*t).*cos(t)-...
    4*beta*mu*k*cos(2*k*t).*sin(t)-4*beta*mu*k^2*sin(2*k*t).*cos(t);
d2ydt2=-sin(t)-mu*sin(2*k*t).*sin(t)+...
    4*mu*k*cos(2*k*t).*cos(t)-4*mu*k^2*sin(2*k*t).*sin(t);

return;
end

