function [xo_p,yo_p,zo_p]=modeuler_tl(tstep,h,j,z,y,x,...
   zval_p,yval_p,xval_p,a,r,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: Tangent linear code of modeuler.m
%
%  %
%  List of main variables
%    a:          sigma coefficient in equations
%    r:          rho coefficient in equations
%    b:          beta coefficient in equations
%    h:          Time step for numerical scheme
%    tstep:      Number of time steps to perform
%    [xval,yval,zval]: Linearization state fields at all times
%    [xval_p,yval_p,zval_p]: Initial perturbations
%    j:          Index to pick up correct initial field
%
%  Output:
%    [xo_p,yo_p,zo_p]: Evolved perturbations
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
x_p=zeros(tstep,1);
y_p=zeros(tstep,1);
z_p=zeros(tstep,1);
%
x_p(1)=xval_p(j+1);
y_p(1)=yval_p(j+1);
z_p(1)=zval_p(j+1);
%
for i = 1:tstep
  kx1 = h*a*(y(i)-x(i));
  ky1 = h*(r*x(i)-y(i)-x(i)*z(i));
  kz1 = h*(x(i)*y(i) - b*z(i));
%
  kx1_p = h*a*(y_p(i)-x_p(i));
  ky1_p = h*(r*x_p(i)-y_p(i)-x_p(i)*z(i)-x(i)*z_p(i));
  kz1_p = h*(x_p(i)*y(i)+x(i)*y_p(i)-b*z_p(i));

  kx2_p = h*a*(y_p(i)+ky1_p-x_p(i)-kx1_p);
  ky2_p = h*(r*(x_p(i)+kx1_p)-y_p(i)-ky1_p ...
          - (x_p(i)+kx1_p)*(z(i)+kz1) - (x(i)+kx1)*(z_p(i)+kz1_p));
  kz2_p = h*((x_p(i)+kx1_p)*(y(i)+ky1)+(x(i)+kx1)*(y_p(i)+ky1_p) ...
          - b*(z_p(i)+kz1_p));
%
  x_p(i+1)=x_p(i)+0.5d0*(kx1_p+kx2_p);
  y_p(i+1)=y_p(i)+0.5d0*(ky1_p+ky2_p);
  z_p(i+1)=z_p(i)+0.5d0*(kz1_p+kz2_p);
end

xo_p=x_p;
yo_p=y_p;
zo_p=z_p;
%
