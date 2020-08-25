function [xval_hat,yval_hat,zval_hat]=modeuler_adj(tstep,h,j,z,y,x,...
   xo_hat,yo_hat,zo_hat,a,r,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: Adjoint code of modeuler.m
%
%   %
%  List of main variables
%    a:          sigma coefficient in equations
%    r:          rho coefficient in equations
%    b:          beta coefficient in equations
%    h:          Time step for numerical scheme
%    tstep:      Number of time steps to perform
%    [xval,yval,zval]: Linearization state fields at all times
%    [xo_hat,yo_hat,zo_hat]: Input adjoint variables
%    j:          Index to pick up correct initial field
%
%  Output:
%    [xval_hat,yval_hat,zval_hat]: Output adjoint variables
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
x_hat=zeros(tstep+1,1);
y_hat=zeros(tstep+1,1);
z_hat=zeros(tstep+1,1);
xval_hat=zeros(tstep+1,1);
yval_hat=zeros(tstep+1,1);
zval_hat=zeros(tstep+1,1);

x_hat(tstep+1) = x_hat(tstep+1) + xo_hat;
y_hat(tstep+1) = y_hat(tstep+1) + yo_hat;
z_hat(tstep+1) = z_hat(tstep+1) + zo_hat;
xo_hat = 0.0d0;
yo_hat = 0.0d0;
zo_hat = 0.0d0;

for i=tstep:-1:1
  kx1 = h*a*(y(i)-x(i));
  ky1 = h*(r*x(i)-y(i)-x(i)*z(i));
  kz1 = h*(x(i)*y(i) - b*z(i));
%
  kz1_hat = 0.5 * z_hat(i+1);
  kz2_hat = 0.5 * z_hat(i+1);
  z_hat(i) = z_hat(i) + z_hat(i+1);
  z_hat(i+1) = 0;
  ky1_hat = 0.5 * y_hat(i+1);
  ky2_hat = 0.5 * y_hat(i+1);
  y_hat(i) = y_hat(i) + y_hat(i+1);
  y_hat(i+1) = 0;
  kx1_hat = 0.5 * x_hat(i+1);
  kx2_hat = 0.5 * x_hat(i+1);
  x_hat(i) = x_hat(i) + x_hat(i+1);
  x_hat(i+1) = 0;
%  
  x_hat(i) = x_hat(i) + h*kz2_hat*(y(i)+ky1);
  kx1_hat = kx1_hat + h*kz2_hat*(y(i)+ky1);
  y_hat(i) = y_hat(i) + h*(x(i)+kx1)*kz2_hat;
  ky1_hat = ky1_hat + h*(x(i)+kx1)*kz2_hat;
  z_hat(i) = z_hat(i) - h*b*kz2_hat;
  kz1_hat = kz1_hat - h*b*kz2_hat;
  kz2_hat = 0.0d0;
%
  x_hat(i) = x_hat(i) + h*r*ky2_hat;
  kx1_hat = kx1_hat + h*r*ky2_hat;
  y_hat(i) = y_hat(i) - h* ky2_hat;
  ky1_hat = ky1_hat - h*ky2_hat;
  x_hat(i) = x_hat(i) - h*ky2_hat*(z(i)+kz1);
  kx1_hat = kx1_hat - h*ky2_hat*(z(i)+kz1);
  z_hat(i) = z_hat(i) - h*(x(i)+kx1)*ky2_hat;
  kz1_hat = kz1_hat - h*(x(i)+kx1)*ky2_hat;
  ky2_hat = 0.0d0;
%
  y_hat(i) = y_hat(i) + h*a*kx2_hat;
  ky1_hat = ky1_hat + h*a*kx2_hat;
  x_hat(i) = x_hat(i) - h*a*kx2_hat;
  kx1_hat = kx1_hat - h*a*kx2_hat;
  kx2_hat = 0.0d0;
%
  x_hat(i) = x_hat(i) + h*kz1_hat*y(i);
  y_hat(i) = y_hat(i) + h*x(i)*kz1_hat;
  z_hat(i) = z_hat(i) - h*b*kz1_hat;
  kz1_hat = 0.0d0;
  x_hat(i) = x_hat(i) + h*r*ky1_hat;
  y_hat(i) = y_hat(i) - h*ky1_hat;
  x_hat(i) = x_hat(i) - h*ky1_hat*z(i);
  z_hat(i) = z_hat(i) - h*x(i)*ky1_hat;
  ky1_hat = 0.0d0;
  y_hat(i) = y_hat(i) + h*a*kx1_hat;
  x_hat(i) = x_hat(i) - h*a*kx1_hat;
  kx1_hat = 0.0d0;
end
xval_hat(j+1) = xval_hat(j+1) + x_hat(1);
x_hat(1) = 0.0d0;
yval_hat(j+1) = yval_hat(j+1) + y_hat(1);
y_hat(1) = 0.0d0;
zval_hat(j+1) = zval_hat(j+1) + z_hat(1);
z_hat(1) = 0.0d0;

