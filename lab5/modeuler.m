function [xo,yo,zo]=modeuler(tstep,h,j,zval,yval,xval,a,r,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: Apply modified Euler scheme to the Lorenz equations
% %
%  6/7/04 Moved all functions inline to make adjoint easier
%
%  List of main variables
%    a:          sigma coefficient in equations
%    r:          rho coefficient in equations
%    b:          beta coefficient in equations
%    h:          Time step for numerical scheme
%    tstep:      Number of time steps to perform
%    [xval,yval,zval]: Initial fields
%    j:          Index to pick up correct initial field
%
%  Output:
%    [xo,yo,zo]: Trajectories of evolved fields
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=zeros(tstep,1);
y=zeros(tstep,1);
z=zeros(tstep,1);
%
x(1)=xval(j+1);
y(1)=yval(j+1);
z(1)=zval(j+1);
%
for i = 1:tstep
  kx1 = h*a*(y(i)-x(i));
  ky1 = h*(r*x(i)-y(i)-x(i)*z(i));
  kz1 = h*(x(i)*y(i) - b*z(i));

  kx2 = h*a*(y(i)+ky1-x(i)-kx1);
  ky2 = h*(r*(x(i)+kx1)-y(i)-ky1-(x(i)+kx1)*(z(i)+kz1));
  kz2 = h*((x(i)+kx1)*(y(i)+ky1)-b*(z(i)+kz1));

  x(i+1)=x(i)+0.5d0*(kx1+kx2);
  y(i+1)=y(i)+0.5d0*(ky1+ky2);
  z(i+1)=z(i)+0.5d0*(kz1+kz2);
end
xo=x;
yo=y;
zo=z;
%

