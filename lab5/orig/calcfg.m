function [f,g] = calcfg(tstep,h,x,y,z,datx,daty,datz,a,r,b,D,freq)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: Calculate cost function and its gradient
%           
%
%  %  List of main variables
%    a:          sigma coefficient in equations
%    r:          rho coefficient in equations
%    b:          beta coefficient in equations
%    D:          Observation weighting matrix
%    h:          Time step for numerical scheme
%    freq:       Frequency of observations
%    tstep:      Number of time steps to perform
%    [x,y,z]:    Forward trajectory
%    [datx,daty,datz]: Observation values
%
%  Output:
%    [f,g]: Cost function and gradient
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=0.0d0;
%
% Calculate cost function
for i=1:freq:tstep
  f = f + 0.5 * ( (x(i)-datx(i))*(x(i)-datx(i))*D(i) ...
        + (y(i)-daty(i))*(y(i)-daty(i))*D(i) ...
        + (z(i)-datz(i))*(z(i)-datz(i))*D(i) );
end
%
% Calculate gradient of cost function
x_hat=zeros(tstep,1);
y_hat=zeros(tstep,1);
z_hat=zeros(tstep,1);
%
x_hat(tstep) = (x(tstep)-datx(tstep))*D(tstep);
y_hat(tstep) = (y(tstep)-daty(tstep))*D(tstep);
z_hat(tstep) = (z(tstep)-datz(tstep))*D(tstep);
%
for i=tstep-1:-1:1
  [x1_hat,y1_hat,z1_hat] = modeuler_adj(1,h,0,z(i),y(i),x(i),...
               x_hat(i+1),y_hat(i+1),z_hat(i+1),a,r,b);
  x_hat(i) = x1_hat(1) + (x(i)-datx(i))*D(i);
  y_hat(i) = y1_hat(1) + (y(i)-daty(i))*D(i);
  z_hat(i) = z1_hat(1) + (z(i)-datz(i))*D(i);
end
g=[x_hat(1),y_hat(1),z_hat(1)];

