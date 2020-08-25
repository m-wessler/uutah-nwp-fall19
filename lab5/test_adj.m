function test_adj
format long
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: To test adjoint code of modeuler.m
%
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1.0d1;       % sigma
b=8.0d0/3.0d0; % rho
r=2.8d1;       % beta
h=0.01d0;        % time step
tstep=10;     % Number of steps
%
x=zeros(tstep+1,1);
y=zeros(tstep+1,1);
z=zeros(tstep+1,1);
x(1)=1.0d0;
y(1)=2.0d0;
z(1)=1.5d0;
%
pert=rand(1,3);

[x,y,z]=modeuler(tstep,h,0,z,y,x,a,r,b);
x_p(1)=pert(1);
y_p(1)=pert(2);
z_p(1)=pert(3);
[x_p,y_p,z_p]=modeuler_tl(tstep,h,0,z,y,x,z_p,y_p,x_p,a,r,b);
%
x_hat(tstep+1,1)=x_p(tstep+1);
y_hat(tstep+1,1)=y_p(tstep+1);
z_hat(tstep+1,1)=z_p(tstep+1);
[x_hat,y_hat,z_hat]=modeuler_adj(tstep,h,0,z,y,x,...
   x_hat(tstep+1,1),y_hat(tstep+1,1),z_hat(tstep+1,1),a,r,b);
%
Forward_prod = x_p(tstep+1)*x_p(tstep+1) ...
  + y_p(tstep+1)*y_p(tstep+1) + z_p(tstep+1)*z_p(tstep+1)
Adj_prod = x_p(1)*x_hat(1) + y_p(1)*y_hat(1) + z_p(1) * z_hat(1)
