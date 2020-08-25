function test_tl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Purpose: To test tangent linear code of modeuler.m
%
%  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1.0d1;       % sigma
b=8.0d0/3.0d0; % rho
r=2.8d1;       % beta
h=0.01d0;        % time step
nsteps=100;     % Number of steps
pert=[1.0d0,-1.0d0,0.5d0];
%
x=zeros(nsteps+1,1);
y=zeros(nsteps+1,1);
z=zeros(nsteps+1,1);
%
gamma = 1.0d3;
for i=1:9
gamma = gamma * 0.1d0;
x(1)=1.0d0;
y(1)=2.0d0;
z(1)=1.5d0;
x_p(1)=gamma*pert(1);
y_p(1)=gamma*pert(2);
z_p(1)=gamma*pert(3);
x2(1)=x(1) + x_p(1);
y2(1)=y(1) + y_p(1);
z2(1)=z(1) + z_p(1);
%
[x,y,z]=modeuler(nsteps,h,0,z,y,x,a,r,b);
[x2,y2,z2]=modeuler(nsteps,h,0,z2,y2,x2,a,r,b);
%
[x_p,y_p,z_p]=modeuler_tl(nsteps,h,0,z,y,x,z_p,y_p,x_p,a,r,b);
%
nl_pert=[x2-x,y2-y,z2-z];
lin_pert=[x_p,y_p,z_p];
%
nl_final = nl_pert(nsteps,1:3);
lin_final = lin_pert(nsteps,1:3);
rel_err(i) = 100d0 * norm(nl_final-lin_final) / norm(lin_final);
gamma_vec(i) = gamma;
end

loglog(gamma_vec,rel_err)
xlabel('gamma');
ylabel('Relative error');
