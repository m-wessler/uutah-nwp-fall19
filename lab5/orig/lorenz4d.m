%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program to perform 4D-Var data assimilation on the
% Lorenz equations
%
%  List of main variables
%    a:          sigma coefficient in equations
%    r:          rho coefficient in equations
%    b:          beta coefficient in equations
%
%    ta:         Length of assimilation window
%    tf:         Length of forecast window
%    h:          Time step for numerical scheme
%    freq:       Frequency of observations in time steps
%    tstep:      Number of time steps for assimilation
%    fcstep:     Total number of time steps for assimilation + forecast
%
%    [xx,yy,zz]: Truth values of trajectory
%    [xf,yf,zf]: First guess values of trajectory
%    [x,y,z]:    Calculated values of trajectory and final analysis
%    [lx,ly,lz]: Gradient values calculated from adjoint
%
%    [datx,daty,datz]: Observation values
%    D:                Observation weighting matrix
%    max_iterations:   Maximum number of minimization iterations
%    tolerance:        Convergence tolerance for minimization
%    s:                Calculated step length in minimization
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------
%   1. Set up
%---------------------------------------------------------------------
% Close any open figure windows
%g=gcf;
% for i=1:g
%   close(i)
% end
% Determine screen size and set up figure positions
clear all
clc
close all
ss=get(0,'ScreenSize');
fig_width=0.55*ss(3);
fig_height=0.4*ss(4);
pos_1=[0.02*ss(3),0.55*ss(4),fig_width,fig_height];
pos_2=[0.02*ss(3),0.05*ss(4),fig_width,fig_height];
pos_3=[0.6*ss(3),0.45*ss(4),0.4*ss(3),0.45*ss(4)];
pl=0.6*ss(3);
pb=0.05*ss(4);
pw=0.3*ss(3);
ph=0.3*ss(4);
pos_4=[pl,pb,pw,ph];
pos_5=[0.02*pw,0.8*ph,0.95*pw,0.2*ph];
pos_6=[0.02*pw,0.02*ph,0.95*pw,0.75*ph];
%
format long
%parameters
a=1.0d1;       % sigma
b=8.0d0/3.0d0; % rho
r=2.8d1;       % beta
%

%---------------------------------------------------------------------
%  2. Input truth parameters and calculate truth
%---------------------------------------------------------------------
% Inputs 1
xstr='True value of x at t=0';
ystr='True value of y at t=0';
zstr='True value of z at t=0';
true_str=inputdlg({xstr,ystr,zstr},'Initial data');
truex=str2num(true_str{1});
truey=str2num(true_str{2});
truez=str2num(true_str{3});

% Inputs 2
xstr='Length of assimilation window';
x2str='Length of subsequent forecast';
ystr='Time step';
zstr='Frequency of observations (in time steps)';
temp_str=inputdlg({xstr,x2str,ystr,zstr},'Run information');
ta=str2num(temp_str{1});
tf=str2num(temp_str{2});
h=str2num(temp_str{3});
freq=str2num(temp_str{4});

tstep = ta/h+1;
fcstep = (ta+tf)/h;

% Run truth
guessx(1)=truex;
guessy(1)=truey;
guessz(1)=truez;
j=0;
[xx,yy,zz]=modeuler(fcstep,h,j,guessz,guessy,guessx,a,r,b);

%---------------------------------------------------------------------
% 3. Calculate observations and plot truth and observations
%---------------------------------------------------------------------
D=zeros(tstep,1);
n_obs = (tstep-1) / freq + 1;
datx=zeros(n_obs,1);
daty=zeros(n_obs,1);
datz=zeros(n_obs,1);
l_noise=menu_asl('Noise on observations?','No','Yes');
l_noise=l_noise-1;
if l_noise==1
  sc_x_noise=randn(n_obs,1);
  sc_y_noise=randn(n_obs,1);
  sc_z_noise=randn(n_obs,1);
  sd_str=inputdlg('Variance of observation error');
  sd=str2num(sd_str{1});
  var=sqrt(sd);
  RX=var*sc_x_noise;
  RY=var*sc_y_noise;
  RZ=var*sc_z_noise;
else
  for i=1:n_obs;
    RX(i)=0.0;
    RY(i)=0.0;
    RZ(i)=0.0;
  end
end
% Set up data and matrix D
j=1;
for i=1:freq:tstep
  datx(i)=xx(i)+RX(j);
  daty(i)=yy(i)+RY(j);
  datz(i)=zz(i)+RZ(j);
  D(i) = 1.;
  j=j+1;
end    
% Plot truth and observations
xvals=0:fcstep;
vec=1:freq:tstep;
v=vec-1;
for i=vec
  x_ob(i)=datx(i);
  z_ob(i)=datz(i);
end
% x_ob and z_ob are temporary arrays for plotting
%
h1=figure('Position',pos_1);
clf;
subplot(2,1,1)
plot(xvals,xx)
hold on
plot(v,x_ob(vec),'om')
hold on
xlabel('Time step')
ylabel('x')
title('Solution for x')
legend('Truth','Observations')
%
subplot(2,1,2)
plot(xvals,zz)
hold on
plot(v,z_ob(vec),'om')
hold on
xlabel('Time step')
ylabel('z')
title('Solution for z')
legend('Truth','Observations')

%---------------------------------------------------------------------
% 4. Input initial guess and run model from this
%---------------------------------------------------------------------
xstr='Initial guess of x at t=0';
ystr='Initial guess of y at t=0';
zstr='Initial guess of z at t=0';
guess_str=inputdlg({xstr,ystr,zstr},'Initial guess');
guessx(2)=str2num(guess_str{1});
guessy(2)=str2num(guess_str{2});
guessz(2)=str2num(guess_str{3});
%
j=1;  
[x,y,z]=modeuler(tstep,h,j,guessz,guessy,guessx,a,r,b);
  [f,g] = calcfg(tstep,h,x,y,z,datx,daty,datz,a,r,b,D,freq);
  lx=g(1);
  ly=g(2);
  lz=g(3);
% Plot first guess
[xf,yf,zf]=modeuler(fcstep,h,j,guessz,guessy,guessx,a,r,b);
subplot(2,1,1)
plot(xvals,xf,':')
legend('Truth','Observations','First guess')
hold on
subplot(2,1,2)
plot(xvals,zf,':')
legend('Truth','Observations','First guess')
hold on

%---------------------------------------------------------------------
% 5. Input minimization info and perform minimization
%---------------------------------------------------------------------
xstr='Maximum number of iterations';
ystr='Tolerance';
zstr={'30','1d-5'};
min_str=inputdlg({xstr,ystr},'Minimization control',1,zstr);
max_iterations=str2num(min_str{1});
tolerance=str2num(min_str{2});

j=2;
lnorm(j)=sqrt(lx(1)*lx(1)+ly(1)*ly(1)+lz(1)*lz(1));
cost(j)=f;
s(j)=0.5d0;

%     minimisation via least squares
l_converged=0;
if lnorm(j)<tolerance
  l_converged=1
end

while lnorm(j)>tolerance & l_converged==0  % 70
  s(j)=0.5;
% normalise gradient vector
  gx=lx(1)/lnorm(j);
  gy=ly(1)/lnorm(j);
  gz=lz(1)/lnorm(j);
       
  guessx(j+1)=guessx(j)-s(j)*gx;
  guessy(j+1)=guessy(j)-s(j)*gy;
  guessz(j+1)=guessz(j)-s(j)*gz;

  [x,y,z]=modeuler(tstep,h,j,guessz,guessy,guessx,a,r,b);
  [f,g] = calcfg(tstep,h,x,y,z,datx,daty,datz,a,r,b,D,freq);
  lx=g(1);
  ly=g(2);
  lz=g(3);

  lnorm(j+1)=sqrt(lx(1)*lx(1)+ly(1)*ly(1)+lz(1)*lz(1));
  cost(j+1)=f;
      
  if (j-1>max_iterations)  % Max number of iterations reached
    l_converged=1;
    disp('Failed to converge in maximum number of iterations')
    Norm_of_gradient=lnorm(j)
  else
    Iteration=j-1
    while (lnorm(j+1) >= lnorm(j) & s(j) > tolerance ...
                       & l_converged == 0) % 75
      s(j)=s(j)*0.5d0;
      guessx(j+1)=guessx(j)-s(j)*gx;
      guessy(j+1)=guessy(j)-s(j)*gy;
      guessz(j+1)=guessz(j)-s(j)*gz;

      [x,y,z]=modeuler(tstep,h,j,guessz,guessy,guessx,a,r,b);
  [f,g] = calcfg(tstep,h,x,y,z,datx,daty,datz,a,r,b,D,freq);
  lx=g(1);
  ly=g(2);
  lz=g(3);

      lnorm(j+1)=sqrt(lx(1)*lx(1)+ly(1)*ly(1)+lz(1)*lz(1));
      cost(j+1)=f;
% Convergence test on norm of gradient
      if (abs((lnorm(j+1)-lnorm(j)))<tolerance)
        disp('Converged on change in gradient')
        number_of_iterations=j-1
        Norm_of_gradient=lnorm(j+1)
        l_converged=1;
      end
    end
    j=j+1;
  end
end

%---------------------------------------------------------------------
% 6. Run forecast from final analysis and plot
%---------------------------------------------------------------------
jtemp=0;
guessx(1)=x(1);
guessy(1)=y(1);
guessz(1)=z(1);
          [x,y,z]=modeuler(fcstep,h,jtemp,guessz,guessy,guessx,a,r,b);

subplot(2,1,1)
plot(xvals,x,'--')
legend('Truth','Observations','First guess','Analysis')
hold on
% Plot vertical line at start of forecast
yminmax=get(gca,'YLim');
yspace=(yminmax(2)-yminmax(1))*0.01;
yvals=yminmax(1):yspace:yminmax(2);
line(tstep-1,yvals,'LineStyle',':','Color','k')
hold off
%
subplot(2,1,2)
plot(xvals,z,'--')
legend('Truth','Observations','First guess','Analysis')
hold on
% Plot vertical line at start of forecast
yminmax=get(gca,'YLim');
yspace=(yminmax(2)-yminmax(1))*0.01;
yvals=yminmax(1):yspace:yminmax(2);
line(tstep-1,yvals,'LineStyle',':','Color','k')
hold off

%---------------------------------------------------------------------
%  7. Calculate and plot error
%---------------------------------------------------------------------
x_error=xx-x;
z_error=zz-z;
h2=figure('Position',pos_2);
clf;
subplot(2,1,1)
plot(xvals,x_error,'k-');
xlabel('Time step');ylabel('Error')
title('Plot of error (Truth-Analysis) in x against time')
hold on
% Plot vertical line at start of forecast
yminmax=get(gca,'YLim');
yspace=(yminmax(2)-yminmax(1))*0.01;
yvals=yminmax(1):yspace:yminmax(2);
line(tstep-1,yvals,'LineStyle',':','Color','k')
%
hold on
subplot(2,1,2)
plot(xvals,z_error,'k-')
xlabel('Time step');ylabel('Error')
title('Plot of error (Truth-Analysis) in z against time')
% Plot vertical line at start of forecast
yminmax=get(gca,'YLim');
yspace=(yminmax(2)-yminmax(1))*0.01;
yvals=yminmax(1):yspace:yminmax(2);
line(tstep-1,yvals,'LineStyle',':','Color','k')

%-----------------------------------------------------------------------
%  8. Plot convergence
%---------------------------------------------------------------------
h3=figure('Position',pos_3);
clf;
xvals=0:j-2;
subplot(2,1,1)
semilogy(xvals,cost(2:j))
title('Convergence of cost function')
xlabel('Iteration')
ylabel('Cost function')
subplot(2,1,2)
semilogy(xvals,lnorm(2:j))
title('Convergence of gradient')
xlabel('Iteration')
ylabel('Norm of gradient')
%-----------------------------------------------------------------------
% 9. Write out options chosen
%---------------------------------------------------------------------
h4=figure('Position',pos_4);
clf;
text1={'List of options chosen'};
text2={['True (x,y,z) at t=0:   (' true_str{1} ',' true_str{2} ',' ...
 true_str{3}  ')']};
text3={['First guess (x,y,z) at t=0:   (' guess_str{1} ',' guess_str{2} ',' ...
 guess_str{3}  ')']};
text4={['Length of assimilation window: ' temp_str{1}]};
text5={['Length of subsequent forecast: ' temp_str{2}]};
text6={['Time step: ' temp_str{3}]};
text7={['Frequency of observations = ' temp_str{4}]};
text9={['Maximum iterations: ' min_str{1}]};
text10={['Tolerance: ' min_str{2}]};
if l_noise==0
  text11={['No noise on observations']};
else
  text11={['Observations have random noise with variance ' num2str(sd)]};
end

str1=[text2;text3;text4;text5;text6;text7;text9;text10;text11];
uicontrol('Style','text','Position',pos_5,'String',text1,...
 'FontSize',14,'FontWeight','bold')
uicontrol('Style','text','Position',pos_6,'String',str1,...
 'FontSize',12,'HorizontalAlignment','left')
