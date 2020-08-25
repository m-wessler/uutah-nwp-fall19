clear all
close all
                                                                                                        
% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
% c
% c     Program to perform a simple Cressman analysis
% c
% c     Z.-X. Pu 09/16/19 for NWP Lab Assignment #3
% c
% c
% cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      n_pts=100;n_obs=5
%       parameter(pi=3.1415926)
%       n_obs_position(1:n_obs)=0;
      xt(1:n_pts)=0;xb(1:n_pts)=0;xa(1:n_pts)=0;w(1:n_pts,1:n_obs)=0;obs(1:n_obs)=0;
%      n_obs_position=[3,15,18,21,43,44,45,74,75,76,77,78,81,82,95];
% n_obs_position=[3,15,21,27,33,39,45,54,60,67,77,83,89,92,96];
     n_obs_position= [3, 5, 7, 9, 11];
      n_obs_sign=1;
      n_obs_big=5;

% c     Generate truth and background
      delta_x=2*pi/(n_pts-1);
        for i=1:n_pts
         xt(i)=sin((i-1)*delta_x);
         xb(i)=0.5*xt(i);
        end

      display('true solution')
      display(xt')
   
      display('background')
      display(xb')
      
% c     Make obs.
      for i=1:n_obs
       obs(i)=xt(n_obs_position(i));
      end 
 
% c     Make some obs incorrect
      obs(n_obs_sign)=-obs(n_obs_sign);
      obs(n_obs_big)=obs(n_obs_big)*1.5;

      display('observation')
      display(obs')

% c     Decide on weighting function
      lw=input('Which weighting function? (1=Square, 2=Exponential)');

     
% c     Set up weighting function
      R=input('Radius of influnce?');
      Rsq=R^2;

% c     Square weighting function
      if lw==1
        for k=1:n_obs
          for i=1:n_pts
           d=abs(n_obs_position(k)-i);
           dsq=d.^2;
           w(i,k)=(Rsq-dsq)./(Rsq+dsq);
             if w(i,k)<0
              w(i,k) = 0.;
             end
          end  
         end  
      end
 
% c     Exponential weighting function
       if lw==2
         for k=1:n_obs
          for i=1:n_pts
            d=abs(n_obs_position(k)-i);
            dsq=d.^2;
            w(i,k)= exp(-dsq./(2*Rsq));
          end
         end  
       end

% c     Perform analysis
      for i=1:n_pts
       d_x = 0.;
       scalar = 0.;
        for k=1:n_obs
          i_ob =n_obs_position(k);
          d_x = d_x +w(i,k) *(obs(k)-xb(i_ob))
          scalar=scalar + w (i,k)
        end
       if scalar==0
        d_x=0;
       else
        d_x=d_x /scalar;
       end
        xa(i)=xb(i)+d_x;
      end

      display('analysis')
      display(xa')

      
%     figure 
      plot(1:n_pts,xa,'LineWidth',2)
      hold on;
      plot(1:n_pts,xt,'o')
      hold on; 
      plot(1:n_pts,xb,'--')
      xlabel('grid points')
      ylabel('x')
            

