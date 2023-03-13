% Plot contour for the policy functions of output gap and inflation
% This is aimed to check whether we can identrify a and b shocks

% clear all
% clc
% 
% warning('off','all'); % åxçêÇÃîÒï\é¶
% 
% disp('Start NK Model with Time Iteration (Linear Interpolation)')
% addpath('./Toolbox')
% addpath('./fun_IRF')
% 
% tstart = tic;   

%% setting conditions 
% zlbflag = 1; policy_flag = 1; 
% script_TL

%%
%{
if strcmp(type,'linear')
    a = 20; b=25;
else
    a =0; b=5;
    if zlbflag == 1
        a = 10; b=15;
    end
end
%}
        
r = floor(mean(size(G.r_grid,2))); % r(-1) is the largest???
%r = 2; % round(min(size(G.r_grid,1))); % at or close to the ZLB???
MP = floor(size(G.MP_grid,2)/2); % monetary shock is zero???
%MP = 1;

Py = squeeze(pf.y(r,:,:,MP));
Ppi = squeeze(pf.pi(r,:,:,MP));
Py_star = squeeze(pf.y_star(r,:,:,MP));
Pr = squeeze(pf.r(r,:,:,MP));

x = G.a_grid;
y = G.b_grid;

[X,Y]=meshgrid(x,y);

 fig_fp_y = figure('Position',[20,20,600,600],'Name','policy function','Color','w');
 subplot(2,2,1)
 %surf(X,Y,Py-Py_star)
 zmin = (min(Py(:)-Py_star(:)));
zmax = (max(Py(:)-Py_star(:)));
zinc = (zmax - zmin) / 20;
zlevs = zmin:zinc:zmax;

 contour(X,Y,Py-Py_star,zlevs,'ShowText','on')
%  contour(X,Y,Py-Py_star,'ShowText','on')
  xlabel('log(A_{t}/A_{t-1})-\gamma_a')
   ylabel('log(Z^b_{t})')     
 %  zlabel('Y_t-Y^*_t')
 colormap jet %copper; %hsv; %parula
 title('policy function of output gap','FontSize',14)
 %  view(3)
 %view([-65,45])
 
 est_date = datestr(date);   
%  name = ['../result/','fig_pf_y_fp_',est_date];
%          saveas(fig_fp_y,name,'fig')
            
 
 subplot(2,2,2)
 %surf(X,Y,Ppi)
  zmin = (min(Ppi(:)));
zmax = (max(Ppi(:)));
zinc = (zmax - zmin) / 20;
zlevs = zmin:zinc:zmax;

contour(X,Y,Ppi,zlevs,'ShowText','on')
xlabel('log(A_{t}/A_{t-1})-\gamma_a')
   ylabel('log(Z^b_{t})')     
%   zlabel('Inflation')
 colormap jet
 title('policy function of inflation','FontSize',14)
%  title('policy function of Y^*_t/A_t ','FontSize',14)
 %  view(3)
% view([-65,45])
  
 est_date = datestr(date);   
%  name = ['../result/','fig_pf_pi_fp_',est_date];
%          saveas(fig_fp_pi,name,'fig')

subplot(2,2,3)
 zmin = (min(Pr(:)));
zmax = (max(Pr(:)));
zinc = (zmax - zmin) / 20;
zlevs = zmin:zinc:zmax;

contour(X,Y,Pr,zlevs,'ShowText','on')
xlabel('log(A_{t}/A_{t-1})-\gamma_a')
   ylabel('log(Z^b_{t})')    
  colormap jet %winter
 title('policy function of interest rate ','FontSize',14)
 
%  view(3)
% view([-65,45])
   
 est_date = datestr(date);   
%  name = ['../result/','fig_pf_r_fp_',est_date];
%          saveas(fig_fp_r,name,'fig')
