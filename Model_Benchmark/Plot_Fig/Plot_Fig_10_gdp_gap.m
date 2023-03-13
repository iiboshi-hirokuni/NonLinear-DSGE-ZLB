%%

%%  Plot Variables
%%  'series_R_lag', 'series_mu_a', 'series_z_b', 'series_r_star', ...
%%       'series_r', 'series_y', 'series_pi', 'series_shock_mp', ...
%%      

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');

addpath('./data_sample')

load('model1_state.mat');
     M1.series_R_lag = series_R_lag;
     M1.series_mu_a  = series_mu_a;
     M1.series_z_b   = series_z_b;
     M1.series_r_star =series_r_star;
     M1.series_r      =series_r   ;
     M1.series_y      =series_y ;
     M1.series_pi     =series_pi ;
     M1.series_shock_mp = series_shock_mp;   
       
              
 s =1983.0;
 nobs = Tobs;
 ti_def = (s+2/4):0.25:s+(nobs+1)/4;
 
 
load('model1_I1__state.mat');
     M2.series_R_lag = series_R_lag;
     M2.series_mu_a  = series_mu_a;
     M2.series_z_b   = series_z_b;
     M2.series_r_star =series_r_star;
     M2.series_r      =series_r   ;
     M2.series_y      =series_y ;
     M2.series_pi     =series_pi ;
     M2.series_shock_mp = series_shock_mp;
     
load('model1_I1_state.mat');
     M5.series_R_lag = series_R_lag;
     M5.series_mu_a  = series_mu_a;
     M5.series_z_b   = series_z_b;
     M5.series_r_star =series_r_star;
     M5.series_r      =series_r   ;
     M5.series_y      =series_y ;
     M5.series_pi     =series_pi ;
%      M5.series_shock_mp = series_shock_mp;   
 
load('model1_V4__state.mat');
     M3.series_R_lag = series_R_lag;
     M3.series_mu_a  = series_mu_a;
     M3.series_z_b   = series_z_b;
     M3.series_r_star =series_r_star;
     M3.series_r      =series_r   ;
     M3.series_y      =series_y ;
     M3.series_pi     =series_pi ;
     M3.series_shock_mp = series_shock_mp;   
load('model1_gap__state.mat');
     M4.series_R_lag = series_R_lag;
     M4series_mu_a  = series_mu_a;
     M4.series_z_b   = series_z_b;
     M4.series_r_star =series_r_star;
     M4.series_r      =series_r   ;
     M4.series_y      =series_y ;
     M4.series_pi     =series_pi ;
     M4.series_shock_mp = series_shock_mp;       

% plot


h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior',...
             'Color','w','File','Natural_Rate_2');
% subplot(1,1,1)
hold on 
plot(ti_def, M1.series_r_star(1:nobs,1),'LineStyle','-','Color','black', 'LineWidth',2.0); 
plot(ti_def, M3.series_r_star(1:nobs,1),'LineStyle',':','Color','r', 'LineWidth',2.5); 
plot(ti_def, M4.series_r_star(2:Tobs,1),'LineStyle','--','Color','b', 'LineWidth',2.0); 
% plot(ti_def, M2.series_r_star(1:nobs,1),'LineStyle',':','Color','r', 'LineWidth',2.0);
% plot(ti_def, M5.series_r_star(1:nobs,1),'LineStyle','-','Color','c', 'LineWidth',2.0);
plot(ti_def, zeros(size(ti_def,2)),'k-'); 
hold off
% ylim([-1.2 1.3])
ylabel(' quarterly, % ')
xlim([1981  2018])
% title('Natural Rate of Interest')
% lgnd = legend('Benchmark', 'Gap & Growth ','Gap ','Growth & Z_b = I(1)'  );
lgnd = legend({'Benchmark', 'Gap & Growth ','Gap '},'FontSize',12);
set(lgnd, 'Box', 'off')

savefig('./output/Fig10_Natural_rate.fig'); 
