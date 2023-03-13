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

load('model0_state.mat');
     M0.series_R_lag = series_R_lag;
     M0.series_mu_a  = series_mu_a;
     M0.series_z_b   = series_z_b;
     M0.series_r_star =series_r_star;
     M0.series_r      =series_r   ;
     M0.series_y      =series_y ;
     M0.series_pi     =series_pi ;
     M0.series_shock_mp = series_shock_mp;
load('model1_state.mat');
     M1.series_R_lag = series_R_lag;
     M1.series_mu_a  = series_mu_a;
     M1.series_z_b   = series_z_b;
     M1.series_r_star =series_r_star;
     M1.series_r      =series_r   ;
     M1.series_y      =series_y ;
     M1.series_pi     =series_pi ;
     M1.series_shock_mp = series_shock_mp;
load('model2_state.mat');
     M2.series_R_lag = series_R_lag;
     M2.series_mu_a  = series_mu_a;
     M2.series_z_b   = series_z_b;
     M2.series_r_star =series_r_star;
     M2.series_r      =series_r   ;
     M2.series_y      =series_y ;
     M2.series_pi     =series_pi ;
     M2.series_shock_mp = series_shock_mp;

% plot

s =1983.0;
 ti = (s+2/4):0.25:s+(Tobs+1)/4; 

h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior',...
      'Color','w','File','Natural_Rate');
% subplot(1,1,1)
hold on 
plot(ti, M1.series_r_star(1:Tobs,1),'LineStyle','-','Color','black', 'LineWidth',2.0); 
plot(ti, M2.series_r_star(1:Tobs,1),'LineStyle','--','Color','r', 'LineWidth',2.0);
plot(ti, M0.series_r_star(1:Tobs,1),'LineStyle',':','Color','b', 'LineWidth',2.0); 
plot(ti, zeros(size(ti,2)),'k-'); 
hold off
% ylim([-0.6 1.3])
ylabel(' quarterly, % ')
xlim([1981  2018])
% title('Natural Rate of Interest')
lgnd = legend({'Model 1', 'Model 2', 'Model w/o ZLB'}, 'FontSize',12);
set(lgnd, 'Box', 'off')

savefig('./output/Fig7_Natural_rate.fig'); 

% h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior','Color','w','File','mu_a');
% % subplot(1,1,1)
% hold on 
% plot(ti, M1.series_mu_a(1:Tobs,1),'LineStyle','-','Color','black', 'LineWidth',2.0); 
% plot(ti, M2.series_mu_a(1:Tobs,1),'LineStyle','--','Color','r', 'LineWidth',2.0);
% plot(ti, M0.series_mu_a(1:Tobs,1),'LineStyle',':','Color','b', 'LineWidth',2.0); 
% plot(ti, zeros(size(ti,2)),'k-'); 
% hold off
% % ylim([-0.6 1.3])
% ylabel(' quarterly, % ')
% xlim([1981  2018])
% title('TFP ')
% lgnd = legend({'Model 1', 'Model 2', 'Model w/o ZLB'},'FontSize',12);
% set(lgnd, 'Box', 'off')
% 
% h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior','Color','w','File','z_b');
% % subplot(1,1,1)
% hold on 
% plot(ti, M1.series_z_b(1:Tobs,1),'LineStyle','-','Color','black', 'LineWidth',2.0); 
% plot(ti, M2.series_z_b(1:Tobs,1),'LineStyle','--','Color','r', 'LineWidth',2.0);
% plot(ti, M0.series_z_b(1:Tobs,1),'LineStyle',':','Color','b', 'LineWidth',2.0); 
% plot(ti, zeros(size(ti,2)),'k-'); 
% hold off
% % ylim([-0.6 1.3])
% ylabel(' quarterly, % ')
% xlim([1981  2018])
% title('Preference ')
% lgnd = legend({'Model 1', 'Model 2', 'Model w/o ZLB'},'FontSize',12);
% set(lgnd, 'Box', 'off')

% %% comparison of Fig 5 and 6
% figure(1)
% plot([M1.series_r_star(1:Tobs,1) 100*log(decomp_t(1:end,6)) ...
%       M1.series_r_star(1:Tobs,1)-100*log(decomp_t(1:end,6))] ,...
%     'LineWidth',2.0   )
% legend('Fig 5', 'Fig 6','Difference')

