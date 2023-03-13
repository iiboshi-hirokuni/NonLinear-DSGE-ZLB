%%  Figure 3  Notional

%%  Plot Variables
%%  'series_R_lag', 'series_mu_a', 'series_z_b', 'series_r_star', ...
%%       'series_r', 'series_y', 'series_pi', 'series_shock_mp', ...
%%      

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');

addpath('../data')
addpath('../fun_smc')
addpath('./data_sample')

load_data_for_graph


% model_indx = 'model1'
% model_indx = 'model2'
% model_indx = 'model0'

% switch model_indx 
%  case 'model0'
load('model0_state.mat');
     M0.series_R_lag = series_R_lag;
     M0.series_mu_a  = series_mu_a;
     M0.series_z_b   = series_z_b;
     M0.series_r_star =series_r_star;
     M0.series_r      =series_r   ;
     M0.series_y      =series_y ;
     M0.series_pi     =series_pi ;
     M0.series_shock_mp = series_shock_mp;
%   case 'model1'
load('model1_state.mat');
     M1.series_R_lag = series_R_lag;
     M1.series_mu_a  = series_mu_a;
     M1.series_z_b   = series_z_b;
     M1.series_r_star =series_r_star;
     M1.series_r      =series_r   ;
     M1.series_y      =series_y ;
     M1.series_pi     =series_pi ;
     M1.series_shock_mp = series_shock_mp;
%  case 'model2'
load('model2_state.mat');
     M2.series_R_lag = series_R_lag;
     M2.series_mu_a  = series_mu_a;
     M2.series_z_b   = series_z_b;
     M2.series_r_star =series_r_star;
     M2.series_r      =series_r   ;
     M2.series_y      =series_y ;
     M2.series_pi     =series_pi ;
     M2.series_shock_mp = series_shock_mp;
% end
     
% plot

s =1983.0;
 ti = (s+2/4):0.25:s+(Tobs+1)/4; 

% 

h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior','Color','w',...
    'File',['Notional_Interest_rate']);
% subplot(1,1,1)
hold on 
plot(ti, yy(1:Tobs,3),'LineStyle','-','Color','b', 'LineWidth',2.0); 
plot(ti, M1.series_R_lag(1:Tobs,1),'LineStyle','--','Color','r', 'LineWidth',2.0);
plot(ti, M2.series_R_lag(1:Tobs,1),'LineStyle',':','Color','k', 'LineWidth',1.5);
plot(ti, zeros(size(ti,2)),'k-'); 
hold off
% ylim([-0.6 1.3])
ylabel(' quarterly, % ')
xlim([1981  2018])
% title('Interest Rate: R_t and Notional Interest Rate: R^*_t')
% lgnd = legend({'Interest Rate: R_t', 'Notional Interest Rate: R^*_t (Model 1)', ...
%                'Notional Interest Rate: R^*_t (Model 2)'},'FontSize',12);
lgnd = legend({'R_t', 'R^*_t (Model 1)', ...
               'R^*_t (Model 2)'},'FontSize',12);           
set(lgnd, 'Box', 'off')

savefig('./output/Fig3_Notional_Interest_rate.fig'); 
