%%

%%  Plot Variables
%%  'series_R_lag', 'series_mu_a', 'series_z_b', 'series_r_star', ...
%%       'series_r', 'series_y', 'series_pi', 'series_shock_mp', ...
%%      

clear all

addpath('./data_sample')

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');

data=csvread('Natural_rate_2.csv',1,1);

load('model1_state.mat');
%      M1.series_R_lag = series_R_lag;
%      M1.series_mu_a  = series_mu_a;
%      M1.series_z_b   = series_z_b;
     M1.series_r_star =series_r_star;
%      M1.series_r      =series_r   ;
%      M1.series_y      =series_y ;
%      M1.series_pi     =series_pi ;
%      M1.series_shock_mp = series_shock_mp;

% plot

s =1983.0;
 ti = (s+7/4+2/4):0.25:s+(Tobs+1)/4; 

h_D = figure('Position',[20,20,900,400],'Name','Filtered Posterior',...
      'Color','w','File','Natural_Rate_3');
% subplot(1,1,1)
hold on 
l1=plot(ti, M1.series_r_star(1+7:Tobs,1),'LineStyle','-','Color','black', 'LineWidth',2.0); 
l2=plot(ti, data(1:126,3)/4,'LineStyle','--','Marker', 'o', 'Color','r', 'LineWidth',2.0);
% l3=plot(ti, data(1:126,2)/4,'LineStyle','--','Marker', 'x','Color','b', 'LineWidth',2.0);
 l4=plot(ti, data(1:126,4)/4,'LineStyle',':','Color','b', 'LineWidth',4.0); 
plot(ti, zeros(size(ti,2)),'k-'); 
hold off
% ylim([-0.6 1.3])
ylabel(' quarterly, % ')
xlim([1985  2018])
% title('Natural Rate of Interest')
% lgnd = legend([l1, l2, l4],{'Benchmark (Model 1)', 'Laubach-Williams (one-sided)',...
%     'Laubach-Williams (two-sided)', 'HP Filter'}, 'FontSize',12);
lgnd = legend([l1, l2, l4],{'Benchmark (Model 1)', 'Laubach-Williams (one-sided)',...
                            'HP Filter'}, 'FontSize',12);
set(lgnd, 'Box', 'off')

savefig('./output/Fig9_Natural_rate.fig'); 

%csvwrite('rate','')



