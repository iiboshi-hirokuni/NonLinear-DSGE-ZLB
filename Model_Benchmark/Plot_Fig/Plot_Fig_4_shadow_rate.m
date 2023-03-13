%%  Figure 4  Notional

%%  Plot Variables
%%  'series_R_lag', 'series_mu_a', 'series_z_b', 'series_r_star', ...
%%       'series_r', 'series_y', 'series_pi', 'series_shock_mp', ...
%%      

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');

% addpath('../data')
% addpath('../fun_smc')
% load_data_for_graph

data=csvread('./data_sample/Natural_rate.csv',1,1);

shadow_rate = 100/4*[NaN(7,1); data(:,3)];

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

% end
     
% plot

s =1983.0;
 ti = (s+2/4):0.25:s+(Tobs+1)/4; 

 model_indx = 'model1' 
h_D = figure('Position',[20,20,900,400],'Name','shadow_rate',...
      'Color','w','File',[ num2str(model_indx) '_shadow_rate']);
% subplot(1,1,1)
hold on 
l1 = plot(ti, M1.series_R_lag(1:Tobs,1),'LineStyle','-','Color','b', 'LineWidth',2.0);
l2 = plot(ti, shadow_rate(1:Tobs,1),'LineStyle',':','Color','r', 'LineWidth',2.0); 
plot(ti, zeros(size(ti,2)),'k-'); 
hold off
% ylim([-0.6 1.3])
ylabel(' quarterly, % ')
xlim([1981  2018])
%title('Notional Rate and Shadow Rate')
lgnd = legend([l1, l2],{'R_t^* (Model 1)','Shadow rate (Ueno, 2017)'},'FontSize',12);
set(lgnd, 'Box', 'off')

savefig('./output/Fig4_Shadow_rate.fig');




