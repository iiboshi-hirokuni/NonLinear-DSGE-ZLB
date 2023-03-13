%  
% 
% 
% Figure 6: Natural Rate of Interest and the Contribution of the Estimated Shocks
% 
% 
% 

clear all

 addpath('../fun_hist_decomp')
 addpath('../fun_smc')
 addpath('../Toolbox')

load_data_for_graph;

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');

 load('model1_state.mat')
 
 d.series_mu_a = series_mu_a;
 d.series_b  = series_z_b;
 d.series_R = series_r;
 d.series_shock_mp = series_shock_mp;
 
 
T = 133; % # of observations
decomp_t =zeros(T,8);
decomp_mp =zeros(T,8);

shock_mu_a = zeros(T,1);
shock_z_b = zeros(T,1);

for i = 2:T
 shock_mu_a(i) = series_mu_a(i) / series_mu_a(i-1)^P.rho_a;
 shock_z_b(i) = series_z_b(i) / series_z_b(i-1)^P.rho_b;
end 


%% Calculate baseline
 for t = 1:Tobs
      state = zeros(8,1);
      state(4,1) = series_R_lag(t)/100;
      state(5,1) = series_R_lag(t)/100;
      state(7,1) = series_mu_a(t)/100;
      state(8,1) = series_z_b(t)/100;
      shock = zeros(3,1);    
      shock(3,1) = series_shock_mp(t)/100;        
      decomp_t(t,:) = gen_hist_decomp(G, S, P, pf, state, shock, zlbflag,policy_flag); 
 end  

 
 %% only TFP shock
 state = zeros(8,1);
  state(4,1) = series_R_lag(1)/100;
  state(5,1) = series_R_lag(1)/100;  
 for t = 1:Tobs    
      state(7,1) = series_mu_a(t)/100;
%       state(8,1) = d.series_b_wos(t)/100;
      state(8,1) = 0;
      shock = zeros(3,1);    
%       shock(3,1) = series_shock_mp(t)/100;           
      decomp_a(t,:) = gen_hist_decomp(G, S, P, pf, state, shock, zlbflag,policy_flag);
      state(4,1) = decomp_a(t,4)/100;
      state(5,1) = decomp_a(t,5)/100;
 end  

 %% only pref shock
 state = zeros(8,1);
  state(4,1) = series_R_lag(1)/100;
  state(5,1) = series_R_lag(1)/100; 
 for t = 1:Tobs
%       state(7,1) = d.series_mu_wos(t)/100;
      state(7,1) = 0;
      state(8,1) = series_z_b(t)/100;
      shock = zeros(3,1);    
%       shock(3,1) = series_shock_mp(t)/100;           
      decomp_b(t,:) = gen_hist_decomp(G, S, P, pf, state, shock, zlbflag,policy_flag); 
      state(4,1) = decomp_b(t,4)/100;
      state(5,1) = decomp_b(t,5)/100;
 end  
 
 %% only MP shock
 state = zeros(8,1);
  state(4,1) = series_R_lag(1)/100;
  state(5,1) = series_R_lag(1)/100; 
 for t = 1:Tobs       
%       state(7,1) = d.series_mu_wos(t)/100;
%       state(8,1) = d.series_b_wos(t)/100;
      state(7,1) = 0;
      state(8,1) = 0;      
      shock = zeros(3,1);    
      shock(3,1) = series_shock_mp(t)/100;           
      decomp_mp(t,:) = gen_hist_decomp(G, S, P, pf, state, shock, zlbflag,policy_flag); 
      state(4,1) = decomp_mp(t,4)/100;
      state(5,1) = decomp_mp(t,5)/100;
 end  
 
 
%%  Plot Graph 
name = { 'Output Growth' 'Y_star' 'Inflation'...
         'Interest Rate' 'Notional Interest Rate'...
         'Natural Rate of Interest'  'Mu_a' 'Z_b'...
         'i.i.d. Shocks'};
     
file_name = char('Fig_y', 'Fig_Ystar', 'Fig_pi', 'Fig_r', ...
                 'Fig_Rstar', 'Fig_NRI', 'Fig_ma','Fig_zb' ); 

s =1983.0;
% s =1985.25;
 ti = (s+1/4):0.25:s+(Tobs)/4; 
 
 
 %% Natural Rate of Interest
i=6;
figure('Position',[20,20,1250,450],...
       'Name','contributions of shocks','File','Natural_rate_Decomp' )
   
 subplot(1,2,1)  
% plot(ti(1:end), 100*log(decomp_t(1:end,i)),'g-o','LineWidth',1.5);
plot(ti(1:end), series_r_star(:),'k-o','LineWidth',1.5);
hold on
  plot(ti(1:end),100*log( decomp_a(1:end,i) ),'k-','LineWidth', 1);
  plot(ti(1:end),100*log( decomp_b(1:end,i) ),'b-x','LineWidth',1);
  plot(ti(1:end),100*log( decomp_mp(1:end,i) ) ,'r--','LineWidth',1.0);
hold off
 title(name(i),'FontSize',16)
xlim([1983 2017])
ylabel('[ quarterly % ]')
if i ==4
   legend('Location','NorthEast')  
else   
 legend('Location','SouthWest')   
end 
 legend({': baseline',...
         ': only technology shock',...
         ': only discount factor shock',...
         ': only monetary policy shock' },'FontSize',12)
 legend('boxoff')
 
 %%
  subplot(1,2,2)
  
  i = 9;
  plot(ti(1:end),series_mu_a(1:end)/100 ,'k-','LineWidth',1);
%  plot(ti(1:end),shock_mu_a(1:end)/100 ,'k-','LineWidth',1);
 hold on 
    plot(ti(1:end),series_z_b(1:end)/100,'b-x','LineWidth',1);
%     plot(ti(1:end),shock_z_b(1:end)/100,'b-x','LineWidth',1);
    plot(ti(1:end), series_shock_mp(1:end)/100,'r--','LineWidth',1.0);
 hold off
      title(name(i),'FontSize',16)
   xlim([1983 2017])
%  ylabel('[ quarterly % ]')
 legend('Location','SouthWest')   
 legend({': technology shock',...
         ': discount factor shock',...
         ': monetary policy shock' },'FontSize',12)
 legend('boxoff')
   
savefig('./output/Fig6_contribution.fig'); 

    

     
 
 