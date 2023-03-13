%
% %   Figure 8: Natural Rate of Interest: Counterfactual Simulation
% 
% 

addpath('../Toolbox')
addpath('../fun_projection')
addpath('../fun_particle_filter')
addpath('../fun_hist_decomp')
addpath('../fun_smc')
addpath('../fun_prior')

addpath('./data_sample')

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');


% %% Options

period = 1;
shock1 = zeros(3,1);

% Endogenous Variable
var_name =  { 'log(Y_t/A_t)';'log(Y^*_t/A_t)';'\pi_t';'R_t';'R^*_{t}';...
              'r^*_t';'mu_a';'z_b';'E_y';'E_pi';'E_y_star' } ;  %%  Modified on Jan 30, 2017 by UEDA
nn = 9;

%%  Case 0:  Benchmark,  non ZLB restriction, non_ZLB parameters and Shocks

zlbflag      =  0;
policy_flag  =  0;

% load('./output/model0_state.mat')
load('model0_state.mat')
pf_model0 = pf;
S_model0 = S;
P_model0 = P;

temp_var = zeros(period,nn,Tobs);

for t = 1:Tobs    
    
    state = zeros(nn,1);
    state(1,1) = exp(series_R_lag(t)/100);
    state(3,1) = series_mu_a(t)/100;
    state(2,1) = series_z_b(t)/100;  
    
    shock1(3,1)=series_shock_mp(t);    

    imp1 = fun_imp(G, S_model0, P_model0, pf_model0, state, shock1, period,zlbflag,policy_flag);
    
    temp_var(:,:,t)= imp1 ;

end   

var_case0 = squeeze(temp_var(1,:,:));




%%  Case 1:  ZLB and Policy type 1 with non_ZLB parameters and Shocks

zlbflag      =  1;
policy_flag  =  1;

% load('./output/model0_state.mat')
load('model0_state.mat')
pf_model0 = pf;
S_model0 = S;
P_model0 = P;

temp_var = zeros(period,nn,Tobs);

for t = 1:Tobs    
    
    state = zeros(nn,1);
    state(1,1) = exp(series_R_lag(t)/100);
    state(3,1) = series_mu_a(t)/100;
    state(2,1) = series_z_b(t)/100;  
    
    shock1(3,1)=series_shock_mp(t);    

    imp1 = fun_imp(G, S_model0, P_model0, pf_model0, state, shock1, period,zlbflag,policy_flag);
    
    temp_var(:,:,t)= imp1 ;

end   

var_case1 = squeeze(temp_var(1,:,:));


%%  Case 2:  ZLB parameters with non ZLB restriction and Shocks

zlbflag      =  0;
policy_flag  =  0;

load('model1_state.mat')
pf_model1 = pf;
S_model1 = S;
P_model1 = P;

load('model0_state.mat')


temp_var = zeros(period,nn,Tobs);

for t = 1:Tobs    
    
    state = zeros(nn,1);
    state(1,1) = exp(series_R_lag(t)/100);
    state(3,1) = series_mu_a(t)/100;
    state(2,1) = series_z_b(t)/100;  
    
    shock1(3,1)=series_shock_mp(t);    

    imp1 = fun_imp(G, S_model1, P_model1, pf_model1, state, shock1, period,zlbflag,policy_flag);
    
    temp_var(:,:,t)= imp1 ;

end   

var_case2 = squeeze(temp_var(1,:,:));


%%  Case 3:  ZLB parameters with non ZLB restriction and Shocks

zlbflag      =  0;
policy_flag  =  0;

load('model0_state.mat')
pf_model3 = pf;
S_model3 = S;
P_model3 = P;

load('model1_state.mat')


temp_var = zeros(period,nn,Tobs);

for t = 1:Tobs    
    
    state = zeros(nn,1);
    state(1,1) = exp(series_R_lag(t)/100);
    state(3,1) = series_mu_a(t)/100;
    state(2,1) = series_z_b(t)/100;  
    
    shock1(3,1)=series_shock_mp(t);    

    imp1 = fun_imp(G, S_model3, P_model3, pf_model3, state, shock1, period,zlbflag,policy_flag);
    
    temp_var(:,:,t)= imp1 ;

end   

var_case3 = squeeze(temp_var(1,:,:));

%% 
s =1983.0;
 ti = (s+2/4):0.25:s+(Tobs+1)/4; 
 
figure('Name','Counter Factual Simulation','File', 'sim_r_star.tif')
hold on
plot(ti,100*(var_case0(6,:)-1),'b-','linewidth',3); 
plot(ti,100*(var_case1(6,:)-1),'r:','linewidth',2);    
plot(ti,100*(var_case2(6,:)-1),'k--','linewidth',2.5); 
plot(ti,100*(var_case3(6,:)-1),'c-','linewidth',1.5);
hold off
xlim([1981  2018])
%title(var_name(6,:));
lgnd = legend('Model w/o ZLB','Model w/o ZLB w/ ZLB','Model w/o ZLB w/ Model 1 Parameters', 'Model w/o ZLB w/ Model 1 Shocks')
set(lgnd, 'Box', 'off')

savefig('./output/Fig8_Counter_factual.fig'); 

% if zlbflag== 1 && policy_flag ==  1
%    save('IRF_model1.mat','save_imp1', 'save_imp2')
% elseif zlbflag== 1 && policy_flag ==  2
%    save('IRF_model2.mat','save_imp1', 'save_imp2')
% elseif zlbflag== 0 && policy_flag ==  0
%    save('IRF_model0.mat','save_imp1', 'save_imp2')
% end   






