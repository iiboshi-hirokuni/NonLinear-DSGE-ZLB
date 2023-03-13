%
%
%   Table 3 
%
%



addpath('./Plot_Fig/data_sample/')

 model = 'model 1';
% model = 'model 2';
% model = 'model 0';

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');

addpath('./fun_hist_decomp')
addpath('./Toolbox')

switch model
    case 'model 1'
         load('model1_state.mat')
         model_idx = 1;
    case 'model 2'
         load('model2_state.mat')
         model_idx = 2;      
   case 'model 0'
         load('model0_state.mat')
         model_idx = 3;
end         

title_model = char( 'Model 1 ',  'Model 2',  'Model w/o ZLB'  );         
         
rng(100) 
shock_idx =3;

period = 4*2;
nsim = 1000;    % sampling number of monte calro
% nsim = 100000;
nn = 9+2;

shock =zeros(3,1);
    shock(1) = P.sigma_b; % preference
    shock(2) = P.sigma_a; % TFP
    shock(3) = P.sigma_r;

% series_r_star;


%% condition
for j = 1:1
    
    if j == 1
         condition = char(' Natural Rate >= 3%');
         l = find( (series_r_star >= 3) );
    elseif j == 2
         condition = char(' 3% >  Natural Rate >= 2%');
         l = find( (series_r_star >= 2) & (series_r_star <3) );
    elseif j == 3
         condition = char(' 2% >  Natural Rate >= 1%');
         l = find( (series_r_star >= 1) & (series_r_star <2) );
    elseif j == 4
         condition = char(' 1% >  Natural Rate >= 0%');      
         l = find( (series_r_star >= 0) & (series_r_star <1) );
    elseif j == 5
         condition = char(' 0% >  Natural Rate >= -1%');      
         l = find( (series_r_star >= -1) & (series_r_star <0) );
    elseif j == 6
         condition = char(' -1% >  Natural Rate' ); 
         l = find( (series_r_star < -1 ) );
    end     

%  disp( l' ) 
 m = size(l,1);
 num2 = 0;
 num0 = 0;

 for i = 1:m
    for n = 1:nsim

    state = zeros(nn,1);
       state(1,1) = exp(series_R_lag(l(i))/100);
       state(3,1) = series_mu_a(l(i))/100;
       state(2,1) = series_z_b(l(i))/100;
    
    yy = fun_data_gen(G, S, P, pf, state, shock, period,zlbflag,policy_flag);
    
%     inflation greater than 2%
    k2 = find( 4*log(yy(:,3))>log(1.02));
%     inflation less than 0 %
    k0 = find( 4*log(yy(:,3))<log(1.00));
    
    if isempty(k2)
       cnt2 = 0;
    else   
       cnt2 = 1;
    end
    
   if isempty(k0)
       cnt0 = 0;
    else   
       cnt0 = 1;
   end
    
       num2 = num2 + cnt2;
       num0 = num0 + cnt0;   
   end  

 end
 
 disp('------------------------------------------------------------------')   
 disp( [ title_model(model_idx,:)   ])
 disp( ['condtion  is '  condition   ])
 disp('------------------------------------------------------------------')     
 disp( [ 'prob ( \pi_{t+h} > 2% | t+1 < t+h < t+8  )= ' num2str(round(num2/(m*nsim)*10000)/100)  ' %' ]);
 disp( [ 'prob ( \pi_{t+h} < 0% | t+1 < t+h < t+8  )= ' num2str(round(num0/(m*nsim)*10000)/100)  ' %' ]);
 disp('------------------------------------------------------------------')   
 
end 


% %% Steady State
% ss_mu_a =0;
% ss_z_b = 0;  % log(Z_b_SS)
% num2 = 0;
% num0 = 0;
% 
% for n = 1:nsim
% 
%     state = zeros(nn,1);
%        state(1,1) = S.R;
%        state(3,1) = ss_mu_a/100;
%        state(2,1) = ss_z_b;
%     
%     yy = fun_data_gen(G, S, P, pf, state, shock, period,zlbflag,policy_flag);
%     
% %     inflation greater than 2%
%     k2 = find( 4*log(yy(:,3))>log(1.02));
% %     inflation less than 0 %
%     k0 = find( 4*log(yy(:,3))<log(1.00));
%     
%     if isempty(k2)
%        cnt2 = 0;
%     else   
%        cnt2 = 1;
%     end
%     
%    if isempty(k0)
%        cnt0 = 0;
%     else   
%        cnt0 = 1;
%    end
%     
%        num2 = num2 + cnt2;
%        num0 = num0 + cnt0;   
% end  
% 
% disp('------------------------------------------------------------------')   
%  disp( [ title_model(model_idx,:)   ])
%  disp( ['Initial Value is  Steady State  '  ])
%  disp('------------------------------------------------------------------')     
%  disp( [ 'prob ( \pi_{t+h} > 2% | t+1 < t+h < t+8  )= ' num2str(round(num2/(1*nsim)*10000)/100)  ' %' ]);
%  disp( [ 'prob ( \pi_{t+h} < 0% | t+1 < t+h < t+8  )= ' num2str(round(num0/(1*nsim)*10000)/100)  ' %' ]);
%  disp('------------------------------------------------------------------')   
 


 