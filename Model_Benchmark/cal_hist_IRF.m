%
%   Plot Impulse Response
%% Canonical New Keynesian Model (Rotemberg Pricing) with Capital 
%   -Imposes the zero lower bound on the interest rate
%   -Endogenous MP rule responds to current inflation (deterministic)
%   -Balanced budget

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');

if zlbflag==0 && policy_flag == 0
       type_model = char('./output/model0');
 elseif zlbflag==1 && policy_flag == 1
       type_model = char('./output/model1');
 elseif zlbflag==1 && policy_flag == 2
       type_model = char('./output/model2');
 end     

% %% Options

shock_idx =3;

period = 10;
num = 1000; % sampling number of generalized impulse response

     shock_b  = 0;              % size of Preference shock
     shock_a  = 0.0;            % size of TFP shock
     shock_MP = 0.01;           % size of MP shock
 
shock1=[shock_b; shock_a; shock_MP];
shock2=[shock_b; shock_a; -1*shock_MP];
shock0 = zeros(3,1);

shock_title ={ '\epsilon^b_t \rightarrow ',  '\epsilon^a_t \rightarrow ', '\epsilon^r_t \rightarrow '  };

% Endogenous Variable
var_name =  { 'log(Y_t/A_t)','log(Y^*_t/A_t)','\pi_t','R_t','R^*_{t}','r^*_t','mu_a','z_b''E_y','E_pi','E_y_star' } ;   %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
nn = 9;
 
load([  type_model  '_state.mat'])

save_imp1 = zeros(period,nn,Tobs);
save_imp2 = zeros(period,nn,Tobs);

for i = 1:Tobs    
    
    state = zeros(nn,1);
    state(1,1) = exp(series_R_lag(i)/100);
    state(3,1) = series_mu_a(i)/100;
    state(2,1) = series_z_b(i)/100;  
    
%     shock0(3,1)=series_shock_mp(i);
%     shock1(3,1)=series_shock_mp(i)+shock_MP;
    
%     state(1)= R ;       % Interest rate state last period
%     state(2)= Z_b ;     % Preference shock             %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
%     state(3)= Mu_a ;

    imp1 = fun_imp(G, S, P, pf, state, shock1, period,zlbflag,policy_flag);

    % without shock
    imp0 = fun_imp(G, S, P, pf, state, shock0, period,zlbflag,policy_flag);
    
    save_imp1(:,:,i)= imp1-imp0 ;

end   

    
for i = 1:Tobs    
    
    state = zeros(nn,1);
    state(1,1) = exp(series_R_lag(i)/100);
    state(3,1) = series_mu_a(i)/100;
    state(2,1) = series_z_b(i)/100;  
    
    shock0(3,1)=series_shock_mp(i)/100;
    shock2(3,1)=series_shock_mp(i)/100+shock_MP*(-1);
    
%     state(1)= R ;       % Interest rate state last period
%     state(2)= Z_b ;     % Preference shock             %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
%     state(3)= Mu_a ;

    imp1 = fun_imp(G, S, P, pf, state, shock2, period,zlbflag,policy_flag);

    % without shock
    imp0 = fun_imp(G, S, P, pf, state, shock0, period,zlbflag,policy_flag);
    
    save_imp2(:,:,i)= imp1-imp0 ;

end   



for i = 1:4
    figure(300+i)
    hist_imp_y = squeeze(save_imp1(:,i,:));
    
    s = 1983;
    [x,y]= meshgrid((s+2/4):0.25:s+(Tobs+1)/4,1:10 );
    surf(x,y,hist_imp_y,'FaceLighting','gouraud','LineWidth',0.3)
    colormap winter
    view(40,20)
    
    if i == 1 
       title('Postive MP Shock \rightarrow Output','FontSize',14)
    elseif i == 2 
       title('Postive MP Shock \rightarrow Potential Output','FontSize',14) 
    elseif i == 3
       title('Postive MP Shock \rightarrow Inflation','FontSize',14)       
    elseif i == 4
       title('Postive MP Shock \rightarrow Interest Rate','FontSize',14)
    end   
       zlabel('[ % ] ','FontSize',12)
    xlabel('period','FontSize',12)
    xlim( [1983 2016])
    ylabel('horizon','FontSize',12)
end   


for i = 1:4
    figure(300+i+4)
    hist_imp_y = squeeze(save_imp2(:,i,:));
    
    s = 1983;
    [x,y]= meshgrid((s+2/4):0.25:s+(Tobs+1)/4,1:10 );
    surf(x,y,hist_imp_y,'FaceLighting','gouraud','LineWidth',0.3)
    colormap winter
    view(40,20)
    
    if i == 1 
       title('Negative MP Shock \rightarrow Output','FontSize',14)
    elseif i == 3
       title('Negative MP Shock \rightarrow Inflation','FontSize',14)       
    elseif i == 4
       title('Negative MP Shock \rightarrow Interest Rate','FontSize',14)
    end   
       zlabel('[ % ] ','FontSize',12)
    xlabel('period','FontSize',12)
    xlim( [1983 2016])
    ylabel('horizon','FontSize',12)    
    
end   


figure(310)
for j = 1:6
   subplot(3,2,j)   
   hold on
     imp1 = squeeze(save_imp1(:,1,(j-1)*20+7));
      plot(imp1, 'k-', 'LineWidth',2 )         
     imp2 = squeeze(save_imp2(:,1,(j-1)*20+7));
      plot(imp2, 'r--', 'LineWidth',2)    
   hold off  
   title([ 'Output  ' num2str((j-1)*5+1985) ' :Q1  ' ],'FontSize',12)
   if j == 6
     legend('Positive MP', 'Negative MP')
   end
   ylim([-0.02 0.02])
end   

figure(311)
for j = 1:6
   subplot(3,2,j)   
   hold on
     imp1 = squeeze(save_imp1(:,3,(j-1)*20+7));
      plot(imp1, 'k-', 'LineWidth',2 )         
     imp2 = squeeze(save_imp2(:,3,(j-1)*20+7));
      plot(imp2, 'r--', 'LineWidth',2)    
   hold off  
   title([ 'Inflation  '  num2str((j-1)*5+1985) ' :Q1  ' ],'FontSize',12)
   if j == 6
     legend('Positive MP', 'Negative MP')
   end
   ylim([-0.0025 0.0025])
end   

figure(312)
for j = 1:6
   subplot(3,2,j)   
   hold on
     imp1 = squeeze(save_imp1(:,4,(j-1)*20+7));
      plot(imp1, 'k-', 'LineWidth',2 )         
     imp2 = squeeze(save_imp2(:,4,(j-1)*20+7));
      plot(imp2, 'r--', 'LineWidth',2)    
   hold off  
   title([ 'Interest Rate '  num2str((j-1)*5+1985) ' :Q1  ' ],'FontSize',12)
   if j == 6
     legend('Positive MP', 'Negative MP')
   end
   ylim([-0.01 0.01])
end   


if zlbflag== 1 && policy_flag ==  1
   save('./output/IRF_model1.mat','save_imp1', 'save_imp2')
elseif zlbflag== 1 && policy_flag ==  2
   save('./output/IRF_model2.mat','save_imp1', 'save_imp2')
elseif zlbflag== 0 && policy_flag ==  0
   save('./output/IRF_model0.mat','save_imp1', 'save_imp2')
end   






