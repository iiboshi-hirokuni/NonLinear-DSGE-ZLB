%
%   Plot Impulse Response
%% Canonical New Keynesian Model (Rotemberg Pricing) with Capital 
%   -Imposes the zero lower bound on the interest rate
%   -Endogenous MP rule responds to current inflation (deterministic)
%   -Balanced budget


% %% Options

period = 20;
num = 1000; % sampling number of generalized impulse response

if shock_idx == 3
     shock_b  = 0;              % size of Preference shock
     shock_a  = 0.0;            % size of TFP shock
     shock_MP = 0.1;           % size of MP shock
elseif shock_idx == 1
     shock_b  = 0.1;              % size of Preference shock
     shock_a  = 0.0;            % size of TFP shock
     shock_MP = 0.0;           % size of MP shock
else
     shock_b  = 0.0;              % size of Preference shock
     shock_a  = 0.1;            % size of TFP shock
     shock_MP = 0.0;           % size of MP shock
end  
 
shock=[shock_b; shock_a; shock_MP];
shock0 = zeros(3,1);

shock_title ={ '\epsilon^b_t \rightarrow ',  '\epsilon^a_t \rightarrow ', '\epsilon^r_t \rightarrow '  };

% Endogenous Variable
var_name =  { 'log(Y_t/A_t)','log(Y^*_t/A_t)','\pi_t','R_t','R^*_{t}','r^*_t','mu_a','z_b''E_y','E_pi','E_y_star' } ;   %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
nn = 9;
 
sum_imp1 = zeros(period,nn);
sum_imp0 = zeros(period,nn);
var_imp = zeros(period,nn);
def_imp = zeros(period,nn);

for j = 1:num  
   
 % State Values       
    R = O.rbound(1)+(O.rbound(2)-O.rbound(1))*rand(1);       %Capital state last period   
    Z_b = O.bbound(1)+(O.bbound(2)-O.bbound(1))*rand(1);      %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
    Mu_a = O.abound(1)+(O.abound(2)-O.abound(1))*rand(1);     
    pie = S.pi_star+randn(1)/100;
    y = S.y_star+randn(1)/100;   
  
   state= zeros(nn,1); 
    state(1)= R ;       % Interest rate state last period
    state(2)= Z_b ;     % Preference shock             %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
    state(3)= Mu_a ;     % Technology shock current period  
%     state(4)= epMP ;    % Monetary policy shock
   
 
% with shock  
imp1 = fun_imp(G, S, P, pf, state, shock, period,zlbflag,policy_flag);

% without shock
imp0 = fun_imp(G, S, P, pf, state, shock0, period,zlbflag,policy_flag);

% save respons to shock
sum_imp1 = sum_imp1 + (imp1);
sum_imp0 = sum_imp0 + (imp0);

var_imp = var_imp + (imp1-imp0).^2;
def_imp = def_imp + (imp1-imp0);

end

% mean of generalized impulse response
mean_imp1 = sum_imp1/num;
mean_imp0 = sum_imp0/num;
def_imp = def_imp/num;
var_imp1 = var_imp/num-def_imp.^2;
std_imp = sqrt(var_imp1); 

ti = 1:1:period;

% plot IRF
h_D = figure('Position',[20,20,450,800],'Name','imp response','Color','w');

%%  IRF 
 for i = 1:6 %size(var_name,2)

subplot(3,2,i)
   y1 = mean_imp1(:,i)-mean_imp0(:,i)-2*std_imp(:,i);
   y2 = 4*std_imp(:,i);
   hh = area(ti, [y1 y2]   ) ;
   
 set(hh(1),'FaceColor',[1 1 1])
   set(hh(2),'FaceColor',[0.5 1 1])   
   set(hh,'LineStyle','none') % Set all to same value
   
  hold on
    plot(ti,mean_imp1(:,i)-mean_imp0(:,i),'LineStyle','-','Color','blue', 'LineWidth',2.0); 

  hold off  
%     ylim([-0.5 2.0]);
    title([ shock_title(shock_idx) var_name(i)], 'Fontsize',10)
    xlabel('Horizon')
%     ylabel('(%)')
    if i == 6
      legend('', '2 \times StDev', 'mean' )
    end
 end   
    
% filename =([ '../result/imp', num2str(zlbflag) ]);
% save(filename, 'ln_imp1','ln_imp0');
    



    
    
    
    
    

