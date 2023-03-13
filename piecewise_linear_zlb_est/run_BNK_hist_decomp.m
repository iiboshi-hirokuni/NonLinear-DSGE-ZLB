%% runsim

format short
   location = 'home_windows';
% % location = 'home_mac';'
setpathdynare4
clear
set(0,'DefaultLineLineWidth',2)
opaths=0;plot_irf=0;colors=char('b','r','k');
close all
tlast=2011.875;


imodel=[1 2 3 ]; opaths=1; maxindi=1; colorpaths='r.';

load PARAM_EXTRA_CALIBRATED
load mle_estimates_fminsearch params_labels params* filtered_errs sample_length err_list
for i=1:numel(params_labels)
    evalc([ cell2mat(params_labels(i)) '= params1(' num2str(i) ')']) ;
end

err_list = char('eps_a','eps_b','eps_r');
% save  PARAM_EXTRA ...
%     BETA BETA1 EC EH ETA JEI M ALPHA PHIK DK LAGP LAGW PIBAR  ...
%     SIGMA TAYLOR_P TAYLOR_R TAYLOR_Y TETAP TETAW XP_SS XW_SS  ...
%     RHO_J RHO_K RHO_P RHO_R RHO_W RHO_Z ...
%     STD_J STD_K STD_P STD_R STD_W STD_Z RHOD

%% -------------------------------
% Load data and Declare observables
%-------------------------------

r_star  = exp(csigma*gamma_a)/cbeta;
data_r_ss1 = 100*(r_star*pi_star-1);

data=csvread('./data/data_jpn_1_def.csv',1,2);
data_y=data(:,1)/100; data_pi=data(:,2)/100; data_r= (data(:,3)-data_r_ss1)/100;
Make_data_zero_rate;

% data_y=data(:,1); data_pi=data(:,2); data_r= (data(:,3)-data_r_ss1);
% data_y=data(:,1); data_pi=data(:,2); data_r=data(:,3)-0.8682;
tstar = 1;

% load data_est_jp
obs_list = char('data_pi','data_r','data_y');

obs=[ data_pi(tstar:end) data_r(tstar:end) data_y(tstar:end) ];
save observables obs
tt_obs = 1983.25:0.25:1983+(size(data_y(tstar:end),1)-1)/4;

for model=[imodel]
    
    
   modnam_00 = 'BNK_model'; % base model (constraint 1 and 2 below don't bind)
   modnam_10 = 'BNK_model_zlb';  % first constraint is true
   % modnam_10 = 'BNK_model';  % first constraint is true
   modnam_01 = 'BNK_model';  % second constraint is true
   modnam_11 = 'BNK_model'; % both constraints bind

%%
constraint1       = 'rnot  < r_zlb'; 
constraint_relax1 = 'r_t > r_zlb';
% constraint1       = 'r_t  < r_zlb'; 
% constraint_relax1 = 'rnot > r_zlb';

constraint2       = 'pi_t  > 10'; 
constraint_relax2 = 'pi_t < 10';
    
    
    curb_retrench = 0;           % if 1 slow down relaxation of constraints
    maxiter = 10;                % number of iterations allowed to look for a solution
    
%      irfshock =char('eps_a','eps_b','eps_r');     
%      sequence1 = [ STD_a STD_b STD_r ];
%      nperiods = 40;
%      sequence=sequence1; sequence(:,setdiff([1:3],[model]))=0;
%     
    sample_length=size(filtered_errs,1);
    yearplot=tlast-(sample_length-1)/4;
    if model==1; model2=[1]; end
    if model==2; model2=[1 2]; end
    if model==3; model2=[1 2 3]; end
%     if model==4; model2=[1 2 3 4 5 6]; end
    irfshock = char(err_list(model2,:));
    sequence = filtered_errs(:,model2);
    nperiods=size(filtered_errs,1);
    
    
    [zdatal zdatap zdatass oobase_ Mbase_] = solve_two_constraints_fast2_temp1(...
        modnam_00,modnam_10,modnam_01,modnam_11,...
        constraint1, constraint2,...
        constraint_relax1, constraint_relax2,...
        sequence,irfshock,nperiods+50,curb_retrench,maxiter);
    
        
    for i=1:Mbase_.endo_nbr
        eval([deblank(Mbase_.endo_names(i,:)),'_l=zdatal(1:nperiods,i);']);
        eval([deblank(Mbase_.endo_names(i,:)),'_p=zdatap(1:nperiods,i);']);
        eval([deblank(Mbase_.endo_names(i,:)),'_ss=zdatass(i);']);
    end
    
    for i = 1:size(Mbase_.param_names,1)
        eval([Mbase_.param_names(i,:),'= Mbase_.params(i);']);
    end
    
    ymod(:,model)=100*data_y_p;
    pmod(:,model)=100*data_pi_p;
    rmod(:,model)=100*data_r_p;
    ydata=100*data_y(tstar:end);
    pdata=100*data_pi(tstar:end);
    rdata=100*data_r(tstar:end)-r_zlb;
    tt=tt_obs(tstar:end);
            
    
end

%%
figure(10)
plot(tt_obs,rdata(2:end)+0.55,'b-');
hold on
%   plot(tt_obs,rmod(2:end,3)-rmod(1,3)+rdata(1),'r:');
  plot(tt_obs,100*data_r_l(2:end)-100*data_r_l(1)+rdata(1)+0.55,'r:');
hold off
% ylim([-0.1 0.5])
 axis tight
 
 plot_filter_err4;


%%
% close all

rdata=100*(data_r(tstar:end)); %-data_r_ss1/2;

bars_color1 = ...
    [    0      0.5000         0
         0      0.4100         0.5500
         0.6900    0.8900    1.0000 ];

sup_title = '';
titlelist = char('Output Growth','Inflation','Interest Rate');
Tvec = tt_obs(tstar);
Tdelta = 1/4;
npers = size(ymod,1);
Tmax = Tvec+npers*Tdelta ;
ylabels = char('%','%');

% opts.legendplace = [.4,0,.2,.1];
opts.legendplace = [0.5,0,0.2,0.1];
% opts.shading = 1;
opts.shading = 0;

tt=tt_obs(tstar:end);

legendlist2 =char('Data','Technology', 'Preference',...
                   'Monetary Policy');

bars_color2=[  1         0.5000     0.5
               0.5       1          0.500
               1.0       1.00       0.0000       ];

figure
makebardec4(bars_color2,titlelist,legendlist2,sup_title,Tvec,Tdelta,Tmax,ylabels,opts,...
    [ydata pdata rdata ],...
    [ymod(:,1) pmod(:,1) rmod(:,1) ],...
    [ymod(:,2)-ymod(:,1) pmod(:,2)-pmod(:,1) rmod(:,2)-rmod(:,1) ],...
    [ymod(:,3)-ymod(:,2) pmod(:,3)-pmod(:,2) rmod(:,3)-rmod(:,2) ])

% for i=1:2
%     hold on
%     subplot(3,1,1)
%     ylim([-40 20])
%     hold on
%     subplot(3,1,2)
%     ylim([-9 6])
% end










