%-------------------------------------------
% Housekeeping
%-------------------------------------------
 clear all
restoredefaultpath
%% setpathdynare4
location = 'home_windows';
% location = 'home_mac';
setpathdynare4
%%
set(0,'DefaultLineLineWidth',2)
format compact

global oo00_  M00_ M10_  M01_  M11_
global params_labels params
global cof cof10 cof01 cof11 ...
    Jbarmat Jbarmat10 Jbarmat01 Jbarmat11 ...
    Dbarmat10 Dbarmat01 Dbarmat11 ...
    decrulea decruleb
global filtered_errs_switch filtered_errs_init model_temp
global datavec irep xstory fstory

%-------------------------------
% Load estimated stuff needed to compile model
%-------------------------------
 load mle_estimates_fminsearch.mat
%  load mle_estimates_temp_test
%  load mle_initial_guess

%-------------------------------
% Overwrite calibrated and estimated parameters if needed
%-------------------------------

% 
% save  PARAM_EXTRA_CALIBRATED ...
%   ALPHA BETA BETA1 DK EC EH ETA JEI LAGP LAGW M PHIK PIBAR 

% XX = [ ];
% save PARAM_EXTRA_BABY XX

%-------------------------------
% Load data and Declare observables
% Create script to speed up filtering
%-------------------------------
csigma=1.5;
gamma_a=0.0025;
cbeta = 0.9975;
pi_star = 1.005 ;
r_star  = exp(csigma*gamma_a)/cbeta;
data_r_ss1 = 100*(r_star*pi_star-1);
r_zlb = -(r_star*pi_star-1);

err_list = char('eps_a','eps_b','eps_r');

data=csvread('./data/data_jpn_1_def.csv',1,2);
data_y=data(:,1)/100; data_pi=data(:,2)/100; data_r= (data(:,3) - data_r_ss1)/100;
Make_data_zero_rate;
tstar = 1;

% load data_est_jp
obs_list = char('data_pi','data_r','data_y');

obs=[ data_pi(tstar:end) data_r(tstar:end) data_y(tstar:end) ];
save observables obs
tt_obs = 1983.25:0.25:1983+(size(data_y(tstar:end),1))/4;
% rzlb = -(r_star*pi_star-1);

ntrain = 20;

call_process_mod_files

%% -----------------------------------
% Filter shocks at the parameter values specified above
%-----------------------------------
filtered_errs_init = zeros(size(obs,1),size(err_list,1));

%  [filtered_errs resids Emat requalzero ] = myfilterzlbrnot(constraint1_difference, constraint2_difference,...
%                       constraint_relax1_difference, constraint_relax2_difference,err_list,obs_list,obs,r_zlb);
   

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
    
%------------------------------------------
% Feed filtered shocks back into model
%------------------------------------------
% filtered_errs(:,[2:6])=0;
  model=3;
  
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

% [zdatal zdatap zdatass oo00_ M00_ ] = solve_two_constraints_fast2_temp1(...
%            modnam_00,modnam_10,modnam_01,modnam_11,...
%            constraint1, constraint2,...
%            constraint_relax1, constraint_relax2,...
%            filtered_errs,err_list,size(obs,1),curb_retrench,maxiter,zeros(M00_.endo_nbr,1));


for i=1:M00_.endo_nbr
  eval([deblank(M00_.endo_names(i,:)),'_l=zdatal(1:nperiods,i);']);
  eval([deblank(M00_.endo_names(i,:)),'_p=zdatap(1:nperiods,i);']);
  eval([deblank(M00_.endo_names(i,:)),'_ss=zdatass(i);']);
end

%%
 obs_name = char('data pi','data r','data y');
 
 co=100*obs(1,2)-100*r_zlb;

 figure('File','natural_rate')
final_sample = [];
for index=1:size(obs_list,1)
  subplot(4,1,index)
  eval([ 'final_sample(:,index) = ' deblank(obs_list(index,:)) '_p;'])
%   l3= plot(tt_obs,eval([deblank(obs_list(index,:)) '_l'])+eval([deblank(obs_list(index,:)) '_ss']),'g'); 
    if index ~= 2 
       l3= plot(tt_obs,100*eval([deblank(obs_list(index,:)) '_l']),'g');
    else 
        l3= plot(tt_obs,100*eval([deblank(obs_list(index,:)) '_l'])...
                  -100*eval([deblank(obs_list(index,:)) '_l(1)']) +co,'g');
    end    
   hold on 
%       l1= plot(tt_obs,obs(:,index)+eval([deblank(obs_list(index,:)) '_ss']),'b');   
%       l2= plot(tt_obs,eval([deblank(obs_list(index,:)) '_p'])+eval([deblank(obs_list(index,:)) '_ss']),'r:'); 
    if index ~= 2 
       l1= plot(tt_obs,100*obs(:,index),'b'); 
       l2= plot(tt_obs,100*eval([deblank(obs_list(index,:)) '_p']),'r:');  
%        axis tight
    else
       l1= plot(tt_obs,100*obs(:,index)-100*r_zlb,'b'); 
       l2= plot(tt_obs,100*eval([deblank(obs_list(index,:)) '_p'])...
               - 100*eval([deblank(obs_list(index,:)) '_p(1)']) +co,'r:'); 
       ylim([-0.05, 2.0]);
    end  
   axis tight
   title(obs_name(index,:),'FontSize',12)
end

  legend([l1,l2,l3],{'data','piece-wise-linear model','linear model'},'FontSize',10)

    subplot(4,1,4)
    l3= plot(tt_obs,r_star_t_l,'g');
    hold on
       l2= plot(tt_obs,r_star_t_p,'r:');
    hold off
    axis tight
    title('natural rate','FontSize',12)
 %%   
%     for i=1:M00_.endo_nbr
%       eval([ 'save result.mat ' deblank(M00_.endo_names(i,:)),'_l, -append;']);
%       eval(['save result.mat ' deblank(M00_.endo_names(i,:)),'_p, -append;']);
%       eval(['save result.mat ' deblank(M00_.endo_names(i,:)),'_ss, -append;']);
%     end
   save( 'pwl_result.mat' , 'r_star_t_l','r_star_t_p','r_star_t_ss'  );
   save( 'pwl_result.mat' , 'y_star_t_l','y_star_t_p','y_star_t_ss'  ,'-append');
   save( 'pwl_result.mat' , 'r_t_l','r_t_p','r_t_ss'  ,'-append');
   save( 'pwl_result.mat' , 'y_t_l','y_t_p','y_t_ss'  ,'-append');
   save( 'pwl_result.mat' , 'pi_t_l','pi_t_p','pi_t_ss'  ,'-append');
   save( 'pwl_result.mat' , 'r_zlb','-append'); 
%%  
  figure('File','natural_rate')
 final_sample = [];
    index= 2;
%   subplot(1,1,index)
  eval([ 'final_sample(:,index) = ' deblank(obs_list(index,:)) '_p;'])
      l3= plot(tt_obs,100*eval([deblank(obs_list(index,:)) '_l'])...
          - 100*eval([deblank(obs_list(index,:)) '_l(1)']) +co,'g');
       
   hold on
       l1= plot(tt_obs,100*obs(:,index)-100*r_zlb,'b'); 
       l2= plot(tt_obs,100*eval([deblank(obs_list(index,:)) '_p'])...
           - 100*eval([deblank(obs_list(index,:)) '_p(1)']) +co,'r:');
   hold off
  axis tight
   title(obs_name(index,:),'FontSize',12) ;
    ylim([-0.3, 0.2]);
   legend([l1,l2,l3],{'data','piece-wise-linear model','linear model'},'FontSize',10)
%    
 %%
%  figure(4000)
%  final_sample = [];
%     index= 1;
% %   subplot(1,1,index)
%   eval([ 'final_sample(:,index) = ' deblank(obs_list(index,:)) '_p;'])
%         l3= plot(tt_obs,100*eval([deblank(obs_list(index,:)) '_l']),'g');
%        
%    hold on 
%        l1= plot(tt_obs,100*obs(:,index),'b'); 
%        l2= plot(tt_obs,100*eval([deblank(obs_list(index,:)) '_p']),'r:');  
%    hold off  
%   axis tight
%    title(obs_name(index,:),'FontSize',12)
%    ylim([-0.8, 0.1])
%    xlim([2000 2016])
%    legend([l1,l2,l3],{'data','piece-wise-linear model','linear model'},'FontSize',10)
%    


