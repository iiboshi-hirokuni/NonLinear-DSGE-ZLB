%-------------------------------------------
% Housekeeping
%-------------------------------------------
clc
clear all 
close all

warning('off','all')
warning


global oo00_  M00_ M10_  M01_  M11_ params_labels params

global cof cof10 cof01 cof11 Jbarmat Jbarmat10 Jbarmat01 Jbarmat11 ...
  Dbarmat10 Dbarmat01 Dbarmat11 decrulea decruleb

global filtered_errs_switch filtered_errs_init model_temp datavec irep xstory fstory

%% setpathdynare4
location = 'home_windows';
% location = 'home_mac';
setpathdynare4
%%

irep=1; datavec=[]; xstory=[]; fstory=[];

set(0,'DefaultLineLineWidth',2)
randn('seed',1);
format compact
filtered_errs_switch=0;

tic
tStart = tic;

%% ----------------------------------------------------------------------
% Invoke calibrated parameters 
%----------------------------------------------------------------------

 paramfile_BNK_zlb


% save  PARAM_EXTRA_CALIBRATED ...
%       pi_star r_star r_zlb phi cbeta omega epsilon csigma gamma_a h chi kappa ...
%       rho_a rho_b rho_r STD_a STD_b STD_r psi_pi psi_y m1 m2 ; 
 
params_matrix = {  ...
    'pi_star '   1.0025   1.00   1.005     1    'NORMAL_PDF' 1.0025  0.0005;
    'h '          0.5    0.45      0.65     1    'BETA_PDF'      0.5   0.2;
    'csigma '     1.5     1.3      2.0      1      'NORMAL_PDF'  1.5     0.3
    'gamma_a '  0.001   -0.0025   0.0025    1      'NORMAL_PDF'  0.00  0.0005;
    'omega '    3.5      2.5       4.5    1      'NORMAL_PDF'  3.5     0.5;
    'kappa '    0.05    0.02       0.07     1      'NORMAL_PDF'  0.05  0.006;     
    'psi_pi '   2.0     1.70      3       1      'GAMMA_PDF'   2.0   0.3 ;
    'psi_y '    0.125   0.05        0.5     1      'GAMMA_PDF' 0.125  0.025 ;
   'rho_a '    0.5    0.25      0.9     1     'BETA_PDF'      0.5   0.2;
   'rho_b '    0.5    0.25      0.9     1     'BETA_PDF'      0.5   0.2;
   'rho_r '    0.5    0.25      0.9     1     'BETA_PDF'      0.5   0.2 ;
    'STD_a '   0.05    0.0001  	0.10     1     'INV_GAMMA_PDF' 0.02   5;
    'STD_b '   0.05    0.0001  	0.10     1     'INV_GAMMA_PDF' 0.02   5;
    'STD_r '   0.05    0.0001  	0.10     1     'INV_GAMMA_PDF' 0.02   5 };


save params_matrix params_matrix
  

% save 

err_list = char('eps_a','eps_b','eps_r'); % alphabet order

H0 = diag(cell2mat(params_matrix(:,5)).^2) ;

params_labels = params_matrix(:,1);
params0 =   cell2mat(params_matrix(:,2));
params_lo = cell2mat(params_matrix(:,3));
params_hi = cell2mat(params_matrix(:,4));
params_mean = cell2mat(params_matrix(:,7));
params_std = cell2mat(params_matrix(:,8));

dist_names = params_matrix(:,6);
codes = dist_names2codes(dist_names);
[p6 p7] = get_dist_inputs(codes,params_mean,params_std);


   [prior] = -priordens(params0, codes, p6, p7, params_lo, params_hi,1);


 
  
for i=1:numel(params_labels)
  evalc([ cell2mat(params_labels(i)) '= params0(' num2str(i) ')']) ;
end


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
obs_list = char('data_pi','data_r','data_y'); % alphabet order

obs=[ data_pi(tstar:end) data_r(tstar:end) data_y(tstar:end) ];
save observables obs
tt_obs = 1983.25:0.25:1983+(size(data_y(tstar:end),1)-1)/4;

r_zlb = -(r_star*pi_star-1);
ntrain = 20;


%% -----------------------------------
% Create script to speed up filtering
%-----------------------------------

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

call_process_mod_files




%   method = 'initial_check'
  method = 'fminsearch'



%%
if strmatch(method,'initial_check')==1
% -----------------------------------
% Check value of the likelihood at the initial guess
%-----------------------------------


  [posterior filtered_errs like prior resids ]=...
   posteriorzlb(params0,params_labels,params_lo,params_hi,...
   modnam_00,modnam_10,modnam_01,modnam_11,...
   constraint1_difference, constraint2_difference,...
   constraint_relax1_difference, constraint_relax2_difference,...
   err_list,obs_list,obs,ntrain,r_zlb, codes, p6, p7,tStart);

  params=params0;
  sample_length = size(obs,1);
      
  save mle_initial_guess

  

end      


%%

if strmatch(method,'fminsearch')==1
  
%    tolerance = 1e-1;      
  tolerance = 1e-3;
  
  options = optimset('Display','Iter','TolFun',tolerance,'TolX',tolerance,'TolFun',tolerance,...
    'MaxFunEvals',10000,'MaxIter',500,...
    'DiffMinChange',tolerance,'Algorithm','Interior-Point'); 
  
  [params,fval]=fminsearchbnd(@(current_params) posteriorzlb...
    (current_params,params_labels,params_lo,params_hi,...
    modnam_00,modnam_10,modnam_01,modnam_11,...
    constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,...
    err_list,obs_list,obs,ntrain,r_zlb,codes,p6,p7,tStart),params0,params_lo,params_hi,options);
  
  params0 = params;
  params1 = params;
  
 %% 
  [posterior filtered_errs like prior] = ...
    posteriorzlb(params1,params_labels,params_lo,params_hi,...
    modnam_00,modnam_10,modnam_01,modnam_11,...
    constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,...
    err_list,obs_list,obs,ntrain,r_zlb,codes,p6,p7,tStart);

  
  params1=params;
 [ hessian_reg stdh_reg hessian_fmin stdh_fmin ] = compute_hessian(xstory,fstory,30);

  save mle_estimates_fminsearch 

  
end



  
