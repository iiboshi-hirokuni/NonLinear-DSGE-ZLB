%-------------------------------------------
% Housekeeping
%-------------------------------------------
clear all
restoredefaultpath

tic
tStart = tic;
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


tstart = tic;                           % Job timer start

%   load mle_estimates_temp_test
%   load mle_initial_guess
  
  load mle_estimates_fminsearch 
datavec=[]; xstory=[]; fstory=[];
irep=1;


target_average=1/3;
init_draws=10;
ndraws=20000;
load_switch = 1; % 1: load file,  0: not load file 
jump= 50.0;

[theta_history, fval_history, accept_average] = ...
   mhalgo('posteriorzlb',params1,hessian_fmin,target_average,init_draws,ndraws,load_switch,jump,...
   params_labels,params_lo,params_hi,...
   modnam_00,modnam_10,modnam_01,modnam_11,...
   constraint1_difference, constraint2_difference,...
   constraint_relax1_difference, constraint_relax2_difference,...
   err_list,obs_list,obs,ntrain,r_zlb, codes, p6, p7,tStart);



save metropolis_estimates


