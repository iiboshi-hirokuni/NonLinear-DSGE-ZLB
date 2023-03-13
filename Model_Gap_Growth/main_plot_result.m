% Time Iteration (Linear Interpolation):
% Canonical New Keynesian Model (Rotemberg Pricing) without Capital 
%   -Imposes the zero lower bound on the interest rate

clear all
% clc

warning('off','all'); % åxçêÇÃîÒï\é¶

disp('Start NK Model with Time Iteration (Linear Interpolation)')
addpath('./Toolbox')
addpath('./fun_projection')
addpath('./fun_particle_filter')
addpath('./fun_smc')
addpath('./fun_prior')

tstart = tic;                           % Job timer start


%%  setting options 
zlbflag     =  1  %1;
policy_flag =  1  % 1:type 1, 2:type 2
nsim      =  1200;           % # of particles of parameters
nstage     =    10;           % # of stages
nparticles =  40000;   % # of particles of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npara      =    20;          % # of parameters

% Load parameters, steady state and grids
 file_name = ['./output/Feb_6_2017/save_para_I0_V4_'  num2str(zlbflag)  num2str(policy_flag) '_' ...
                                      num2str(nsim) '_' num2str(nparticles) ]; 
                                  
 load(file_name);
 P = parameters;      
 P = parameters_update(mean(para_Resamp,1),P);
 P.tol     = 1e-2; %10;     % Convergence criterion
%   P.omega   = 10;
%   P.kappa   = 0.10;
 
 S = steadystate(P);

 Pr = fun_prior_setting;

% Specify grid options
% O.loadpf = 'guess';

% Load discretized state space
G = grids_even(O,P);

%%--------------------------------------------------------------------------
%% Update Policy Functions
%%--------------------------------------------------------------------------

 P.disp_it=1;
  [ pf, ~ ] = solve_model(P,S,G,zlbflag,policy_flag,tstart);
  

%%  conduct post estimations
stats_sample_para;

  plot_dist(parasim, npara,Pr,nstage)
% 
  plot_trace(parasim, npara,Pr,nstage)

%% conduct particle filter  

%     cal_particle;


%%  Plot non-linear version of polcy functions
  type = 'nonlinear'
  shock_type = 'a_shock';
%     plot_pf2;
    
   type = 'nonlinear'
 shock_type = 'MP_shock';
%     plot_pf; 
 
 
% % %% Plot inpulse response functions
%   % preference shock  
%    shock_idx = 1;
%   impulse_response;
%   
%   % TFP shock 
%    shock_idx = 2;
%   impulse_response;
%   
%   % MP shock
%     shock_idx = 3;
%   impulse_response;
  
  
  
  
  
  
  

