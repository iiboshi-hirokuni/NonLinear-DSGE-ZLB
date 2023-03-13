% Time Iteration (Linear Interpolation):
% Canonical New Keynesian Model (Rotemberg Pricing) without Capital 
%   -Imposes the zero lower bound on the interest rate

clear all
clc


set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');

warning('off','all'); % åxçêÇÃîÒï\é¶

disp('Start NK Model with Time Iteration (Linear Interpolation)')
addpath('./Toolbox')
addpath('./fun_projection')
addpath('./fun_particle_filter')
addpath('./fun_hist_decomp')
addpath('./fun_smc')
addpath('./fun_prior')

tstart = tic;                           % Job timer start


%%  setting options 
zlbflag     =  0  %1;
policy_flag =  0  % 1:type 1, 2:type 2
nsim      =  1200;           % # of particles of parameters
nstage     =    10;           % # of stages
nparticles =  40000;   % # of particles of variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npara      =    20;          % # of parameters

% data_type = 'artificial' 
% data_s = 1

data_type = 'actual'

switch data_type
    case 'actual'
      file_name = ['./output/Feb_17_2017/save_para_I0_def_'...
                    num2str(zlbflag)  num2str(policy_flag) '_' ...
                    num2str(nsim) '_' num2str(nparticles)  ];  
      datapath = './data/';
      datafilename = 'data_jpn_1_def';         % actual data
      yy1 = csvread(strcat(datapath, [ datafilename '.csv' ] ), 1, 1);
      yy = [ yy1(:,2) yy1(:,3) yy1(:,4)];
      
      %%  after 2016Q3
      yy2 = csvread(strcat(datapath, 'data-2017.csv'  ), 1, 1);
      yy2017 = [ yy2(:,1) yy2(:,2)/4 yy2(:,3)/4 ];
      
      yy = [ yy; yy2017]; 
      
      
    case 'artificial'
        file_name = ['./output/Jun052017/save_para_I0_def_' ...
                     num2str(zlbflag)  num2str(policy_flag) '_' ...
                     num2str(nsim) '_' num2str(nparticles) '_1-20-20' ]; 
     load( [ 'art_data.mat' ] );      
      s0010110 = 201*data_s 
      e0010110 = s0010110 + 120
      yy = [ y(s0010110:e0010110) pi(s0010110:e0010110) r(s0010110:e0010110)];
      shock_art = [ s_mp(s0010110:e0010110) s_a(s0010110:e0010110) s_b(s0010110:e0010110) ];
end

ZZ = [100; 100; 100]; % coefficients between observables and endogenous variables of the measurement equation

                                  
 load(file_name);
 P = parameters;      
 P = parameters_update(mean(para_Resamp,1),P);
 P.tol     = 1e-3; %10;     % Convergence criterion
 P.disp_it = 1;
 P.ite      = 100; 
%   P.omega   = 10;
%   P.kappa   = 0.10;
 
 S = steadystate(P);

 Pr = fun_prior_setting;

% Specify grid options
O.loadpf = 'guess';
O.rbound = [0.90  1.1];
O.MPbound = [-0.03  0.03];
O.bbound = [-0.1  0.1 ];   % mu^b_t = log(Z_t+1/Z_t) - gamma_b
O.abound = [-0.1  0.1 ];   % mu^a_t = log(A_t+1/A_t) - gamma_a

O.r_pts = 21;  % 11
O.MP_pts = 11; % 5
O.b_pts = 9; 
O.a_pts = 9; 

% Load discretized state space
G = grids_even(O,P);

%%--------------------------------------------------------------------------
%% Update Policy Functions
%%--------------------------------------------------------------------------

  disp('Start NK Model with Time Iteration (Linear Interpolation)')
  P.disp_it=1;
   [ pf, ~ ] = solve_model_par(P,S,G,zlbflag,policy_flag,tstart);
  
  if zlbflag==0 && policy_flag == 0 
       type_model = char('model0');
  elseif zlbflag==1 && policy_flag == 1 
       type_model = char('model1');     
  elseif zlbflag==1 && policy_flag == 2
       type_model = char('model2');
  end      
%      
% Å@ save( [ 'save_pf/' type_model '_pf.mat'], 'pf','O', 'G','P','S');
%   
% 
% %%  conduct post estimations
  stats_sample_para;

  plot_dist(parasim, npara,Pr,nstage)
% % 
  plot_trace(parasim, npara,Pr,nstage)

%% conduct particle filter  

%%
switch data_type
    case 'actual'
     cal_particle;   %% it is necessary 
    case 'artificial'
     cal_particle_art;
end 
%   
     cal_hist_IRF_2;   %% option
%     
% 
% %%  Plot non-linear version of polcy functions
  type = 'nonlinear'
  shock_type = 'a_shock';
    plot_pf2;
%     
   type = 'nonlinear'
%  shock_type = 'MP_shock';
 shock_type = 'a_shock';
    plot_pf; 
%  
 
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
  
  
  
  
  
  
  

