% 
%  Sequential Monte Carlo Methods with Particle Filter
%  Canonical New Keynesian Model (Rotemberg Pricing)    
%   -Imposes the zero lower bound on the interest rate
%

clear all
clc

addpath('./Toolbox')
addpath('./fun_projection')
addpath('./fun_particle_filter')
addpath('./fun_smc')
addpath('./fun_prior')
% addpath()

tstart = tic;                           % Job timer start

disp('Start SMC^2 ')
ncores  =  16

%% setting 
zlbflag = 0 %1;
policy_flag = 0  % 1:type 1, 2:type 2
nsim      = ncores*50           % # of particles of parameters
nstage     =    5           % # of stages
nparticles =  10000   % # of particles of variables

m_err_flag = 1 
 
% npara      =    20;          % # of parameters
% cc1        =   0.035 ;       % adjustment coefficient of MH

%% parallel  
 delete(gcp('nocreate')) %   
 parpool('local', ncores)

main_smc2_I0_def



%% end


