% Time Iteration (Linear Interpolation):
% Canonical New Keynesian Model (Rotemberg Pricing) without Capital   
%   -Imposes the zero lower bound on the interest rate

% clear all
% clc
% 
% addpath('./Toolbox')
% addpath('./fun_projection')
% addpath('./fun_particle_filter')

tstart = tic;                           % Job timer start

disp('Start Particle Filter ')

%% option of parallel computing 
parallel_flag = 1; % 1-> on,  0-> off

%%  load policy functions
% load('pf_TL.mat','pf')
% load('pf_TL.mat')
 
%% data
datapath = './data/';

%% Data (year on year)
datafilename = 'data_jpn_1';         % actual data
yy1 = csvread(strcat(datapath, [ datafilename '.csv' ] ), 1, 1);

yy = [ yy1(:,1) yy1(:,3) yy1(:,4)];

ZZ = [100; 100; 100]; % coefficients between observables and endogenous variables of the measurement equation

%%
Tobs = size(yy,1);

%% Covariance Matrices
 %  Measurement errors of output, inflation, interest rate
%    HH = 1.0*diag([ Pr.pmean(18)^2  1/1*Pr.pmean(19)^2  Pr.pmean(20)^2 ]);   
 
 % structural shocks
 QQ = eye(3);
     QQ(1,1)= P.sigma_b^2; % preference
     QQ(2,2)= P.sigma_a^2; % TFP
     QQ(3,3)= P.sigma_r^2; % monetary policy
   
     
%% initial setting of particles filter
   nshocks=3;          % # of shocks
%    nparticles = 2000;   % # of particles
npar = 4;
  disp( ['# of particles is ' num2str(nparticles) ]);
  
%    file_name = ['./data/particle_' num2str(nparticles) ];
%   load(file_name,'stock_shock','stock_state');
  stock_shock = zeros(nshocks,ceil(nparticles),Tobs); % matrix of 3 dimensions ( nshocks X nparticles X Tobs )
  
   rng(100)  % generator of random numbers is fixed
    
   for i =1:Tobs
        stock_shock(:,:,i) = mvnrnd(zeros(nshocks,1), eye(3), ceil(nparticles) )'; 
   end
   
   stock_state= [ rand(nshocks,nparticles); randn(1,nparticles) ]; % initial values of state variables (t=0)
   
   
 %% conduct particle filter
  switch parallel_flag 
      case  1    
    [ log_lik, ypf] = ...
                fun_ParticleFilter_parallel(yy,ZZ,HH,QQ, O, G, S, P, pf,...
                      nparticles, zlbflag, policy_flag, npar, stock_shock,stock_state) ; 
                  
      case  0                
    [ log_lik, ypf] = ...
                fun_ParticleFilter(yy,ZZ,HH,QQ, O, G, S, P, pf,...
                      nparticles, zlbflag, policy_flag, npar, stock_shock,stock_state) ;             
  end
                  
     disp( ['log likelihood is ' num2str(log_lik) ]);
                  
                  
%% time of computing
 dec = 10^1;
T = toc(tstart);
hh = floor(T/3600);
mm = floor(T/60)-60*hh;
ss = round((T-60*mm-3600*hh)*dec)/dec;

hh = num2str(hh);
mm = num2str(mm);
ss = num2str(ss);

display(['Total Computing time: ' hh 'h' mm 'm' ss 's']);      
         
                  
 %% plot graph               
    m = 1;     % # of layers of interval 
    a = 0.90;  % percentage of interval 
    
 V = variables; 
 
 i = V.y;
 graph_title={'Output: Observable vs. Particle Filtered variables'};
 series_y = fun_plot(yy,ypf,i,a,m,Tobs,graph_title,V,P,ZZ);
 
 i = V.pi;
 graph_title={'Inflation: Observable vs. Particle Filtered variables'};
 series_pi = fun_plot(yy,ypf,i,a,m,Tobs,graph_title,V,P,ZZ);
 
  i = V.r;
 graph_title={'Interest Rate: Observable vs. Particle Filtered variables'};
 series_r =fun_plot(yy,ypf,i,a,m,Tobs,graph_title,V,P,ZZ);
 
 i = V.r_star;
 graph_title={'Natural Rate of Interest, r^*_t : Particle Filter '};
 series_r_star =fun_plot(yy,ypf,i,a,m,Tobs,graph_title,V,P,ZZ);
 
 i = V.r_lg;
 graph_title={'Target Rate, R^*_t : Particle Filter '};
 series_R_lag = fun_plot(yy,ypf,i,a,m,Tobs,graph_title,V,P,ZZ);
 
 i = V.mu_a;
 graph_title={'TFP, \mu^a_t : Particle Filter '};
 series_mu_a = fun_plot(yy,ypf,i,a,m,Tobs,graph_title,V,P,ZZ);
 
 i = V.z_b;
 graph_title={'Preference, Z^b_t : Particle Filter '};
 series_z_b = fun_plot(yy,ypf,i,a,m,Tobs,graph_title,V,P,ZZ);
 
 i = 9;
 graph_title={'Monetary Policy Shock : Particle Filter '};
 series_shock_mp = fun_plot(yy,ypf,i,a,m,Tobs,graph_title,V,P,ZZ);
 
 
 if zlbflag==0 && policy_flag == 0
       type_model = char('./output/model0_gap_');
 elseif zlbflag==1 && policy_flag == 1
       type_model = char('./output/model1_gap_');
 elseif zlbflag==1 && policy_flag == 2
       type_model = char('./output/model2_gap_');
 end      
     
 save( [ type_model '_state.mat'], 'series_R_lag', 'series_mu_a', 'series_z_b', 'series_r_star', ...
       'series_r', 'series_y', 'series_pi', 'series_shock_mp', ...
       'pf', 'Tobs','G', 'S', 'P', 'zlbflag', 'policy_flag' );
   
% save(  'state.mat', 'series_R_lag', 'series_mu_a', 'series_z_b', 'series_r_star', ...
%        'series_r', 'series_y', 'series_pi', 'series_shock_mp', ...
%        'pf', 'Tobs','G', 'S', 'P', 'zlbflag', 'policy_flag' ); 
   
                  
                  