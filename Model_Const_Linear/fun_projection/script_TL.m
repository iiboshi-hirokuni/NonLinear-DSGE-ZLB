% Time Iteration (Linear Interpolation):
% Canonical New Keynesian Model (Rotemberg Pricing) with Capital 
%   -Imposes the zero lower bound on the interest rate

clear all
clc

warning('off','all'); % Œx‚Ì”ñ•\Ž¦

disp('Start NK Model with Time Iteration (Linear Interpolation)')
addpath('./Toolbox')

tstart = tic;                           % Job timer start
%--------------------------------------------------------------------------
% Initialize Policy Functions
%--------------------------------------------------------------------------
% Specify initial conjecture
%  guess: linear policies
%  current: current nonlinear policy functions
O.loadpf = 'guess';

% %%  constraint of ZLB 
zlbflag = 1 %1;
policy_flag = 1;  % 1:type 1, 2:type 2

% Load parameters, steady state and grids
P = parameters;
S = steadystate(P);

% Specify grid options
O.rbound = [0.99  1.1];
O.MPbound = [-0.01  0.01];
O.bbound = [-0.025  0.025];   % mu^b_t = log(Z_t+1/Z_t) - gamma_b
O.abound = [-0.025  0.025];   % mu^a_t = log(A_t+1/A_t) - gamma_a

O.r_pts = 5; 
O.MP_pts = 5; 
O.b_pts = 11; 
O.a_pts = 9; 

O.e_pts = 1; % number of node of shock b for stochastic foresight
             %  1-> perfect foresight
O.u_pts = 1; % shock of a
O.v_pts = 1; % shock of  MP

% Load discretized state space
G = grids_even(O,P);

% Retrieve initial policy functions
if strcmp(O.loadpf,'guess')
    pf = guess_TL(P,S,G);    %% calculate linear version of policy functions
elseif strcmp(O.loadpf,'current')  
    load('pf_TL.mat','pf')
end  

%%  Plot linear version of polcy functions
%    type = 'linear';
%    shock_type = 'a_shock';
%    plot_pf;
%    
%    shock_type = 'MP_shock';
%    plot_pf;
   
   
% %--------------------------------------------------------------------------
% % Update Policy Functions
% %--------------------------------------------------------------------------

it = 1;                                 % Iteration Counter
converged = 0;                          % Convergence Flag
dist_max = 0;                           % Max distance vector
 perbind = 0;                           % percentage of ZLB 
 no_solution = 0;                       % when number of iteration is over 50, then there is no solution. 
 
% 
% % Preallocate arrays to store policy function updates
pf_y_up = zeros(G.griddim);
pf_y_star_up = zeros(G.griddim);
pf_pi_up = zeros(G.griddim);
pf_r_up = zeros(G.griddim);
 
% % Time iteration/linear interpolation algorithm
% 
while converged == 0 && perbind < 50 && no_solution == 0
    istart = tic;                       % Iteration timer start
   
    %     parfor i = 1:G.nodes  %% for pararell computing
    for i = 1:G.nodes
%     
%         % Find optimal policy functions on each node.
%         % csolve finds the zeros of 'eqm'
%         % Start csolve with the current policy function  
          start = [pf.y(i) pf.pi(i) pf.y_star(i) ]';
          state = [G.r_gr(i) G.b_gr(i) G.a_gr(i) G.MP_gr(i)];
       
%         RR =eqm_TL(start,state,pf,P,S,G,zlbflag);
%         RR
        
         %  csolve(FUN,x,gradfun,crit,itmax,varargin)
         argzero = csolve('eqm_TL',start,[],1e-3,5,state,pf,P,S,G,zlbflag,policy_flag);
%  
%       % Store updated policy functions       
          pf_y_up(i) = argzero(1);
          pf_pi_up(i) = argzero(2);
          pf_y_star_up(i) = argzero(3);
   
    %% 
     Y = argzero(1);
     pie = argzero(2);
     Y_star = argzero(3);
%      Z_b_pf = P.rho_b * state(2) ;   %%  2017/Jan/27 C³
%      Bpf = exp(Z_b_pf - state(2) );   %%  2017/Jan/27 C³
     Mu_a_pf = P.rho_a * state(3) ;
     r_star_p = S.r_star*(1+P.sigma*Mu_a_pf- (1-P.rho_b)*state(2) ); %%  2017/Jan/27 C³
     R_star = state(1)^P.rho_r * (r_star_p*S.pi_star...
            * (pie/S.pi_star)^P.psi_pi*(Y/Y_star)^P.psi_y)^(1-P.rho_r)*exp(state(4));
     Rp = R_star;
     if zlbflag == 1
         Rp = max(1,Rp);
     end        
     pf_r_up(i) = Rp;
    
%     display(i)
    end
 
     % Find where ZLB binds   
     %   Interest rate rule  
      locs = find(pf_r_up == 1);
     %   Percent nodes binding
      perbind = 100*numel(locs)/G.nodes;
 
    % Policy function distances
    dist_y = abs(pf_y_up - pf.y);
    dist_pi = abs(pf_pi_up - pf.pi);
    dist_y_star = abs(pf_y_star_up - pf.y_star);
    

%     % Maximum distance   
    dist_max(it) = max([dist_y(:)' dist_pi(:)' dist_y_star(:)']);
% 
    % Update policy functions
    pf.y = pf_y_up;
    pf.pi = pf_pi_up;
    pf.y_star = pf_y_star_up;
    pf.r = pf_r_up;
% 
%     % Check convergence criterion
    if it > 5 && all(dist_max(end-2:end) < P.tol)
        converged = 1;
    elseif it > 50
        no_solution = 1;
        disp('no_solution since the number of iteration is over criterion ');
    end
% 
    % Iteration Information and save policy functions to disk
    if mod(it,2) == 1 || converged == 1 
        it = itinfo(istart,tstart,1,it,dist_max(it));
        disp(['nodes binding ZLB are : ' num2str(perbind) '%']);
        save('pf_TL.mat','pf','O','P','S','G');
    else
        it = it+1;
    end  
   
end
  
% 

%%  Plot non-linear version of polcy functions
  type = 'nonlinear'
  shock_type = 'a_shock';
    plot_pf;
    
   type = 'nonlinear'
 shock_type = 'MP_shock';
    plot_pf; 
%     
% %% Plot inpulse response functions
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
  
  
  
  
  
  
  

