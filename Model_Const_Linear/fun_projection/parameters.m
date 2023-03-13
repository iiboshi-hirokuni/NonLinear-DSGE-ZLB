function P = parameters

P.pi_star = 1 ;

P.sigma = 1.5;
P.beta  = 0.995;
P.chi   = 1.2;
P.gamma_a = 0.00;
% P.gamma_b = 0.00;  %%  2017/Jan/27 C³
P.omega = 1;
P.epsilon = 6;
P.kappa = 0.024;    %%% a la Woodford (2003) %%% MODIFIED BY UEDA %%%
P.phi   = (P.epsilon-1)*(P.omega+P.sigma)/P.kappa/P.pi_star;  %%% MODIFIED BY UEDA %%%

% Monetary Policy 
P.psi_pi = 1.5; 
P.psi_y  = 0.5/4; 

% Stochastic Processes
P.rho_a = 0.5; 
P.rho_b = 0.5;
P.rho_r = 0.5;
P.sigma_a = 0.5;
P.sigma_b = 0.5;
P.sigma_r = 0.5;

% Steady states
% P.pi_star = 1 ;
 P.r_star  = exp(P.sigma*P.gamma_a)/P.beta; %%  2017/Jan/27 C³

% Algorithm
P.tol     = 1e-2; %10;     % Convergence criterion
P.ite      = 30;   % max number of iteration  
P.disp_it  =  0;            % display number of iterations, on ->1, off->0
