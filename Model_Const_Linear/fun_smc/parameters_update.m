function P = parameters_update(para,P)

P.sigma = para(1);
P.beta  = para(2);
P.chi   = para(3);
P.gamma_a =para(4)/100;
% P.gamma_b = para(5)/100;
P.omega   = para(6);
P.epsilon = para(7);
P.kappa =   para(8);    %%% a la Woodford (2003) %%% MODIFIED BY UEDA %%%

% Steady states
P.pi_star = 1 + para(9)/100 ;
P.phi   = (P.epsilon-1)*(P.omega+P.sigma)/P.kappa/P.pi_star;  %%% MODIFIED BY UEDA %%%
 
% % Monetary Policy 
 P.psi_pi = para(10); 
 P.psi_y  = para(11);
% 
% Stochastic Processes
 P.rho_a = para(12); 
 P.rho_b = para(13);
 P.rho_r = para(14);
 P.sigma_a = para(15);
 P.sigma_b = para(16);
 P.sigma_r = para(17);
 
 P.sigma_y = para(18);
 P.sigma_pi = para(19);
 P.sigma_R = para(20); 
