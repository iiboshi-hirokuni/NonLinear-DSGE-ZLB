function [T, M, eu] = linmodel(P,S,V)

% [T, M, eu] = linmodel(P,S,V)
%   Solves the log-linear model with GENSYS
% Inputs:
%     P     :   Structure of parameters
%     S     :   Structure of steady state values
%     V     :   Structure of variable locations and names
% Output:
%     T     :   Transition matrix
%     M     :   Impact matrix
%     eu    :   [existence, uniqueness]

%---------------------------------------------
%   Initialize GENSYS components
%---------------------------------------------
GAM0j = zeros(V.nvar);
GAM1j = zeros(V.nvar);
PSI0j = zeros(V.nvar,V.nshock);
PPIj = zeros(V.nvar,V.nfore);
CC = zeros(V.nvar,1);

 % steady state
%  c_t=0;
%  r_t=0;
%  pi_t=0;
%  y_t =0;
%  g_t =0;
%  z_t =0;
%  r_t_lg = 0;
%  Et_c_t1 =0;
%  Et_pi_t1 =0;
%  Et_y_t1 =0; 
%  Et_z_t1 =0;

GAM0j(1,1) = -((P.omega+P.sigma)*(P.epsilon-1.0))/(P.phi*S.pi_star);
GAM0j(1,2) = ((P.omega+P.sigma)*(P.epsilon-1.0))/(P.phi*S.pi_star);
      GAM0j(1,3) = 1.0;
      GAM0j(1,11) = -P.beta*exp(-P.gamma_a*(P.sigma-1.0));  %%  2017/Jan/27 èCê≥
      GAM0j(2,1) = 1.0;
      GAM0j(2,2) = -1.0;
      GAM0j(2,4) = 1.0/(S.pi_star*S.r_star*P.sigma);
      GAM0j(2,6) = -1.0/(S.r_star*P.sigma);
      GAM0j(2,9) = -1.0;
      GAM0j(2,10) = 1.0;
      GAM0j(2,11) = -1.0/(S.pi_star*P.sigma);
      GAM0j(3,1) = P.psi_y*(P.rho_r-1.0);
      GAM0j(3,2) = -P.psi_y*(P.rho_r-1.0);
      GAM0j(3,3) = (P.psi_pi*(P.rho_r-1.0))/S.pi_star;
      GAM0j(3,4) = 1.0/(S.pi_star*S.r_star);
      GAM0j(3,5) = -P.rho_r/(S.pi_star*S.r_star);
      GAM0j(4,6) = 1.0;
      GAM0j(4,7) = -S.r_star*P.rho_a*P.sigma;
      GAM0j(4,8) = S.r_star*(P.rho_b-1);   %%  2017/Jan/27 èCê≥
      GAM0j(5,2) = -1.0;
      GAM0j(5,7) = -P.rho_a;
      GAM0j(5,10) = 1.0;
      GAM0j(6,1) = 1.0;
      GAM0j(7,3) = 1.0;
      GAM0j(8,7) = 1.0;
      GAM0j(9,8) = 1.0;
      GAM0j(10,2) = 1.0;
      GAM0j(10,7) = -1.0;
      GAM0j(11,5) = 1.0;

 GAM1j(6,9) = 1.0;
      GAM1j(7,11) = 1.0;
      GAM1j(8,7) = P.rho_a;
      GAM1j(9,8) = P.rho_b;
      GAM1j(10,2) = 1.0;
      GAM1j(11,4) = 1.0;

      PSI0j(3,3) = 1.0;
      PSI0j(8,1) = 1.0;
      PSI0j(9,2) = 1.0;

      PPIj(6,1) = 1.0;
      PPIj(7,2) = 1.0;

% CC = zeros(V.nvar,1);
%---------------------------------------------
%   Solve Linear Model
%---------------------------------------------
[T,~,M,~,~,~,~,eu] = gensys(GAM0j,GAM1j,CC,PSI0j,PPIj,1);

% T
% GAM0j
% GAM1j
% PPIj
% PSI0j

