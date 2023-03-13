
m1 = 1;
m2 = 1;

h =  0.50;
%pi_star_100 = 0 ;
csigma = 1.5;
cbeta  = 0.99875;
chi   = 1.0;
%gamma_a_100 = 0.1;
gamma_a = 0.0025 ;
omega = 3;
epsilon = 6;
kappa = 0.05;    

% Monetary Policy 
psi_pi = 1.5; 
psi_y  = 0.125; 
pi_star = 1.005 ;
% Stochastic Processes
rho_a = 0.5; 
rho_b = 0.5;
rho_r = 0.5;
STD_a = 0.005;
STD_b = 0.005;
STD_r = 0.005;

r_star  = exp(csigma*gamma_a)/cbeta; 
r_zlb = -(r_star*pi_star-1);
%  pi_star = 1;
%   load PARAM_EXTRA_CALIBRATED
%   


  