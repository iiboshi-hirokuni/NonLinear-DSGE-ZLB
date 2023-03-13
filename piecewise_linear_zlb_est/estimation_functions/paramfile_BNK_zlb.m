
m1 = 1;
m2 = 1;
h =  0.40;
pi_star = 1.01 ;
sigma = 1.5;% P.sigma = 1.5;
beta  = 0.9975;
chi   = 1.0;
gamma_a = 0.01;
omega = 2;
epsilon = 6;
kappa = 0.05;    

% Monetary Policy 
psi_pi = 1.5; 
psi_y  = 0.4; 

% Stochastic Processes
 rho_a = 0.5; 
 rho_b = 0.5;
 rho_r = 0.5;
 sigma_a = 0.05;
 sigma_b = 0.05;
 sigma_r = 0.05;

phi   = (epsilon-1)*(omega+sigma)/kappa/pi_star;  

% Steady states
pi_star = 1 ;
r_star  = exp(sigma*gamma_a)/beta; 

r_zlb = -(r_star*pi_star-1);

load PARAM_EXTRA_CALIBRATED

% load PARAM_EXTRA_BABY