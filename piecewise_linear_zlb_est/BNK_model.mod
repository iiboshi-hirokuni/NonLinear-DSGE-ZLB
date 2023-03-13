
%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

%clear all;
close all;
%clear all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var pi_t y_t y_star_t r_t rnot r_star_t mu_t zb_t
    data_y data_pi data_r data_rnot;

varexo  eps_a eps_b eps_r;


// 2. パラメータの宣言
parameters r_zlb r_star phi
           pi_star gamma_a cbeta omega epsilon csigma  h chi kappa
           psi_pi psi_y  rho_a rho_b rho_r STD_a STD_b STD_r
           m1 m2; 

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
m1 = 1;
m2 = 1;

% h =  0.40;
% pi_star = 1.01 ;
% csigma = 1.5;  
% cbeta  = 0.9975;
% chi   = 1.0;
% gamma_a = 0.01;
%omega = 2;
%epsilon = 6;
%kappa = 0.05;    

% Monetary Policy 
%psi_pi = 1.5; 
%psi_y  = 0.4; 

% Stochastic Processes
%rho_a = 0.5; 
%rho_b = 0.5;
%rho_r = 0.5;
%STD_a = 0.05;
%STD_b = 0.05;
%STD_r = 0.05;

%phi   = (epsilon-1)*(omega+csigma)/kappa/pi_star;  

% Steady states
%pi_star = 1 ;
%r_star  = exp(csigma*gamma_a)/cbeta;
%r_zlb = -(r_star*pi_star-1);

%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------


model(linear);


 // Eq (22) NKPC, page 4
pi_t  = m2*cbeta*exp((1-csigma)*gamma_a)*(pi_t(+1))
   +(epsilon-1)/(phi*pi_star)*(omega+csigma/(1-h*exp(-gamma_a)))*(y_t - y_star_t );

// Eq (23), page 4
(y_t - y_star_t) =  m1*( y_t(+1) - y_star_t(+1) -h*exp(-gamma_a)*(y_t -y_star_t)
-(1-h*exp(-gamma_a))/csigma*( r_t/(r_star*pi_star) -(pi_t(+1))/pi_star - (r_star_t)/r_star )); 
    
// Eq (24), Taylor Rule page 4
rnot/(r_star*pi_star)  =  rho_r*(r_t(-1)/(r_star*pi_star) )
 +(1-rho_r)*( psi_pi*(pi_t)/pi_star + psi_y*(y_t - y_star_t ) ) +eps_r ;

//  ZLB constraint
r_t = rnot;


// Eq (25), page 4
(r_star_t)/r_star =  
(1/m1-(omega + csigma/(1-h*exp(-gamma_a)))^(-1)* 
csigma*h*exp(-gamma_a)/(1-h*exp(-gamma_a)))*omega*y_star_t 
+ (1-m1*rho_b)/m1*zb_t + ( (omega + csigma/(1-h*exp(-gamma_a)))^(-1)
* omega*csigma*h*exp(-gamma_a)/(1-h*exp(-gamma_a))+csigma)*rho_a*mu_t;
    
// Eq (26), page 4
y_star_t  = (omega + csigma/(1-h*exp(-gamma_a)))^(-1)
 *csigma*h*exp(-gamma_a)/(1-h*exp(-gamma_a))*(y_t(-1) - mu_t);


// Eq (27), page 4
mu_t = rho_a*mu_t(-1)+  eps_a;

// Eq (30), page 4
zb_t =  rho_b*zb_t(-1) +  eps_b;    

// Eq (26), page 4
//y_star_t = (omega + csigma/(1-h*exp(-gamma_a)))^(-1)*csigma*h*exp(-gamma_a)/(1-h*exp(-gamma_a))*(y_t(-1)-mu_t);

//% observation equations
data_y   = 1*(mu_t + y_t - y_t(-1));
data_pi  = 1*(pi_t+pi_star -1);
data_r     =1* r_t;
data_rnot  =1* rnot;
%data_r   = 1*(r_t + r_star*pi_star-1);
%data_rnot  = 1*(rnot + r_star*pi_star-1 );

end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

// 全てゼロになることを念のためチェック
%steady;

// モデルのチェック
%check;

shocks ;
var eps_a =STD_a^2;
var eps_b =STD_b^2;
var eps_r =STD_r^2;
end ;

stoch_simul(order=1,noprint,nomoments,irf=0) data_y data_pi data_r;

//estimated_params;
// PARAM NAME, INITVAL, LB, UB, PRIOR_SHAPE, PRIOR_P1, PRIOR_P2, PRIOR_P3, PRIOR_P4, JSCALE
// PRIOR_SHAPE: BETA_PDF, GAMMA_PDF, NORMAL_PDF, INV_GAMMA_PDF
% stderr e_a       ,0.5  ,,,inv_gamma_pdf, 0.02,5;
% stderr e_b       ,0.5  ,,,inv_gamma_pdf, 0.02,5;
% stderr e_r       ,0.5  ,,,inv_gamma_pdf, 0.02,5;
% csigma , ,     ,     ,NORMAL_PDF, 1.5, 0.3;
% gamma_a_100 , ,     ,     ,NORMAL_PDF, 0.0, 0.05;
% omega , ,       ,     ,NORMAL_PDF, 3.0, 0.5;
% kappa , ,       ,     ,NORMAL_PDF, 0.05, 0.006; 
% pi_star_100 , ,     ,     ,NORMAL_PDF, 0, 0.5;
% psi_pi  , ,     ,     ,gamma_pdf,1.500 ,0.15 ;
% psi_y  , ,     ,      ,gamma_pdf,0.125 ,0.025 ;
% rho_a   , ,     ,     ,beta_pdf ,0.500,0.200;
% rho_b   , ,     ,     ,beta_pdf ,0.500,0.200;
% rho_r   , ,     ,     ,beta_pdf ,0.500,0.200;
% end;

% varobs data_y data_pi data_r;

% estimation(datafile='./data_jpn.mat', 
%        mode_compute=4,
%        %load_mh_file,
%        mode_check,
 %       graph_format = fig,
 %       plot_priors=0,
%        first_obs=1,
%        presample=0,
%        mcmc_jumping_covariance = prior_variance,   % identity_matrix|hessian|
%        mh_replic=1000,
%        mh_nblocks=2,
 %       mh_drop=0.45,
  %      mh_jscale=0.8,
   %     bayesian_irf,
%        forecast=100,
     %   irf=40
%        ) data_y data_pi data_r;

//shock_decomposition data_y data_pi data_r ;






