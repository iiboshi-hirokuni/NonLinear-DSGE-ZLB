

warning off;


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


// 2. ÉpÉâÉÅÅ[É^ÇÃêÈåæ
parameters r_zlb r_star phi
           pi_star gamma_a cbeta omega epsilon csigma  h chi kappa
           psi_pi psi_y  rho_a rho_b rho_r STD_a STD_b STD_r
           m1 m2; 

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------


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
r_t = r_zlb;


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


end;




