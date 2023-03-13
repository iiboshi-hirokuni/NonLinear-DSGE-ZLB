
function [ys,check]=BNK_model_steadystate(junk,ys)

global M_

 paramfile_BNK_zlb
 
 load PARAM_EXTRA_EST;

% Steady states
 phi   = (epsilon-1)*(omega+csigma)/kappa/pi_star;  
 r_star  = exp(csigma*gamma_a)/cbeta; 
%  r_zlb = -(r_star*pi_star-1);

% pi_star_100  pi_star gamma_a  gamma_a_100  cbeta omega epsilon csigma  h chi kappa
% psi_pi psi_y rho_a rho_b rho_r sigma_a sigma_b sigma_r m1 m2;

nparams = size(M_.param_names,1);
for icount = 1:nparams
    eval(['M_.params(icount) = ',M_.param_names(icount,:),';'])
%     eval(['M_.params(icount) = ',M_.param_names(icount,:)])
end


% Steady states

check=0;

 y_t =0; % log((epsilon - 1)/(epsilon*chi)*(1-h*exp(-gamma_a))^(-csigma))^(1/(csigma+omega));
 y_star_t =0; %log((epsilon - 1)/(epsilon*chi)*(1-h*exp(-gamma_a))^(-csigma))^(1/(csigma+omega));
 
 pi_t = 0; %pi_star;
 r_star_t = 0; %r_star;
 r_t = 0; %(r_star_t*pi_star) ;
 rnot = 0 ;

 mu_t = 0;
 zb_t = 0;

data_y = 0;
data_r = 0;
data_rnot = 0;

% data_pi = 100*(pi_star-1);
% data_r = 100*(r_star*pi_star-1);
% data_rnot = 100*(r_star*pi_star-1);

data_pi = 1*(pi_star-1);
% data_r = 1*(r_star*pi_star-1);
% data_rnot = 1*(r_star*pi_star-1);

fprintf('\n h \t \t pi_star r_star r_zlb  sigma gamma_a  omega ');
fprintf('\n %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f \n', h,pi_star, r_star, r_zlb, csigma, gamma_a, omega )
fprintf('\n kappa \t psi_pi  psi_y  rho_a  rho_b  rho_r ');
fprintf('\n %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f \n', kappa, psi_pi, psi_y, rho_a, rho_b, rho_r )

% var b bnot c ec lb maxlev y ;
ys = [ pi_t y_t y_star_t r_t rnot r_star_t mu_t zb_t ...
       data_y data_pi data_r data_rnot ]' ;



