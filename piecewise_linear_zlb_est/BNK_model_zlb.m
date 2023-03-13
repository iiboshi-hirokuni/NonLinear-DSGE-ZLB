%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic;
global M_ oo_ options_ ys0_ ex0_ estimation_info
options_ = [];
M_.fname = 'BNK_model_zlb';
%
% Some global variables initialization
%
global_initialization;
diary off;
M_.exo_names = 'eps_a';
M_.exo_names_tex = 'eps\_a';
M_.exo_names = char(M_.exo_names, 'eps_b');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_b');
M_.exo_names = char(M_.exo_names, 'eps_r');
M_.exo_names_tex = char(M_.exo_names_tex, 'eps\_r');
M_.endo_names = 'pi_t';
M_.endo_names_tex = 'pi\_t';
M_.endo_names = char(M_.endo_names, 'y_t');
M_.endo_names_tex = char(M_.endo_names_tex, 'y\_t');
M_.endo_names = char(M_.endo_names, 'y_star_t');
M_.endo_names_tex = char(M_.endo_names_tex, 'y\_star\_t');
M_.endo_names = char(M_.endo_names, 'r_t');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_t');
M_.endo_names = char(M_.endo_names, 'rnot');
M_.endo_names_tex = char(M_.endo_names_tex, 'rnot');
M_.endo_names = char(M_.endo_names, 'r_star_t');
M_.endo_names_tex = char(M_.endo_names_tex, 'r\_star\_t');
M_.endo_names = char(M_.endo_names, 'mu_t');
M_.endo_names_tex = char(M_.endo_names_tex, 'mu\_t');
M_.endo_names = char(M_.endo_names, 'zb_t');
M_.endo_names_tex = char(M_.endo_names_tex, 'zb\_t');
M_.endo_names = char(M_.endo_names, 'data_y');
M_.endo_names_tex = char(M_.endo_names_tex, 'data\_y');
M_.endo_names = char(M_.endo_names, 'data_pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'data\_pi');
M_.endo_names = char(M_.endo_names, 'data_r');
M_.endo_names_tex = char(M_.endo_names_tex, 'data\_r');
M_.endo_names = char(M_.endo_names, 'data_rnot');
M_.endo_names_tex = char(M_.endo_names_tex, 'data\_rnot');
M_.param_names = 'r_zlb';
M_.param_names_tex = 'r\_zlb';
M_.param_names = char(M_.param_names, 'r_star');
M_.param_names_tex = char(M_.param_names_tex, 'r\_star');
M_.param_names = char(M_.param_names, 'phi');
M_.param_names_tex = char(M_.param_names_tex, 'phi');
M_.param_names = char(M_.param_names, 'pi_star');
M_.param_names_tex = char(M_.param_names_tex, 'pi\_star');
M_.param_names = char(M_.param_names, 'gamma_a');
M_.param_names_tex = char(M_.param_names_tex, 'gamma\_a');
M_.param_names = char(M_.param_names, 'cbeta');
M_.param_names_tex = char(M_.param_names_tex, 'cbeta');
M_.param_names = char(M_.param_names, 'omega');
M_.param_names_tex = char(M_.param_names_tex, 'omega');
M_.param_names = char(M_.param_names, 'epsilon');
M_.param_names_tex = char(M_.param_names_tex, 'epsilon');
M_.param_names = char(M_.param_names, 'csigma');
M_.param_names_tex = char(M_.param_names_tex, 'csigma');
M_.param_names = char(M_.param_names, 'h');
M_.param_names_tex = char(M_.param_names_tex, 'h');
M_.param_names = char(M_.param_names, 'chi');
M_.param_names_tex = char(M_.param_names_tex, 'chi');
M_.param_names = char(M_.param_names, 'kappa');
M_.param_names_tex = char(M_.param_names_tex, 'kappa');
M_.param_names = char(M_.param_names, 'psi_pi');
M_.param_names_tex = char(M_.param_names_tex, 'psi\_pi');
M_.param_names = char(M_.param_names, 'psi_y');
M_.param_names_tex = char(M_.param_names_tex, 'psi\_y');
M_.param_names = char(M_.param_names, 'rho_a');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_a');
M_.param_names = char(M_.param_names, 'rho_b');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_b');
M_.param_names = char(M_.param_names, 'rho_r');
M_.param_names_tex = char(M_.param_names_tex, 'rho\_r');
M_.param_names = char(M_.param_names, 'STD_a');
M_.param_names_tex = char(M_.param_names_tex, 'STD\_a');
M_.param_names = char(M_.param_names, 'STD_b');
M_.param_names_tex = char(M_.param_names_tex, 'STD\_b');
M_.param_names = char(M_.param_names, 'STD_r');
M_.param_names_tex = char(M_.param_names_tex, 'STD\_r');
M_.param_names = char(M_.param_names, 'm1');
M_.param_names_tex = char(M_.param_names_tex, 'm1');
M_.param_names = char(M_.param_names, 'm2');
M_.param_names_tex = char(M_.param_names_tex, 'm2');
M_.exo_det_nbr = 0;
M_.exo_nbr = 3;
M_.endo_nbr = 12;
M_.param_nbr = 22;
M_.orig_endo_nbr = 12;
M_.aux_vars = [];
M_.Sigma_e = zeros(3, 3);
M_.H = 0;
options_.linear = 1;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('BNK_model_zlb_static');
erase_compiled_function('BNK_model_zlb_dynamic');
M_.lead_lag_incidence = [
 0 5 17;
 1 6 18;
 0 7 19;
 2 8 0;
 0 9 0;
 0 10 0;
 3 11 0;
 4 12 0;
 0 13 0;
 0 14 0;
 0 15 0;
 0 16 0;]';
M_.nstatic = 6;
M_.nfwrd   = 2;
M_.npred   = 3;
M_.nboth   = 1;
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:3];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(12, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(3, 1);
M_.params = NaN(22, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 41;
M_.NNZDerivatives(2) = 0;
M_.NNZDerivatives(3) = -1;
warning off;
close all;
save('BNK_model_zlb_results.mat', 'oo_', 'M_', 'options_');


disp(['Total computing time : ' dynsec2hms(toc) ]);
