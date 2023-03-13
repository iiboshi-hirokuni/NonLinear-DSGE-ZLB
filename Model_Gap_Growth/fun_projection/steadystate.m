function S = steadystate(P)

% S = steadystate(P)
%   Computes the deterministic steady state
% Input:
%   P : Structure of parameters
% Output:
%   S : structure of steady state values

% inflation
S.pi_star = 1 ;

% potential output
S.y_star  = ((P.epsilon - 1)/(P.epsilon*P.chi))^(1/(P.sigma+P.omega)) ;

% real interest rate
S.r_star  = exp(P.sigma*P.gamma_a)/P.beta; %%  2017/Jan/27 C³

S.R = S.pi_star*S.r_star;



