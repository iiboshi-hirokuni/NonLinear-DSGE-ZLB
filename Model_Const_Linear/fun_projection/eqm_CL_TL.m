function RR =eqm_CL_TL(x,state,pf,P,S,G,zlbflag,policy_flag)

% R = eqm_TL(x,state,pf,P,S,G)
%   Outputs residuals of the equilibrium system of equations for time
%   iteration/linear interpolation method.
% Inputs:
%   x       :   Policy function values on node i
%   state   :   State variable values on node i
%   pf      :   Structure of policy functions
%   P       :   Structure of parameters
%   S       :   Structure of steady state values
%   G       :   Structure of grids
% Output:
%   R       :   Residuals
%

% Preallocate function output
RR = zeros(size(x));
Rdim = size(RR,2);

% State Values
R    = state(1);       % Interest rate state last period
% Mu_b = state(2);     % Preference shock
Z_b = state(2); 
Mu_a = state(3);     % Technology shock current period  
epMP = state(4);    % Monetary policy shock

%  
% pf.y =  S.y_star*(1 + linpf_y);
% pf.y_star =  S.y_star*(1 + linpf_y_star);
% pf.pi = S.pi_star + linpf_pi;
% pf.r  = S.R  + linpf_r;


for j = 1:Rdim   
    % Policy Function Guesses
    Y = x(1,j);         % Output current period    
    pie = x(2,j);       % Inflation current period 
    Y_star = x(3,j);    %  
    
    % logarithm % 
    y = log(Y/S.y_star);
    y_star = log(Y_star/S.y_star);
    
    %----------------------------------------------------------------------
    % Solve for variables
    %----------------------------------------------------------------------
        
    % Technology law of motion  % Eq(15)
    Mu_a_pf = P.rho_a * Mu_a + G.u_nodes;
    
    % Preference law of motion  % Eq(16)
%     Mu_b_pf = P.rho_b * Mu_b + G.e_nodes;
    Z_b_pf = P.rho_b * Z_b + G.e_nodes;
          
    % real rate Eq.(14)
    r_star_p = S.r_star*(1+P.sigma*Mu_a_pf+(1-P.rho_b)*Z_b);
    
    % Monetary Policy Rule (13)
    R_star = S.r_star*S.pi_star*(1+P.rho_r*(R/(S.r_star*S.pi_star)-1) ...
             + (1-P.rho_r)*(P.psi_pi*(pie-S.pi_star)/S.pi_star +P.psi_y*(y-y_star))+epMP);
    
    % Monetary Policy Rule (4)
    Rp = R_star;
    if zlbflag == 1
        Rp = max(1,Rp);
    end        

    if policy_flag==1
        R_lag = R_star;
    else 
        R_lag = Rp;
    end
  
    %----------------------------------------------------------------------
    % Linear interpolation of the policy variables 
    %----------------------------------------------------------------------              
    
      [Ypf, piepf, Y_s_pf] = ...
         allterp430(G.r_grid,G.b_grid,G.a_grid, G.MP_grid,...
                R_lag, Z_b_pf, Mu_a_pf,  G.w_nodes,  pf.y, pf.pi, pf.y_star); 
%                 R_lag, Mu_b_pf, Mu_a_pf,  G.w_nodes,  pf.y, pf.pi, pf.y_star); 
           
    % logarithm %        
    ypf = log(Ypf/S.y_star);
    y_star_pf = log(Y_s_pf/S.y_star);    
         
    %----------------------------------------------------------------------
    % First-order conditions
    %----------------------------------------------------------------------
    
    % Consumption Euler Equation  Eq(11)
    RR(1,j) = -(pie-S.pi_star)+P.beta*exp((1-P.sigma)*P.gamma_a)*(piepf-S.pi_star)...
            + (P.epsilon-1)*(P.omega+P.sigma)/(P.phi*S.pi_star)*(y-y_star);
            
  
    % Firm Pricing Equation  Eq(12) 
    RR(2,j) = -(y-y_star) + ( ypf - y_star_pf) ...
        - 1/P.sigma*( (Rp-S.r_star*S.pi_star)/(S.r_star*S.pi_star) ...
        - (piepf-S.pi_star)/S.pi_star - (r_star_p-S.r_star)/S.r_star ) ;
       

end
