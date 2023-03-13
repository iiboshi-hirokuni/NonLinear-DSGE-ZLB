function RR =eqm_TL(x,state,pf,P,S,G,zlbflag,policy_flag)

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
Z_b = state(2);     % Preference shock
Mu_a = state(3);     % Technology shock current period  
epMP = state(4);    % Monetary policy shock

%  
for j = 1:Rdim   
    % Policy Function Guesses
    Y = x(1,j);         % Output current period    
    pie = x(2,j);       % Inflation current period 
    Y_star = x(3,j);    %  
    
    %----------------------------------------------------------------------
    % Solve for variables
    %----------------------------------------------------------------------
    
    % consumption Eq.(5)
    Cp = Y - P.phi/2*(pie-S.pi_star)^2*Y;
    
    % Lambda, Eq.(1)
    Lp= Cp^(-1*P.sigma);
    
    % A*_t/A_t Eq.(6)
    L_star = Y_star^(-1*P.sigma);
    
    % Preference law of motion  % Eq(30)
    Z_b_pf =real( P.rho_b * Z_b + G.e_nodes);     %%  2017/Jan/27 C³
    Bpf = exp(Z_b_pf - Z_b);                  %%  2017/Jan/27 C³
        
    % Technology law of motion  % Eq(10)
    Mu_a_pf = real(P.rho_a * Mu_a + G.u_nodes);
    Apf = exp(Mu_a_pf+P.gamma_a); 
        
    % real rate Eq.(31)
    r_star_p = S.r_star*(1+P.sigma*Mu_a_pf+ (1-P.rho_b)*Z_b ); %%  2017/Jan/27 C³
    
    % Monetary Policy Rule (4)
    R_star = R^P.rho_r * (S.r_star*S.pi_star...
        * (pie/S.pi_star)^P.psi_pi*(Y/Y_star)^P.psi_y)^(1-P.rho_r)*exp(epMP);
    
    % Monetary Policy Rule (4)
    Rp = R_star;
    if zlbflag == 1
        Rp = max(1,Rp);
    end        

    if policy_flag==1
        R_lag = real(R_star);
    else 
        R_lag = real(Rp);
    end
  
    %----------------------------------------------------------------------
    % Linear interpolation of the policy variables 
    %----------------------------------------------------------------------              
    
      [Ypf, piepf, Y_s_pf] = ...
         allterp430(G.r_grid,G.b_grid,G.a_grid, G.MP_grid,...
                R_lag, Z_b_pf, Mu_a_pf,  G.w_nodes,  pf.y, pf.pi, pf.y_star); 
                          
          
    %----------------------------------------------------------------------        
    % Solve for variables inside expectations
    %----------------------------------------------------------------------    
    
     % consumption at t+1, Eq.(5)
    Cpf = Ypf - P.phi/2*(piepf-S.pi_star)^2*Ypf;
    
%     % Lambda, at t+1,  Eq.(1)
    Lpf= Cpf^(-1*P.sigma);
    
    % A*_t+1/A_t+1 Eq.(6)
    L_star_pf = Y_s_pf^(-1*P.sigma);
        
    %----------------------------------------------------------------------
    % Numerical integration using Gauss-Hermite Quadrature
    %----------------------------------------------------------------------
    
	 % Expectations;   RHS of Eq(2)
        Econs_1 = (Lpf/Lp)*(Apf^(-1*P.sigma) )*Bpf*Rp/piepf;
                    
     % Expectations;   RHS of Eq(6)
        Econs_2 = (L_star_pf/L_star)*(Apf^(-1*P.sigma) )*Bpf*r_star_p;               

     % Expectations; 4th term of LHS of Eq(3)
        Efp =  (Lpf/Lp)*(Apf^-P.sigma)*Bpf*P.phi...
             *(piepf-S.pi_star)*piepf*Ypf/Y*Apf;        
     
    %----------------------------------------------------------------------
    % First-order conditions
    %----------------------------------------------------------------------
    
    % Consumption Euler Equation  Eq(2)
    RR(1,j) = ( 1-P.beta*Econs_1 ); 
        
   % Consumption Euler Equation  Eq(6)
    RR(2,j) = 1-P.beta*Econs_2;
      
    % Firm Pricing Equation  Eq(3) --> move all RHS to LHS
    RR(3,j) = 1-P.phi*(pie-S.pi_star)*pie...
           -P.epsilon*(1-P.chi*Y^P.omega/Lp-P.phi/2*(pie-S.pi_star)^2)...
           +P.beta*Efp; 
             
%       disp(  [ Y_star^(-1*(P.sigma+P.omega)) Cp Y^P.omega Lp Y^P.omega/Lp] );
       

end
