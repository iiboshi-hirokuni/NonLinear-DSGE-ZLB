function imp = fun_imp(G, S, P, pf, state, shock, period,zlbflag,policy_flag)

 
% State Values
R    = state(1);       % Interest rate state last period
Z_b = state(2);     % Preference shock   %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
Mu_a = state(3);     % Technology shock current period  
epMP = state(4);    % Monetary policy shock

shock_b = shock(1);
shock_a = shock(2);
shock_MP = shock(3);

imp =zeros(period,size(state,1));

for j = 1:period     
     
    % Preference law of motion  % Eq(9)
    Z_bp = P.rho_b * Z_b + shock_b;  %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
%     Bp = exp(Mu_bp+P.gamma_b);
        
    % Technology law of motion  % Eq(10)
    Mu_ap = P.rho_a * Mu_a + shock_a;
%     Ap = exp(Mu_ap+P.gamma_a); 
    
    % Monetary Policy Rule 
    if (policy_flag==2) && (zlbflag == 1)
        R = max(1,R);
    end  
      
  %----------------------------------------------------------------------
    % Linear interpolation of the policy variables 
    %----------------------------------------------------------------------              
    
      [Yp, piep, Y_star_p] = ...
         allterp430(G.r_grid,G.b_grid,G.a_grid, G.MP_grid,...
                R, Z_bp, Mu_ap, shock_MP,  pf.y, pf.pi, pf.y_star);  %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
                          
          
    %----------------------------------------------------------------------        
    % Solve for variables inside expectations
    %----------------------------------------------------------------------    
    
%      % consumption at t+1, Eq.(5)
%     Cpf = Ypf - P.phi/2*(piepf-S.pi_star)^2*Ypf;
%     
% %     % Lambda, at t+1,  Eq.(1)
%     Lpf= Cpf^(-1*P.sigma);
%     
%     % A*_t+1/A_t+1 Eq.(6)
%     L_star_pf = Y_s_pf^(-1*P.sigma);

    Z_b_pf = P.rho_b * Z_bp;     %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
    Mu_a_pf = P.rho_a * Mu_ap;    
     
    % consumption Eq.(5)
    Cp = Yp - P.phi/2*(piep-S.pi_star)^2*Yp;
    
    % Lambda, Eq.(1)
    Lp= Cp^(-1*P.sigma);
    
    % A*_t/A_t Eq.(6)
    L_star = Y_star_p^(-1*P.sigma);    
        
    % real rate Eq.(34)
    r_star_p = S.r_star*(1+P.sigma*Mu_a_pf+(1-P.rho_b)*Z_b);    %  r_star_p = S.r_star*(1+P.sigma*Mu_a_pf-Mu_b_pf); %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
    
    % Monetary Policy Rule (4)
    R_star = R^P.rho_r * (S.r_star*S.pi_star...
        * (piep/S.pi_star)^P.psi_pi*(Yp/Y_star_p)^P.psi_y)^(1-P.rho_r)*exp(shock_MP);   % S.r_star instead of r_star_p %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
    
    % Monetary Policy Rule (4)
    Rp = R_star;
    if zlbflag == 1
        Rp = max(1,Rp);
    end 
    
   
    
 %  { 'y','y_star','pi','r','r_lg','r_star','mu_a','mu_b''E_y','E_pi','E_y_star' } ; 
  imp(j,:) = real( [ log(Yp) log(Y_star_p) piep Rp R_star r_star_p shock_b shock_a  shock_MP]);

    if policy_flag==1 % interest rate state last period
        R = R_star;
    else 
        R = Rp;
   end
   
      
%      pie = piep; 
%      Y   = Yp;
  
     Mu_a = Mu_ap;
     Z_b = Z_bp;   %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
     
     shock_a = 0;
     shock_b = 0;
     shock_MP = 0;

end

