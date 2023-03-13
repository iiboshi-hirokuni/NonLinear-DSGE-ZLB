function [ pf, converged ] = solve_model_CL(P,S,G,zlbflag,policy_flag,tstart,disp_it)


 pf = guess_TL(P,S,G);    %% calculate linear version of policy functions

% %--------------------------------------------------------------------------
% % Update Policy Functions
% %--------------------------------------------------------------------------

it = 1;                                 % Iteration Counter
converged = 0;                          % Convergence Flag
dist_max = 0;                           % Max distance vector
 perbind = 0;                           % percentage of ZLB 
 no_solution = 0;                       % when number of iteration is over 50, then there is no solution. 
 
% 
% % Preallocate arrays to store policy function updates
pf_y_up = zeros(G.griddim);
pf_y_star_up = zeros(G.griddim);
pf_pi_up = zeros(G.griddim);
pf_r_up = zeros(G.griddim);
 
% % Time iteration/linear interpolation algorithm
% 
while converged == 0 && perbind < 50 && no_solution == 0
    istart = tic;                       % Iteration timer start
   
    %     parfor i = 1:G.nodes  %% for pararell computing
    for i = 1:G.nodes
%     
%         % Find optimal policy functions on each node.
%         % csolve finds the zeros of 'eqm'
%         % Start csolve with the current policy function  
          start = [pf.y(i) pf.pi(i) pf.y_star(i) ]';
          state = [G.r_gr(i) G.b_gr(i) G.a_gr(i) G.MP_gr(i)];
          
         %  csolve(FUN,x,gradfun,crit,itmax,varargin)
%          argzero = csolve('eqm_TL',start,[],1e-3,5,state,pf,P,S,G,zlbflag,policy_flag);
         argzero = csolve('eqm_CL_TL',start,[],1e-3,5,state,pf,P,S,G,zlbflag,policy_flag);
%  
%       % Store updated policy functions       
          pf_y_up(i) = real(argzero(1));
          pf_pi_up(i) =real( argzero(2));
          pf_y_star_up(i) = real(argzero(3));
   
    %% calculate policy function of interest rate
     Y = argzero(1);
     pie = argzero(2);
     Y_star = argzero(3);
     epMP = G.MP_gr(i);
     R    = G.r_gr(i);
%      Mu_b_pf = P.rho_b * state(2) ;
       Z_b = state(2);
%       Z_b_pf = P.rho_b * state(2) ;
     Mu_a_pf = P.rho_a * state(3) ;
    
    r_star_p = S.r_star*(1+P.sigma*Mu_a_pf+(1-P.rho_b )*Z_b);
    y = log(Y/S.y_star);
    y_star = log(Y_star/S.y_star);
    R_star = S.r_star*S.pi_star*(1+P.rho_r*(R/(S.r_star*S.pi_star)-1) ...
             + (1-P.rho_r)*(P.psi_pi*(pie-S.pi_star)/S.pi_star +P.psi_y*(y-y_star))+epMP);
     Rp = R_star;
     if zlbflag == 1
         Rp = max(1,Rp);
     end        
     pf_r_up(i) = Rp;
     
    end
 
     % Find where ZLB binds   
     %   Interest rate rule  
      locs = find(pf_r_up == 1);
     %   Percent nodes binding
      perbind = 100*numel(locs)/G.nodes;
 
    % Policy function distances
    dist_y = abs(pf_y_up - pf.y);
    dist_pi = abs(pf_pi_up - pf.pi);
    dist_y_star = abs(pf_y_star_up - pf.y_star);
    

%     % Maximum distance   
    dist_max(it) = max([dist_y(:)' dist_pi(:)' dist_y_star(:)']);
% 
    % Update policy functions
    pf.y = pf_y_up;
    pf.pi = pf_pi_up;
    pf.y_star = pf_y_star_up;
    pf.r = pf_r_up;
% 
%     % Check convergence criterion
    if it > 5 && all(dist_max(end-1:end) < P.tol)
        converged = 1;
    elseif it > P.ite
        no_solution = 1;
        converged = 2;
        disp('no_solution since the number of iteration is over criterion ');
    elseif perbind >50
        no_solution = 1;
        converged = 2;
        disp('no_solution since area of the ZLB is over 50% ');
    end
% 
    % Iteration Information and save policy functions to disk
     if converged == 1     
        it = itinfo(istart,tstart,1,it,dist_max(it));
        disp(['nodes binding ZLB are : ' num2str(perbind) '%']);
        %         disp('');
        %         disp([ 'model is solved']);
     elseif mod(it,10) == 1 && P.disp_it == 1 
        it = itinfo(istart,tstart,1,it,dist_max(it));        
        disp(['nodes binding ZLB are : ' num2str(perbind) '%']);
    else
        it = it+1;
    end  
   
end
  
% 