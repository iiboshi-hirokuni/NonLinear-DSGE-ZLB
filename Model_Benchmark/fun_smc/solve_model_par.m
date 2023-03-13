function [ pf, converged ] = solve_model_par(P,S,G,zlbflag,policy_flag,tstart)


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
    parfor i = 1:G.nodes
%     
%         % Find optimal policy functions on each node.
%         % csolve finds the zeros of 'eqm'
%         % Start csolve with the current policy function  
          start = [pf.y(i) pf.pi(i) pf.y_star(i) ]';
          state = [G.r_gr(i) G.b_gr(i) G.a_gr(i) G.MP_gr(i)];
          
         %  csolve(FUN,x,gradfun,crit,itmax,varargin)
         argzero = csolve('eqm_TL',start,[],1e-3,5,state,pf,P,S,G,zlbflag,policy_flag);
%  
%       % Store updated policy functions       
          pf_y_up(i) = real(argzero(1));
          pf_pi_up(i) =real( argzero(2));
          pf_y_star_up(i) = real(argzero(3));
   
    %% calculate policy function of interest rate
     Y = argzero(1);
     pie = argzero(2);
     Y_star = argzero(3);
%      Z_b_pf = P.rho_b * state(2) ;   %%  2017/Jan/27 C³
%      Bpf = exp(Z_b_pf - state(2) );   %%  2017/Jan/27 C³
     Mu_a_pf = P.rho_a * state(3) ;
     r_star_p = S.r_star*(1+P.sigma*Mu_a_pf- (1-P.rho_b)*state(2) ); %%  2017/Jan/27 C³
     R_star = state(1)^P.rho_r * (r_star_p*S.pi_star...
            * (pie/S.pi_star)^P.psi_pi*(Y/Y_star)^P.psi_y)^(1-P.rho_r)*exp(state(4));
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
    if it > 5 && all(dist_max(end-2:end) < P.tol)
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
     elseif mod(it,4) == 1 && P.disp_it == 1 
        it = itinfo(istart,tstart,1,it,dist_max(it));        
        disp(['nodes binding ZLB are : ' num2str(perbind) '%']);
    else
        it = it+1;
    end  
   
end
  
% 