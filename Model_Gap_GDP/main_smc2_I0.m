% 
%  Sequential Monte Carlo Methods with Particle Filter
%  Canonical New Keynesian Model (Rotemberg Pricing)    
%   -Imposes the zero lower bound on the interest rate
%


%% setting 
npara      =    20;          % # of parameters
cc1        =   0.25 ;       % adjustment coefficient of MH
set_itr_disp = 1;

%%  sampling parameters from prior 
    Pr = fun_prior_setting;
    cc = cc1*diag(Pr.pstdd );
    cc(16,16)=cc1*0.01; % stdev of shock a
    cc(17,17)=cc1*0.01; % stdev of shock b
    cc(18,18)=cc1*0.01; % stdev of shock r
    
%%  Measurement errors of output, inflation, interest rate
if m_err_flag == 1  
    m_error = [0.005 0.005 0.0025]; % measurement errors are of the sizes of 1%, 1%, 0.5% of SD, respectively
elseif m_err_flag == 2  
    m_error = [0.001 0.001 0.0005];
end    
    
HH = 100^2*diag([ (m_error(1)*Pr.pmean(18))^2 (m_error(2)*Pr.pmean(19))^2 (m_error(3)*Pr.pmean(20))^2 ]);
  

%% Load Data
  load_data

%% Load parameters, steady state and grids
P = parameters;
S = steadystate(P);

% Specify grid options
O.rbound = [0.90  1.1];
O.MPbound = [-0.03  0.03];
O.bbound = [-0.1  0.1 ];   % mu^b_t = log(Z_t+1/Z_t) - gamma_b
O.abound = [-0.1  0.1 ];   % mu^a_t = log(A_t+1/A_t) - gamma_a

O.r_pts = 9; 
O.MP_pts = 5; 
O.b_pts = 9; 
O.a_pts = 9; 

O.e_pts = 1; % number of node of shock b for stochastic foresight
             %  1-> perfect foresight
O.u_pts = 1; % shock of a
O.v_pts = 1; % shock of  MP

% Load discretized state space
G = grids_even(O,P);
     
%% initial setting of particles filter
  nshocks=3;          % # of shocks
  npar = 1;           % for parallel computing
  disp( ['# of particles is ' num2str(nparticles) ]);
  
%   file_name = ['./data/particle_' num2str(nparticles) ];
%   load(file_name,'stock_shock','stock_state');
  
  stock_shock = zeros(nshocks,ceil(nparticles),Tobs); % matrix of 3 dimensions ( nshocks X nparticles X Tobs )
  rng(100)  % generator of random numbers is fixed
  
   for i =1:Tobs
        stock_shock(:,:,i) = mvnrnd(zeros(nshocks,1), eye(3), ceil(nparticles) )'; 
   end
   
   stock_state= [ rand(nshocks,nparticles); randn(1,nparticles) ]; % initial values of state variables (t=0)
      
   
  %% Covariance Matrices  
    parasim =  sample_pri(Pr,nsim,npara); 
    lik_stock = zeros(nsim,3);
    
    stock_accept_rate = zeros(nstage,nsim);
    psi = 0;
    weight = ones(1,nsim);
    priornew = zeros(nsim,1);
 
 %% start SMC^2   
for i = 1:nstage 
    disp( ' ' );   
  disp([ ' ' num2str(i), ' th-stage' ]) ;
 
%% Step 1. Correction

for  k = 1:nsim
   P_nsim(k) = parameters_update(parasim(k,:),P);
end

   %% update of psi of likelihood function 
        if i > 1
          psi = (i/nstage)^2-((i-1)/nstage)^2;
        else
          psi = (i/nstage)^2;
        end    
         disp( [ 'psi = ' num2str(psi) ]  )

  parfor k = 1:nsim  
% for k = 1:nsim    
  if mod(i,set_itr_disp)==0
         disp( ' ' );   
         disp([ ' Step 1: Correction ' num2str(k), ' th-particle of ' num2str(i), ' th-stage ']) ;
  end       
  
       P=P_nsim(k) ;    
  
    % Projection method for solution of DSGE model  
     [check_para] = fun_check_para(parasim(k,:),Pr);
     
     if check_para == 1
          likenew = -1E6;
           if mod(i,set_itr_disp)==0
                   disp([ 'out of bound of parameters' ]);
           end        
     elseif i == 1
             [ likenew ]  =   step_solve_particle(P,O,tstart,yy,ZZ,HH, ...
                                 nparticles, zlbflag, policy_flag, npar, ...
                                 stock_shock,stock_state) ;  
     else
             dummy =  lik_stock(k,:);
             likenew = dummy(1,2);
     end        
                        
%       % calculate prior density  
        if i ==1
              priornew(k) =priodens(parasim(k,:), Pr.pmean, Pr.pstdd, Pr.pshape);   
              postnew = psi*likenew + priornew(k) ;
        else
             dummy =  lik_stock(k,:);
             post_old = dummy(1,1);
             postnew  = post_old + psi*likenew;
        end     
            
        % incremental weights
        incwt    =  psi*likenew;
     
       % save lik and posterior
        if mod(i,set_itr_disp)==0
               disp( [ 'post = ' num2str(postnew) ',   like = '  num2str(likenew) ...
                 ]);
        end     
%                '  para ='  num2str(parasim(k,:))    ] );
       lik_stock(k,:)=[postnew likenew incwt];
     
  end    
% end 
 
 %% Step 2. Selection
  disp( ' ' );   
  disp([ ' Step 2: Selection of ' num2str(i), ' th-stage ']);
      [ para_Resamp , lik_Resamp, weight, ESS ] = resample_para(parasim', lik_stock', npara,nsim, weight);
        para_Resamp = para_Resamp';
        lik_Resamp = lik_Resamp';    
        
        post_Resamp = zeros(nsim,3);
     
        
 
 %% Step 3. Mutation 
      para_new =zeros(nsim,npara);
      
       % variances of parameters of Algorithm 10 of Herbst and Schorfheide (2016, Ch5, p113)    
      if i >1
          V=cov(para_Resamp);
          V = V + diag( 1E-7*ones(npara,1));
          cc = cc1*chol(V,'upper');
      end
          
 for  j = 1:nsim
    
       % RW sampling of candidates of parameters    
       para_new(j,:) = para_Resamp(j,:) + randn(1,npara)*cc  ;
          para_new(j,:) = para_new(j,:).*(1-Pr.pmask)'+para_Resamp(j,:).*Pr.pmask';

       % Projection method for solution of DSGE model
       P_nsim(j) = parameters_update(para_new(j,:),P);
end  
 
parfor j = 1:nsim 
%  for j = 1:nsim  
      if mod(j,set_itr_disp)==0
             disp( ' ' );  
            disp([ ' Step 3: Mutation ' num2str(j), ' th-particle of ' num2str(i), '/' ...
                     num2str(nstage), ' th-stage ' ]) ;  
      end           
    
         P=P_nsim(j) ;
         
        [check_para] = fun_check_para(parasim(j,:),Pr);
     
        if check_para == 1
               likenew = -1E6;
                 if mod(j,set_itr_disp)==0
                       disp([ 'out of bound of parameters' ]);
                 end      
        else
            [ likenew] = step_solve_particle(P,O,tstart,yy,ZZ,HH, ...
                            nparticles, zlbflag, policy_flag, npar,...
                            stock_shock,stock_state);  
        end                
               
        %  calculate prior density   
        priornew(j) =priodens(para_new(j,:), Pr.pmean, Pr.pstdd, Pr.pshape);   
        postnew =  (i/nstage)^2*likenew + priornew(j); 
        
        % MH_step
        r = min(1,exp(postnew-lik_Resamp(j,1) ));   
        
        if (rand < r)   
           weight(j)=1; 
           stock_accept_rate(i,j) = 1;
           para_Resamp(j,:) = para_new(j,:) ;           
           post_Resamp(j,:) = [postnew  likenew   0  ];
            if mod(j,set_itr_disp)==0
                 disp( [ 'post = ' num2str(postnew) ',   like = '  num2str(likenew) ...
                         ]);
            end         
%                  ', para = '  num2str(para_new(j,:)) ] );
        else
            post_Resamp(j,:) = [ lik_Resamp(j,1) lik_Resamp(j,2) 0 ];
             if mod(j,set_itr_disp)==0
                    disp( 'no change of parameters')
             end       
        end 
        
end   
         R =   mean(stock_accept_rate(i,:))*100;
         % scaling factor of Algorithm 10 of Herbst and Schorfheide (2016, Ch5, p113)        
         x=6; 
         cc1 = cc1*(0.95+0.1*exp(x*(R-0.25))/(1+exp(x*(R-0.25))) );

         disp( ' ' ); 
         disp( [  num2str(i), ' th-stage_s  Accept Rate = ' num2str(R) ' %' ]);   
         disp( [  num2str(i), ' th-stage_s  scaling factor = ' num2str(cc1)  ]); 
         
         parasim = para_Resamp;
         lik_stock = post_Resamp;
         
        %% save file 
          file_name = ['./output/save_temp_para_' num2str(zlbflag)  num2str(policy_flag) ...
                       '_' num2str(nparticles) '_' num2str(nsim) '_'  num2str(nstage) ];                      
          save(file_name,'para_Resamp','post_Resamp', 'stock_accept_rate');
          
        %% plot graph
          plot_dist(parasim, npara,Pr,1)
          plot_trace(parasim, npara,Pr,1)
          
        
          
%%  end of Iterations           
end 



% save parameters

  file_name = ['./output/save_para_I0_'  num2str(zlbflag)  num2str(policy_flag) '_' ...
                                      num2str(nsim) '_' num2str(nparticles) ];
                      
     save(file_name,'para_Resamp','post_Resamp', 'stock_accept_rate',...
           'O','Pr','HH');  %,'stock_shock','stock_state');
     
%  load(file_name,'para_Resamp','post_Resamp');
     
  stats_sample_para

  plot_dist(parasim, npara,Pr,i)
  
  plot_trace(parasim, npara,Pr,i)

                  
%% time of computing
 dec = 10^1;
T = toc(tstart);
hh = floor(T/3600);
mm = floor(T/60)-60*hh;
ss = round((T-60*mm-3600*hh)*dec)/dec;

hh = num2str(hh);
mm = num2str(mm);
ss = num2str(ss);

display(['Total Computing time: ' hh 'h' mm 'm' ss 's']);      
         
                  
 