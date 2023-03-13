function [ log_lik, xt]=...
    fun_ParticleFilter_parallel(yy,ZZ,HH,QQ, O, G, S, P, pf,...
                      particles, zlbflag, policy_flag, ncores, shock_save,rand_save)  
                  
% (yy,ZZ,HH,QQ, O,G, S, P, pf, particles, zlbflag,...
%                           ncores, ACflag, shock_save,rand_save) 
%  input
%
%  output
%   State Space Model
%   (1) Measurement Eq.
%              
%   (2) Transition Eq.
   nn = 8+3;  % number of state variables;  
%    nshock =3;
   
  t           = size(yy,1);                         % number of period
  xt         = zeros(t, nn, particles );     % particles of factors (3*1) at period t+
  prob_x = zeros(1,particles );     % weight of each particles
  xt_i      = zeros(nn,particles);
      
  % Generation of perticles of xt,  at period 1   
       n_particles = particles/ncores;
       wt=zeros(n_particles,ncores); 
       x1=zeros(nn,n_particles,ncores); 
  
       shock = squeeze(shock_save(:,:,1));
       n_par= n_particles;
       
       xt_Resamp = resample_ini(O,S,particles,nn,rand_save);
       
        log_lik=0;
       
  % Generation of perticles of xt,  at period 1 through T
  for i=1:t                   % t is period from 1 to T 
      
        x_temp      =  zeros(nn,n_particles,4); 
        
        shock = squeeze(shock_save(:,:,i));
        
     parfor k =1:ncores          
        
          [wt(:,k), a]= sub_parallel_particle(yy,ZZ,HH,QQ,G, S, P, pf, ...
                                            n_particles, zlbflag,policy_flag, nn,i,  ...
                                            shock(:,1+(k-1)*n_par:k*n_par),...
                                            xt_Resamp(:,1+(k-1)*n_par:k*n_par));       
                                                                           
           x_temp(:,:,k)     =  a;      %  sampling of x_t
     end
       
    for k = 1:ncores;
         prob_x(1,1+(k-1)*n_particles:k*n_particles)  = wt(:,k);  
         xt_i(:,1+(k-1)*n_particles:k*n_particles)        = squeeze(x_temp(:,:,k)); 
    end  
    
%        xt(i-1,:,:) = re_xt_i;       % resampling of x_t-1                         
                                 
       %% Step 4:  calculating likelihood
        lik_t = mean(prob_x);         
        log_lik = log_lik + real(log(lik_t));  % likelihood of full period  
   
    
    %% Step 3 : Resampling xt(j) with Prob (prob_x(j))
     xt_Resamp = resample(xt_i, prob_x, particles);
   
    xt(i,:,:) = xt_Resamp;
   
  end 
    
 %% End of function        
end
 
  
  

  
  
  
  
  