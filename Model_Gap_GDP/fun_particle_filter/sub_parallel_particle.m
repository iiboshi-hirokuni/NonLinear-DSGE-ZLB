function [wt,xt_t ]= ...
               sub_parallel_particle(yy,ZZ,HH,QQ, G, S, P, pf, ...
                       particles, zlbflag,policy_flag,  nn,t, shock_save,xt_Resamp)
 
%        xt_Resamp = zeros(nn,particles);               
%                          
%         for j=1:particles    
%             xt_Resamp(:,j) = mnrnd(1, prob_x/sum(prob_x)) * xt';  
%         end             
     

  wt = zeros(1,particles); 
  xt_t = zeros(nn,particles); 
  
   num = particles; 
   AA=chol(QQ);
    
 for j=1:num 

  %% Step 1:  generating particle x_t        
         %   Transition Eq.    
%             xt(:,j)
           state = squeeze(xt_Resamp(:,j)) ;                               % factors at period t-1         
           shock = AA*shock_save(:,j);
           xt_t(:,j) = gen_particle(G, S, P, pf, state, shock, zlbflag,policy_flag);
                 
  %% Step 2 : Calculating Probabilities corresponding to particles 
          %   Measurement Eq.        
         xx = [ ZZ(1)*( xt_t(1,j)-  xt_t(2,j) ); ...
                ZZ(2)* ( log(xt_t(3,j))  ) ; ...                          % inflation
                ZZ(3)* ( log( xt_t(4,j) )  )    ];                       % interest rate

          wt(1,j) = mvnpdf(yy(t,:)',  xx, HH );  
          
          if  isnan(wt(1,j))   %% if wt(1,j) ==NaN
                wt(1,j) = 0;
          end               
           
 end  
 
 %% End of function        
end
 
 