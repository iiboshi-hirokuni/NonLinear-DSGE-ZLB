function  [ mean_decomp, all_factor  ] = ...
         fun_decomp(Tobs,d,G, S, P, pf, zlbflag, policy_flag, shock_idex, p,q,a )  

sum_decomp  = zeros(Tobs,8);

mat = zeros(2^4,4);
 n = 1;
 for i = 0:1
   for j = 0:1
      for k = 0:1
        for l = 0:1   
           mat(n,1) = i;
           mat(n,2) = j;
           mat(n,3) = k;
           mat(n,4) = l;
           n = n+1;
        end   
      end
   end
 end  

 for t = 1:Tobs

for i = 1:2^4     
        state = zeros(8,1);
        
        %% setting state variables
        if  mat(i,1)==1
           state(7,1) = d.series_mu_a(t)/100*(p/q)^a(1);  
       else    
           state(7,1) = d.series_mu_a(t)/100*((p-1)/q)^a(1);
        end
       
       if  mat(i,2)==1
           state(8,1) = d.series_z_b(t)/100*(p/q)^a(2);
       else    
           state(8,1) =  d.series_z_b(t)/100*((p-1)/q)^a(2);
       end
       
        
       if  mat(i,3)==1
            state(4,1) = d.series_R_lag(t)/100*(p/q)^a(3); 
            state(5,1) = d.series_R_lag(t)/100*(p/q)^a(3); 
       else    
             state(4,1) = d.series_R_lag(t)/100*((p-1)/q)^a(3); 
             state(5,1) = d.series_R_lag(t)/100*((p-1)/q)^a(3); 
       end
       
      
        if  mat(i,4)==1
           shock = zeros(3,1);    
           shock(3,1) = d.series_shock_mp(t)/100*(p/q)^a(4);
        else    
           shock = zeros(3,1);    
           shock(3,1) = d.series_shock_mp(t)/100*((p-1)/q)^a(4);
        end
       
        
      %% with shock 
       if  shock_idex == 1
           state(7,1) = d.series_mu_a(t)/100*(p/q)^a(1);    
       elseif shock_idex == 2
          state(8,1) = d.series_z_b(t)/100*(p/q)^a(2);
       elseif shock_idex == 3
          state(4,1) = d.series_R_lag(t)/100*(p/q)^a(3);
          state(5,1) = d.series_R_lag(t)/100*(p/q)^a(3);
       else    
          shock = zeros(3,1);    
          shock(3,1) = d.series_shock_mp(t)/100*(p/q)^a(4);  
       end     
           
           decomp_1 = gen_hist_decomp(G, S, P, pf, state, shock, zlbflag, policy_flag); 
     
      %% without shock   
      if  shock_idex == 1
           state(7,1) = d.series_mu_a(t)/100*((p-1)/q)^a(1);    
       elseif shock_idex == 2
          state(8,1) = d.series_z_b(t)/100*((p-1)/q)^a(2);
       elseif shock_idex == 3
          state(4,1) = d.series_R_lag(t)/100*((p-1)/q)^a(3);
          state(5,1) = d.series_R_lag(t)/100*((p-1)/q)^a(3);
       else    
          shock = zeros(3,1);    
          shock(3,1) = d.series_shock_mp(t)/100*((p-1)/q)^a(4);
       end     
       decomp_0 = gen_hist_decomp(G, S, P, pf, state, shock, zlbflag, policy_flag); 
     
      %%  sum of decompositions 
       sum_decomp(t,:) = sum_decomp(t,:) + decomp_1 - decomp_0;
  
   end  

 end
    
     %%  mean of decompositions 
       mean_decomp = 100*sum_decomp/2^4;
       
       
 for t = 1:Tobs   
 
    state = zeros(8,1);
    state(4,1) = d.series_R_lag(t)/100;
    state(5,1) = d.series_R_lag(t)/100;
    state(7,1) = d.series_mu_a(t)/100;
    state(8,1) = d.series_z_b(t)/100;
    
    shock = zeros(3,1);    
    shock(3,1) = d.series_shock_mp(t)/100;
    
     all_factor(t,:) =gen_hist_decomp(G, S, P, pf, state, shock, zlbflag, policy_flag); 
     
 end    
      