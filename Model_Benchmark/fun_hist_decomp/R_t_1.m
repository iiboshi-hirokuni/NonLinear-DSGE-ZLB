function  decomp_R1 = fun_recursive(decomp_R,decomp_a,decomp_b,decomp_mp,series_R_lag,Tobs)

% series_R_lag;

con_R = decomp_R(:,5) ;

% decomp_b(:,5) ;
% 
% decomp_a(:,5);
% 
% decomp_mp(:,5);

ratio_R = decomp_R(:,5)./series_R_lag;

ratio_a = decomp_a(:,5)./series_R_lag;

ratio_b = decomp_b(:,5)./series_R_lag;

ratio_mp = decomp_mp(:,5)./series_R_lag;


R = zeros(Tobs,4);


for t = Tobs:-1:1


  for i = 1:t-1
    
      R1 = 1;
      
      if i == 1
          R1 = 1;
      else    
         for j = 1:i
             R1 = R1 * ratio_R(Tobs-j);
         end
      end   
   
       R(t,1) = R(t,1) + con_R(t)*R1*ratio_a(t-i);
       R(t,2) = R(t,2) + con_R(t)*R1*ratio_b(t-i);
       R(t,3) = R(t,3) + con_R(t)*R1*ratio_mp(t-i); 
  end  
  
end
  
     R(:,4)= decomp_R(:,5)-R(:,1)- R(:,2) - R(:,3);

%       figure(2)
%       plot(R);
%       legend('a','b','mp','R_{t-1}' )
      
  for t = 1:Tobs
    for i = 1:4  
        if R(t,i) > 0 
           decomp_R1_p(t,i) = R(t,i);
        else
           decomp_R1_n(t,i) = R(t,i);  
        end   
    end
 end      
      
      
      