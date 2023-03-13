function re_xt = resample_ini(O,S,n_par,nn,rand_save)

re_xt = zeros(nn,n_par);

% rng('default');

for i = 1:n_par
       R = O.rbound(1)+(O.rbound(2)-O.rbound(1))*rand_save(1,i);       %Capital state last period   
       Mu_b = O.bbound(1)+(O.bbound(2)-O.bbound(1))*rand_save(2,i);        %
       Mu_a = O.abound(1)+(O.abound(2)-O.abound(1))*rand_save(3,i);    
%        pie = S.pi_star+randn(1)/100;
%        y = S.y_star+randn(1)/100;   
       
%       r = O.rbound(1)+(O.rbound(2)-O.rbound(1))*rand_save(1,i);       %Capital state last period   
%       g = O.gbound(1)+(O.gbound(2)-O.gbound(1))*rand_save(2,i);        %
%       z = O.zbound(1)+(O.zbound(2)-O.zbound(1))*rand_save(3,i);    
% %       pie = P.pi+randn(1)/100;
%        y = S.y+rand_save(4,i)/25;
      
      % Policy Function       
     state= zeros(nn,1); 
     state(4)= R; 
%      state(3)= pie; 
     state(5)= R;
     state(7)= Mu_a;
     state(8)= Mu_b; 
     
     re_xt(:,i) = state;


end
