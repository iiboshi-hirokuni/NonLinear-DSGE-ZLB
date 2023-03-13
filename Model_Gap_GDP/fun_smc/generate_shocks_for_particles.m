
 nshocks = 3;
 nparticles = 2000;
 Tobs = 134;
 

stock_shock = zeros(nshocks,ceil(nparticles),Tobs); % matrix of 3 dimensions ( nshocks X nparticles X Tobs )

   for i =1:Tobs
        stock_shock(:,:,i) = mvnrnd(zeros(nshocks,1), eye(3), ceil(nparticles) )'; 
   end
   
   stock_state= [ rand(nshocks,nparticles); randn(1,nparticles) ]; % initial values of state variables (t=0)
   
   file_name = ['./data/particle_' num2str(nparticles) ];
   save(file_name,'stock_shock','stock_state');
      