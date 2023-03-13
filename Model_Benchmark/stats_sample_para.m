
para_names_p = char('sigma','beta','chi','gamma_a', 'gamma_b', ...
                    'omega', 'epsilon', 'kappa', 'pi_star', ...
                    'psi_pi', 'psi_y', 'rho_a', 'rho_b', 'rho_r', ...
                    'sigma_a', 'sigma_b', 'sigma_r', ...
                    'sigma_Y', 'sigma_Pi', 'sigma_R','phi', 'r^*',...
                    'post', 'lik',...
                    'accept_rate(1)', 'accept_rate(2)','accept_rate(3)',...
                    'accept_rate(4)', 'accept_rate(5)','accept_rate(6)',...
                    'accept_rate(7)', 'accept_rate(8)','accept_rate(9)',...
                    'accept_rate(10)', 'accept_rate(11)','accept_rate(12)'...
                    ,'accept_rate(13)',...
                    'accept_rate(14)', 'accept_rate(15)','accept_rate(16)',...
                    'accept_rate(17)', 'accept_rate(18)','accept_rate(19)',...
                    'accept_rate(20)'   );
                 
%   psi = 1 - ((nstage-1)/nstage)^2;               
%                  
%   post_Resamp(:,1) = post_Resamp(:,1:2);  

%% P.phi   = (P.epsilon-1)*(P.omega+P.sigma)/P.kappa/P.pi_star;
para_phi = ( para_Resamp(:,7)-1 ).*( para_Resamp(:,6)+para_Resamp(:,1) )./...
            para_Resamp(:,8)./(1+para_Resamp(:,9)/100);
        

r_star   = 100*log(exp(para_Resamp(:,1).*para_Resamp(:,4)/100)./ (para_Resamp(:,2).*exp(para_Resamp(:,5))));        

             
parasim = [ para_Resamp para_phi r_star post_Resamp(:,1:2) stock_accept_rate' ] ;

parasim(:,16:18) = 100*parasim(:,16:18);



save( [ './output/para_' num2str(zlbflag) '-' num2str(policy_flag) '.mat' ] , 'parasim');

%% Output posterior file

   % Calculating of Posterior estimates 

a  =0.90;   % (1-2*rate) percentage of Credibile interval of Parameters  
rate = (1-a)/2;


% calculation of posterior estimates of parameters
%

sort_para=zeros( nsim, npara+3+nstage );


  for i=1:1:(npara+3+nstage)
     sort_para(:,i) = sort(parasim(:,i),1);
  end


%  The Fisrt equation
  para_low = sort_para(ceil(( nsim )*rate),:); 
  para_up  = sort_para(floor(( nsim )*(1-rate)),:);

iBm = min([500, nsim/2]); 

fprintf( [ '\n\n  zlbflag= ' num2str(zlbflag) ', policy_flag= ' num2str(policy_flag) ] );
fprintf( [ '\n\n  nstage = ' num2str(nstage) ', particle (para) = ' num2str(nsim) ...
                  ', particle (var) = ' num2str(nparticles) ] ); 
fprintf( [ '\n\n  Measurement error = ' num2str(HH(1,1)) ', ' num2str(HH(2,2)) ...
                  ', ' num2str(HH(3,3)) ] );              
              
fprintf('\n--------------------------------------------');
fprintf('------------------------------------');
fprintf('\n\n                        [ESTIMATION RESULT]');
fprintf('\n--------------------------------------------');
fprintf('------------------------------------');
fprintf('\nParameter    Mean        Stdev     ');
fprintf('95%%Low     95%%Up    Geweke     Inef.');
fprintf('\n----------------------------------');
fprintf('--------------------------------------------\n');

if nsim > 25
for i = 1:(npara)
 if (Pr.pmask(i)==0) 
     fprintf('%s %10.4f  %10.4f %9.3f %9.3f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      fGeweke(parasim(:,i), iBm), ...
      ftsvar(parasim(:,i), iBm)/var(parasim(:,i)) ]  );
 end
end

for i = (npara+1):(npara+3+nstage)
% fprintf( '%s ',  para_names_p(i)  );
fprintf('%s %10.4f  %10.4f %9.3f %9.3f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      fGeweke(parasim(:,i), iBm), ...
      ftsvar(parasim(:,i), iBm)/var(parasim(:,i)) ]  );
  if (i == npara+1)||(i== npara+3)
      fprintf('----------------------------------');
      fprintf('-------------------------------------------- \n');
  end    
end

else
  for i = 1:(npara)
      if (Pr.pmask(i)==0)   
% fprintf( '%s ',  para_names_p(i)  );
fprintf('%s %10.4f  %10.4f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      ]  );
      end
  end
  fprintf('\n----------------------------------');
fprintf('--------------------------------------------\n');

  for i = (npara+1):(npara+4+nstage)
     % fprintf( '%s ',  para_names_p(i)  );
fprintf('%s %10.4f  %10.4f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      ]  );
  end
  
end    

fprintf('-----------------------------------');
fprintf('---------------------------------------------\n');

 
est_date = datestr(date);         
result_name = ['./output/ESTIMATE_para_',num2str(zlbflag),'-', num2str(policy_flag),...
               '-',num2str(nsim), '-',num2str(nparticles), est_date,'.txt'];          

fileID = fopen(result_name,'w');

fprintf(fileID, [ '\n\n  zlbflag = ' num2str(zlbflag) ', policy_flag = ' num2str(policy_flag) ] ); 
fprintf(fileID, [ '\n\n  nstage = ' num2str(nstage) ', particle (para) = ' num2str(nsim) ...
                  ', particle (var) = ' num2str(nparticles) ] ); 

fprintf(fileID, [ '\n\n  Measurement error = ' num2str(HH(1,1)) ', ' num2str(HH(2,2)) ... 
    ', ' num2str(HH(3,3)) ] );                

fprintf(fileID,'\n\n                        [ESTIMATION RESULT]');
fprintf(fileID,'\n----------------------------------');
fprintf(fileID,'------------------------------------');
fprintf(fileID,'\nParameter         Mean        Stdev     ');
fprintf(fileID,'95%%Low     95%%Up    Geweke     Inef.');
fprintf(fileID,'\n----------------------------------');
fprintf(fileID,'--------------------------------------------\n');

if nsim > 25
for i = 1:(npara)
 if (Pr.pmask(i)==0)   
     fprintf(fileID,'%s %10.4f  %10.4f %9.3f %9.3f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      fGeweke(parasim(:,i), iBm), ...
      ftsvar(parasim(:,i), iBm)/var(parasim(:,i)) ]  );
 end
end
for i = (npara+1):(npara+3+nstage)
% fprintf( '%s ',  para_names_p(i)  );
fprintf(fileID,'%s %10.4f  %10.4f %9.3f %9.3f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      fGeweke(parasim(:,i), iBm), ...
      ftsvar(parasim(:,i), iBm)/var(parasim(:,i)) ]  );
end

else
  for i = 1:(npara)
      if (Pr.pmask(i)==0) 
% fprintf( '%s ',  para_names_p(i)  );
fprintf(fileID,'%s %10.4f  %10.4f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      ]  );
      end
  end
  for i = (npara+1):(npara+3+nstage)
     % fprintf( '%s ',  para_names_p(i)  );
fprintf(fileID,'%s %10.4f  %10.4f %9.3f %9.3f \n' ,...
         para_names_p(i,:), ...   
    [ mean(parasim(:,i)) std(parasim(:,i)) para_low(i) para_up(i) ...
      ]  );
  end
  
end    

fprintf(fileID,'-----------------------------------');
fprintf(fileID,'-----------------------------------');
fclose(fileID);
