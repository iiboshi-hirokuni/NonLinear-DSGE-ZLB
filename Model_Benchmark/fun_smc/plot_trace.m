% 
% 
%  Plot traces of sampling parameters
% 
% 
% 

function plot_trace(parasim, npara, P,ii)


para_names_p = {'sigma','beta','chi','gamma_a', 'gamma_b', ...
                    'omega', 'epsilon', 'kappa', '\pi_{star}', ...
                    'psi_{\pi}', 'psi_y', 'rho_a', 'rho_b', 'rho_r', ...
                    'sigma_a', 'sigma_b', 'sigma_r', ...
                    'sigma_Y', 'sigma_{pi}', 'sigma_R', ...
                    'post', 'lik',...
                    'accept_rate(1)', 'accept_rate(2)','accept_rate(3)',...
                    'accept_rate(4)', 'accept_rate(5)','accept_rate(6)',...
                    'accept_rate(7)', 'accept_rate(8)','accept_rate(9)',...
                    'accept_rate(10)', 'accept_rate(11)','accept_rate(12)'...
                     };
 
%% draw graph both of Prior and posterior densities
    k=0;
    j = 1;
    
 for i =1:1:npara
     
  if P.pmask(i)==0 
    
  if mod(i,4*5)==1
     k = k +1 ;
     figure('Name',[num2str(k*100+ii)],'Position',[750,100,750,500])  
%       figure('Position',[250,250,1000,250], 'Name', 'Fitting of Macro Factors','NumberTitle','off', 'FileName','Fig3001')  
  end
 
 subplot(4,4,j) % 4(çs)Å~4(óÒ)
 
%      [density0,x0]  = ksdensity(parasim_pri(1:end,i));
%      [density,x]  = ksdensity(parasim_post(1:end,i));
     
           plot(parasim(:,i)','LineStyle','-','Color','b',...
          'LineWidth',1.0);
      
%      if both_plot == 1 
%       hold on
%           plot(x0,density0,'LineStyle','--','Color','r',...
%          'LineWidth',2.5);
%        hold off 
%      end
      
     title( char(para_names_p(i)),'FontSize',12) ; 
     
      j= j+1; 
  end
   
  
 end
 
 pause(0.05)





