% 
% 
%  Plot distributions of parameters
% 
% 
% 
% 

function plot_dist(parasim, npara, P,ii)

% parasim = [ para_Resamp  post_Resamp stock_accept_rate' ] ;

both_plot = 1; % on =1, off =0 

para_names_p = {'sigma','beta','chi','gamma_a', 'gamma_b', ...
                    'omega', 'epsilon', 'kappa', '\pi_{star}', ...
                    'psi_{\pi}', 'psi_y', 'rho_a', 'rho_b', 'rho_r', ...
                    'sigma_a', 'sigma_b', 'sigma_r', ...
                    'sigma_Y', 'sigma_{\pi}', 'sigma_R', ...
                    'post', 'lik',...
                    'accept_rate(1)', 'accept_rate(2)','accept_rate(3)',...
                    'accept_rate(4)', 'accept_rate(5)','accept_rate(6)',...
                    'accept_rate(7)', 'accept_rate(8)','accept_rate(9)',...
                    'accept_rate(10)', 'accept_rate(11)','accept_rate(12)'...
                     };
 

%% Output posterior file
%

%% draw graph only Prior density
% if output_type == 1
 k = 0;
% 
%  for i =1:1:npara
%     
%  if mod(i,3*4)==1
%      k = k +1 ;
%      figure(1000+k*10)     
%  end
%  
%  subplot(3,4,i-(k-1)*3*4) % 3(çs)Å~4(óÒ)
%  
%      [density,x1]  = ksdensity(parasim_pri(1:end,i));
%      
%     plot(x1,density,'LineStyle','--','Color','r',...
%         'LineWidth',2.5);
%      title( para_names_p(i),'FontSize',12 );
%  
%  
%  end 
% end

%% draw graph both of Prior and posterior densities


 j=1 ;
    
 for i =1:1:npara
     
  if P.pmask(i)==0   
    
  if mod(i,4*5)==1
     k = k +1 ;
     figure('Name',[num2str(k*1000+ii)],'Position',[25,50,800,600] )    
%       figure('Position',[250,250,1000,250], 'Name', 'Fitting of Macro Factors','NumberTitle','off', 'FileName','Fig3001')  
  end
 
 subplot(4,4,j) % 4(çs)Å~4(óÒ)
 
%      [density0,x0]  = ksdensity(parasim_pri(1:end,i));
     [density,x]  = ksdensity(parasim(1:end,i));
     
           plot(x,density,'LineStyle','-','Color','b',...
          'LineWidth',2.5);
           y=linspace(min(density),max(density),size(x,2));

    
%      if both_plot == 1 
%       hold on
%           plot(x0,density0,'LineStyle','--','Color','r',...
%          'LineWidth',2.5);
%        hold off 
%      end
      
     title( char(para_names_p(i)),'FontSize',12) ; 
%     if i == 1 && mode_flag==1
%      legend('density','max', 'second', 'third');
%     end  
     j= j+1; 
  end
   
 end
 
 pause(0.05)




