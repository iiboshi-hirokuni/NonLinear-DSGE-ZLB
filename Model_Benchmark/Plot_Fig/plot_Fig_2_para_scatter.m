%
%  Figure 2: Particles for Parameter sigma and their Likelihood
%
%
%

addpath('./data_sample')

set(0,'defaultAxesFontSize',12);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',12);
set(0,'defaultTextFontName','century');


para_names_p = char('sigma','beta','chi','gamma_a', 'gamma_b', ...
                    'omega', 'epsilon', 'kappa', 'pi_star', ...
                    'psi_pi', 'psi_y', 'rho_a', 'rho_b', 'rho_r', ...
                    'sigma_a', 'sigma_b', 'sigma_r', ...
                    'sigma_Y', 'sigma_Pi', 'sigma_R','phi', ...
                    'post', 'lik',...
                    'accept_rate(1)', 'accept_rate(2)','accept_rate(3)',...
                    'accept_rate(4)', 'accept_rate(5)','accept_rate(6)',...
                    'accept_rate(7)', 'accept_rate(8)','accept_rate(9)',...
                    'accept_rate(10)', 'accept_rate(11)','accept_rate(12)'...
                     );
             

load('para_1-1.mat')

x1 = parasim(:,1);  %sigma
x4 = parasim(:,4);  %gamma_a
% x5 = parasim(:,5);  %gamma_b
x6 = parasim(:,6);  %omega
x8 = parasim(:,8);  %kappa
x9 = parasim(:,9);  %pi\star
x10 = parasim(:,10);  %psi_pi
x11 = parasim(:,11);  %psi_y
x12 = parasim(:,12);  %rho_a

lik = parasim(:,23); 
[y, i] = max(parasim(:,23))
y1 = parasim(i,:);

m1 = median(parasim(:,1));
m4 = median(parasim(:,4));
m6 = median(parasim(:,6));
m8 = median(parasim(:,8));
mlik = median(parasim(:,23));



figure('Name','Scatter of Parameters','File','Fig2_scatter');

subplot(1,1,1) 
    
%    scatter(x10,lik0,2,'c','x');
%    scatter(x1,lik,5,'b','o','filled');
    scatter(x1,lik,5,'b','x');
    hold on
%    scatter(x12,lik2,2,'r','+'); 
   scatter(y1(1),y1(23),100,'+','k','LineWidth',2.5)
   scatter(m1,mlik, 100,'o','filled','r')
%    scatter(y2(1),y2(23),50,'o','filled','k')
   hold off
%    title('\sigma vs. likelihood','FontSize',12);
   xlabel('\sigma','FontSize',16);
   ylabel('Log Likelihood','FontSize',14); 
   xlim([1.1 1.7])
      lgnd = legend({ 'Particles', 'Particle with maximum log lik',...
             'Particle for median'},'FontSize',12);
      set(lgnd, 'Box', 'off')
      set(lgnd,'Location','southwest')
      
savefig('output\Fig2_scatter.fig');    
      

   