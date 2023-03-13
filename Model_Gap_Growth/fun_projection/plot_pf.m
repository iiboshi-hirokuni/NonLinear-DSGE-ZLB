%  shock_type = 'a_shock';
%  shock_type = 'MP_shock';

if strcmp(type,'linear')
    a = 20; b=25;
else
    a =0; b=5;
    if zlbflag == 1
        a = 10; b=15;
    end
end


switch shock_type
    case 'MP_shock'
        
b = round(size(G.b_grid,2)/2);
% MP = round(size(G.MP_grid,1)/2);
r_lag = round(size(G.r_grid,2)/2)-2;


Py = squeeze(pf.y(r_lag,b,:,:));
Ppi = squeeze(pf.pi(r_lag,b,:,:));
Py_star = squeeze(pf.y_star(r_lag,b,:,:));
Pr = squeeze(pf.r(r_lag,b,:,:));


y = G.a_grid;
x = G.MP_grid;

[X,Y]=meshgrid(x,y);

 fig_mp_y = figure('Position',[20,20,600,600],'Name',...
             'policy function of Output by Monetary Policy','Color','w');
 surf(X,Y,Py)
   xlabel('MP shock')
   ylabel('log(A_{t+1}/A_t) - gamma_a')  
   zlabel('Output')
 colormap jet %copper; %hsv; %parula
 title('policy function of Output','FontSize',14)
%  view(3)
  view([55,25])

 est_date = datestr(date);   
%  name = ['../output/','fig_pf_y_mp_',est_date];
%          saveas( fig_mp_y,name,'fig')
        
 %%       
 
 fig_mp_pi = figure('Position',[20,20,600,600],'Name',...
              'policy function of Inflation by Monetary Policy','Color','w');
 surf(X,Y,Ppi)
  xlabel('MP shock')
   ylabel('log(A_{t+1}/A_t) - gamma_a') 
   zlabel('Inflation')
 colormap jet
 title('policy function of Inflation','FontSize',14)
%  view(3)
  view([55,25])

 est_date = datestr(date);   
%  name = ['../output/','fig_pf_pi_mp_',est_date];
%   saveas(fig_mp_pi,name,'fig')
        
 %%       
 
 fig_mp_y_star = figure('Position',[20,20,600,600],'Name',...
                   'policy function of Y^*_t/A_t by Monetary Policy','Color','w');
 surf(X,Y,Py_star)
 xlabel('MP shock')
   ylabel('log(A_{t+1}/A_t) - gamma_a')   
   zlabel('Y^*_t/A_t')
  colormap jet %winter
%  title('policy function of Interest Rate','FontSize',14)
 title('policy function of Y^*_t/A_t ','FontSize',14)
%  view(3)
 view([55,45])
 
 est_date = datestr(date);   
%  name = ['../output/','fig_pf_r_mp_',est_date];
%          saveas(fig_mp_y_star,name,'fig')
         
         
%%       
 
 fig_mp_r = figure('Position',[20,20,600,600],'Name',...
                   'policy function of interest rate by Monetary Policy','Color','w');
 surf(X,Y,Pr)
 xlabel('MP shock')
   ylabel('log(A_{t+1}/A_t) - gamma_a')   
   zlabel('interest rate')
  colormap jet %winter
%  title('policy function of Interest Rate','FontSize',14)
 title('policy function of interest rate ','FontSize',14)
%  view(3)
 view([25,35])
 
 est_date = datestr(date);   
%  name = ['../output/','fig_pf_r_mp_',est_date];
%          saveas(fig_mp_r,name,'fig')


        
 %%
 
    case 'a_shock'       
        
r = floor(mean(size(G.r_grid,2)));
r = 2; % round(min(size(G.r_grid,1)));
MP = floor(size(G.MP_grid,2)/2);
MP = 1;

Py = squeeze(pf.y(r,:,:,MP));
Ppi = squeeze(pf.pi(r,:,:,MP));
Py_star = squeeze(pf.y_star(r,:,:,MP));
Pr = squeeze(pf.r(r,:,:,MP));

 x = G.a_grid;
y = G.b_grid;

[X,Y]=meshgrid(x,y);

 fig_fp_y = figure('Position',[20,20,600,600],'Name','policy function of Output ','Color','w');
 surf(X,Y,Py)
   xlabel('log(A_{t+1}/A_t)-gamma_a')
   ylabel('log(Z^b_{t}')   
   zlabel('Output')
 colormap jet %copper; %hsv; %parula
 title('policy function of Output','FontSize',14)
 %  view(3)
 view([-65,45])
 
 est_date = datestr(date);   
%  name = ['../output/','fig_pf_y_fp_',est_date];
%          saveas(fig_fp_y,name,'fig')
        
 %%     
 
fig_fp_pi = figure('Position',[20,20,600,600],'Name','policy function of Inflation ','Color','w');
 surf(X,Y,Ppi)
 xlabel('log(A_{t+1}/A_t)-gamma_a')
   ylabel('log(Z^b_{t}')    
   zlabel('Inflation')
 colormap jet
 title('policy function of Inflation','FontSize',14)
%  title('policy function of Y^*_t/A_t ','FontSize',14)
 %  view(3)
 view([-65,45])
  
 est_date = datestr(date);   
%  name = ['../output/','fig_pf_pi_fp_',est_date];
%          saveas(fig_fp_pi,name,'fig')
        
 %%     
 
 fig_fp_y_star = figure('Position',[20,20,600,600],'Name','policy function ofY^*_t/A_t   ','Color','w');
 surf(X,Y,Py_star)
 xlabel('log(A_{t+1}/A_t)-gamma_a')
   ylabel('log(Z^b_{t}')   
   zlabel('Y^*_t/A_t')
  
  colormap jet %winter
%  title('policy function of Interest Rate','FontSize',14)
 title('policy function of Y^*_t/A_t ','FontSize',14)
 
 %  view(3)
 view([-65,45])
  
 est_date = datestr(date);   
%  name = ['../output/','fig_pf_r_fp_',est_date];
%          saveas(fig_fp_y_star,name,'fig')
        
 %%     
 
 fig_fp_r = figure('Position',[20,20,600,600],'Name','policy function of  interest rate ','Color','w');
 surf(X,Y,Pr)
 xlabel('log(A_{t+1}/A_t)-gamma_a')
    ylabel('log(Z^b_{t}')   
    zlabel('interest rate')
  colormap jet %winter
%  title('policy function of Interest Rate','FontSize',14)
 title('policy function of interest rate ','FontSize',14)
 
%  view(3)
 view([-65,45])
   
 est_date = datestr(date);   
%  name = ['../output/','fig_pf_r_fp_',est_date];
%          saveas(fig_fp_r,name,'fig')




end