function [ YY_median ] = fun_plot_shock(yy,ypf,i,a,m,Tobs, graph_title,V,P,ZZ,shock_art);

V = variables; 

rate    =   (1-a)/2;
b       =   a/2/m;      %  
col1    =   0.55;   
col2    =   0.45;
col3    =   0.1;

aa = 1;

yy1_filter = squeeze(ypf(:,i,:));
nsim = size(yy1_filter,2);

if i == V.y
    z_y = squeeze(ypf(2:end,1,:));
    z_y_l = squeeze(ypf(1:end-1,1,:));
    z_mu_a = squeeze(ypf(2:end,7,:));
    z_filter = ZZ(1)*( z_y-z_y_l+z_mu_a+P.gamma_a);
    yy_filter =  z_filter;
    
%     z_filter = squeeze(ypf(:,6,:));
%     yy_filter = ZZ(1)*(log(yy_filter(2:end,:))-log(yy_filter(1:end-1,:))+log(z_filter(2:end,:)) +log(P.gamma));
%     yy_filter = ZZ(1)*(log(yy_filter(2:end,:))-log(yy_filter(1:end-1,:)) +log(P.gamma)); 

 elseif i == V.r
     yy_filter = ZZ(3)*(log(yy1_filter(1:end,:)) ) ; 
     
 elseif i == V.pi 
      yy_filter = ZZ(2)*( log(yy1_filter(1:end,:)) ) ; 
 elseif i == V.r_lg      
       yy_filter = ZZ(3)*log(yy1_filter(1:end,:))  ; 
 elseif i == V.r_star      
       yy_filter = ZZ(3)*log(yy1_filter(1:end,:))  ;  
  elseif i == V.mu_a      
       yy_filter = ZZ(3)*(yy1_filter(1:end,:))  ;        
   elseif i == V.z_b      
       yy_filter = ZZ(3)*(yy1_filter(1:end,:))  ;   
    elseif i >= 9      
       yy_filter = ZZ(3)*(yy1_filter(1:end,:))  ;          
 end    
    
% size(yy_filter)
if (i ==V.y)
    nobs = size(yy_filter,1);
else
    nobs = size(yy,1);
end
    
    
 sort_filter = zeros( nsim, nobs );

%  Sort of sampling forecast
for k = 1:nobs
    sort_filter(:,k) = sort(yy_filter(k,:)',1);
end 

% size(sort_filter)

% %  Graph of Forecast 
% YY_mean = mean(yy_filter(:,:)',1)';
YY_band = zeros(nobs,2*m);
%   Sort_forecast_1 = sort_forecast(:,:);
YY_median = sort_filter(ceil((nsim)*0.5),:)'; 

for j = 1:m 
    YY_band(:,j)        = sort_filter(ceil((nsim)*(rate+b*(j-1))),:)'; 
    YY_band(:,2*m-j+1)  = sort_filter(ceil((nsim)*(1-(rate+b*(j-1)))),:)';
end

YY_f_save = YY_band(:,1);
for j = 1:m-1
    YY_f_save = [ YY_f_save YY_band(:,j+1)-YY_band(:,j) ];
end  
YY_f_save = [ YY_f_save YY_median-YY_band(:,m) YY_band(:,m+1)-YY_median ];
for j = m:-1:2
    YY_f_save = [ YY_f_save YY_band(:,2*m+2-j)-YY_band(:,2*m+1-j) ];
end

s =1983.0;
% s =1985.25;

if i == V.y
 ti = (s+2/4):0.25:s+(nobs+1)/4; 

h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior','Color','w');
title_name = graph_title; %{'Particle Filter'};
% subplot(1,1,1)
hold on 
hh = area(ti,YY_f_save(1:nobs,:)  ) ; 

%     ti = (s+2/4):0.25:s+(nobs)/4;  
for j = 1:m    
    set(hh(j),'FaceColor',[col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3)])    
    set(hh(2*(m+1)-j),'FaceColor',[col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3)])
end
set(hh(m+1),'FaceColor',[1-col1 1-col2 1-col3])
set(hh(1),'FaceColor',[1 1 1])    
set(hh,'LineStyle','none') % Set all to same value

plot(ti, YY_median(1:nobs,1),'LineStyle','-','Color','black', 'LineWidth',2.0); 
xlim([1981  2018])

% out = zeros(Tobs,1);

else
    
    ti = (s+1/4):0.25:s+(nobs)/4;

h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior','Color','w');
title_name = graph_title; %{'Particle Filter'};
% subplot(1,1,1)
hold on 
hh = area(ti,YY_f_save(1:nobs,:)  ) ; 
  
for j = 1:m    
    set(hh(j),'FaceColor',[col1*(1-(j-1)*1/m)+(1-col1) col2*(1-(j-1)*1/m)+(1-col2) col3*(1-(j-1)*1/m)+(1-col3)])    
    set(hh(2*(m+1)-j),'FaceColor',[col1*(1-(j)*1/m)+(1-col1) col2*(1-(j)*1/m)+(1-col2) col3*(1-(j)*1/m)+(1-col3)])
end
set(hh(m+1),'FaceColor',[1-col1 1-col2 1-col3])
set(hh(1),'FaceColor',[1 1 1])    
set(hh,'LineStyle','none') % Set all to same value
l3 = hh(1+m);

l2= plot(ti, YY_median(1:nobs,1),'LineStyle','-','Color','black', 'LineWidth',2.0); 
xlim([1981  2018])

out = zeros(Tobs,1);
% ti = (s+1/4):0.25:s+(nobs)/4;  
end

if   i == 9
      plot(ti, shock_art(1:Tobs,1),'LineStyle','-.','Color','red', 'LineWidth',1.0);
      legend('', '', [num2str(a*100) '% interval'], 'mean', 'actual'); 
      out = YY_median(1:nobs,1);
      shock_i = shock_art(1:Tobs,1);
elseif  i == 10
      plot(ti, shock_art(1:Tobs,2),'LineStyle','-.','Color','red', 'LineWidth',1.0);
      legend('', '', [num2str(a*100) '% interval'], 'mean', 'actual');   
      out = YY_median(1:nobs,1);
      shock_i = shock_art(1:Tobs,2);
elseif  i == 11
      l1= plot(ti, shock_art(1:Tobs,3),'LineStyle','-.','Color','red', 'LineWidth',1.0);
      legend([l1, l2, l3], { 'actual', 'mean', [num2str(a*100) '% interval']});   
      out = YY_median(1:nobs,1);
      shock_i = shock_art(1:Tobs,3);
end

title( title_name(1),'FontSize',14 ) 

hold off


 cusum=zeros(Tobs,2); 
 c1 = 0;
 c2 = 0;
 mean1 = mean(shock_i);
 mean2 = mean(YY_median);
 
 for i = 1:Tobs
   c1 = c1 + sqrt((shock_i(i)-mean1)^2);
   c2 = c2 + sqrt((YY_median(i)-mean1)^2);
   cusum(i,1) = c1; 
   cusum(i,2) = c2;
 end
 
 figure('Name' ,'cusum')
   plot(cusum(:,1), cusum(:,2),'b-')
   hold on
   plot( [0 cusum(end,1)], [0 cusum(end,1)],'r-.'  )
   hold off
   title( title_name(1),'FontSize',14 ) 
   xlabel('Data')
   ylabel('Mean')
   
   
