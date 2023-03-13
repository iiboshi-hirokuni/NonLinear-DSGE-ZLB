function [ YY_median ] = fun_plot(yy,ypf,i,a,m,Tobs, graph_title,V,P,ZZ);

V = variables; 

rate    =   (1-a)/2;
b       =   a/2/m;      %  
col1    =   0.75;   
col2    =   0.65;
col3    =   0.2;

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
    elseif i == 9      
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

plot(ti, YY_median(1:nobs,1),'LineStyle','-','Color','black', 'LineWidth',2.0); 
xlim([1981  2018])

out = zeros(Tobs,1);
% ti = (s+1/4):0.25:s+(nobs)/4;  
end

if i ==V.y       % output
%  size(ti,2)  
%  size(yy(2:end,1))
      plot(ti, yy(2:end,1),'LineStyle','--','Color','red', 'LineWidth',2.0);
      legend('', '', [num2str(a*100) '% interval'], 'mean', 'actual');
%      est_y = yy(1:Tobs,1);
%      save('./result/est_y.mat', 'est_y');
%      ylim([-4.0 5.0])

elseif i == V.pi  % inlation
     plot(ti, yy(1:Tobs,2),'LineStyle','--','Color','red', 'LineWidth',2.0);
      legend('', '', [num2str(a*100) '% interval'], 'mean', 'actual');
%      est_i = yy(1:Tobs,2);
%      save('./result/est_i.mat', 'est_i');
%      ylim([-4 7])  
elseif i == V.r   % interest rate
     plot(ti, yy(1:Tobs,3),'LineStyle','--','Color','red', 'LineWidth',2.0);
      legend('', '', [num2str(a*100) '% interval'], 'mean', 'actual');
%      est_r = yy(1:Tobs,3);
%      save('./result/est_r.mat', 'est_r');
%      ylim([-1 10])     
elseif i ==10    %% 
     out = YY_median(1:nobs,1);
%     plot(ti, yy(1:Tobs,7),'LineStyle','--','Color','red', 'LineWidth',1.5);
%      ylim([0.5  2])   
elseif i ==V.r_star    %%
     out = YY_median(1:nobs,1);
end

title( title_name(1),'FontSize',14 ) 

hold off
