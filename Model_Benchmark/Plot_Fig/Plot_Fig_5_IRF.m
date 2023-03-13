% 
% 
%  Figure 5: Impulse Responses to a Monetary Policy Shock
% 
%
%


clear all

addpath('./data_sample')

set(0,'defaultAxesFontSize',11);
set(0,'defaultAxesFontName','century');
set(0,'defaultTextFontSize',11);
set(0,'defaultTextFontName','century');


 addpath('../fun_projection')
 
V=variables;


 load('IRF_model1.mat')
     IRF_model1_pos =  save_imp1;
     IRF_model1_neg =  save_imp2;
     
 load('IRF_model2.mat')
     IRF_model2_pos =  save_imp1;
     IRF_model2_neg =  save_imp2;
     
 load('IRF_model0.mat')
     IRF_model0_pos =  save_imp1;
     IRF_model0_neg =  save_imp2; 
     
     
%    [ V.y V.
     
% h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior','Color','w','File',['IRF_y.tif']);
% for j = 1:2
%    subplot(2,2,j)   
%    i = V.y;
%    hold on
%      imp1 = squeeze(IRF_model1_pos(:,i,(j-1)*100+7));
%         plot(imp1, 'k-o', 'LineWidth',1.5 )
%      imp1 = squeeze(IRF_model2_pos(:,i,(j-1)*100+7));
%         plot(imp1, 'r-x', 'LineWidth',2.5 ) 
%      imp1 = squeeze(IRF_model0_pos(:,i,(j-1)*100+7));
%          plot(imp1, 'b:*', 'LineWidth',3 ) 
%      imp2 = squeeze(IRF_model1_neg(:,i,(j-1)*100+7));
%          plot(imp2, 'k-', 'LineWidth',1.5)
%      imp2 = squeeze(IRF_model2_neg(:,i,(j-1)*100+7));
%          plot(imp2, 'r-', 'LineWidth',2.5) 
%      imp2 = squeeze(IRF_model0_neg(:,i,(j-1)*100+7));
%          plot(imp2, 'b:', 'LineWidth',3)    
%     
%    hold off  
%    title([ 'Output  ' num2str((j-1)*25+1985) ' :Q1  ' ],'FontSize',12)
%    if j == 2
%      legend('Pos: Model 1', 'Pos: Model2', 'Pos: Model w/o ZLB', 'Neg: Model 1', 'Neg: Model2', 'Neg: Model w/o ZLB')
%    end
%    ylim([-0.02 0.02])
%    xlabel('horizon','FontSize',12);
% end   

h_D = figure('Position',[20,20,750,750],'Name','Filtered Posterior','Color','w',...
      'File',['IRF_Fig5']);
for j = 1:2
   subplot(3,2,j)   
  i = V.pi;
   hold on
    imp1 = squeeze(IRF_model1_pos(:,i,(j-1)*100+7));
        l1 = plot(100*imp1, 'k-', 'LineWidth',1.5 );
     imp1 = squeeze(IRF_model2_pos(:,i,(j-1)*100+7));
        l2= plot(100*imp1, 'r-.', 'LineWidth',2.0 ) ;
     imp1 = squeeze(IRF_model0_pos(:,i,(j-1)*100+7));
        l3= plot(100*imp1, 'b:', 'LineWidth',2 ) ;
     imp2 = squeeze(IRF_model1_neg(:,i,(j-1)*100+7));
        l4= plot(100*imp2, 'k-o', 'LineWidth',1.5);
     imp2 = squeeze(IRF_model2_neg(:,i,(j-1)*100+7));
        l5= plot(100*imp2, 'r-x', 'LineWidth',1.5) ;
     imp2 = squeeze(IRF_model0_neg(:,i,(j-1)*100+7));
        l6 = plot(100*imp2, 'b:*', 'LineWidth',2)    ;
   hold off  
   title({'\pi_t  ' [ num2str((j-1)*25+1985) ': Q1  ' ]},'FontSize',12)
   
   ylim(100*[-0.0025 0.0025]);
   if j ==1
   ylabel('[  %  ]')
   end
   
   
   if j == 1
     lgnd=legend([l1, l2, l3 ],...
                 {'Pos: Model 1', 'Pos: Model 2', 'Pos: Model w/o ZLB',...
                  'Neg: Model 1', 'Neg: Model 2', 'Neg: Model w/o ZLB'},... 
                 'FontSize',10);    
     set(lgnd, 'Box', 'off')
     set(lgnd, 'Location','southeast')
   end
   
    if j == 2
     lgnd=legend([l4, l5, l6 ],...
                 {'Neg: Model 1', 'Neg: Model 2', 'Neg: Model w/o ZLB'},... 
                 'FontSize',10);    
     set(lgnd, 'Box', 'off')
     set(lgnd, 'Location','northeast')
   end
   
%     xlabel('horizon','FontSize',10);
end   

% h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior','Color','w','File',['IRF_r.tif']);
for j = 1:2
   subplot(3,2,j+2)     
  i = V.r;
   hold on
   imp1 = squeeze(IRF_model1_pos(:,i,(j-1)*100+7));
        plot(100*imp1, 'k-', 'LineWidth',1.5 );
     imp1 = squeeze(IRF_model2_pos(:,i,(j-1)*100+7));
        plot(100*imp1, 'r-.', 'LineWidth',2.0 ) ;
     imp1 = squeeze(IRF_model0_pos(:,i,(j-1)*100+7));
         plot(100*imp1, 'b:', 'LineWidth',2 ) ;
     imp2 = squeeze(IRF_model1_neg(:,i,(j-1)*100+7));
         l4 = plot(100*imp2, 'k-o', 'LineWidth',1.5);
     imp2 = squeeze(IRF_model2_neg(:,i,(j-1)*100+7));
         l5 = plot(100*imp2, 'r-x', 'LineWidth',1.5) ;
     imp2 = squeeze(IRF_model0_neg(:,i,(j-1)*100+7));
         l6 = plot(100*imp2, 'b:*', 'LineWidth',2)  ;  
   hold off 
   
   title({'R_t  '  [ num2str((j-1)*25+1985) ': Q1  '] },'FontSize',12)
   
   
   
   ylim(100*[-0.01 0.01])
   if j ==1
      ylabel('[  %  ]'  )
   end   
%     xlabel('horizon','FontSize',12);
end   

% h_D = figure('Position',[20,20,900,300],'Name','Filtered Posterior','Color','w','File',['IRF_r.tif']);
for j = 1:2
   subplot(3,2,j+4)     
  i = V.r_lg;
   hold on
   imp1 = squeeze(IRF_model1_pos(:,i,(j-1)*100+7));
        l1= plot(100*imp1, 'k-', 'LineWidth',1.5 );
     imp1 = squeeze(IRF_model2_pos(:,i,(j-1)*100+7));
        l2= plot(100*imp1, 'r-.', 'LineWidth',2.0 ) ;
     imp1 = squeeze(IRF_model0_pos(:,i,(j-1)*100+7));
        l3=  plot(100*imp1, 'b:', 'LineWidth',2 ) ;
     imp2 = squeeze(IRF_model1_neg(:,i,(j-1)*100+7));
         l4 = plot(100*imp2, 'k-o', 'LineWidth',1.5);
     imp2 = squeeze(IRF_model2_neg(:,i,(j-1)*100+7));
         l5 = plot(100*imp2, 'r-x', 'LineWidth',1.5) ;
     imp2 = squeeze(IRF_model0_neg(:,i,(j-1)*100+7));
         l6 = plot(100*imp2, 'b:*', 'LineWidth',2)  ;  
   hold off 
   
   title({'R^*_t  '  [ num2str((j-1)*25+1985) ': Q1  '] },'FontSize',12)
   
   
   ylim(100*[-0.01 0.01])
   if j ==1
      ylabel('[  %  ]'  )
   end   
    xlabel('horizon','FontSize',12);
end   

savefig('./output/Fig5_IRF.fig'); 
 
         