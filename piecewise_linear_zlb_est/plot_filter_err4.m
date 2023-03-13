
err_name={'shock a','shock b','shock MP'};

filtered_errs_3 =filtered_errs;
filtered_errs_3(103:106,3) =-1*filtered_errs(103:106,3);
j=2;

figure(101)
for i=1:3
subplot(3,1,i)
 plot(tt_obs,0*filtered_errs(j:end,i),'k:');
 hold on
 plot(tt_obs,filtered_errs(j:end,i));
 
 if i==3  
%      plot(tt_obs,filtered_errs_3(j:end,i),'r:.'); 
 end
 
 hold off 
 title([ char(err_name(i)) ],'FontSize',12);
 axis tight
end 

 