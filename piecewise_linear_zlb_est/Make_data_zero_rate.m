
data_r = data(:,3);

for i=1:length(data_r)
   if data_r(i) < 0.03
      data_r_1(i) = 0;
   else
      data_r_1(i) = data_r(i); 
   end   
end    

figure(1)
plot(data_r)
hold on
  plot(data_r_1)
hold off
 ylim([0 0.2])

data_r= (data_r_1-data_r_ss1)'/100;
