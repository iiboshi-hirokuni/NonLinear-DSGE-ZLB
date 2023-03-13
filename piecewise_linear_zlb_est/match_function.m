function [resids grad init_out E] = match_function(...
      err_vals,err_list,obs_list,current_obs,init_val,...
      constraint1,constraint2,...
      constraint_relax1,constraint_relax2,i1,i2,i_obs_list)

    
    
global M00_

global filtered_errs_switch filtered_errs_init this_period sample_length obs_temp model_temp


% %  [~, i1, ~]=intersect(M00_.endo_names,obs_list,'rows');
% % [i1]=find(ismember(M00_.endo_names,obs_list,'rows'));
% 
%  if size(obs_list,1)==6
%     i1=[9 10 11 12 13 15]';
%  else
%      i1=[9 10 11 12 13]';
%  end
% 
% %  [~, i2, ~]=intersect(M00_.exo_names,err_list,'rows');
% %  [i2]=find(ismember(M00_.exo_names,err_list,'rows'));

%  if size(err_list,1)==6
%      i2=[1 2 3 4 5 6]';
%  else
%      i2=[1 2 3 5 6]';
%  end


nper = 60;
curb_retrench = 0;

maxiter = 10;

    
[ zdatal zdatap zdatass init_out error_flag E] = ...
    solve_two_constraints_nextcall(constraint1, constraint2,...
    constraint_relax1, constraint_relax2,...
    err_vals',err_list,nper,curb_retrench,maxiter,init_val) ;

grad = E(i1,i2);


  
  
  
nobs = size(obs_list,1);
resids = zeros(nobs,1);
% error_flag = 0;

if error_flag == 0

  eval_zdata_transformed_script
  
% -- add observation block in model ---%  
  
%   % put in model file
   

%    for this_obs = 1:nobs
%         eval(['resids(this_obs) = ' deblank(obs_list(this_obs,:)) '_p(1)-current_obs(this_obs);']);
%    end
  resids=zdatap(1,i_obs_list)'-current_obs';
%         keyboard

  
else
  
   resids = resids+100;
   
   
end


plot_temp=0;

if plot_temp == 1
    
  if error_flag==0

    obs_temp(this_period,1:nobs)=current_obs;

    for this_obs = 1:nobs
      eval([ 'model_temp(this_period,this_obs) = ' deblank(obs_list(this_obs,:)) '_p(1);'])
    end

    if this_period==sample_length
    if plot_temp==1
      for index=1:size(obs_list,1)
        subplot(3,2,index)
        plot(obs_temp(:,index)); hold on
        plot(model_temp(:,index),'r')
        axis tight
        title(obs_list(index,:))
      end
    end
    end
    
  end

end



end