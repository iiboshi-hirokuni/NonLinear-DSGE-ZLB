%------------------------------------------
% Feed filtered shocks back into model
% -- add observation block in model, check that final_samples matches obs
%------------------------------------------

init0=zeros(M00_.endo_nbr,1);
%   init0(strmatch('rnot',M00_.endo_names))=obs(1,strmatch('r',obs_list));
%   init0(strmatch('dp',M00_.endo_names))=obs(1,strmatch('dp',obs_list));
%   init0(strmatch('k',M00_.endo_names))=k_ss*obs(1,strmatch('ik',obs_list));


sample_length = size(obs,1);
load PARAM_EXTRA_CALIBRATED
for i=1:numel(params_labels)
  evalc([ cell2mat(params_labels(i)) '= params1(' num2str(i) ')']) ;
end
eval([ 'save PARAM_EXTRA_BABY ' cell2mat(params_labels') ]);




filtered_errs2=filtered_errs;
[amax imax]=min(obs(:,strmatch('ctot',obs_list)));


%   imax=10
% imax=find(tt_obs==2009.375);
%     filtered_errs2(imax:end,:)=0;

extra_t=0;


[zdatal zdatap zdatass oo00_ M00_ ] = solve_two_constraints_fast2_temp1(...
  modnam_00,modnam_10,modnam_01,modnam_11,...
  constraint1, constraint2,...
  constraint_relax1, constraint_relax2,...
  filtered_errs2,err_list,sample_length+extra_t,curb_retrench,maxiter,init0);
tt_obs2=[tt_obs
  tt_obs(end)+(1:extra_t)'/4];

for i=1:M00_.endo_nbr
  eval([deblank(M00_.endo_names(i,:)),'_l=zdatal(1:sample_length+extra_t,i);']);
  eval([deblank(M00_.endo_names(i,:)),'_p=zdatap(1:sample_length+extra_t,i);']);
  eval([deblank(M00_.endo_names(i,:)),'_ss=zdatass(i);']);
end



figure
final_sample = [];
for index=1:size(obs_list,1)
  subplot(3,3,index)
  eval([ 'final_sample(:,index) = ' deblank(obs_list(index,:)) '_p;'])
  plot(tt_obs,obs(:,index)); hold on
  plot(tt_obs2,eval([deblank(obs_list(index,:)) '_p']),'r')
  axis tight
  title(obs_list(index,:))
end
legend('data','model')



%   subplot(4,2,7)
%   plot(tt_obs,rnot_m,'r'); axis tight



%   display('The magnitude of the largest difference between the original sample and the synthetic sample is:')
%   max(max(abs(final_sample-obs)))

