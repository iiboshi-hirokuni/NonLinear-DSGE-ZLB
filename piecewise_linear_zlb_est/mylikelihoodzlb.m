function [like filtered_errs] = mylikelihoodzlb(current_params,params_labels,lowerbound,upperbound,...
    modnam_00_,modnam_10_,modnam_01_,modnam_11_,...
    constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,...
    err_list,obs_list,obs,ntrain,rzlb)


global M_ oo_

global oo00_  M00_ M10_  M01_  M11_

global cof cof10 cof01 cof11 ...
    Jbarmat Jbarmat10 Jbarmat01 Jbarmat11 ...
    Dbarmat10 Dbarmat01 Dbarmat11 ...
    decrulea decruleb

global filtered_errs_switch filtered_errs_init model_temp


if size(current_params,1)<size(current_params,2)
    current_params=current_params';
    params=current_params;
else
    params = min(max(lowerbound,current_params),upperbound);
end


% save to disk -- will be read by the parameter file invoked
% when the model is solved again.
for i=1:numel(params_labels)
    x=cell2mat(params_labels(i));
    evalc([ x '= params(' num2str(i) ')']) ;
end
eval([ 'save PARAM_EXTRA_BABY ' cell2mat(params_labels') ]);

paramfile_baby00


eval(modnam_00_)
oo00_ = oo_;
M00_ = M_;





zdatass = oo00_.dr.ys;

[hm1,h,hl1,Jbarmat] = get_deriv(M00_,zdatass);
cof = [hm1,h,hl1];


M10_.params = M00_.params;
[hm1,h,hl1,Jbarmat10,resid] = get_deriv(M10_,zdatass);
cof10 = [hm1,h,hl1];
Dbarmat10 = resid;

M01_.params = M00_.params;
[hm1,h,hl1,Jbarmat01,resid] = get_deriv(M01_,zdatass);
cof01 = [hm1,h,hl1];
Dbarmat01 = resid;

M11_.params = M00_.params;
[hm1,h,hl1,Jbarmat11,resid] = get_deriv(M11_,zdatass);
cof11 = [hm1,h,hl1];
Dbarmat11 = resid;

[decrulea,decruleb]=get_pq(oo00_.dr);

try

sample_length = size(obs,1);
nerrs = size(err_list,1);

if ~filtered_errs_switch
    filtered_errs_init = zeros(sample_length,nerrs);
end

[filtered_errs resids Emat requalzero ] = myfilterzlbrnot(constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,err_list,obs_list,obs,rzlb);
nobs=size(filtered_errs,1);

filtered_errs_switch = 1;
filtered_errs_init = filtered_errs;

err_list_zlb = setdiff(err_list,'eps_r','rows');

[~, ~, ishocksfe ]=intersect(err_list,err_list,'rows');
[~, ~, ishocksfe_zlb ]=intersect(err_list_zlb,err_list,'rows');


% 1 - Use this if you DO NOT want to concentrate the likelihood
% Assign to COVMAT the values of standard deviation declared in paramfile

for i=1:size(err_list,1)
    eval( [ 'COVMAT1(i,i) = STD_' upper(err_list(i,5)) '^2;'] )
end
for i=1:size(err_list_zlb,1)
    eval( [ 'COVMAT0(i,i) = STD_' upper(err_list_zlb(i,5)) '^2;'] )
end


% 2 - Use this if you WANT to concentrate the likelihood
% Estimate COVMAT based on filtered errors
% for i=1:size(err_list,1)
%     COVMAT1(i,i) = std(filtered_errs(:,i))^2 ;
%     evalc([ 'STD_' upper(err_list(i,5)) '= COVMAT1(i,i)^0.5']) ;
% end
% COVMAT0=COVMAT1(ishocksfe_zlb,ishocksfe_zlb');




% Technical issues ZLB related. If R @ ZLB, Emat has a row of zeros corresponding
% to R and Gmat is not full rank. In that case we drop one row from the selector matrix
%
% This corresponds to picking obs_list_rnot instead of obs_list below


obs_list_zlb = setdiff(obs_list,'data_r','rows');
selector_matrix0=zeros(size(obs_list_zlb,1),M_.endo_nbr);
for iobs=1:size(obs_list_zlb,1)
    [~, ~, iobscols0]=intersect(obs_list_zlb,M_.endo_names,'rows');
    selector_matrix0(iobs,iobscols0(iobs))=1;
end

selector_matrix1=zeros(size(obs_list,1),M_.endo_nbr);
for iobs=1:size(obs_list,1)
    [~, ~, iobscols1]=intersect(obs_list,M_.endo_names,'rows');
    selector_matrix1(iobs,iobscols1(iobs))=1;
end


[~, ~, ishocks ]=intersect(err_list,M_.exo_names,'rows');
[~, ~, ishocks_zlb ]=intersect(err_list_zlb,M_.exo_names,'rows');




likeall=0;

for t = 1:nobs
    
    if requalzero(t)==1
        % When r=0, remove r from list of observables and monetary shock
        % from shocks
        Gmat0 = selector_matrix0*Emat(:,ishocks_zlb,t);
        log_det_jacobian(t) = log(det(COVMAT0)) + 2*log(abs(det(Gmat0)));
        trace_term(t) = filtered_errs(t,ishocksfe_zlb)*inv(COVMAT0)*filtered_errs(t,ishocksfe_zlb)';
    else
        Gmat1  = selector_matrix1*Emat(:,ishocks,t);
        log_det_jacobian(t) = log(det(COVMAT1)) + 2*log(abs(det(Gmat1)));
        trace_term(t) = filtered_errs(t,ishocksfe)*inv(COVMAT1)*filtered_errs(t,ishocksfe)';
    end
    
    likei(t,1) = log_det_jacobian(t)/2 + trace_term(t)/2;
    likeall = likeall + likei(t);
    
end

like = sum(likei(ntrain+1:end));


if max(abs(params-current_params))>1e-8
    like = like + 1e6*max(abs(params-current_params)) ;
end

maxresid = max(abs(resids(:)));
if maxresid>1e-3
    disp('Penalize failure of residuals to be zero')
    like = like + sum(resids(:).^2)*1e7;
end

catch myerror

      filtered_errs = obs*0+nan;
      like = 1e8;
      resids = obs*0+nan;

end


