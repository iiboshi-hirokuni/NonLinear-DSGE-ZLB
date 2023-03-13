function [filtered_errs resids Emat requalzero ] = myfilterzlbrnot(constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,err_list,obs_list,obs,rzlb)

global M00_

global filtered_errs_switch filtered_errs_init this_period sample_length obs_temp model_temp


%-------------------------------------
% Filter shocks
%-------------------------------------

options_fsolve = optimset('Display','None','MaxFunEvals',1e10,'MaxIter',1e5,'TolFun',1e-4,...
    'Algorithm','trust-region-dogleg');
verbose = 0;

sample_length = size(obs,1);
nerrs = size(err_list,1);
init_val = zeros(M00_.endo_nbr,1);
err_vals = zeros(nerrs,1);
err_vals_zlb = zeros(nerrs-1,1);

resids = zeros(sample_length,nerrs);

pos_r = find(strcmp('data_r',cellstr(obs_list)));
pos_eps_r = find(strcmp('e_r',cellstr(err_list)));  %% 2019/07/11

my_list = 1:nerrs;

err_list_zlb = setdiff(err_list,'e_r','rows'); %% 2019/07/11
obs_list_zlb = setdiff(obs_list,'data_r','rows');
obs_list_rnot = union(obs_list_zlb,'data_rnot','rows');

maxiters = 10;
requalzero=zeros(1,sample_length);
tolresidr = 1e-5;
tolsolve = 1e-6;
tolzlb = 1e-5;


for this_period=1:sample_length
    
    if verbose
        tic
    end
    
    current_obs = obs(this_period,:);
    init_val_old = init_val;
    
       
    if current_obs(pos_r)<rzlb+tolzlb
        
        requalzero(this_period)=1;
        
        current_obs_zlb = current_obs(find(my_list~=pos_r));
        err0 = filtered_errs_init(this_period,1:numel(err_vals_zlb));
        
        
        [ err_vals_out_zlb em ]= csolve_grad('match_function',...
            err0',tolsolve,maxiters,...
            err_list_zlb,obs_list_zlb,current_obs_zlb,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference);
        
        err_vals_out = zeros(nerrs,1);
        err_vals_out(find(my_list~=pos_eps_r)) = err_vals_out_zlb ;
        filtered_errs(this_period,:)=err_vals_out';
        
        [resids(this_period,:), ~, init_val, Emat(:,:,this_period) ] = match_function(...
            err_vals_out,err_list,obs_list,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,constraint_relax1_difference,constraint_relax2_difference);
                
    end
        
    
    if current_obs(pos_r)<rzlb+tolzlb && max(abs(resids(this_period,:)))>tolresidr
        
        if verbose==1
            disp([ this_period NaN 100*resids(this_period,:)])
        end
        
        err0 = filtered_errs_init(this_period,1:numel(err_vals_zlb));
        
        [ err_vals_out_zlb ] = fsolve(@(err_vals) match_function(...
            err_vals,err_list_zlb,obs_list_zlb,current_obs_zlb,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,...
            constraint_relax2_difference), err0',options_fsolve);
        
        err_vals_out = zeros(nerrs,1);
        err_vals_out(find(my_list~=pos_eps_r)) = err_vals_out_zlb ;
        filtered_errs(this_period,:)=err_vals_out';
        
        [resids(this_period,:), ~, init_val, Emat(:,:,this_period) ] = match_function(...
            err_vals_out,err_list,obs_list,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,constraint_relax1_difference,constraint_relax2_difference);
        
        if verbose==1
            disp([ this_period NaN 100*resids(this_period,:)])
            disp('I just called fsolve at ZLB, compare residuals before and after ')
        end
        
    end
    
       
    if current_obs(pos_r)>=rzlb+tolzlb
        
        err0 = filtered_errs_init(this_period,1:numel(err_vals));
        
        
        [ err_vals_out em ] = csolve_grad('match_function',...
            err0',tolsolve,maxiters,...
            err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference);
        
        filtered_errs(this_period,:)=err_vals_out';
        
        [ resids(this_period,:), ~, init_val, Emat(:,:,this_period)] = match_function(...
            err_vals_out,err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference);
        
    end
        
    
    if (resids(this_period,pos_r))>tolresidr
        
        if verbose==1
            disp([ this_period NaN resids(this_period,:)])
        end
        
        err0 = filtered_errs_init(this_period,1:numel(err_vals));
        
        [ err_vals_out em ] = csolve_grad('match_function',...
            err0',tolsolve,maxiters,...
            err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference);
        
        filtered_errs(this_period,:)=err_vals_out';
        
        [ resids(this_period,:), ~, init_val, Emat(:,:,this_period)] = match_function(...
            err_vals_out,err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference);
        
        if verbose==1
            disp([ this_period NaN resids(this_period,:)])
            disp(err_vals_out(4))
            disp('I just added monetary shocks because notional rate was above zero')
            disp(' ')
            keyboard
        end
        
    end
    
       
    if max(abs(resids(this_period,:)))>tolresidr && current_obs(pos_r)>=rzlb+tolzlb
        
        if verbose==1
            disp([ this_period NaN 100*resids(this_period,:)])
        end
        
        err0 = filtered_errs_init(this_period,1:numel(err_vals));
        
        [ err_vals_out ] = fsolve(@(err_vals) match_function(...
            err_vals,err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,...
            constraint_relax2_difference), err0',options_fsolve);
        
        filtered_errs(this_period,:)=err_vals_out';
        
        [ resids(this_period,:), ~, init_val, Emat(:,:,this_period)] = match_function(...
            err_vals_out,err_list,obs_list_rnot,current_obs,init_val_old,...
            constraint1_difference,constraint2_difference,...
            constraint_relax1_difference,constraint_relax2_difference);
        
        if verbose==1
            disp([ this_period NaN 100*resids(this_period,:)])
            disp('I just called fsolve, compare residuals before and after ')
            keyboard
        end
        
    end
    
        
    if max(abs(resids(this_period,:)))>0.05
        init_val_old=0*init_val_old;
        error('huge resids, give up')
    end
    
    
    if verbose
        toc
        elapsed_time(this_period,:) = toc;
    end
    
    
    
end


end


