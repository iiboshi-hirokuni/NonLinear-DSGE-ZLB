function [posterior filtered_errs like prior resids ] =posteriorzlb(current_params,params_labels,lowerbound,upperbound,...
    modnam_00_,modnam_10_,modnam_01_,modnam_11_,...
    constraint1_difference, constraint2_difference,...
    constraint_relax1_difference, constraint_relax2_difference,...
    err_list,obs_list,obs,ntrain,rzlb, codes, p6, p7,tStart)



global M_ oo_ oo00_  M00_ M10_  M01_  M11_
global cof cof10 cof01 cof11 Jbarmat Jbarmat10 Jbarmat01 Jbarmat11 Dbarmat10 Dbarmat01 Dbarmat11 decrulea decruleb
global filtered_errs_switch filtered_errs_init model_temp
global datavec irep xstory fstory

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
eval([ 'save PARAM_EXTRA_EST  ' cell2mat(params_labels') ]);

paramfile_BNK_zlb
% load PARAM_EXTRA_CALIBRATED

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


%---------------------------------------------
% Calculate likelihood
%---------------------------------------------

%try
    
    sample_length = size(obs,1);
    nerrs = size(err_list,1);
    
    filtered_errs_init = zeros(sample_length,nerrs);
    
    [filtered_errs resids Emat requalzero ] = myfilterzlbrnot(constraint1_difference, constraint2_difference,...
                      constraint_relax1_difference, constraint_relax2_difference,err_list,obs_list,obs,rzlb);
   
    nobs=size(filtered_errs,1);
    
    err_list_zlb = setdiff(err_list,'eps_r','rows');
    
    [~, ~, ishocksfe ]=intersect(err_list,err_list,'rows');
    [~, ~, ishocksfe_zlb ]=intersect(err_list_zlb,err_list,'rows');
    
    % If shocks are among the params_labels, estimate them
    if numel(findstr('STD',cell2mat(params_labels(:)')))==nerrs
        
        
        disp('Will estimate shocks without concentrating')
        for i=1:size(err_list,1)
%             eval( [ 'COVMAT1(i,i) = STD_' upper(err_list(i,5)) '^2;'] )
             eval( [ 'COVMAT1(i,i) = STD_' (err_list(i,5)) '^2;'] )
        end
        for i=1:size(err_list_zlb,1)
%             eval( [ 'COVMAT0(i,i) = STD_' upper(err_list_zlb(i,5)) '^2;'] )
            eval( [ 'COVMAT0(i,i) = STD_' (err_list_zlb(i,5)) '^2;'] )
        end
        
    else
        
        % Otherwise concentrate likelihood
        disp('Will concentrate likelihood')
        for i=1:size(err_list,1)
            COVMAT1(i,i) = std(filtered_errs(:,i))^2 ;
            evalc([ 'STD_' upper(err_list(i,5)) '= COVMAT1(i,i)^0.5']) ;
        end
        COVMAT0=COVMAT1(ishocksfe_zlb,ishocksfe_zlb');
        
    end
    
    
    
    %-------------------------------------
    % Calculate the selection matrix
    %-------------------------------------
    
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
            % When r=0, remove r from list of observables and monetary shock from shocks
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
        disp('Penalize params outside bound')
        like = like + 1e6*max(abs(params-current_params)) ;
    end
    
    maxresid = max(abs(resids(:)));
    if maxresid>1e-3
        disp('Penalize failure of residuals to be zero')
        like = like + sum(resids(:).^2)*1e7;
    end
    
   
    [prior] = -priordens(params, codes, p6, p7, lowerbound, upperbound,1);
    if prior == Inf
       % If parameters outside prior bound, minus prior density is very large
       prior= 1e8;
    end

    
% catch myerror
%     
%     disp(myerror)
%     filtered_errs = obs*0+nan;
%     like = 1e8;
%     prior = 1e8;
%     resids = obs*0+nan;
%     
% end

if isinf(like)==1
    like=1e8;
end
    

%---------------------------------------------
% Calculate posterior
%---------------------------------------------

% remember that the likelihood has already been multiplied by -1
% hence, posterior is -1 times the log of the prior
posterior = like+prior;

if isnan(posterior)==1
    keyboard
end

%-------------------------------------
% Display info on screen
%-------------------------------------

disp(['Current posterior is ' num2str(-1*round(posterior*10)/10), ' ( = ', num2str(-1*round(like*10)/10),...
     ' + ', num2str(-1*round(prior*10)/10) ' )  , iteration # ' num2str(irep) ])
% datavec(irep,:) = [ posterior SIGMA TETAP TETAW EC RHOD TAYLOR_P TAYLOR_R TAYLOR_Y BETA1];
datavec(irep,:) = [ posterior   gamma_a omega  kappa pi_star psi_pi psi_y rho_a rho_b rho_r ];

fstory(irep,:) = [ posterior like prior ];
xstory(irep,:) = [ current_params ];
save datavec datavec fstory xstory params*

tElapsed = toc(tStart);
disp(['Total computing time: ', num2str(floor(tElapsed/3600)), 'h',num2str(floor(mod(tElapsed,3600)/60)), 'm'  ]);


if datavec(irep,1)==min(datavec(:,1))
    disp('minimum found, save into mle_estimates_temp_test')
    disp(' gamma_a omega  kappa pi_star psi_pi psi_y rho_a rho_b rho_r')
% %     disp('     SIGMA      TETAP    TETAW     EC     RHOD     TAYLOR_P    TAYLOR_R  TAYLOR_Y    BETA1')
    disp(datavec(irep,2:end))
    
    params1=params;
    save mle_estimates_temp_test like posterior prior datavec filtered_errs err_list obs_list obs params*
end

irep = irep+1;

disp(' ')

%  save temp
% params_matrix = {  ...
%   'BETA1 '    0.9925    0.9000      0.9940     1     'BETA_PDF'      0.8743 0.025^2
%   'EC '       0.6       0.0001      0.9999     1     'BETA_PDF'      0.7    0.10
%   'EH '       0.8       0.0001      0.9999     1     'BETA_PDF'      0.7    0.10
%   'PHIK '     4.0       0.0001     50.0000     1     'GAMMA_PDF'     5      2
%   'SIGMA '    0.5       0.0001      0.9900     1     'BETA_PDF'      0.333  0.20
%   'TAYLOR_P ' 1.7       1.0000      5.0000     1     'NORMAL_PDF'    1.5    0.25
%   'TAYLOR_R ' 0.5       0.0000      0.9500     1     'BETA_PDF'      0.75   0.1
%   'TAYLOR_Y ' 0.1       0.0000      0.5000     1     'NORMAL_PDF'    0.125  0.025
%   'TETAP '    0.92  	0.1000      0.9900     1     'BETA_PDF'      0.5    0.075
%   'TETAW '    0.92      0.1000   	0.9900     1     'BETA_PDF'      0.5    0.075
%   'RHOD '     0.40      0.0000  	0.9990     1     'BETA_PDF'      0.75   0.10
%   'RHO_J '    0.98      0.0000  	0.9990     1     'BETA_PDF'      0.75   0.10
%   'RHO_K '    0.8       0.0000  	0.9990     1     'BETA_PDF'      0.75   0.10 
%   'RHO_R '    0.6       0.0000      0.9990     1     'BETA_PDF'      0.5    0.10 
%   'RHO_Z '    0.7       0.0000  	0.9990     1     'BETA_PDF'      0.75   0.10
%   'STD_J '    0.08      0.0001  	0.9999     1     'INV_GAMMA_PDF' 0.01   1
%   'STD_K '    0.04      0.0001   	0.9999     1     'INV_GAMMA_PDF' 0.01   1
%   'STD_P '    0.003     0.0001  	0.9999     1     'INV_GAMMA_PDF' 0.01   1
%   'STD_R '    0.001     0.0001  	0.9999     1     'INV_GAMMA_PDF' 0.01   1
%   'STD_W '    0.010     0.0001  	0.9999     1     'INV_GAMMA_PDF' 0.01   1
%   'STD_Z '    0.016     0.0001  	0.9999     1     'INV_GAMMA_PDF' 0.01   1     } ;   

