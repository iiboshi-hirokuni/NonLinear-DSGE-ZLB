

%% setpathdynare4
location = 'home_windows';
% location = 'home_mac';
setpathdynare4


warning off
do_plots=1;

% folder='C:\Users\iiboshi\Dropbox\matlab_code\DSGE_piese-wise\BNK_zlb_est';
% folder='C:\Users\Hirokuni Iiboshi\Dropbox\matlab_code\Behav_DSGE_Habit\BNK_zlb_est';
folder='../BNK_zlb_est'; 

eval([ 'load ' folder '\mle_estimates_fminsearch params1 params_labels '])
eval([ 'load ' folder '\metropolis_chain theta_history fval_history indx'])
eval([ 'load ' folder '\BNK_model_zlb_results M_ '])
eval([ 'load ' folder '\params_matrix params_matrix '])

if do_plots==1
    
    figure(100)
    for ii=1:numel(params_labels)
        subplot(4,4,ii)
        ksdensity(theta_history(ii,2:indx)); hold on
        
        params_labels = (params_matrix(:,1));
        
        plab = cell2mat(params_labels(ii,:));
        if plab(1:2)=='ST'; axis tight; end
        
        params0 =   cell2mat(params_matrix(:,2));
        params_lo = cell2mat(params_matrix(:,3));
        params_hi = cell2mat(params_matrix(:,4));
        params_mean = cell2mat(params_matrix(:,7));
        params_std = cell2mat(params_matrix(:,8));
        dist_names = params_matrix(ii,6);
        codes = dist_names2codes(params_matrix(ii,6));
        [p6 p7] = get_dist_inputs(codes,params_mean(ii),params_std(ii));
        hold on
        myplot_priors2(codes,params_lo(ii),params_hi(ii),p6,p7,params_labels(ii,:))
        hold on
        plot_lines_green(params1(ii));
        
        v = axis;
        title(cell2mat(params_labels(ii,:)))
        xlim([v(1) max(0.02,v(2))])
    end
    
    
    figure(1010)
    legend('Posterior','Prior','mode')
    subplot(4,4,ii+1)
    plot(-fval_history(1:indx))
    title([ 'The posterior draws ' num2str(max(-fval_history),'%0.4f')])
    
    for ii=1:numel(params_labels)
        subplot(4,4,ii)
        plot(theta_history(ii,2:indx)); hold on; axis tight
        plab = cell2mat(params_labels(ii,:));
        if plab(1:2)=='ST'; axis tight; end
        title(cell2mat(params_labels(ii,:)))
    end
    
end


[ amax b ]=max(-fval_history);

mode_metropolis = theta_history(:,b);

disp('-----------')
disp(folder)
disp(['Number of runs is ' num2str(indx)])
disp(['Objective @ mode  ' num2str(amax,6)])
disp(' ')
disp('PARAMETER   MODE MODE-METRO    PRIOR-DIS    PRIOR-MEAN PRIOR-STD      5%      median    95%    post.mean')
for ii=1:numel(params_labels)
    trspaces=blanks(10-size(char(params_labels(ii,:)),2));
    trspaces2=blanks(13-size(char(params_matrix(ii,6)),2));
    disp([ char(params_labels(ii,:)) trspaces ' ' ...
        num2str(params1(ii),'%0.5f') '   '  ...
        num2str(mode_metropolis(ii),'%0.4f') '      '  ...
        char(params_matrix(ii,6)) trspaces2 '  ' ...
        num2str(cell2mat(params_matrix(ii,7)),'%0.5f')   '    ' ...
        num2str(cell2mat(params_matrix(ii,8)),'%0.5f')   '      ' ...
        num2str(prctile(theta_history(ii,1:indx),5),'%0.5f')   '   ' ...
        num2str(prctile(theta_history(ii,1:indx),50),'%0.5f')  '   ' ...
        num2str(prctile(theta_history(ii,1:indx),95),'%0.5f')  '   ' ...
        num2str(mean(theta_history(ii,1:indx)),'%0.5f') '   '...
        num2str(std(theta_history(ii,1:indx)),'%0.5f'  )   ])
end


disp(' ')
disp('Calibrated parameters')
[ calibrated_names ical]=setdiff(M_.param_names,char(params_labels),'rows');
disp('-----------------')
for i=1:numel(ical)
    trspaces3=blanks(11-size(char(M_.param_names(ical(i),:)),2)) ;
    disp([  M_.param_names(ical(i),:)  '=' trspaces3  num2str(M_.params(ical(i),:),'%0.6f') ';'   ])
    i=i+1;
end






ndraws=numel(fval_history);

cut_percent_obs=0.2;

theta_history_cut = theta_history(:,(cut_percent_obs*indx+1:indx));
fval_history_cut = fval_history(cut_percent_obs*indx+1:indx);

trunc=0.1:0.1:0.9;

md = mhm_marginal_density(theta_history_cut,fval_history_cut,trunc,fval_history);

disp(['Log Marginal data density is ' num2str(mean(md))])

