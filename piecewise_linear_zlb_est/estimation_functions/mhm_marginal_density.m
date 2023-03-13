function md = mhm_marginal_density(thetaHistorycut,fvalHistorycut,truncation_points,fvalHistory)

ni = numel(truncation_points);


lpost_mode=fvalHistory(1);

chain_length=size(thetaHistorycut,2);
npar = size(thetaHistorycut,1);


MU = mean(thetaHistorycut,2);

varcov = (cov(thetaHistorycut'));

logdetSIGMA = log(det(varcov));

invSIGMA = inv(varcov);

logpi2 = log(2*pi);

for ii=1:ni
    
    critval = chi2inv(truncation_points(ii),npar);
    
    log_truncation_points = log(truncation_points(ii));
    
    tmp = 0;
    
    for i = 1:chain_length
        
        this_theta = thetaHistorycut(:,i);
        
        deviation = (this_theta-MU)'*invSIGMA*(this_theta-MU);
        
        if deviation < critval
            
            lftheta = -log_truncation_points...
                -(npar/2)*logpi2...
                -0.5*logdetSIGMA...
                -0.5*deviation ;
            tmp =tmp+exp(lftheta + fvalHistorycut(i) - fvalHistorycut(1));
            
        end
        
        %                     deviation  = (x2(i,:)-MU)*invSIGMA*(x2(i,:)-MU)';
        %                     if deviation <= critval
        %                         lftheta = -log(p)-(npar*log(2*pi)+log(detSIGMA)+deviation)/2;
        %                         tmp = tmp + exp(lftheta - logpo2(i) + lpost_mode);
        %                     end
        %         marginal(linee,:) = [p, lpost_mode-log(tmp/((TotalNumberOfMhDraws-TODROP)*nblck))];
        
    end
    
    md(ii,1)=-fvalHistorycut(1)-log(tmp/chain_length);
    
end
