function [ re_xt,  re_prob  ] = resample(xt, prob_xt, nparticles)

w = prob_xt(1,:)/sum(prob_xt(1,:));
total = 0;

re_xt = xt;
re_prob = prob_xt;

for i = 1:nparticles
    
    n = round(w(i)*nparticles);
    
    if (n>0)&&(total<nparticles)
        
        if (total + n) > nparticles
            total_last = nparticles;
        else
            total_last = total + n;
        end
        
        for j = total+1:total_last
            re_xt(:,j)   = xt(:,i);
            re_prob(:,j) = prob_xt(:,j);
        end
    end
    
    total = total + n;
end