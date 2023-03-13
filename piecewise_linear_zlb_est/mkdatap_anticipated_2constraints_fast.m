function [zdata Ecurrent ]=mkdatap_anticipated_2constraints_fast(nperiods,decrulea,decruleb,...
    cof,Jbarmat,...
    cof10,Jbarmat10,Dbarmat10,...
    cof01,Jbarmat01,Dbarmat01,...
    cof11,Jbarmat11,Dbarmat11,...
    regime1,regimestart1,...
    regime2,regimestart2,...
    violvecbool,endog_,exog_,...
    irfshock,scalefactormod,init)


nvars = size(endog_,1);


if nargin<16
    init=zeros(nvars,1);
end

if nargin<15
    scalefactormod=1;
end


nshocks = size(irfshock,1);

exog_cell = cellstr(exog_);

for i = 1:nshocks
    shockpos = find(strcmp(irfshock(i,:),exog_cell));
    
    if ~isempty(shockpos)
        irfshockpos(i) = shockpos;
    else
        error(['Shock ',irfshock(i,:),' is not in the model']);
    end
end



Cbarmat = cof(:,1:nvars);
Bbarmat = cof(:,nvars+1:2*nvars);
Abarmat = cof(:,2*nvars+1:3*nvars);


% cofstar contains the system for the model when the constraint binds


Cbarmat10 = cof10(:,1:nvars);
Bbarmat10 = cof10(:,nvars+1:2*nvars);
Abarmat10 = cof10(:,2*nvars+1:3*nvars);

Cbarmat01 = cof01(:,1:nvars);
Bbarmat01 = cof01(:,nvars+1:2*nvars);
Abarmat01 = cof01(:,2*nvars+1:3*nvars);

Cbarmat11 = cof11(:,1:nvars);
Bbarmat11 = cof11(:,nvars+1:2*nvars);
Abarmat11 = cof11(:,2*nvars+1:3*nvars);

% get the time-dependent decision rules
nregimes1 = length(regime1);
nregimes2 = length(regime2);

Tmax = max([regimestart1(nregimes1) regimestart2(nregimes2)])-1;  % Tmax is the position of the last period
% when the constraint binds


if Tmax > 0
    P = zeros(nvars,nvars,Tmax);
    D = zeros(nvars,Tmax);
    
    
    
    if (violvecbool(Tmax,1) & ~violvecbool(Tmax,2))
    noinvmat = ((Abarmat10*decrulea+Bbarmat10));
    P(:,:,Tmax) = -noinvmat\Cbarmat10;
    D(:,Tmax) = -noinvmat\Dbarmat10;  
    elseif (violvecbool(Tmax,1) & violvecbool(Tmax,2))
    noinvmat = ((Abarmat11*decrulea+Bbarmat11));
    P(:,:,Tmax) = -noinvmat\Cbarmat11;
    D(:,Tmax) = -noinvmat\Dbarmat11;
    else
    noinvmat = ((Abarmat01*decrulea+Bbarmat01));
    P(:,:,Tmax) = -noinvmat\Cbarmat01;
    D(:,Tmax) = -noinvmat\Dbarmat01;  
    end
    
    
    
    
    for i = Tmax-1:-1:1        
        
        if (violvecbool(i,1) & ~violvecbool(i,2))
            noinvmat = (Bbarmat10+Abarmat10*P(:,:,i+1));
            P(:,:,i)=-noinvmat\Cbarmat10;
            D(:,i) = -noinvmat\(Abarmat10*D(:,i+1)+Dbarmat10);            
        elseif (~violvecbool(i,1) & violvecbool(i,2))
            noinvmat = (Bbarmat01+Abarmat01*P(:,:,i+1));
            P(:,:,i)=-noinvmat\Cbarmat01;
            D(:,i) = -noinvmat\(Abarmat01*D(:,i+1)+Dbarmat01);
        elseif (violvecbool(i,1) & violvecbool(i,2))
            noinvmat = (Bbarmat11+Abarmat11*P(:,:,i+1));
            P(:,:,i)=-noinvmat\Cbarmat11;
            D(:,i) = -noinvmat\(Abarmat11*D(:,i+1)+Dbarmat11);
        else
            noinvmat = (Bbarmat+Abarmat*P(:,:,i+1));
            P(:,:,i)=-noinvmat\Cbarmat;
            D(:,i) = -noinvmat\(Abarmat*D(:,i+1));
        end
        
    end

    
% Double check the appropriate invmat in each case
% right now -- inherited from previous loop
if Tmax > 1
    
    if ( ~violvecbool(1,1) & violvecbool(1,2) )
        E = -noinvmat\Jbarmat01;
    elseif ( violvecbool(1,1) & ~violvecbool(1,2) )
        E = -noinvmat\Jbarmat10;
    elseif ( violvecbool(1,1) & violvecbool(1,2) )
        E = -noinvmat\Jbarmat11;
    else
        E = -noinvmat\Jbarmat;
    end
    
else  % Tmax is equal to 1    
    if ( ~violvecbool(1,1) & violvecbool(1,2) )
       noinvmat = ((Abarmat01*decrulea+Bbarmat01));
       E = -noinvmat\Jbarmat01;  
    elseif ( violvecbool(1,1) & violvecbool(1,2) )
        noinvmat = ((Abarmat11*decrulea+Bbarmat11));
        E = -noinvmat\Jbarmat11;     
    else
        noinvmat = ((Abarmat10*decrulea+Bbarmat10));
        E = -noinvmat\Jbarmat10;
            
    end
    
end

Ecurrent = E;

else
  
  Ecurrent = decruleb;

end

% generate data
% history will contain data, the state vector at each period in time will
% be stored columnwise.
history = zeros(nvars,nperiods+1);
history(:,1) = init;
errvec = zeros(size(exog_,1),1);

for i = 1:nshocks
    errvec(irfshockpos(i)) = scalefactormod(i);
end

% deal with shocks
irfpos =1;
if irfpos <=Tmax
    history(:,irfpos+1) = P(:,:,irfpos)* history(:,irfpos)+...
        D(:,irfpos) + E*errvec;
else
    history(:,irfpos+1) = decrulea*history(:,irfpos)+decruleb*errvec;
end

% all other periods
for irfpos=2:nperiods+1
    if irfpos <=Tmax
        history(:,irfpos+1) = P(:,:,irfpos)* history(:,irfpos)+...
            D(:,irfpos);
    else
        history(:,irfpos+1) = decrulea*history(:,irfpos);
    end
end


history=history';
zdata = history(2:end,:);