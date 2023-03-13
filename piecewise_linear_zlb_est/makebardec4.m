function makebardec4(cmp1,titlelist,legendlist,figlabel,yearshock,step,maxperiod,ylabels,opts,...
    zdata1,zdata2,zdata3,zdata4,zdata5,zdata6,zdata7,zdata8)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

nobs = size(zdata1,1);

%if maxperiod is defined, the nobs is calculated as such
if maxperiod > 0
    nobs = (maxperiod-yearshock)/step;
end

ndsets=7;       % default, changed below as applicable
if nargin==10
    zdata2=nan*zdata1;
    zdata3=nan*zdata1;
    zdata4=nan*zdata1;
    zdata5=nan*zdata1;
    zdata6=nan*zdata1;
    zdata7=nan*zdata1;
    ndsets =1;
    zdata(:,:,1) = zdata1;
elseif nargin==11
    zdata3=nan*zdata1;
    zdata4=nan*zdata1;
    zdata5=nan*zdata1;
    zdata6=nan*zdata1;
    zdata7=nan*zdata1;
    ndsets =2;
    zdata(:,:,1) = zdata1;
    zdata(:,:,2) = zdata2;
elseif nargin == 12
    zdata4 =nan*zdata1;
    zdata5 =nan*zdata1;
    zdata6=nan*zdata1;
    zdata7=nan*zdata1;
    ndsets=3;
    zdata(:,:,1) = zdata1;
    zdata(:,:,2) = zdata2;
    zdata(:,:,3) = zdata3;
elseif nargin == 13
    zdata5 =nan*zdata1;
    zdata6=nan*zdata1;
    zdata7=nan*zdata1;
    ndsets=4;
    zdata(:,:,1) = zdata1;
    zdata(:,:,2) = zdata2;
    zdata(:,:,3) = zdata3;
    zdata(:,:,4) = zdata4;
elseif nargin == 14
    zdata6 =nan*zdata1;
    zdata7=nan*zdata1;
    ndsets=5;
    zdata(:,:,1) = zdata1;
    zdata(:,:,2) = zdata2;
    zdata(:,:,3) = zdata3;
    zdata(:,:,4) = zdata4;
    zdata(:,:,5) = zdata5;
elseif nargin == 15
    zdata7=nan*zdata1;
    ndsets=6;
    zdata(:,:,1) = zdata1;
    zdata(:,:,2) = zdata2;
    zdata(:,:,3) = zdata3;
    zdata(:,:,4) = zdata4;
    zdata(:,:,5) = zdata5;
    zdata(:,:,6) = zdata6;
elseif nargin == 16
    zdata8=nan*zdata1;
    ndsets=7;
    zdata(:,:,1) = zdata1;
    zdata(:,:,2) = zdata2;
    zdata(:,:,3) = zdata3;
    zdata(:,:,4) = zdata4;
    zdata(:,:,5) = zdata5;
    zdata(:,:,6) = zdata6;
    zdata(:,:,7) = zdata7;
elseif ((nargin>=16) || (nargin <=4))
    error ('makechart takes 5 to 11 arguments')
end


if yearshock>-100
    xvalues = yearshock+(0:nobs-1)'*step; % Matteo plot year on x axis
else
    xvalues = (1:nobs)'; % Matteo plot year on x axis
end



nvars = size(titlelist,1);
nshocks = size(legendlist,1)-1;

if nvars==1
    nrows=1;
    ncols = 1;
elseif nvars==2
    nrows =3;
    ncols = 1;
elseif nvars == 3
    nrows = 3;
    ncols = 1;
elseif nvars==4
    nrows = 2;
    ncols = 2;
elseif (nvars==5 || nvars ==6)
    nrows = 3;
    ncols = 2;
elseif (nvars==7 || nvars==8)
    nrows = 4;
    ncols = 2;
elseif nvars>8 && nvars<=12
    nrows = 3;
    ncols = 4;
elseif nvars>12 && nvars<=15
    nrows = 5;
    ncols = 3;
else
    error('too many variables (makechart)')
end


for i = 1:nvars
    var(i).title = titlelist(i,:);
    for t  = 1:nobs
        for indi = 2 : nshocks+1
            var(i).min(t,indi-1) = min(zdata(t, i, indi),0);
            var(i).max(t,indi-1) = max(zdata(t, i, indi),0);
            var(i).total(t) = (zdata(t, i, 1));
        end
    end
end

colormap default
colormap(cmp1)

% year       = floor(xvalues);
% month      = floor(12*(xvalues - year)+1);
% smonth     = size(month,1);
% 
% NumTicks = smonth/3;
% 
% monthlist  = {'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',...
%             'Oct','Nov','Dec'};
%         
% tstep      = nobs*step/NumTicks;
% mstep      = floor(nobs/NumTicks);
% tick_d     = NumTicks+1;
% x_tick     = NaN(1,tick_d);
% x_label    = cell(1,tick_d);
% M_0        = month(1,1);
% M_0        = 12*(yearshock - floor(yearshock));
% Y_0        = year(1,1);
% w_c        = 0;

% for i = 1:tick_d;
%     t_index = i - 1; 
%     x_tick(1,i) = yearshock + tstep * t_index;
%     m_index = M_0 + mstep * t_index;
%     tw_c       = 0;
%     
%     while m_index > 12;
%         tw_c    = tw_c + 1;
%         m_index = m_index - 12;
%         w_c     = tw_c;
%     end;
%     
%     Y_index      =  floor(yearshock) + w_c;
%     x_year       =  num2str(Y_index);
%     x_year       =  x_year(1,3:4);
% 
%     % CHECK
%     if m_index==0
%     x_month      =  cellstr(monthlist(1+round(m_index)));
%     else
%     x_month      =  cellstr(monthlist(round(m_index)));
%     end
%     
%     x_l          =  strcat(x_month,x_year);
%     x_label(1,i) =  x_l;
%     
% end;

% 
map = [0 0 0.3
    0 0 0.4
    0 0 0.5
    0 0 0.6
    0 0 0.8
    0 0 1.0];
           
colormap default
% colormap(map)

for i = 1:nvars

    var(i).subplot=subplot(nrows,ncols,i);
    
    plot(xvalues,var(i).total,'k','Linewidth',2);
    
    hold on
    bmin = bar(xvalues,var(i).min,'stack','FaceColor','flat'); colormap(cmp1);
    set(bmin,'ShowBaseLine','off')
    for k = 1:3
       bmin(k).CData=k; 
    end
    colormap default; colormap(cmp1);
    bmax = bar(xvalues,var(i).max,'stack','FaceColor','flat'); colormap(cmp1);
    set(bmax,'ShowBaseLine','off')
    for k = 1:3
       bmax(k).CData=k; 
    end
    hold on
    line_p=plot(xvalues,var(i).total,'k','Linewidth',2) ;
    
%     if opts.shading==1
%     shade(1990+7/12,1991+3/12,[0.8 0.8 0.8]);
%     shade(2001+3/12,2001+11/12,[0.8 0.8 0.8]);
%     shade(2007+12/12,2009+6/12,[0.8 0.8 0.8]);
%     end    
    
    
    axis('tight')

    hold on
    
    title(var(i).subplot,[var(i).title],'FontSize',12)
    
    if i==nvars
        if numel(strvcat(legendlist(1,:)))
            h=legend(legendlist,'Location','Northeast','Orientation','Vertical');
%             legend('boxoff')
            set(h,'Fontsize',10)
        end
    end


    
    if i==nvars
        if nvars>3
            text('String',figlabel,'Units','normalized','Position',[1.2 1.21],...
                'FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
        else
            text('String',figlabel,'Units','normalized','Position',[0.4 1.24],...
                'FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
        end
    end
    
    
end

newUnits = 'normalized';
set(h,'Position',opts.legendplace,'Units',newUnits)
