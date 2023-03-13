%% data
datapath = './data/';

%% Data (year on year)
datafilename = 'data_jpn_1_def';         % actual data
% datafilename = 'data_jpn_1';         % old version

yy1 = csvread(strcat(datapath, [ datafilename '.csv' ] ), 1, 1);

yy = [ yy1(:,2) yy1(:,3) yy1(:,4)];

% yy = [ yy1(:,2) yy1(:,3) yy1(:,4)];

% yy = [ 100*yy1(:,1) 100*yy1(:,3) yy1(:,4)]; % old version

ZZ = [100; 100; 100]; % coefficients between observables and endogenous variables of the measurement equation

% ZZ = [100; 100; 400];  % old version

%%
Tobs = size(yy,1);