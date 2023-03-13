% uncomment a location below -- adapt to local installation

%   location = 'home_windows';
% % location = 'home_mac';

restoredefaultpath

if strmatch(location,'home_windows','exact')    
%     dir1='C:\dynare\4.3.3\matlab';
    dir1='C:\dynare\4.3.1\matlab';
    dir2='.\toolkit_files';
    dir3='.\estimation_functions';
elseif strmatch(location,'home_mac','exact')
   dir1='/Applications/dynare/4.3.1/matlab';
    dir2='./toolkit_files';
    dir3='./estimation_functions';
    %% 

else 
    error('Specify path to Dynare installation')
end

path(dir1,path);
path(dir2,path);
 path(dir3,path);


dynare_config

