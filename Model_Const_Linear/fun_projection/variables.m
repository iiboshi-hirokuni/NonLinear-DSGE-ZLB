function V = variables

% V = variables
%   Names variables and assigns locations
% Output:
%   V : Structure of variable names and locations

% Variables names
V.names = { 'y','y_star','pi','r','r_lg','r_star','mu_a','z_b', ...     %%%%%%%%%%  Modified on Jan 30, 2017 by UEDA
            'E_y','E_pi','E_y_star' } ;

% Variables titles       
V.desc = {  'Output (%)',
            'Potential',
            'Inflation (% point)',
            'Nominal Interest Rate',
            'lag Nominal Interest Rate',
            'Real Rate',            
            'Productivity (%)',
            'Preference'
         };

%Shocktypes: 
V.shocktypes = {'a','b','MP'}; 

% Forecast errors
V.foretypes = {'y','pi'}; 
             
% Number of variables
V.nvar = length(V.names);
% Number of shocks
V.nshock = length(V.shocktypes);
% Number of forecast errors
V.nfore = length(V.foretypes);

% Establish variable index
for j = 1:size(V.names,2)
   eval(['V.' V.names{j} ' = j;']);      
end
% Establish shock index
for j = 1:V.nshock
    eval(['V.eps_' V.shocktypes{j} ' = j;'])
end
% Establish forecast error index
for j = 1:V.nfore
    eval(['V.fe' V.foretypes{j} ' = j;'])
end