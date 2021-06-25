% Set variables
% fcnHandle = @GeneTherapySystemElimination;
tspan = [0, 720];

% Read in a table of fit values
T = readtable('terms-creAndDna.csv');

% Pick which fit you want to use
fit = 'fit1';

% Get the variables 
FitTable = T(:,{'varName', [fit, '_fit'], [fit, '_std']});

% List all fit values and standard deviations for variables
vars = table2cell(FitTable);

varsToUse = zeros(1, length(vars));
for i = 1:length(vars)
    varsToUse(i) = vars{i, 2};
end

[tint, yint] = creAndDna(tspan, varsToUse);

figure;
plot(tint, yint)