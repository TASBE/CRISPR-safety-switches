% Set seed for random number generator for reproducibility in testing
rng('default')

% Select model to use (use exact name of function)
modelName = 'perturbableFitCas9AssumeDelay';

% Read in a table of fit values
T = readtable(['terms-', modelName, '.csv']);

% Pick which fit you want to use
fit = 'name';

% Get the variables 
FitTable = T(:,{'varName', [fit, '_fit'], [fit, '_std']});

% List all fit values and standard deviations for variables
vars = table2cell(FitTable);

% Which variable do you want to vary
perturbedVars = ["lag"];

% How many times do you want to run it
nRuns = 10;

% Set time span
tspan = 336;

% Pick perturbation type
perturbType = 'random';

% Pass to perturbation script
[Dist, results] = feval([perturbType, 'perturbation'], modelName, vars, perturbedVars, nRuns, tspan);

% Save results
save('20210404-rand-lag.mat', 'Dist', 'results');

