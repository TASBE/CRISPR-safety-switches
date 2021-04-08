% Set seed for random number generator for reproducibility in testing
rng('default')

%% Set variables for run
% Select model to use (use exact name of function)
fcnHandle = @GeneTherapySystemElimination;

% Read in a table of fit values
T = readtable('terms-GeneTherapySystemElimination.csv');

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
tspan = [0, 336];

% Pick perturbation type
random = false;

% % Set path for output
% outpath = '20210407-test/';
% 
% %% Make outptut folder
% if ~isfolder(outpath)
%     mkdir(outpath);
% end

%% Pass to perturbation script
[whichVars, Dist, results] = logNormalPerturbation(fcnHandle, vars, perturbedVars, nRuns, tspan, random);

%% Save results
save('202104046-lag.mat', 'Dist', 'results', 'random', 'perturbedVars',...
    'tspan', 'whichVars');

