% Set seed for random number generator for reproducibility in testing
rng('default')

%% Set variables for run
% Select model to use (use exact name of function)
fcnHandle = @linear;

% Read in a table of fit values
T = readtable('terms-linear.csv');

% Pick which fit you want to use
fit = 'fit1';

% Get the variables 
FitTable = T(:,{'varName', [fit, '_fit'], [fit, '_std']});

% List all fit values and standard deviations for variables
vars = table2cell(FitTable);

% Which variable do you want to vary
perturbedVars = ["Tc"];

% How many times do you want to run it
nRuns = 3;

% Set time span
tspan = [0, 720];

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
[whichVars, results] = logNormalPerturbation(fcnHandle, vars, perturbedVars, nRuns, tspan, random);

%% Save results
save('20210411-test.mat', 'results', 'random', 'perturbedVars',...
    'tspan', 'whichVars');

