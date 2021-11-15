% all_perturbations.m
% Last updated: 2021-11-14, by Helen Scott
%
% REVISION HISTORY:
%   11/15/2021
%       * 
%   11/14/2021 - Helen Scott
%       * Initial implementation
%       * Goal: Have a one-click method to (re)generate info for all plots
%       * I wanted the perturbation part to be seperate from the plotting
%       * Not sure which perturbation combinations to do
%       * Not sure about the log normal perturbation calculation now
%

%% Set-Up
% Load in the model catalog
load('model-catalog.mat')

% Select a subset of models
% Have to overwrite n_models and models

% Load in the paramters
load('base-parameters.mat')
% We want to perturb every parameter individually
parameterNames = keys(parameters);

% Set the initial
initial = containers.Map();
initial('AAV') = 50; % Dash used 50 nM
initial('genome') = 10; % Dash used 10 nM
initial('Cre_regulated_region') = 50; % Assume same as AAV
initial('CreH_regulated_region') = 50;

% Set the time span
tspan = [0 312];

% Set the number of runs to do
% Note: An odd number gives you a 50% percentile
nRuns = 51; 

% Set the output path
resultsPath = './results/';

%% Loop
% For all models
for i = 1:n_models
    
    % Collect model info
    modelName = models{i, MODEL_NAME};
    modelFun = models{i, MODEL_FUN};
    
    % Clean up name (for making directory/files)
    cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
        ' ', '_'), '/', '-');
    
    % Make an output folder for this circuit model
    outpath = [resultsPath, cleanModelName, '/'];
    
    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        fprintf('Directory does not exist, attempting to make it: %s', ...
            outpath);
        mkdir(outpath);
    end
    
    % Print to screen
    disp(['Currently on: ', modelName]);
    
    % For every perturbation combination listed
    for j = 1:length(parameterNames)
        
        % TODO: Skip the irrelevant variables
        
        % Print to screen
        disp(['Perturbed Parameter(s): ', parameterNames{j}]);
        
        % Do the perturbation
        results = logNormalPerturbation(modelFun, parameters, ...
            perturbationCombinations{j}, nRuns, tspan, initial, false);
        
        % Extract just what is needed for plotting
        plotData = {};
        
        % Save results as a .mat file
        save([outpath, cleanModelName, '-perturb-', parameterNames{j}, ...
            '.mat'], 'plotData')
    end
    
end
