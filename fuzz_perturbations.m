% fuzz_perturbations
% Last updated: 2021-12-06, by Helen Scott
%
% REVISION HISTORY:
%   12/06/2021 - Helen Scott
%       * Initial implementation
%       * Goal: One push to generate results of 1000 random perturbations 
%           of tweaked parameters for the good circuit 
%       * This does NOT use the old "fuzz" function in the matlab util
%           folder, but will use the densityheatmap.m for plotting

%% Set-Up
% Name of the model I am interested in
modelName = 'Sequential Activator \rightarrow Cre-ON';
modelFun = @Chain_Activator_Cre_on_Kill_Switch;

% Clean up name (for making directory/files)
cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
    ' ', '_'), '/', '-');

% Set the initial
initial = containers.Map();
initial('AAV') = 10;
initial('Cre_regulated_region') = initial('AAV'); % if present, starts unmodified
initial('genome') = 1;

% Set the time span
tspan = [0 720];

% Set the number of runs to do
nRuns = 10000; 

% Set the output path
resultsPath = './fuzz-perturbation-results/';

%% Loop
% For all parameter sets
for i = [5, 10, 20]
    % Set the parameter name
    parameterSet = ['tweaked_parameters_', int2str(i)];
    
    % Load the parameters
    eval(parameterSet);
    
    % We want to perturb every parameter individually
    parameterNames = keys(parameters);
    
    % Make an output folder for this sets of tweaked parameters
    outpath = [resultsPath, parameterSet, '/'];
    
    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        disp(['Directory does not exist, attempting to make it: ', ...
            outpath]);
        mkdir(outpath);
    end
    
    % Print to screen
    disp(['Currently on: ', parameterSet]);
    
        % Do the perturbation
    results = logNormalPerturbation(modelFun, parameters, ...
        parameterNames, nRuns, tspan, initial, true);

    % Extract just what is needed for plotting
    percentiles = results{:, 1};
    perturbedValues = results{:, 2};
    x = results{1, 4}{1};
    ys = cell(length(results), 1);
    for k = 1:length(results)
        ys{k} = results{k, 4}{2};
    end

    % Save results as a .mat file
    save([outpath, 'random-perturb-all.mat'], 'x', 'ys')
end