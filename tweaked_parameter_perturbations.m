% tweaked_parameter_perturbations.m
% Last updated: 2021-12-06, by Helen Scott
%
% REVISION HISTORY:
%   12/06/2021 - Helen Scott
%       * Initial implementation, copied all_base_parameter_perturbations.m
%

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
% Note: An odd number gives you a 50% percentile
nRuns = 101; 

% Set the output path
resultsPath = './tweaked-parameter-perturbation-results/';

%% Loop
% For all models
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
        fprintf('Directory does not exist, attempting to make it: %s', ...
            outpath);
        mkdir(outpath);
    end
    
    % Print to screen
    disp(['Currently on: ', parameterSet]);
    
    % For every perturbation combination listed
    for j = 1:length(parameterNames)
        
        % TODO: Skip the irrelevant variables
        
        % Print to screen
        disp(['Perturbed Parameter(s): ', parameterNames{j}]);
        
        % Do the perturbation
        results = logNormalPerturbation(modelFun, parameters, ...
            string(parameterNames{j}), nRuns, tspan, initial, false);
        
        % Extract just what is needed for plotting
        percentiles = results(:, 1);
        perturbedValues = results(:, 2);
        x = results{1, 4}{1};
        ys = cell(length(results), 1);
        for k = 1:length(results)
            ys{k} = results{k, 4}{2};
        end
        
        % Save results as a .mat file
        save([outpath, cleanModelName, '-perturb-', parameterNames{j}, ...
            '.mat'], 'percentiles', 'perturbedValues', 'x', 'ys')
    end
    
end
