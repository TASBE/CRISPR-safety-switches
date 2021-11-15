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

% Load in the paramters
load('base-parameters.mat')

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

% List all of the perturbation sets
perturbationCombinations = {...
    {'alpha_p_Cas9', 'alpha_p_TF', 'alpha_p_Cre', 'alpha_p_TF2'},
    {'delta_Cas9', 'delta_TF', 'delta_TF2', 'delta_Cre', 'delta_CreH'},
    {'alpha_r_sgRNA2', 'alpha_r_sgRNA1'},
    {'delta_g'},
    {'Cas_degradation'},
    {'Cas_gRNA_binding'},
    {'k_cat'},
    {'K_A'},
    {'K_R'},
    {'k_cre'}
    };

%% Loop
% For all models
for i = 31:n_models % CHANGE BACK TO 1!
    
    % Collect model info
    modelName = models{i, MODEL_NAME};
    modelFun = models{i, MODEL_FUN};
    
    % Print to screen
    disp(['Currently on: ', modelName]);
    
    % For every perturbation combination listed
    for j = 1:length(perturbationCombinations)
        
        % Print to screen
        disp(['Perturbed Parameter(s): ', perturbationCombinations{j}]);
        
        % Do the perturbation
        results = logNormalPerturbation(modelFun, parameters, ...
            perturbationCombinations{j}, nRuns, tspan, initial, false);

        % Save results as a .mat file
        save([strrep(strrep(strrep(modelName, '\rightarrow ', ''), ' ', '_'), '/', '-'), ...
            '-perturb-', strjoin(perturbationCombinations{j}, '-'), ...
            '.mat'], 'results')
    end
    
end
