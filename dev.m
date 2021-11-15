% Load results
load('Activator-perturb-alpha_p_Cas9-alpha_p_TF-alpha_p_Cre-alpha_p_TF2.mat')

% Hardocde for now
variable = {'alpha_p_Cas9', 'alpha_p_TF', 'alpha_p_Cre', 'alpha_p_TF2'};
variablePos = NaN;
tspan = [0 312]; % Remove from plotting function
name = 'Activator-perturb-alpha';
units = 'AAV (nM)';

% Call plotting function
rainbowplot(results, variable, variablePos, tspan, name, units)

%% Debug the variable perturbation
% load('model-catalog.mat')
% load('base-parameters.mat')
% 
% % Set the initial
% initial = containers.Map();
% initial('AAV') = 50; % Dash used 50 nM
% initial('genome') = 10; % Dash used 10 nM
% initial('Cre_regulated_region') = 50; % Assume same as AAV
% initial('CreH_regulated_region') = 50;
% 
% % Set the time span
% tspan = [0 312];
% 
% % Set the number of runs to do
% % Note: An odd number gives you a 50% percentile
% nRuns = 51; 
% 
% % List all of the perturbation sets
% perturbationCombinations = {...
%     {'alpha_p_Cas9', 'alpha_p_TF', 'alpha_p_Cre', 'alpha_p_TF2'},
%     {'delta_Cas9', 'delta_TF', 'delta_TF2', 'delta_Cre', 'delta_CreH'},
%     {'alpha_r_sgRNA2', 'alpha_r_sgRNA1'},
%     {'delta_g'},
%     {'Cas_degradation'},
%     {'Cas_gRNA_binding'},
%     {'k_cat'},
%     {'K_A'},
%     {'K_R'},
%     {'k_cre'}
%     };
% % Set to look at Cas_degradation
% j = 5;
%     
% % Collect model info
% % Set to look at no-delay
% i = 1;
% modelName = models{i, MODEL_NAME};
% modelFun = models{i, MODEL_FUN};
% 
% results = logNormalPerturbation(modelFun, parameters, ...
%             perturbationCombinations{j}, nRuns, tspan, initial, false);