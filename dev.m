% Load in the model catalog
load('model-catalog.mat')

% Load in the paramters
load('base-parameters.mat')

    
% Collect model info
modelName = models{1, MODEL_NAME};
modelFun = models{1, MODEL_FUN};
    
% Do the perturbation
results = logNormalPerturbation(...
    modelFun, parameters, {'Cas_degradation'}, 51, [0 312], initial, false);

[strrep(modelName, ' ', '_'), '-perturb-', strjoin(perturbationCombinations{1}, '-'), '.mat']