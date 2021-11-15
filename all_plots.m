% all_plots.m
% Last updated: 2021-11-14, by Helen Scott
%
% REVISION HISTORY:
%   11/14/2021 - Helen Scott
%       * Initial implementation
%       * Goal: Have a one-click method to (re)generate info for all plots
%       * Want the plotting to be separrate from the actual perturbation

% Load in the model catalog
load('model-catalog.mat')

% Load in the paramters % Do I need this for the plotting?
load('base-parameters.mat')

% For all models
for i = 1:n_models
    
    % Collect model info
    modelName = models{i, MODEL_NAME};
    modelFun = models{i, MODEL_FUN};
    
    % Load the perturbation results
    
    
    % Generatee XXX plot
    % Generate XXX plot
end
