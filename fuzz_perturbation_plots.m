% fuzz_perturbation_plots.m
% Last updated: 2021-12-06, by Helen Scott
%
% REVISION HISTORY:
%   12/06/2021 - Helen Scott
%       * Initial implementation
%
%% Set-Up
% Name of the model I am interested in
modelName = 'Sequential Activator \rightarrow Cre-ON';
modelFun = @Chain_Activator_Cre_on_Kill_Switch;

% Clean up name (for making directory/files)
cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
    ' ', '_'), '/', '-');

% Set the output path
resultsPath = './fuzz-perturbation-results/';

% Set the output path for all of the plots
outpath = [resultsPath, 'plots/'];

% If directory doesn't exist, try to create it
if ~isfolder(outpath)
    disp(['Directory does not exist, attempting to make it: ', ...
        outpath]);
    mkdir(outpath);
end

% Set color bar limits
cLimits = [0 500];

% Set number of bins
nBins = [200 50];

%% Loop
% For all parameter sets
for i = [5, 10, 20]    
    % Set the parameter set name
    parameterSet = ['tweaked_parameters_', int2str(i)];
    
    % Load the parameters
    eval(parameterSet);
    
    % We want to perturb every parameter individually
    parameterNames = keys(parameters);
    
    % Make an output folder for this sets of tweaked parameters
    folderName = [resultsPath, parameterSet, '/'];
    
    % Load the perturbation results
    load([folderName, 'random-perturb-all.mat'])
        
    % Plot
    densityheatmap(x, ys, nBins, cLimits, i, outpath);
end