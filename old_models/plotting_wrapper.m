%% Load results
load('20210411-test.mat');

%% Set variables
% outpath = 'a/b/c/'; % Must end in /
experimentName = '210411 test';

%% Plot
if random
    % Set color bar limits
    cLimits = [0 inf];

    % Set number of bins
    nBins = [336 100];
    
    % Plot
    densityheatmap(results, nBins, cLimits, experimentName);
else
    variablePos = find(whichVars);
    rainbowplot(results, perturbedVars, variablePos, tspan, experimentName)
end
