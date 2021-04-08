%% Load results
load('dev-results/20210322-lag.mat');

%% Set variables
outpath = 'a/b/c/'; % Must end in /
experimentName = '210406 test';

%% Plot
if random
    % Set color bar limits
    cLimits = [0 inf];

    % Set number of bins
    nBins = [336 100];
    
    % Plot
    densityheatmap(results, nBins, cLimits, experimentName, outpath);
else
    variablePos = find(whichVars);
    rainbowplot(results, perturbedVars, variablePos, tspan, experimentName)
end
