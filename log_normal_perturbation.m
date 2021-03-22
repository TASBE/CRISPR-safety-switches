% List all fit values and standard deviations for variables
vars = {"lag" 44.2974 2; "Tc" 2.6514 2; "T1" 0.6807 2};

% Which variable do you want to vary
perturbedVars = ["Tc", "T1"];

% How many times do you want to run it
nRuns = 100;

% Set time span
tspan = 336;

% Boolean vector, True if the variable is being varied 
whichVars = zeros(1, length(vars));
for i = 1:length(perturbedVars)
    for j = 1:length(vars)
        if whichVars(j) == 0
            whichVars(j) = vars{j} == perturbedVars{i};
        end
    end
end

% Create vector to hold results.
results = cell(nRuns, 3);

% Header for ticker
tic
if nRuns > 10
    disp('PERTURBATIONS RUN:')
end
% Pass to the single perturbation that many times
for i = 1:nRuns
    % Ticker for progress, every 10 loops
    if mod(i, 10) == 0
        fprintf('%d ', i); 
    end
    
    % Define which percentiles I want for the perturbations
    % Lower and uperlimit set to capture center 95% (2 standard deviations)
    percentiles = linspace(0.023, 0.977, nRuns);
    
    % Normal distribution values at percentiles
    pd = makedist('Normal');
    normalDist = zeros(size(percentiles));
    for j = 1:length(percentiles)
        normalDist(j) = icdf(pd, percentiles(j));
    end
    
    % Generate all variable perturbations
    varsToUse = zeros(nRuns, length(vars));
    for k = 1:length(whichVars)
        if whichVars(k)
            mean = vars{k, 2};
            stddev = vars{k, 3};
            logDist = zeros(size(percentiles));
            for l = 1:length(percentiles)
                logDist(l) = 10^(normalDist(l)*log10(stddev)+log10(mean));
            end
            varsToUse(:, k) = logDist;
        else
            varsToUse(:, k) = vars{k, 2};
        end
    end
    
    % Save the percentiles to the results array
    results(i, 1) = {percentiles(i)};
    % Save the variables to the results array
    results(i, 2) = {varsToUse(i, :)};
    % Run the model, save to results array
    results(i, 3) = {perturbableFitCas9AssumeDelay([0, tspan], varsToUse(i, 1), varsToUse(i, 2), varsToUse(i, 3))};
    
    % Carriage return at the end
    if i == nRuns
        fprintf('\n');
    end
end
toc

save('20210322-TcT1.mat', 'results');