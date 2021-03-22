% Set seed of rnadom number generator 
rng('default');

% List all fit values for variables
vars = ["lag" 44.2974; "Tc" 2.6514; "T1" 0.6807]; % Change to cell array

% Which variable do you want to varry
perturbedVar = "Tc";

% Create a random value (from the log normal distribution) for the variable
% to varry, use the means for the others
whichVars = perturbedVar == vars(:, 1);

% How many times do you want to run it
nRuns = 1000;

% Create vector to hold pertubed values.
perturbedValues = zeros(1, nRuns);

% Create vector to hold results.
% make variable
results = zeros(14, 2, nRuns);
% results = zeros(336, 2, nRuns);

% Header for ticker
disp(datetime('now'));
if nRuns > 10
    disp('PERTURBATIONS RUN:')
end
% Pass to the single perturbation that many times
for i = 1:nRuns
    % Ticker for progress, every 10 loops
    if mod(i, 10) == 0
        fprintf('%d ', i); 
    end
    % Carriage return at the end
    if i == nRuns
        fprintf('\n');
    end
    
    % Get a new value for the perturbed var
    varsToUse = zeros(1, length(vars));
    for k = 1:length(whichVars)
        if whichVars(k)
            perturbedValue = lognrnd(log(str2double(vars(k, 2))), log(2)); % Change to 10^(randn()*log10(stddev)+log10(mean))
            varsToUse(k) = perturbedValue;
        else
            varsToUse(k) = str2double(vars(k, 2));
        end
    end
    
    % Add the perturbed variable value to a vector
    % Run the model with the perturbed variable, add to a new dimension of
    % an array
    perturbedValues(i) = perturbedValue;
    results(:,:,i) = perturbableFitCas9AssumeDelay(varsToUse(1), varsToUse(2), varsToUse(3));
end

save('20210311-lag03.mat', 'perturbedValues', 'results');
disp(datetime('now'));