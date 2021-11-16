function results = logNormalPerturbation(...
        fcnHandle, vars, perturbedVars, nRuns, tspan, initial, random)
    % FUNCTION NAME:
    %   logNormalPerturbation
    %
    % DESCRIPTION:
    %   Generate a series of values for one or more variables perturbed on
    %   a log normal scale. The variables are then passes to a specified
    %   model differential equation.
    %
    % INPUT:
    %   fcnHandle - (function handle) Function handle for differential 
    %       equation function to be used
    %   vars - (map) Map object listing the parameter names and values as 
    %       determined from parameter fitting and literature
    %   perturbedVars - (cell array) Name(s) of variables to be perturbed
    %   nRuns - (double) Number of perturbations to run
    %   tspan - (double) Time span to model ([start, end])
    %   initial - (map) Map object listing initial value names and values
    %   random - (logical) Is the perturbation to be random on the log 
    %       normal scale?
    %
    % OUTPUT:
    %   results - (cell) Column 1 = percentile of the perturbed value in
    %       the log normal distribution (NaN for random distributions), 
    %       column 2 = variables used for the differential equation, column
    %       3 = results from the model function
    %
    % ASSUMPTIONS AND LIMITATIONS:
    %   Currently asusming std devs for all parameters are 2
    %
    % REVISION HISTORY:
    %   11/14/2021 - Helen Scott
    %       * Change to use with new models
    %           * Got rid of boolean vector (whichVars) because it won't 
    %               work with the parameter map
    %   04/08/2021 - Helen Scott
    %       * Fix for loop (Was recreating varsToUse every loop)
    %
    
    %% Set-up
    % Header for ticker
    tic
    if nRuns > 10
        disp('PERTURBATIONS RUN:')
    end

    % Create cell array to hold results
    results = cell(nRuns, 4);
    
    %% Generate Variables to Use in the DE
    % If making a regularly spaced log normal distribution, find
    % percentiles of a normal distribution
    if ~random
        % Define which percentiles I want for the perturbations
        % Lower and uperlimit set to capture center 95%
        percentiles = linspace(0.023, 0.977, nRuns);

        % Normal distribution values at percentiles
        pd = makedist('Normal');
        normalDist = zeros(size(percentiles));
        for j = 1:length(percentiles)
            normalDist(j) = icdf(pd, percentiles(j));
        end
        
    end
    
    % For every percentile
    for l = 1:length(normalDist)
        % Copy the parameter map
        varsToUse = containers.Map(vars.keys, vars.values);
        % Assume that all parameters have a standard deviation of 2
        stddev = 2;
        
        % For each perturbed var
        for k = 1:length(perturbedVars)
            % Generate value to use
            perturbedValue = 10^(normalDist(l) * ...
                log10(stddev) + log10(vars(perturbedVars{k})));
            varsToUse(perturbedVars{k}) = perturbedValue;
            results{l, 2} = perturbedValue; % Probelematic if we have more than 1 perturbed variable...
        end
        
        % Save the percentile and the perturbed vars to the results array
        if random
            results{l, 1} ='NaN';
        else
            results{l, 1} = percentiles(l);
        end
        results{l, 3} = varsToUse;
        
    end
    
    %% Run the DE model
    % Pass to the single perturbation that many times
    for i = 1:nRuns
        % Ticker for progress, every 10 loops
        if mod(i, 10) == 0
            fprintf('%d ', i); 
        end
        
        % Run the model, save to results array
        [x, y_out, y] = fcnHandle(tspan, results{i, 3}, initial);
        results{i, 4} = {x, y_out, y};

        % Carriage return at the end
        if i == nRuns
            fprintf('\n');
        end
        
    end
    toc