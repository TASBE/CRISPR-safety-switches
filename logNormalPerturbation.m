function [whichVars, results] = logNormalPerturbation(...
        fcnHandle, vars, perturbedVars, nRuns, tspan, random)
    
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
    %   vars - (cell array) Cell array holding all variables needed for
    %       differential equation, column 1 = variable 1, column 2 = 
    %       fit value (used as mean of the log normal distribution), 
    %       column 3 = standard deviation of fit value
    %   perturbedVars - (string array) Name(s) of variables to be perturbed
    %   nRuns - (double) Number of perturbations to run
    %   tspan - (double) Time span to model ([start, end])
    %   random - (logical) Is the pertuebation to be random on the log 
    %       normal scale?
    %
    % OUTPUT:
    %   whichVars - (logical) Logical if the variable in that row of vars
    %       is being varied in the perturbation
    %   Dist - (double) List of perturbed values of variable of interest 
    %   results - (cell) Column 1 = percentile of the perturbed value in
    %       the log normal distribution (NaN for random distributions), 
    %       column 2 = varaibles used for the differential equation, column
    %       3 = results from the differential equation (sol)
    %
    % ASSUMPTIONS AND LIMITATIONS:
    %   Order of variables in the vars array, must match the order the 
    %   variables are "unpacked" in the DE function.
    %
    % REVISION HISTORY:
    %   04/08/2021 - hscott
    %       * Fix for loop (Was recreating varsToUse every loop)
    %
    
    %% Set-up
    % Header for ticker
    tic
    if nRuns > 10
        disp('PERTURBATIONS RUN:')
    end

    % Create vector to hold results
    results = cell(nRuns, 3);
    
    %% Generate Variables to Use in the DE
    % Boolean vector, True if the variable is being varied 
    whichVars = zeros(1, length(vars));
    for i = 1:length(perturbedVars)
        for j = 1:length(vars)
            if whichVars(j) == 0
                whichVars(j) = (vars{j,1} == perturbedVars(i));
            end
        end
    end
    
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
    
    % Generate variables to use
    varsToUse = zeros(nRuns, length(vars));
    for k = 1:length(whichVars)
        if whichVars(k)
            mean = vars{k, 2};
            stddev = vars{k, 3};
            Dist = zeros(1, nRuns);
            for l = 1:nRuns
                if random
                    Dist(l) = 10^(randn()*log10(stddev)+log10(mean));
                else
                    Dist(l) = 10^(normalDist(l)*log10(stddev)+log10(mean));
                end
            end
            varsToUse(:, k) = Dist;
        else
            varsToUse(:, k) = vars{k, 2};
        end
    end
    
    %% Run the DE model
    % Pass to the single perturbation that many times
    for i = 1:nRuns
        % Ticker for progress, every 10 loops
        if mod(i, 10) == 0
            fprintf('%d ', i); 
        end
        
        % Save percentile values
        if random
            % Random values don't have percentiles- save as NaN
            results(i,1) = {NaN};
        else
            % Save the percentiles to the results array
            results(i,1) = {percentiles(i)};
        end
        
        % Save the variables to the results array
        results(i, 2) = {varsToUse(i, :)};
        
        % Run the model, save to results array
        results(i, 3) = {fcnHandle(tspan, varsToUse(i, :))};

        % Carriage return at the end
        if i == nRuns
            fprintf('\n');
        end
    end
    toc