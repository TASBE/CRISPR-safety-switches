function [Dist, results] = randomperturbation(modelName, vars, perturbedVars, nRuns, tspan)

    % Boolean vector, True if the variable is being varied 
    whichVars = zeros(1, length(vars));
    for i = 1:length(perturbedVars)
        for j = 1:length(vars)
            if whichVars(j) == 0
                whichVars(j) = (vars{j,1} == perturbedVars(i));
            end
        end
    end

    % Create vector to hold results.
    results = cell(nRuns, 2);

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

        % Generate all variable perturbations
        varsToUse = zeros(nRuns, length(vars));
        for k = 1:length(whichVars)
            if whichVars(k)
                mean = vars{k, 2};
                stddev = vars{k, 3};
                Dist = zeros(1, nRuns);
                for l = 1:nRuns
                    Dist(l) = 10^(randn()*log10(stddev)+log10(mean));
                end
                varsToUse(:, k) = Dist;
            else
                varsToUse(:, k) = vars{k, 2};
            end
        end

        % Save the variables to the results array
        results(i, 1) = {varsToUse(i, :)};
        % Run the model, save to results array
        results(i, 2) = {feval(modelName, [0, tspan], varsToUse)};

        % Carriage return at the end
        if i == nRuns
            fprintf('\n');
        end
    end
    toc