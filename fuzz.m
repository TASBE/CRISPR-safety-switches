function fuzz(fcnHandle, tspan, varsArray, units, verbose, save)
    % Default path is current path
    if nargin < 5
       verbose = true;
        save = true;
    end 
    
    % Extract values of variables into an array
    vars = zeros(1, length(varsArray));
    for i = 1:length(varsArray)
        vars(i) = varsArray{i, 2};
    end
    
    % Solve for all values being normal
    [tint, yint] = fcnHandle(tspan, vars);

    % For each variable
    for i = 1:length(vars)
        % For tracking
        if verbose
            fprintf('Currently on variable %i/%i\n', i, length(vars));
        end
        % Extract variable name
        varName = varsArray{i};
        % Plot normal value
        figure(i);
        plot(tint, yint, 'Color', 'black');
        hold on;
        % Calculate +20%
        varsToUse = vars;
        varsToUse(i) = varsToUse(i) * 0.8;
        % Get plotting values
        [tint, ylow] = fcnHandle(tspan, varsToUse);
        % Plot
        plot(tint, ylow, 'Color', 'blue');
        % Calculate -20%
        varsToUse = vars;
        varsToUse(i) = varsToUse(i) * 1.2;
        % Get plotting values
        [tint, yhigh] = fcnHandle(tspan, varsToUse);
        % Plot
        plot(tint, yhigh, 'Color', 'red');
        % Style
        title(sprintf('Uncertainty Analysis for %s', varName));
        grid on;
        xlabel('Time (Hours)');
        ylabel('Cas9 Plasmid (nM)');
        % Save plot
        if save
            saveas(gcf, [varName '.eps'], 'epsc');
            saveas(gcf, [varName '.png'], 'png');
            saveas(gcf, [varName '.fig']);
        end
    end