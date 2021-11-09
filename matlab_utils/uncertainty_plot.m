function uncertainty_plot(fcnHandle, tspan, varsArray, units)
    
    % Extract values of variables into an array
    vars = zeros(1, length(varsArray));
    for i = 1:length(varsArray)
        vars(i) = varsArray{i, 2};
    end
    
    % Solve for all values being normal
    [tint, yint] = fcnHandle(tspan, vars);
    
    figure; 
    % Plot normal Array
    plot(tint, yint, 'Color', 'black', 'LineWidth', 2);
    hold on;
    
    % Get unique colors for the number of variables
    jetcustom = hsv(length(vars));
%     custom = ["#e6194B", "#3cb44b", "#ffe119", "#4363d8", "#f58231", ...
%         "#42d4f4", "#f032e6", "#fabed4", "#469990", "#dcbeff", ...
%         "#9A6324", "#fffac8", "#800000", "#aaffc3", "#000075"];
    
    h = zeros(1, length(vars) * 3);
    % For each variable
    for i = 1:length(vars)
        % For tracking
        fprintf('Currently on variable %i/%i\n', i, length(vars))
        % Extract variable name
        varName = varsArray{i};
        % Plot normal value
        h(((i-1)*3)+1) = plot(tint, yint, 'Color', jetcustom(i, :), 'DisplayName', varName);
        hold on;
        % Calculate -20%
        varsToUse = vars;
        varsToUse(i) = varsToUse(i) * 0.8;
        % Get plotting values
        [tint, ylow] = fcnHandle(tspan, varsToUse);
        % Plot
        h(((i-1)*3)+2) = plot(tint, ylow, '--', 'Color', jetcustom(i, :));
        % Calculate +20%
        varsToUse = vars;
        varsToUse(i) = varsToUse(i) * 1.2;
        % Get plotting values
        [tint, yhigh] = fcnHandle(tspan, varsToUse);
        % Plot
        h(((i-1)*3)+3) = plot(tint, yhigh, '-.', 'Color', jetcustom(i, :));
    end
    
    hold off;
    % Style
    title(sprintf('Uncertainty Analysis for %s Circuit', func2str(fcnHandle)));
    grid on;
    xlabel('Time (Hours)');
    ylabel(units);
    
    legend(h(1:3:end));
    lgd = legend;
    lgd.FontSize = 12;
    lgd.Title.String = 'Variables';
    lgd.Location = 'northeastoutside';
    % Save plot
    saveas(gcf, [func2str(fcnHandle) '_uncertainty.eps'], 'epsc');
    saveas(gcf, [func2str(fcnHandle) '_uncertainty.png']);
    saveas(gcf, [func2str(fcnHandle) '_uncertainty.fig']);