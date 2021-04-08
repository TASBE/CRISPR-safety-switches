function rainbowplot(results, variable, variablePos, tspan, fcnHandle, name, path)
    %% Make outout path path
    % Default path is current path
    if nargin < 7
        path = './plots/';  % Note: frontslash works for both Windows and Mac/Unix
    else
        path = [path, 'plots/'];
    end 
    
    disp(path)
    
    % If directory doesn't exist, try to create it
    if ~isfolder(path)
        sanitized_path = strrep(path, '/', '&#47;');
        sanitized_path = strrep(sanitized_path, '\', '&#92;');
        sanitized_path = strrep(sanitized_path, ':', '&#58;');
        fprintf('Directory does not exist, attempting to create it: %s',sanitized_path);
        mkdir(path);
    end
    
    %% Plot
    % Generate colors for all runs
    jetcustom = jet(length(results));
    % Define center (mean value for the distribution)
    center = length(results)/2;

    % Plot
    figure;
    for i = 1:length(results)
        sol = results{i, 3};
        tint = 1:tspan(end);
        y = deval(sol, tint);
        yint = fcnHandle(y);
        if i == center
            plot(tint, yint, 'Color', 'black', 'LineWidth', 2);
        else
            plot(tint, yint, 'Color',  jetcustom(i, :));
        end
        hold on;
    end
    colormap(jet(length(results)));
    
    if length(variable) == 1
        % Variable values for each 10 percentile
        % FIND A CLEANER WAY TO DO THIS
        percentiles = [10:10:100];
        c = zeros(1, length(percentiles));
        for i = 1:length(percentiles)
            c(i) = results{percentiles(i),2}(variablePos);
        end
        % Color bar for actual variable values
        caxis([results{1,1} results{end,1}]);
        cb = colorbar('YTickLabel',c);
        ylabel(cb, sprintf('Perturbed Value of %s', variable));
        title(sprintf('Perturbation of %s', variable));
    else
        % Color bar for percentiles
        cb = colorbar;
        caxis([results{1, 1} results{end, 1}]);
        ylabel(cb,'Percentile of Perturbed Varriable(s)');
        title(sprintf('Perturbation of %s and %s together', variable(1), variable(2))); % Technically only good for two variables
    end
    grid on;
    xlabel('Time (Hours)');
    ylabel('Cas9 Plasmid (nM)');

    %%  Save figure
    % Make a clean name
    cleanName = sanitize_filename(name);
    % Output
    saveas(gcf, [path cleanName '.eps'], 'epsc');
    saveas(gcf, [path cleanName '.png'], 'png');
    saveas(gcf, [path cleanName '.fig']);
end