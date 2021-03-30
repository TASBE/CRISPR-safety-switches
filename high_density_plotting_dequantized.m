% Load results
load('20210322-TcT1.mat');

nRuns = length(results);
variable = ["Tc", "T1"];
variablePos = [2, 3]; % GET RID OF MAGIC NUMBERS!
stddev = 2;

% Generate colors for all runs
jetcustom = jet(length(results));

center = length(results)/2;

if length(variable) == 1
    % Variable values for each 10 percentile
    % FIND A CLEANER WAY TO DO THIS
    percentiles = [10:10:100];
    c = zeros(1, length(percentiles));
    for i = 1:length(percentiles)
        c(i) = results{percentiles(i),2}(variablePos);
    end

    % Plot
    figure;
    for i = 1:length(results)
        sol = results{i, 3};
        tint = 1:336;
        yint = zeros(size(tint));
        for j = tint
            ddeRes = deval(sol, j);
            yint(j) = 100* ddeRes(5) / (ddeRes(4) + ddeRes(5));
        end
        if i == center
            plot(tint, yint, 'Color', 'black', 'LineWidth', 2);
        else
            plot(tint, yint, 'Color',  jetcustom(i, :));
        end
        hold on;
    end

    colormap(jet(length(results)));
%     cb = colorbar;
    caxis([results{1,1} results{end,1}]);
    cb = colorbar('YTickLabel',c);
    ylabel(cb, sprintf('Perturbed Value of %s', variable));
    title(sprintf('Perturbation of %s', variable));
else
    % Plot
    figure;
    for i = 1:length(results)
        sol = results{i, 3};
        tint = 1:336;
        yint = zeros(size(tint));
        for j = tint
            ddeRes = deval(sol, j);
            yint(j) = 100* ddeRes(5) / (ddeRes(4) + ddeRes(5));
        end
        if i == center
            plot(tint, yint, 'Color', 'black', 'LineWidth', 2);
        else
            plot(tint, yint, 'Color',  jetcustom(i, :));
        end
        hold on;
    end
    colormap(jet(length(results)));
    cb = colorbar;
    caxis([results{1, 1} results{end, 1}]);
    ylabel(cb,'Percentile of Perturbed Varriable(s)');
    title(sprintf('Perturbation of %s and %s together', variable(1), variable(2))); % Technically only good for two variables
end

% General style
grid on;
xlabel('Time (Hours)');
ylabel('Expression (%)');

% Save file
if length(variable) == 1
    saveas(gcf, variable + '-RainbowPlot.png');
else
    saveas(gcf, variable(1) + '-' + variable(2) + '-RainbowPlot.png');
end