% all_plots.m
% Last updated: 2021-12-03, by Helen Scott
%
% REVISION HISTORY:
%   12/03/2021 - Helen Scott
%       * Turned off the visibility of the figures because a full run ends
%           up making ~600 plots
%   12/02/2021 - Helen Scott
%       * Plot all perturbations for all models in the results folder
%       * Colorbar just shows the percents because perturbation isn't
%           reporting a list for the percentiles or the perturbed values
%   11/15/2021 - Helen Scott
%       * Create skeleton file
%       * Goal: Have a one-click method to (re)generate all plots
%       * Want the plotting to be separrate from the actual perturbation
%      

% Load in the model catalog
model_catalog

% Load in the paramters
base_parameters
% We want to perturb every parameter individually
parameterNames = keys(parameters);

% For loop for all of the models
for i = 1:n_models
    % Collect model info
    modelName = models{i, MODEL_NAME};
    modelFun = models{i, MODEL_FUN};

    % Clean up name (for making directory/files)
    cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
        ' ', '_'), '/', '-');

    % Folder name for results
    resultsPath = './single-parameter-perturbation-results/'; % TODO; change to './single-base-parameter-perturbation-results/'
    folderName = [resultsPath, cleanModelName, '/'];

    % Set the output path for the plots for this circuit model
    outpath = [folderName, 'plots/'];

    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        fprintf('Directory does not exist, attempting to make it: %s', ...
            outpath);
        mkdir(outpath);
    end

    % Make this a for loop- for all of the variables
    for j = 1:length(parameterNames)
        % Collect variable info
        varName = parameterNames{j};

        % Load the perturbation results
        load([folderName, cleanModelName, '-perturb-', varName, ...
                    '.mat'])

        %% Plot
        % Generate colors for all runs
        jetcustom = jet(length(ys));
        % Define center (mean value for the distribution)
        center = length(ys)/2;

        % Plot
        figure('visible', 'off'); % Make but don't show the figure
        for k = 1:length(ys)
            if i == center
                plot(x, ys{k}, 'Color', 'black', 'LineWidth', 2);
            else
                plot(x, ys{k}, 'Color',  jetcustom(k, :));
            end
            hold on;
        end

        % Add color bar
        colormap(jet(length(ys)));
        % Variable values for each 10 percentile
        % Saving the perturbed values is broken in the perturbation code
        % percentiles = [10:10:100];
        % c = zeros(1, length(percentiles));
        % for l = 1:length(percentiles)
        %     c(l) = results{percentiles(i),2}(variablePos);
        % end
        % Color bar for actual variable values
        % caxis([results{1,1} results{end,1}]);
        % cb = colorbar('YTickLabel',c);
        % Placholder color bar
        cb = colorbar;
        ylabel(cb, sprintf('Perturbed Value of %s', varName));
        title(sprintf('Perturbation of %s', varName));

        % Style
        grid on;
        xlabel('Time (Hours)');
        ylabel('AAV');

        % Save figure
        saveas(gcf, [outpath, cleanModelName, '-perturb-', varName, '.png'], 'png');
    end
end

