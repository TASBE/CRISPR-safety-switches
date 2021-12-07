% tweaked_parameter_perturbation_plots.m
% Last updated: 2021-12-06, by Helen Scott
%
% REVISION HISTORY:
%   12/06/2021 - Helen Scott
%       * Initial implementation, copied all_plots.m
%

% Name of the model I am interested in
modelName = 'Sequential Activator \rightarrow Cre-ON';
modelFun = @Chain_Activator_Cre_on_Kill_Switch;

% Clean up name (for making directory/files)
cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
    ' ', '_'), '/', '-');

% Load in the paramters, just care about the names, not the values so base
% is fine
base_parameters
% We want to perturb every parameter individually
parameterNames = keys(parameters);

% Folder with the results for all of the parameter sets
resultsPath = './tweaked-parameter-perturbation-results/';

% For loop for all of the models
for i = [5, 10, 20]
    % Set the parameter name
    parameterSet = ['tweaked_parameters_', int2str(i)];
    
    % Load the parameters
    eval(parameterSet)    
    
    % Make an output folder for this sets of tweaked parameters
    folderName = [resultsPath, parameterSet, '/'];

    % Set the output path for the plots for this circuit model
    outpath = [folderName, 'plots/'];

    % If directory doesn't exist, try to create it
    if ~isfolder(outpath)
        disp(['Directory does not exist, attempting to make it: ', ...
            outpath]);
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
                plot(x/24, ys{k}, 'Color', 'black', 'LineWidth', 2);
            else
                plot(x/24, ys{k}, 'Color',  jetcustom(k, :));
            end
            hold on;
        end

        % Add color bar
        colormap(jet(length(ys)));
        % Variable values for each 10 percentile
%         percentiles = [10:10:100];
%         c = zeros(1, length(percentiles));
%         for l = 1:length(percentiles)
%             c(l) = results{percentiles(i),2}(variablePos);
%         end
% %         Color bar for actual variable values
%         caxis([results{1,1} results{end,1}]);
%         cb = colorbar('YTickLabel',c);
        % Placholder color bar
        cb = colorbar;
        ylabel(cb, sprintf('Perturbed Value of %s', varName));
        title(sprintf('Perturbation of %s', varName));

        % Style
        grid on;
        xlabel('Time (Days)');
        ylabel('AAV');

        % Save figure
        saveas(gcf, [outpath, cleanModelName, '-perturb-', varName, '.png'], 'png');
    end
end