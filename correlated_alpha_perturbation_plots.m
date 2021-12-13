% alpha_perturbation_plots.m
% Last updated: 2021-12-06, by Helen Scott
%
% REVISION HISTORY:
%   12/06/2021 - Helen Scott
%       * Initial implementation, copied all_plots.m
%
%% Set-Up
% Name of the model I am interested in
modelName = 'Sequential Activator \rightarrow Cre-ON';
modelFun = @Chain_Activator_Cre_on_Kill_Switch;

% Clean up name (for making directory/files)
cleanModelName = strrep(strrep(strrep(modelName, '\rightarrow ', ''),...
    ' ', '_'), '/', '-');

% Set the output path
resultsPath = './correlated-alpha-perturbation-results/';

 % Set the output path for all the plots
outpath = [resultsPath, 'plots/'];

% If directory doesn't exist, try to create it
if ~isfolder(outpath)
    disp(['Directory does not exist, attempting to make it: ', outpath]);
    mkdir(outpath);
end

%% Loop
% For all parameter sets
for i = [5, 10, 20]
    
    % Set the parameter set name
    parameterSet = ['tweaked_parameters_', int2str(i)];
    
    % Make an output folder for this sets of tweaked parameters
    folderName = [resultsPath, parameterSet, '/'];
    
    % Load the results
    load([folderName, cleanModelName, '-perturb-alphas.mat']);
    
    %% Plot
    % Generate colors for all runs
    jetcustom = jet(length(ys));
    % Define center (mean value for the distribution)
    center = median(1:length(ys));

    % Plot
    figure('visible', 'off'); % Make but don't show the figure
    for k = 1:length(ys)
        if k == center
            plot(x/24, ys{k}, 'Color', 'black', 'LineWidth', 2);
        else
            plot(x/24, ys{k}, 'Color',  jetcustom(k, :));
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
    ylabel(cb, sprintf('Percentile of Perturbed Values'));
    title(sprintf('Correlated Perturbation of Production Rates (50%% off at %d days)', i));

    % Style
    grid on;
    xlabel('Time (Days)');
    ylabel('AAV');
    ylim([-0.5 10.5]);

    % Save figure
    saveas(gcf, [outpath, parameterSet, '-perturb-alphas.png'], 'png');
end