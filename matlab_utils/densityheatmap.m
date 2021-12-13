function densityheatmap(x, ys, nBins, cLimits, name, path)
    % FUNCTION NAME:
    %   densityheatmap
    %
    % DESCRIPTION:
    %   Create ans save a plot of the 
    %
    % INPUT:
    %   results - (cell) The output from logNormalPerturbation (Column 1 = 
    %       percentile of the perturbed value in the log normal 
    %       distribution (NaN for random distributions), column 2 = 
    %       varaibles used for the differential equation, column 3 = 
    %       results from the differential equation (sol))
    %   nBins - (double) The number of bins in each dimension of the 
    %       histogram
    %   cLimits - (double) The colormap limits: a two-element vector of the
    %       form [cmin cmax]
    %   name - (char) Title to be used in plot title, file names
    %   path - (char) output path, will have 'plots/' added to it, default
    %       = './'
    %
    % OUTPUT:
    %   .eps, .fig, .png figure saved in the 'plots/' subfolder of outpath
    %
    % ASSUMPTIONS AND LIMITATIONS:
    %   Path is assumed to end in '/' if supplied
    %
    % REVISION HISTORY:
    %   12/06/2021 - Helen Scott
    %       * Changes to work with new results structure
    %   04/08/2021 - hscott
    %       * Initial implement
    %

    %% Add points to arrays
    xPoints = [];
    yPoints = [];
    for i = 1:length(ys)
        xPoints = [xPoints, x/24];
        yPoints = [yPoints, ys{i}];
    end

    %% Plot
    figure('visible', 'off'); % Make but don't show the figure
    hist3([xPoints', yPoints'], nBins, 'CdataMode','auto', 'EdgeColor', 'none');
    xlabel('Time (Days)');
    ylabel('AAV');
    cb = colorbar;
    caxis(cLimits);
    ylabel(cb,'Number of Trajectories Per Cell');
    view(2);
    title(sprintf('Trajectory Density of Random Perturbation (50%% off at %d days)', name));
    
    %%  Save figure
    saveas(gcf, [path, 'random-perturb-all-tweaked-parameters-', int2str(name), '.png'], 'png');