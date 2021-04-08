function densityheatmap(results, nBins, cLimits, name, path)
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
    %   04/08/2021 - hscott
    %       * Initial implement
    %
    
    %% Make outout path path
    % Default path is current path
    if nargin < 5
        path = './plots/';  % Note: frontslash works for both Windows and Mac/Unix
    else
        path = [path, 'plots/'];
    end 
    
    % If directory doesn't exist, try to create it
    if ~isfolder(path)
        sanitized_path = strrep(path, '/', '&#47;');
        sanitized_path = strrep(sanitized_path, '\', '&#92;');
        sanitized_path = strrep(sanitized_path, ':', '&#58;');
        disp('Directory does not exist, attempting to create it: %s',sanitized_path);
        mkdir(path);
    end

    %% Add points to arrays
    x = [];
    y = [];
    for i = 1:length(results)
        sol = results{i, end};
        tint = 1:336;
        yint = zeros(size(tint));
        for j = tint
            ddeRes = deval(sol, j);
            yint(j) = 100* ddeRes(5) / (ddeRes(4) + ddeRes(5));
        end
        x = [x, tint];
        y = [y, yint];
    end

    %% Plot
    figure;
    hist3([x', y'], nBins, 'CdataMode','auto', 'EdgeColor', 'none');
    xlabel('Time (hours)');
    ylabel('% EGFP neg');
    cb = colorbar;
    caxis(cLimits);
    ylabel(cb,'Number of Trajectories Per Cell');
    view(2);
    title(['Trajectory Density of ', name, ' perturbation']);
    
    %%  Save figure
    % Make a clean name
    cleanName = sanitize_filename(name);
    % Output
    saveas(gcf, [path cleanName '.eps'], 'epsc');
    saveas(gcf, [path cleanName '.png'], 'png');
    saveas(gcf, [path cleanName '.fig']);