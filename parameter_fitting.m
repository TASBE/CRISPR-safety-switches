% Hardcode the data
xFig1 = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 10 10 10 13 13 13];
yFig1 = [0.53 0.54, 0.59 5.62 5.83 6.26 36.9 36.2 38.7 67.8 65.6 66.9 79 76 78.7 96 97.2 96.7 95.2 96.8 96];
% Get the average for each unique time point
[G, days] = findgroups(xFig1);
yMeansFig1 = splitapply(@mean, yFig1, G);
% Convert to hours
hours = days .* 24;

% Plot just to check the data (i.e. for typos)
plot(xFig1.*24, yFig1, 'ko')
hold on;
plot(hours, yMeansFig1, 'r-')

% Create the model, see the end of the file

% Use a starting point of everything equal to 1
% parameters0 = ones(1, 7);
% Fitting was not changing the starting values
% Try initial values on different scales to see if that helps
parameters0 = [4 -2 4 0 -6 -4 -4];
% parameters(1) = alpha_P (molecules/hour)
% parameters(2) = delta_P (fraction/hour): for a several-day halflife, 10^-2
%   Note: max_P = alpha_P/delta_P - should be in range of 10^1 - 10^8
%   With a several-day halflife, this implies alpha_P in 10^-1 - 10^6
% parameters(3) = alpha_g (molecules/hour)
% parameters(4) = delta_g (fraction/hour): note that 10^-1 => ~6 hour halflife
% parameters(5) = k_cre
% parameters(6) = k_cas
% parameters(7) = k_cat

% Set biologically releveant upper and lower bounds % TODO
lb = [-1 -4 -1 -2 -10 -10 -10];
ub = [ 6  0  6  0  -4  -2  -2];

% Fit
% x = lsqcurvefit(@fit_Cre_on_Kill_Switch, parameters0, hours, yMeansFig1, lb, ub);

x = fminsearch(@(x)(sum((fit_Cre_on_Kill_Switch(x,hours) - yMeansFig1).^2)),parameters0)
%  x = 3.0415   -1.6977    3.3090    0.0003   -7.1535   -4.2577   -4.7518

% Run the model with the fit parameters
res = fit_Cre_on_Kill_Switch(x, 1:350);

% Plot the model results
h = figure('PaperPosition',[1 1 5 4]);
plot(xFig1.*24, yFig1, 'ko'); hold on;
plot(xFig1.*24, yFig1, 'r-');
plot(1:350, res, 'b--');
xlabel('Hours');
ylabel('% negative GFP');
legend('Location','SouthEast','Data','Mean','Fit');
outputfig(h,'parameter_fit','plots');

