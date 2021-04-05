% Load results
load('dev-results/20210330-rand-lag.mat');

tspan = 336; % Change to pass this on from the data generation script

x = [];
y = [];

% Add points to arrays
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

% Plot
figure;
hist3([x', y'], [336 100], 'CdataMode','auto')
xlabel('Time (hours)')
ylabel('% EGFP neg')
colorbar
caxis([0 20])
view(2)