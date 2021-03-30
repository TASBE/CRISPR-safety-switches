% Load results
load('20210322-lag.mat');

% Generate colors for all runs
nColors = (ceil(range(logDist)) + 1) * 1000;
jetcustom = jet(nColors);

offset = ceil(logDist(1)-1);
center = length(results)/2;

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
        plot(tint, yint, 'Color',  jetcustom(((ceil(logDist(i)) - offset) * 1000), :));
    end
    hold on;
end
set(gca,'ColorScale','log')
colormap(jetcustom);
cb = colorbar;
caxis([0 ceil(logDist(end))]);
ylabel(cb,'Value of Perturbed Varriable');
grid on;
title('Perturbation of T1 Alone');
xlabel('Time (Hours)');
ylabel('Expression (%)');
