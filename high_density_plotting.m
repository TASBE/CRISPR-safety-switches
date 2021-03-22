% Load results
load('20210322-TcT1.mat');

% Generate colors for all runs
jetcustom = jet(length(results));

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
    plot(tint, yint, 'Color',  jetcustom(i, :));
    hold on;
end
colormap(jet(length(results)));
cb = colorbar;
caxis([results{1, 1} results{end, 1}]);
ylabel(cb,'Percentile of Perturbed Varriable(s)');
grid on;
title('Perturbation of Tc and T1 Together');
xlabel('Time (Hours)');
ylabel('Expression (%)');
