x = randn(1,1000);
y = randn(1,1000);
n = hist3([x', y']);
pcolor(n)

scatter(x, y)

load carbig
X = [MPG,Weight];
hist3(X,'CdataMode','auto')
xlabel('MPG')
ylabel('Weight')
colorbar
view(2)


x2 = []; % Change to be pre-determined size
y2 =[]; % Change to be pre-determined size

% Save points rather than plot
for i = 1:length(results)
    sol = results{i, 3};
    tint = 1:336;
    yint = zeros(size(tint));
    for j = tint
        ddeRes = deval(sol, j);
        yint(j) = 100* ddeRes(5) / (ddeRes(4) + ddeRes(5));
    end
    x2 = [x2, tint];
    y2 = [y2, yint];
end

% Scatterplot
figure(1);
scatter(x2, y2)
% 3D Histogram
figure(2);
hist3([x2', y2']);
% Heat map
figure(3);
n = hist3([x2', y2']);
pcolor(n)
% 3D histogram view 2
figure(4);
hist3([x2', y2'],'CdataMode','auto')
xlabel('Time (hours)')
ylabel('% EGFP neg')
colorbar
view(2)
% Try at a higher resolution
figure(5);
hist3([x2', y2'], [336 100], 'CdataMode','auto')
xlabel('Time (hours)')
ylabel('% EGFP neg')
colorbar
view(2)