[x,y] = GeneTherapySystemElimination(tspan, varsToUse)

figure;
plot(x,y(:,10))

sol = results{2,3};
x = linspace(0,720,720);
y = deval(sol, x)

figure;
yint = y(10, :)
plot(x,yint)

% For safe keeping
%         yint = zeros(size(tint));
%         for j = tint
%             ddeRes = deval(sol, j);
%             yint(j) = 100* ddeRes(5) / (ddeRes(4) + ddeRes(5));
%         end

%% Functions to extract the value to plot from the results of the DE
function yint = solvePerturbableFitCas9AssumeDelay(y)
    yint = 100 * y(5) / (y(4) + y(5));