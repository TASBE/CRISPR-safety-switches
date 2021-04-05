% sol = GeneTherapySystemElimination(720, 0.6);

figure;
tint = 1:720;
yint = zeros(size(tint));
for j = tint
    OdeRes = deval(sol, j);
    yint(j) = y(10);
end
plot(tint, yint);