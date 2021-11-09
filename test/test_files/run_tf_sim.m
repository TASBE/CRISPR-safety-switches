parameters = containers.Map();
parameters('K_R') = 100;
parameters('alpha_p_Cas9') = 1;
parameters('alpha_p_TF') = 1;
parameters('delta_Cas9') = 1;
parameters('delta_TF') = 1;
parameters('n') = 2;

initial = containers.Map();
initial('AAV') = 10;
initial('genome') = 1;

range = 10.^(0:0.1:2);
y_out = zeros(numel(range),73);
for i = 1:numel(range),
    parameters('alpha_p_TF') = range(i);
    [time_interval, y_out(i,:), y] = simple_repression([0 72],parameters,initial);
end

% uncomment to see figures
% figure;
% plot(time_interval, y_out);
% xlabel('time'); ylabel('[Cas9]');
%
% figure;
% semilogx(range,y_out(:,end),'*-')
% xlabel('\alpha_{p,TF}'); ylabel('Final [Cas9]');

test_end_points = y_out([1 6 11 15 21],end)'
expected = [9.9010    9.0909    5.0000    1.3681    0.0990]
