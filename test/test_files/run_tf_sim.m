parameters = containers.Map();
parameters('K_R') = 100;
parameters('alpha_p_Cas9') = 1;
parameters('alpha_p_TF') = 1;
parameters('delta_Cas9') = 1;
parameters('delta_TF') = 1;
parameters('n') = 2;

initial = containers.Map();
initial('AAV') = 10;

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

test_end_points = y_out([1 6 11 15 21],end)';
expected = [9.8951    9.0897    4.9997    1.3675    0.0990];
