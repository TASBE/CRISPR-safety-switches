parameters = containers.Map();
parameters('K_R') = 100;
parameters('alpha_p_Cas9') = 1;
parameters('alpha_p_TF') = 1;
parameters('delta_Cas9') = 1;
parameters('delta_TF') = 1;
parameters('n') = 2;

range = 10.^(0:0.1:2);
final = zeros(numel(range),1);
for i = 1:numel(range),
    parameters('alpha_p_TF') = range(i);
    [time_interval, y_out(i,:), y] = simple_repression([0 72],parameters,10);
end

if false % turn on to see figures
    figure; 
    plot(time_interval, y_out);
    xlabel('time'); ylabel('[Cas9]');
    
    figure; 
    semilogx(range,y_out(:,end),'*-')
    xlabel('\alpha_{p,TF}'); ylabel('Final [Cas9]');
end

test_end_points = y_out([1 6 11 15 21],end)';
expected = [9.8951    9.0897    4.9997    1.3675    0.0990];
