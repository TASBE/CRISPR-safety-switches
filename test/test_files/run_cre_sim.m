parameters = containers.Map();
parameters('alpha_p_Cas9') = 1;
parameters('delta_Cas9') = 1;
parameters('k_cre') = 0.001;
parameters('alpha_p_Cre') = 0.1;
parameters('delta_Cre') = 0.1;

initial = containers.Map();
initial('AAV') = 10;
initial('genome') = 1;
initial('Cre_regulated_region') = 10;

range = 10.^(-2:0.1:-1);
y_out = zeros(numel(range),73);
for i = 1:numel(range),
    parameters('alpha_p_Cre') = range(i);
    [time_interval, y_out(i,:), y] = simple_recombinase([0 72],parameters,initial);
end

% uncomment to see figures
% figure;
% plot(time_interval, y_out);
% xlabel('time'); ylabel('[Cas9]');
% 
% figure;
% plot(time_interval, y);
% xlabel('time'); ylabel('Molecule');
% legend({'AAV', 'Cas9', 'Cre', 'Cre_regulated_region', 'edited_Cre_regulated_region'});
% 
% figure;
% semilogx(range,y_out(:,end),'*-')
% xlabel('\alpha_{p,TF}'); ylabel('Final [Cas9]');

test_end_points = y_out([1 3 5 7 9 11],end)'
expected = [0.2471    1.3687    4.3349    5.0009    5.0008    5.2009]
