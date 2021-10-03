addpath('generated_models')

parameters = containers.Map();
parameters('Cas_degradation') = 0.01;
parameters('Cas_gRNA_binding') = 1;
parameters('alpha_r_sgRNA1') = 0.1;
parameters('alpha_r_sgRNA2') = 0.1;
parameters('delta_g') = 0.2;
parameters('k_cat') = 0.1;

parameters('K_A') = 10;
parameters('K_R') = 1;
parameters('alpha_p_Cas9') = 0.1;
parameters('alpha_p_TF') = 0.1;
parameters('delta_Cas9') = 0.1;
parameters('delta_TF') = 1;
parameters('n') = 2;

parameters('k_cre') = 0.001;
parameters('alpha_p_Cre') = 0.01;
parameters('delta_Cre') = 0.1;


initial = containers.Map();
initial('AAV') = 10;
initial('Cre_regulated_region') = initial('AAV'); % if present, starts unmodified
initial('genome') = 1;

time = [0 720];
y_out = zeros(3,73);
y = cell(3,1);
[~, y_out(1,:), y{1}] = Basic_kill_switch(time,parameters,initial,10);
[~, y_out(2,:), y{2}] = TF_delayed_kill_switch(time,parameters,initial,10);
[~, y_out(3,:), y{3}] = Cre_delayed_kill_switch(time,parameters,initial,10);
[~, y_out(4,:), y{4}] = Joint_delayed_kill_switch(time,parameters,initial,10);
[time_interval, y_out(5,:), y{5}] = Chain_delayed_kill_switch(time,parameters,initial,10);

figure; 
plot(time_interval, y_out);
legend('Base','TF','Cre','Joint','Chain');
xlabel('time'); ylabel('[AAV]');
