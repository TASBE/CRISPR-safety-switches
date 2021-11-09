addpath('generated_models')

parameters = containers.Map();
% Therapy parameters - cannot be modified
parameters('Cas_degradation') = 0.01;  % (assuming bound, unbound, and post-edit are equal)
parameters('Cas_gRNA_binding') = 1; % assuming identical for all gRNAs
parameters('alpha_r_sgRNA2') = 0.1;
parameters('alpha_p_Cas9') = 0.1;
parameters('delta_Cas9') = parameters('Cas_degradation');
parameters('delta_g') = 0.2;
parameters('k_cat') = 0.1;

% Base kill switch parameters
parameters('alpha_r_sgRNA1') = parameters('alpha_r_sgRNA2'); % assume not readily modified 

% Delay mechanism parameters
parameters('alpha_p_TF') = 0.01; % Can be modulated - explore values
parameters('delta_TF') = parameters('alpha_p_Cas9'); % Assume stable proteins
parameters('K_A') = 10;         % Cannot be modulated; determine by fit
parameters('K_R') = 1;          % Cannot be modulated; determine by fit
parameters('n') = 2;            % Cannot be modulated; determine by fit

parameters('alpha_p_Cre') = 0.01; % Can be modulated - explore values
parameters('delta_Cre') = parameters('alpha_p_Cas9'); % Assume stable proteins
parameters('k_cre') = 0.001;    % Cannot be modulated; determine by fit


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

h = figure('PaperPosition',[1 1 6 4]); 
plot(time_interval, y_out);
legend('Base','Activator','Cre-on','Joint Cre-on/Activator','Chain Repressor Cre-on');
xlabel('Hours'); ylabel('[AAV]');
title('Sampler of models (non-fit parameters)');
outputfig(h,'sampler','generated_models/plots');