% This file generates paper figure 2: comparison of all circuits with basal parameters

base_parameters;
model_catalog;

initial = containers.Map();
initial('AAV') = 10;
initial('Cre_regulated_region') = initial('AAV'); % if present, starts unmodified
initial('genome') = 1;

FOCUS = 10; 
%models{FOCUS,MODEL_NAME} % confirm we've got 'Chain Activator \rightarrow Cre-ON'

% Base parameters code for ~10^5 proteins max
% Therefore, it's quite reasonable to adjust the alpha values up and down by 100x
tuning = -3:0.1:1;
n_tunings = numel(tuning);

time = [0 720];
y_out = zeros(n_tunings,721);
y_complete = cell(n_tunings,1);

fprintf('Tuning model %s',models{FOCUS,MODEL_NAME});
count = 0;
for t1=1:n_tunings
    tuned_param = containers.Map(parameters.keys, parameters.values);
    tuned_param('alpha_p_Cre') = tuned_param('alpha_p_Cre')*0.03;
    tuned_param('delta_TF') = tuned_param('delta_TF')*10^tuning(t1);
    tuned_param('alpha_p_TF') = tuned_param('alpha_p_TF')*0.01*10^tuning(t1);
    try
        [time_interval, y_out(t1,:), y_complete{t1}] = models{FOCUS,MODEL_FUN}(time,tuned_param,initial,1);
        count = count+1;
        if mod(count,10)==0, fprintf('.'); end;
    catch
        fprintf('!');
    end
end
fprintf('\n');

% Not saving this because it ends up being too big, at 68 MB
%save('parameter_exploration.mat','tuning','time_interval','y_out','y_complete');

h = figure('PaperPosition',[1 1 6 4]); 
for t1=1:n_tunings
    % Red = low alpha_p_TF and delta_TF
    color = [2*max(0,0.5-t1/n_tunings), 0, 2*max(0,t1/n_tunings-0.5)];
    plot(time_interval/24, squeeze(y_out(t1,:)), 'Color', color); hold on;
end
xlabel('Days'); ylabel('[AAV]');
xlim([0 30]); ylim([-0.5 10.5]);
outputfig(h,'tuning_of_chain_activator_cre','paper_figures');
