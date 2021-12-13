tweaked_parameters_20
model_catalog;

initial = containers.Map();
initial('AAV') = 10;
initial('Cre_regulated_region') = initial('AAV'); % if present, starts unmodified
initial('genome') = 1;

FOCUS = find(cellfun(@(x)(strcmp(x,'Sequential Activator \rightarrow Cre-ON')),models(:,MODEL_NAME)));

time = [0 720];
tuned_param = containers.Map(parameters.keys, parameters.values);

[time_interval, y_out, y_complete] = models{FOCUS,MODEL_FUN}(time,tuned_param,initial,1);

h = figure('PaperPosition',[1 1 6 4]); 
plot(time_interval/24, squeeze(y_out)); hold on;
xlabel('Days'); ylabel('[AAV]');
xlim([0 30]); ylim([-0.5 10.5]);