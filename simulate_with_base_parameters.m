% This file generates paper figure 2: comparison of all circuits with basal parameters

base_parameters;
model_catalog;

initial = containers.Map();
initial('AAV') = 10;
initial('Cre_regulated_region') = initial('AAV'); % if present, starts unmodified
initial('genome') = 1;

time = [0 720];
y_out = zeros(3,73);
fprintf('Simulating with base parameters');
for i=1:n_models,
    try
        [time_interval, y_out(i,:), ~] = models{i,MODEL_FUN}(time,parameters,initial,10);
        fprintf('.');
    catch
        fprintf('!');
    end
end
fprintf('\n');

h = figure('PaperPosition',[1 1 6 4]); 
plot(time_interval, y_out);
legend(models{:,MODEL_NAME});
xlabel('Hours'); ylabel('[AAV]');
%title('Models with Base Parameters');
outputfig(h,'base_parameters','paper_figures');
