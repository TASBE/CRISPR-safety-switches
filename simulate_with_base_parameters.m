% This file generates paper figure 2: comparison of all circuits with basal parameters

base_parameters;
model_catalog;

initial = containers.Map();
initial('AAV') = 10;
initial('Cre_regulated_region') = initial('AAV'); % if present, starts unmodified
initial('genome') = 1;

time = [0 720];
y_out = zeros(n_models,721);
y_complete = cell(n_models,1);
fprintf('Simulating with base parameters');
for i=1:n_models,
    try
        [time_interval, y_out(i,:), y_complete{i}] = models{i,MODEL_FUN}(time,parameters,initial,1);
        fprintf('.');
    catch
        fprintf('!');
    end
end
fprintf('\n');

representative = [1, 2, 5, 7, 11, 15];
clusters = {[1 3 4 8 9 10 12 16] [2 6 14] [5 13 17] [7] [11] [15]};
colors = {'k', 'r', 'g', 'b', 'c', 'm'};

% clusters based on: 
% y_out(:,[50 200 end])

legend_names = cell(numel(clusters),1);
for i=1:numel(clusters)
    cluster = clusters{i};
    legend_names{i} = models{cluster(1),MODEL_NAME};
    for j=2:numel(cluster) % plot other cluster members
        legend_names{i} = sprintf('%s, %s',legend_names{i},models{cluster(j),MODEL_NAME});
    end    
end
% No Delay, Repressor, Cre-ON, Chain Cre-ON \rightarrow Repressor, Chain Cre-OFF \rightarrow Repressor, Chain Activator \rightarrow Cre-ON, Chain Repressor \rightarrow Cre-ON, Joint Repressor/Cre-ON
legend_names{1} = 'No Delay, (Chain \rightarrow) Repressor, (Chain \rightarrow) Cre-ON, Parallel Repressor/Cre-ON';

h = figure('PaperPosition',[1 1 6 6]); 
for i=1:numel(representative)
    plot(time_interval/24, y_out(representative(i),:), colors{i}); hold on;
end
for i=1:numel(clusters)
    cluster = clusters{i};
    for j=2:numel(cluster) % plot other cluster members
        plot(time_interval/24, y_out(cluster(j),:), colors{i});
    end
end
legend('Location','SouthOutside',legend_names);
xlabel('Days'); ylabel('[AAV]');
xlim([0 14]); ylim([-0.5 10.5]);
%title('Models with Base Parameters');
outputfig(h,'base_parameters','paper_figures');
