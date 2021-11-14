% This file generates paper figure 2: comparison of all circuits with basal parameters

base_parameters;
model_catalog;

initial = containers.Map();
initial('AAV') = 10;
initial('Cre_regulated_region') = initial('AAV'); % if present, starts unmodified
initial('CreH_regulated_region') = initial('AAV'); % if present, starts unmodified
initial('genome') = 1;

time = [0 720];
y_out = nan(n_models,721);
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

% raw plot for clustering:
h = figure('PaperPosition',[1 1 6 6]); 
for i=1:n_models
    plot(time_interval, y_out(i,:)); hold on;
end
% 8 clusters
% clusters based on: 
% y_out(:,[20 50 100 end])

clusters = {
    [1 3 4 7 8 9 12 13 14 16 20 21 24 27 29 30 31] % no significant delay
    [5 11 15 19 28] % 1<X<3 at 20: (Repressor)->Cre-OFF
    [2 18 25]       % 4<X<6 at 20: (Repressor ->/Parallel) Activator
    [6 22]          % 1<X<8 at 50, X<1 at 100: Dual Activator
    [26]            % 0.5<2<2 at 50
    [10]            % 1<X<2 at 100 - 'Sequential Cre-OFF \rightarrow Activator'
    [17]            % 1<X<8 at end
    [23]            % 8<x at end - Parallel Activator/Cre-OFF
    };
n_clusters = numel(clusters);
representative = zeros(n_clusters,1);
for i=1:n_clusters, representative(i) = clusters{i}(1); end;
colors = {'k-', 'r-', 'g-', 'b-', 'c-', 'm-','k--','b--'};


legend_names = cell(numel(clusters),1);
for i=1:numel(clusters)
    cluster = clusters{i};
    legend_names{i} = models{cluster(1),MODEL_NAME};
    for j=2:numel(cluster) % plot other cluster members
        legend_names{i} = sprintf('%s, %s',legend_names{i},models{cluster(j),MODEL_NAME});
    end    
end
legend_names{1} = 'No Delay, (Any) Cre-ON, (Any \rightarrow/Parallel) Repressor, Sequential Activator \rightarrow Cre-OFF';
legend_names{2} = '(Cre-OFF/Cre-ON/Repressor \rightarrow) Cre-OFF, Parallel Cre-OFF/Repressor';
legend_names{3} = '(Cre-ON/Repressor \rightarrow/Parallel) Activator';
legend_names{7} = '(Cre-Off/Repressor \rightarrow/Parallel) Cre-OFF';

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
