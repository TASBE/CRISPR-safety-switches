% This file generates paper figure 2: comparison of all circuits with basal parameters

base_parameters;
model_catalog;

initial = containers.Map();
initial('AAV') = 10;
initial('Cre_regulated_region') = initial('AAV'); % if present, starts unmodified
initial('genome') = 1;

% Base parameters code for ~10^5 proteins max
% Therefore, it's quite reasonable to adjust the alpha values up and down by 100x
tuning = -2:0.5:2;
n_tunings = numel(tuning);

time = [0 720];
y_out = zeros(n_models,n_tunings,n_tunings,721);
y_complete = cell(n_models,n_tunings,n_tunings);
for i=2:n_models,
    fprintf('Tuning model %s',models{i,MODEL_NAME});
    count = 0;
    for t1=1:n_tunings
        for t2=1:n_tunings
            tuned_param = containers.Map(parameters.keys, parameters.values);
            tuned_param('alpha_p_TF') = tuned_param('alpha_p_TF')*10^tuning(t1);
            tuned_param('alpha_p_Cre') = tuned_param('alpha_p_Cre')*10^tuning(t2);
            try
                [time_interval, y_out(i,t1,t2,:), y_complete{i,t1,t2}] = models{i,MODEL_FUN}(time,tuned_param,initial,1);
                count = count+1;
                if mod(count,50)==0, fprintf('.'); end;
            catch
                fprintf('!');
            end
        end
    end
    fprintf('\n');
end

% Not saving this because it ends up being too big, at 68 MB
%save('parameter_exploration.mat','tuning','time_interval','y_out','y_complete');

greymap = colormap;
greymap(end,:) = [0.5 0.5 0.5];
for i=2:n_models,
    halfway_time = inf(n_tunings,n_tunings);
    for t1=1:n_tunings
        for t2=1:n_tunings
            idx = find(y_out(i,t1,t2,:) < 5,1);
            if numel(idx), halfway_time(t1,t2) = idx; end;
        end
    end
    
    h = figure('PaperPosition',[1 1 6 4]); 
    imagesc(tuning,tuning,halfway_time/24);
    ylabel('alpha_{p,TF}'); xlabel('alpha_{p,Cre}');
    c = colorbar;
    c.Label.String = 'Days to 50% elimination';
    colormap(greymap);
    caxis([0 30]);
    title(models{i,MODEL_NAME});
    outputfig(h,sprintf('heatmap_%s',models{i,MODEL_NAME}),'paper_figures/exploration');
end

for i=2:n_models,
    h = figure('PaperPosition',[1 1 6 4]); 
    for t1=1:n_tunings
        for t2=1:n_tunings
            %color = [1-(t1/n_tunings), 0, t2/n_tunings];
            color = [max(0,-tuning(t1)/2), 1-t2/n_tunings, max(0,tuning(t1)/2)];
            plot(time_interval/24, squeeze(y_out(i,t1,t2,:)), 'Color', color); hold on;
        end
    end
    xlabel('Days'); ylabel('[AAV]');
    xlim([0 30]); ylim([-0.5 10.5]);
    title(models{i,MODEL_NAME});
    outputfig(h,sprintf('parameter_exploration_%s',models{i,MODEL_NAME}),'paper_figures/exploration');
end
