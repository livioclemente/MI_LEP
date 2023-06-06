% Creazione dei dati
Ivalues_pth = cell(num_subjects_pth, 1);
for s = 1:num_subjects_pth
     Ivalues_pth{s} = I_values_pth{s};
end
data_I_pth = cell2mat(Ivalues_pth);

Ivalues_ct = cell(num_subjects_ct, 1);
for s = 1:num_subjects_ct
     Ivalues_ct{s} = I_values_ct{s};
end
data_I_ct = cell2mat(Ivalues_ct);

% Passaggio 1: Calcola il valore t per confrontare i due gruppi per ogni campione
t_values = zeros(1, 282);

for sample = 1:282
    [~, ~, ~, stats] = ttest2(data_I_ct(:, sample), data_I_pth(:, sample));
    t_values(sample) = stats.tstat;
end

% Step 2: Calcolare il threshold
alpha = 0.05; % livello di significatività
df = num_subjects_ct + num_subjects_pth - 2; % gradi di libertà
t_critical = tinv(1 - alpha/2, df); % valore t critico per un test a due code
threshold = t_critical; % soglia basata sulla distribuzione del t-value

% Passaggio 3: Seleziona i campioni con valori t superiori a threshold
positive_samples = t_values > threshold;

% Passaggio 4: Seleziona i campioni con valori t inferiori a -threshold
negative_samples = t_values < -threshold;

% Clustering dei campioni significativi positivi
positive_clusters = bwconncomp(positive_samples);
for i = 1:numel(positive_clusters.PixelIdxList)
    cluster = positive_clusters.PixelIdxList{i};
    % Fai qualcosa con il cluster dei campioni positivi
end

% Clustering dei campioni significativi negativi
negative_clusters = bwconncomp(negative_samples);
for i = 1:numel(negative_clusters.PixelIdxList)
    cluster = negative_clusters.PixelIdxList{i};
    % Fai qualcosa con il cluster dei campioni negativi
end

% Calcola le statistiche del cluster per i campioni positivi
positive_cluster_stats = cellfun(@(cluster) sum(t_values(cluster)), positive_clusters.PixelIdxList);

% Calcola le statistiche del cluster per i campioni negativi
negative_cluster_stats = cellfun(@(cluster) sum(t_values(cluster)), negative_clusters.PixelIdxList);

% Seleziona la statistica del cluster più grande per i campioni positivi
max_positive_cluster_stat = max(positive_cluster_stats);

% Seleziona la statistica del cluster più grande per i campioni negativi
max_negative_cluster_stat = max(negative_cluster_stats);

figure;
hold on;
plot(time_axis_corrected_ct, t_values, 'k', 'LineWidth', 1);
line([time_axis_corrected_ct(1), time_axis_corrected_ct(end)], [threshold, threshold], 'Color', 'k', 'LineStyle', '--');
line([time_axis_corrected_ct(1), time_axis_corrected_ct(end)], [-threshold, -threshold], 'Color', 'k', 'LineStyle', '--');
line([time_axis_corrected_ct(1), time_axis_corrected_ct(end)], [0, 0], 'Color', 'k', 'LineWidth', 1);

% Plot del cluster positivo più grande
if ~isempty(positive_clusters.PixelIdxList)
    positive_cluster_sizes = cellfun(@numel, positive_clusters.PixelIdxList);
    [~, max_positive_cluster_idx] = max(positive_cluster_sizes);
    max_positive_cluster = positive_clusters.PixelIdxList{max_positive_cluster_idx};
    x = time_axis_corrected_ct(max_positive_cluster);
    y = t_values(max_positive_cluster);
    patch([x, fliplr(x)], [y, zeros(size(y))], 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Plot del cluster negativo più grande
if ~isempty(negative_clusters.PixelIdxList)
    negative_cluster_sizes = cellfun(@numel, negative_clusters.PixelIdxList);
    [~, max_negative_cluster_idx] = max(negative_cluster_sizes);
    max_negative_cluster = negative_clusters.PixelIdxList{max_negative_cluster_idx};
    x = time_axis_corrected_ct(max_negative_cluster);
    y = t_values(max_negative_cluster);
    patch([x, fliplr(x)], [y, zeros(size(y))], 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% Imposta limiti dell'asse y
ylim([-3, 8]);
xlim([-100, 1000]);

xlabel('Tempo');
ylabel('t-value');
title('Cluster Positivo e Negativo più Grande');
hold off;
