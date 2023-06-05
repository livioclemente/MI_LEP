% Step 1: Calcolare le differenze tra le medie dei due gruppi ad ogni timepoint
diff_means = mean(data_ct, 1) - mean(data_pth, 1);

% Step 2: Calcolare la statistica t ad ogni timepoint
[~, ~, ~, stats] = ttest2(data_ct, data_pth);
t_stats = stats.tstat;

% Step 3: Trovare i cluster di timepoints contigui in cui t_stats supera una certa soglia
df = size(data_ct, 1) - 1 + size(data_pth, 1) - 1;  % gradi di libertà
threshold = tinv(1 - 0.05/2, df);
clusters = bwconncomp(abs(t_stats) > threshold);

% Step 4: Calcolare la somma dei valori di t_stats all'interno di ogni cluster
cluster_stats = cellfun(@(x) sum(t_stats(x)), clusters.PixelIdxList);

% Step 5: Permutare i labels dei gruppi e ripetere i passaggi 1-4 per ogni permutazione
n_permutations = 1000;
max_cluster_stats = zeros(n_permutations, 1);
combined_data = [data_ct; data_pth];

for i_permutation = 1:n_permutations
    permuted_labels = randperm(size(combined_data, 1));
    permuted_data_ct = combined_data(permuted_labels(1:size(data_ct, 1)), :);
    permuted_data_pth = combined_data(permuted_labels((size(data_ct, 1) + 1):end), :);

    permuted_diff_means = mean(permuted_data_ct, 1) - mean(permuted_data_pth, 1);
    [~, ~, ~, permuted_stats] = ttest2(permuted_data_ct, permuted_data_pth);
    permuted_t_stats = permuted_stats.tstat;
    permuted_clusters = bwconncomp(abs(permuted_t_stats) > threshold);
    permuted_cluster_stats = cellfun(@(x) sum(permuted_t_stats(x)), permuted_clusters.PixelIdxList);

    if ~isempty(permuted_cluster_stats)
        max_cluster_stats(i_permutation) = max(abs(permuted_cluster_stats));
    end
end

% Step 6: Confrontare ciascuna statistica del cluster originale con la distribuzione delle statistiche del cluster massimo permutato per calcolare i valori p
p_values = arrayfun(@(x) mean(max_cluster_stats >= abs(x)), cluster_stats);
disp('P-values for each cluster:');
disp(p_values);

% Draw the t-values line
plot(t_stats, 'LineWidth', 1.5);

% Set y-axis limits
ylim([-4 4]);

% Draw a horizontal line at zero
line([0 length(t_stats)], [0 0], 'Color', 'k', 'LineWidth', 1);

% Draw the threshold lines
hold on;
plot([0 length(t_stats)], [threshold threshold], '--k');
plot([0 length(t_stats)], [-threshold -threshold], '--k');

% Highlight significant clusters
for i_cluster = 1:length(clusters.PixelIdxList)
    if p_values(i_cluster) < 0.05
        timepoints = clusters.PixelIdxList{i_cluster};
        start_time = min(timepoints);
        end_time = max(timepoints);
        area(start_time:end_time, [threshold threshold], 'FaceColor', 'r');
    end
end

% Axis labels
xlabel('Timepoint');
ylabel('t-value');

hold off;

