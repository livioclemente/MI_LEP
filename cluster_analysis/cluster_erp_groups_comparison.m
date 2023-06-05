erp_pth = cell(num_subjects_pth, 1);
for s = 1:num_subjects_pth
     erp_pth{s} = erp_mean_baseline_corrected_pth{s}(11, :);
end
data_pth = cell2mat(erp_pth);

erp_ct = cell(num_subjects_ct, 1);
for s = 1:num_subjects_ct
     erp_ct{s} = erp_mean_baseline_corrected_ct{s}(11, :);
end
data_ct = cell2mat(erp_ct);

%%
% Calcola il numero di soggetti e il numero di timepoints per ogni gruppo
[num_subjects_ct, num_timepoints] = size(data_ct);
[num_subjects_pth, ~] = size(data_pth);

% Controlla che i gruppi abbiano lo stesso numero di timepoints
if size(data_ct, 2) ~= size(data_pth, 2)
    error('Il numero di timepoints non coincide tra i due gruppi.');
end

% Inizializza una matrice per memorizzare la statistica del test
t_values = zeros(num_timepoints, 1);

% Calcola la statistica del test per ogni timepoint
for timepoint = 1:num_timepoints
    [~, p, ~, stats] = ttest2(data_ct(:, timepoint), data_pth(:, timepoint));
    t_values(timepoint) = stats.tstat;
end

alpha = 0.05;  % livello di significatività
df = (num_subjects_pth+num_subjects_ct)-2;  % gradi di libertà
t_critical = tinv(1-alpha/2, df);  % valore t critico per un test a due code
% Applica la soglia per definire i clusters
threshold = t_critical; % questo dipende dal tuo numero di soggetti e dal tuo livello di significatività

clusters_positive = t_values > threshold;
clusters_negative = t_values < -threshold;

% Calcola la dimensione di ciascun cluster
cluster_sizes_positive = sum(clusters_positive, 1);
cluster_sizes_negative = sum(clusters_negative, 1);

% Concatena i dati dei due gruppi lungo la prima dimensione
data_all = [data_ct; data_pth];  % dimensione: [(num_soggetti_pth + num_soggetti_ct), num_timepoints]

% Calcola la dimensione di ciascun cluster permutato
num_permutations = 1000;
permuted_cluster_sizes_positive = zeros(num_permutations, 1);
permuted_cluster_sizes_negative = zeros(num_permutations, 1);
for permutation = 1:num_permutations
    % Permuta le etichette dei gruppi
    permuted_data = data_all(randperm(num_subjects_ct + num_subjects_pth), :);
    
    % Calcola la statistica del test per ogni timepoint
    permuted_t_values = zeros(num_timepoints, 1);
    for timepoint = 1:num_timepoints
        [~, ~, ~, stats] = ttest2(permuted_data(1:num_subjects_ct, timepoint), permuted_data((num_subjects_ct+1):end, timepoint));
        permuted_t_values(timepoint) = stats.tstat;
    end

    % Applica la soglia per definire i clusters
    permuted_clusters_positive = permuted_t_values > threshold;
    permuted_clusters_negative = permuted_t_values < -threshold;
    
    % Dimensione del cluster più grande
    permuted_cluster_sizes_positive(permutation) = max(sum(permuted_clusters_positive, 1));
    permuted_cluster_sizes_negative(permutation) = max(sum(permuted_clusters_negative, 1));
end

% p-values dei clusters
p_values_positive = zeros(size(cluster_sizes_positive));
p_values_negative = zeros(size(cluster_sizes_negative));
for cluster = 1:length(cluster_sizes_positive)
    p_values_positive(cluster) = mean(permuted_cluster_sizes_positive >= cluster_sizes_positive(cluster));
end
for cluster = 1:length(cluster_sizes_negative)
    p_values_negative(cluster) = mean(permuted_cluster_sizes_negative >= cluster_sizes_negative(cluster));
end

% p-values
disp('P-values for positive clusters:');
disp(p_values_positive);
disp('P-values for negative clusters:');
disp(p_values_negative);

figure;
hold on;

% Draw the t-values line
plot(time_axis_corrected_pth, t_values, 'LineWidth', 1.5);

% Set x-axis limits
xlim([-100 1000]);

% Set y-axis limits
ylim([-3 5]);

% Draw a horizontal line at zero
line([time_axis_corrected_pth(1) time_axis_corrected_pth(end)], [0 0], 'Color', 'k', 'LineWidth', 1);

% Draw the threshold lines
plot(time_axis_corrected_pth, threshold * ones(size(time_axis_corrected_pth)), '--k');
plot(time_axis_corrected_pth, -threshold * ones(size(time_axis_corrected_pth)), '--k');

significant_clusters_positive = find(cluster_sizes_positive > 0);
significant_clusters_negative = find(cluster_sizes_negative > 0);

% Draw areas for significant clusters
for cluster = significant_clusters_positive
    start_time = find(clusters_positive(cluster,:), 1, 'first');
    end_time = find(clusters_positive(cluster,:), 1, 'last');
    x = [time_axis_corrected_pth(start_time:end_time), fliplr(time_axis_corrected_pth(start_time:end_time))];
    y = [t_values(start_time:end_time), fliplr(t_values(start_time:end_time))];
    fill(x, y, 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

for cluster = significant_clusters_negative
    start_time = find(clusters_negative(cluster,:), 1, 'first');
    end_time = find(clusters_negative(cluster,:), 1, 'last');
    x = [time_axis_corrected_pth(start_time:end_time), fliplr(time_axis_corrected_pth(start_time:end_time))];
    y = [t_values(start_time:end_time), fliplr(t_values(start_time:end_time))];
    fill(x, y, 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end
 % la parte per colorare le aree non funziona

% Axis labels
xlabel('Time (ms)');
ylabel('Average t-value');
title('Time-course of the t-test statistic');


hold off;
