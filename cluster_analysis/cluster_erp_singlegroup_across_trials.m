erp_dlpthpiede = cell(num_subjects, 1);
for s = 1:num_subjects
     erp_dlpthpiede{s} = erp_mean_baseline_corrected{s}(11, :);
end
data_foot = cell2mat(erp_dlpthpiede);

erp_dlpthmano = cell(num_subjects, 1);
for s = 1:num_subjects
     erp_dlpthmano{s} = erp_mean_baseline_corrected_dlpthmano{s}(11, :);
end
data_hand = cell2mat(erp_dlpthmano);

erp_dlpthgin = cell(num_subjects, 1);
for s = 1:num_subjects
     erp_dlpthgin{s} = erp_mean_baseline_corrected_dlpthgin{s}(11, :);
end
data_knee = cell2mat(erp_dlpthgin);

%%

% Concatena i dati delle tre condizioni lungo la terza dimensione
data_all = cat(3, data_foot, data_knee, data_hand);

% Calcola il numero di soggetti e il numero di timepoints
[num_subjects, num_timepoints, num_conditions] = size(data_all);

% Inizializza una matrice per memorizzare la statistica del test
t_values = zeros(num_timepoints, num_conditions * (num_conditions - 1) / 2);

% Calcola la statistica del test per ogni coppia di condizioni
counter = 1;
for condition1 = 1:(num_conditions - 1)
    for condition2 = (condition1 + 1):num_conditions
        [~, p, ~, stats] = ttest(data_all(:, :, condition1) - data_all(:, :, condition2));
        t_values(:, counter) = stats.tstat;
        counter = counter + 1;
    end
end

alpha = 0.05;  % livello di significatività
df = num_subjects-1;  % gradi di libertà
t_critical = tinv(1-alpha/2, df);  % valore t critico per un test a due code
% Applica la soglia per definire i clusters
threshold = t_critical; % questo dipende dal tuo numero di soggetti e dal tuo livello di significatività

clusters_positive = t_values > threshold;
clusters_negative = t_values < -threshold;

% Calcola la dimensione di ciascun cluster
cluster_sizes_positive = sum(clusters_positive, 1);
cluster_sizes_negative = sum(clusters_negative, 1);

% Calcola la dimensione di ciascun cluster permutato
num_permutations = 10000;
permuted_cluster_sizes_positive = zeros(num_permutations, 1);
permuted_cluster_sizes_negative = zeros(num_permutations, 1);
for permutation = 1:num_permutations
    % Permuta le etichette delle condizioni
    permuted_data = data_all;
    for subject = 1:num_subjects
        permuted_data(subject, :, :) = permuted_data(subject, :, randperm(num_conditions));
    end
    
    % statistica del test per ogni coppia di condizioni
    counter = 1;
    permuted_t_values = zeros(num_timepoints, num_conditions * (num_conditions - 1) / 2);
    for condition1 = 1:(num_conditions - 1)
        for condition2 = (condition1 + 1):num_conditions
            [~, ~, ~, stats] = ttest(permuted_data(:, :, condition1) - permuted_data(:, :, condition2));
            permuted_t_values(:, counter) = stats.tstat;
            counter = counter + 1;
        end
    end

    % Applica la soglia per definire i clusters
    permuted_clusters_positive = permuted_t_values > threshold;
    permuted_clusters_negative = permuted_t_values < -threshold;
    
    % dimensione del cluster più grande
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

%%

% p-values
disp('P-values for positive clusters:');
disp(p_values_positive);
disp('P-values for negative clusters:');
disp(p_values_negative);

% statistica media del test attraverso le coppie di condizioni
mean_t_values = mean(t_values, 2);

% Crea un nuovo grafico
figure;

% statistica media del test nel tempo
plot(time_axis_corrected, mean_t_values, 'LineWidth', 1.5);
hold on;

% linea orizzontale alla soglia
plot(time_axis_corrected, threshold * ones(size(time_axis_corrected)), '--k');
plot(time_axis_corrected, -threshold * ones(size(time_axis_corrected)), '--k');

% indici dei clusters significativi
significant_clusters_positive = find(p_values_positive < 0.05);
significant_clusters_negative = find(p_values_negative < 0.05);

% Per ogni cluster significativo colora l'area tra il grafico e la soglia
for cluster = significant_clusters_positive
    start_time = time_axis_corrected(clusters_positive(cluster, 1));
    end_time = time_axis_corrected(clusters_positive(cluster, end));
    fill([start_time end_time end_time start_time], [threshold threshold mean_t_values(cluster) mean_t_values(cluster)], 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

for cluster = significant_clusters_negative
    start_time = time_axis_corrected(clusters_negative(cluster, 1));
    end_time = time_axis_corrected(clusters_negative(cluster, end));
    fill([start_time end_time end_time start_time], [-threshold -threshold mean_t_values(cluster) mean_t_values(cluster)], 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none');
end

% assegnamenti degli assi
xlabel('Time (s)');
ylabel('Average t-value');
title('Time-course of the t-test statistic');

