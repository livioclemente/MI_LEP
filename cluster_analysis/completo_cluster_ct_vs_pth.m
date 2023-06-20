% Load groups info
load("group_ct.mat");
load("group_pth.mat");
load("ct_vs_pth_results.mat");

% Data creation
erp_ct = cell(num_subjects_ct, 1);
for s = 1:num_subjects_ct
    erp_ct{s} = erp_mean_ct{s}(11, :);
end
erp_ct = cell2mat(erp_ct);

erp_pth = cell(num_subjects_pth, 1);
for s = 1:num_subjects_pth
    erp_pth{s} = erp_mean_pth{s}(11, :);
end
erp_pth = cell2mat(erp_pth);

% Number of permutations
num_perms = 1000;
% Alpha value
alpha = 0.05;

% Calcola i valori t per ogni timepoint tra data_ct e data_pth
[~, ~, ~, stats] = ttest2(erp_ct, erp_pth);
t_values_erp_ctvspth = stats.tstat;

% Permuta i dati e calcola i valori t per ogni permutazione
perm_t_values_erp_ctvspth = zeros(num_perms, Nt);
for p = 1:num_perms
    perm_data = [erp_ct; erp_pth];
    perm_data = perm_data(randperm(size(perm_data, 1)), :);
    perm_data_ct = perm_data(1:size(erp_ct, 1), :);
    perm_data_pth = perm_data(size(erp_ct, 1)+1:end, :);
    
    for sample = 1:Nt
        [~, ~, ~, stats] = ttest2(perm_data_ct(:, sample), perm_data_pth(:, sample));
        perm_t_values_erp_ctvspth(p, sample) = stats.tstat;
    end
end

% Calcola i p-value
p_values_erp_ctvspth = sum(perm_t_values_erp_ctvspth > repmat(t_values_erp_ctvspth, num_perms, 1)) / num_perms;

% Confronta i p-value con il livello alpha
significantly_different_erp_ctvspth = p_values_erp_ctvspth < alpha;

% Stampa i risultati
for sample = 1:Nt
    if significantly_different_erp_ctvspth(sample)
        fprintf('Il campione %d è significativamente diverso tra le due condizioni (p = %.4f).\n', sample, p_values_erp_ctvspth(sample));
    else
        fprintf('Il campione %d non è significativamente diverso tra le due condizioni (p = %.4f).\n', sample, p_values_erp_ctvspth(sample));
    end
end

% Soglia per l'identificazione dei cluster - qui sto usando la soglia che hai calcolato precedentemente
threshold_erp_ctvspth = quantile(perm_t_values_erp_ctvspth(:), 0.975);

% Identifica i cluster di campioni adiacenti per i quali il valore t supera la soglia
clusters_pos_erp_ctvspth = bwconncomp(t_values_erp_ctvspth > threshold_erp_ctvspth);
clusters_neg_erp_ctvspth = bwconncomp(t_values_erp_ctvspth < -threshold_erp_ctvspth);

% Calcola la statistica del cluster per ogni cluster come la somma dei valori t nel cluster
masses_pos_erp_ctvspth = cellfun(@(I) sum(t_values_erp_ctvspth(I)), clusters_pos_erp_ctvspth.PixelIdxList);
masses_neg_erp_ctvspth = cellfun(@(I) sum(abs(t_values_erp_ctvspth(I))), clusters_neg_erp_ctvspth.PixelIdxList);

% Trova il cluster positivo/negativo più grande
[max_mass_pos_erp_ctvspth, ~] = max(masses_pos_erp_ctvspth);
[max_mass_neg_erp_ctvspth, ~] = max(masses_neg_erp_ctvspth);

% Assicurati che time_axis_corrected_ct abbia la stessa lunghezza di t_values
assert(length(time_axis_corrected) == length(t_values_erp_ctvspth), 'Time axis and t values must have the same length');

% Trova l'indice del cluster positivo/negativo più grande
[~, idx_max_mass_pos_erp_ctvspth] = max(masses_pos_erp_ctvspth);
[~, idx_max_mass_neg_erp_ctvspth] = max(masses_neg_erp_ctvspth);

% Save variables for plot
save('ERP_cluster_ct_vs_pth.mat', 'time_axis_corrected', 't_values_erp_ctvspth', 'threshold_erp_ctvspth', 'clusters_pos_erp_ctvspth', 'clusters_neg_erp_ctvspth', 'idx_max_mass_pos_erp_ctvspth', 'idx_max_mass_neg_erp_ctvspth');

% Data creation
MI_ct = cell(num_subjects_ct, 1);
for s = 1:num_subjects_ct
    MI_ct{s} = I_values_ct{s};
end
MI_ct = cell2mat(MI_ct);

MI_pth = cell(num_subjects_pth, 1);
for s = 1:num_subjects_pth
    MI_pth{s} = I_values_pth{s};
end
MI_pth = cell2mat(MI_pth);

% Number of permutations
num_perms = 1000;
% Alpha value
alpha = 0.05;

% Calcola i valori t per ogni timepoint tra data_ct e data_pth
[~, ~, ~, stats] = ttest2(MI_ct, MI_pth);
t_values_MI_ctvspth = stats.tstat;

% Permuta i dati e calcola i valori t per ogni permutazione
perm_t_values_MI_ctvspth = zeros(num_perms, Nt);
for p = 1:num_perms
    perm_data = [MI_ct; MI_pth];
    perm_data = perm_data(randperm(size(perm_data, 1)), :);
    perm_data_ct = perm_data(1:size(MI_ct, 1), :);
    perm_data_pth = perm_data(size(MI_ct, 1)+1:end, :);
    
    for sample = 1:Nt
        [~, ~, ~, stats] = ttest2(perm_data_ct(:, sample), perm_data_pth(:, sample));
        perm_t_values_MI_ctvspth(p, sample) = stats.tstat;
    end
end

% Calcola i p-value
p_values_MI_ctvspth = sum(perm_t_values_MI_ctvspth > repmat(t_values_MI_ctvspth, num_perms, 1)) / num_perms;

% Confronta i p-value con il livello alpha
significantly_different_MI_ctvspth = p_values_MI_ctvspth < alpha;

% Stampa i risultati
for sample = 1:Nt
    if significantly_different_MI_ctvspth(sample)
        fprintf('Il campione %d è significativamente diverso tra le due condizioni (p = %.4f).\n', sample, p_values_MI_ctvspth(sample));
    else
        fprintf('Il campione %d non è significativamente diverso tra le due condizioni (p = %.4f).\n', sample, p_values_MI_ctvspth(sample));
    end
end

% Soglia per l'identificazione dei cluster - qui sto usando la soglia che hai calcolato precedentemente
threshold_MI_ctvspth = quantile(perm_t_values_MI_ctvspth(:), 0.975);

% Identifica i cluster di campioni adiacenti per i quali il valore t supera la soglia
clusters_pos_MI_ctvspth = bwconncomp(t_values_MI_ctvspth > threshold_MI_ctvspth);
clusters_neg_MI_ctvspth = bwconncomp(t_values_MI_ctvspth < -threshold_MI_ctvspth);

% Calcola la statistica del cluster per ogni cluster come la somma dei valori t nel cluster
masses_pos_MI_ctvspth = cellfun(@(I) sum(t_values_MI_ctvspth(I)), clusters_pos_MI_ctvspth.PixelIdxList);
masses_neg_MI_ctvspth = cellfun(@(I) sum(abs(t_values_MI_ctvspth(I))), clusters_neg_MI_ctvspth.PixelIdxList);

% Trova il cluster positivo/negativo più grande
[max_mass_pos_MI_ctvspth, ~] = max(masses_pos_MI_ctvspth);
[max_mass_neg_MI_ctvspth, ~] = max(masses_neg_MI_ctvspth);

% Assicurati che time_axis_corrected_ct abbia la stessa lunghezza di t_values
assert(length(time_axis_corrected) == length(t_values_MI_ctvspth), 'Time axis and t values must have the same length');

% Trova l'indice del cluster positivo/negativo più grande
[~, idx_max_mass_pos_MI_ctvspth] = max(masses_pos_MI_ctvspth);
[~, idx_max_mass_neg_MI_ctvspth] = max(masses_neg_MI_ctvspth);

% Save variables for plot
save('MI_cluster_ct_vs_pth.mat', 'time_axis_corrected', 't_values_MI_ctvspth', 'threshold_MI_ctvspth', 'clusters_pos_MI_ctvspth', 'clusters_neg_MI_ctvspth', 'idx_max_mass_pos_MI_ctvspth', 'idx_max_mass_neg_MI_ctvspth');

% General info
t_start = -0.1;
t_end = 1;
Fs = 256;
time_plot = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, Nt);

% Create a figure with three subplots
figure('Position', [0, 0, 1200, 800])

% Main subplot: Mutual Information
axm = subplot(5, 5, [2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
imagesc(time_plot, time_plot, diff_h0);
hold on;
contour(time_plot, time_plot, zmap, 1, 'linecolor', 'k', 'LineWidth', 1);
colormap(brewermap([], '*RdBu'));
colorbar;
axis square;
lim = max(abs(min(min(diff_h0))), max(max(diff_h0)));
clim([-lim, lim]);
xlabel('Time (ms)');
ylabel('Time (ms)');
title('Difference in Interaction Information between Groups with Cluster-based Correction');

% Second subplot: Average ERP of all conditions
subplot(5, 5, [22.45 23 24 24.34]);
plot(time_plot, t_values_erp_ctvspth, 'k');
hold on;
title('Valori t per ogni punto temporale');
xlabel('Tempo');
ylabel('Valore t');

line([min(time_axis_corrected), max(time_axis_corrected)], [0, 0], 'Color', 'k');
line([min(time_axis_corrected), max(time_axis_corrected)], [threshold_erp_ctvspth, threshold_erp_ctvspth], 'Color', 'k', 'LineStyle', '--');
line([min(time_axis_corrected), max(time_axis_corrected)], [-threshold_erp_ctvspth, -threshold_erp_ctvspth], 'Color', 'k', 'LineStyle', '--');

try
    [~, idx_max_mass_pos_erp_ctvspth] = max(masses_pos_erp_ctvspth);
    for idx = clusters_pos_erp_ctvspth.PixelIdxList{idx_max_mass_pos_erp_ctvspth}
        area(time_axis_corrected(idx), t_values_erp_ctvspth(idx), 'FaceColor', 'r', 'EdgeColor', 'none');
    end
catch
    disp('No positive clusters found');
end

try
    [~, idx_max_mass_neg_erp_ctvspth] = max(masses_neg_erp_ctvspth);
    for idx = clusters_neg_erp_ctvspth.PixelIdxList{idx_max_mass_neg_erp_ctvspth}
        area(time_axis_corrected(idx), t_values_erp_ctvspth(idx), 'FaceColor', 'b', 'EdgeColor', 'none');
    end
catch
    disp('No negative clusters found');
end

% Imposta i limiti dell'asse x
xlim([-100, 1000]);

hold off;

% Third subplot: Mutual Information as a function of time
subplot(5, 5, [1 6 11 16]);
plot(time_axis_corrected, t_values_MI_ctvspth, 'k');
hold on;
xlabel('Tempo');
ylabel('Valore t');

line([min(time_axis_corrected), max(time_axis_corrected)], [0, 0], 'Color', 'k');
line([min(time_axis_corrected), max(time_axis_corrected)], [threshold_MI_ctvspth, threshold_MI_ctvspth], 'Color', 'k', 'LineStyle', '--');
line([min(time_axis_corrected), max(time_axis_corrected)], [-threshold_MI_ctvspth, -threshold_MI_ctvspth], 'Color', 'k', 'LineStyle', '--');

try
    [~, idx_max_mass_pos_MI_ctvspth] = max(masses_pos_MI_ctvspth);
    for idx = clusters_pos_MI_ctvspth.PixelIdxList{idx_max_mass_pos_MI_ctvspth}
        area(time_axis_corrected(idx), t_values_MI_ctvspth(idx), 'FaceColor', 'r', 'EdgeColor', 'none');
    end
catch
    disp('No positive clusters found');
end

try
    [~, idx_max_mass_neg_MI_ctvspth] = max(masses_neg_MI_ctvspth);
    for idx = clusters_neg_MI_ctvspth.PixelIdxList{idx_max_mass_neg_MI_ctvspth}
        area(time_axis_corrected(idx), t_values_MI_ctvspth(idx), 'FaceColor', 'b', 'EdgeColor', 'none');
    end
catch
    disp('No negative clusters found');
end

% Imposta i limiti dell'asse x
xlim([-100, 1000]);
set(gca, 'CameraUpVector', [-1 0 0]);
hold off;
