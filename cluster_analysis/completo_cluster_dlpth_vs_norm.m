% Load groups info
load("group_dlpth.mat");
load("group_norm.mat");

% Initialization
Nt1 = size(x_dlpth, 2);
Nt2 = size(x_norm, 2);
Nt_cluster = min(Nt1, Nt2);  % If x_pth1 and x_pth2 have different numbers of time points, use the smaller number

nperm = 1000;  % Adjust this to your needs

% Get data
x1 = x_dlpth;
x2 = x_norm;
stim1 = stim_dlpth;
stim2 = stim_norm;
n1 = size(x1, 1);
n2 = size(x2, 1);

% Calculate interaction information for each group
h0_dlpth = zeros(Nt_cluster, Nt_cluster);
h0_norm = zeros(Nt_cluster, Nt_cluster);

for t1 = 1:Nt_cluster
    for t2 = (t1+1):Nt_cluster
        JMI_tot1 = mi_gg([x1(:, t1) x1(:, t2)], stim1, false);
        II_tot1 = JMI_tot1 - mi_gg(x1(:, t1), stim1, false) - mi_gg(x1(:, t2), stim1, false);
        h0_dlpth(t1, t2) = II_tot1;

        JMI_tot2 = mi_gg([x2(:, t1) x2(:, t2)], stim2, false);
        II_tot2 = JMI_tot2 - mi_gg(x2(:, t1), stim2, false) - mi_gg(x2(:, t2), stim2, false);
        h0_norm(t1, t2) = II_tot2;
    end
end

% Calculate the difference between the groups
diff_h0 = h0_dlpth - h0_norm;

% Perform permutation test
diff_perm_t_values = zeros(Nt_cluster, Nt_cluster, nperm);
for iperm = 1:nperm
    % Combine and permute the data
    combined_x = [x1; x2];
    combined_stim = [stim1; stim2];
    perm_indices = randperm(n1 + n2);
    
    % Divide the permuted data back into two groups
    perm_x1 = combined_x(perm_indices(1:n1), :);
    perm_stim1 = combined_stim(perm_indices(1:n1));
    perm_x2 = combined_x(perm_indices(n1+1:end), :);
    perm_stim2 = combined_stim(perm_indices(n1+1:end));
    
    % Calculate interaction information for permuted groups
    h0_perm1 = zeros(Nt_cluster, Nt_cluster);
    h0_perm2 = zeros(Nt_cluster, Nt_cluster);
    for t1 = 1:Nt_cluster
        for t2 = (t1+1):Nt_cluster
            JMI_perm1 = mi_gg([perm_x1(:, t1) perm_x1(:, t2)], perm_stim1, false);
            II_perm1 = JMI_perm1 - mi_gg(perm_x1(:, t1), perm_stim1, false) - mi_gg(perm_x1(:, t2), perm_stim1, false);
            h0_perm1(t1, t2) = II_perm1;

            JMI_perm2 = mi_gg([perm_x2(:, t1) perm_x2(:, t2)], perm_stim2, false);
            II_perm2 = JMI_perm2 - mi_gg(perm_x2(:, t1), perm_stim2, false) - mi_gg(perm_x2(:, t2), perm_stim2, false);
            h0_perm2(t1, t2) = II_perm2;
        end
    end

    % Calculate the difference between permuted groups
    diff_perm = h0_perm1 - h0_perm2;
    diff_perm_t_values(:,:,iperm) = diff_perm;
end

% Fill in the other half of the matrices
diff_h0 = diff_h0 + diff_h0';
for iperm = 1:nperm 
    diff_perm_t_values(:,:,iperm) = diff_perm_t_values(:,:,iperm) + diff_perm_t_values(:,:,iperm)';
end 

% Cluster-based correction
preclust_pval = 0.05;  % Set this as needed
clust_pval = 0.05;  % Set this as needed

% Collect the largest suprathreshold clusters
clust_max = zeros(nperm, 1);
for iperm = 1:nperm
    perms = true(1,nperm);
    perms(iperm) = 0;
    zvals = squeeze((diff_perm_t_values(:,:,iperm) - mean(diff_perm_t_values(:,:,perms),3)) ./ std(diff_perm_t_values(:,:,perms),[],3));
    zvals(abs(zvals) < norminv(1 - preclust_pval)) = 0; 
    
    clust_info = bwconncomp(zvals);
    clust_max(iperm) = max([0 cellfun(@numel, clust_info.PixelIdxList)]);
end

% Identify significant clusters in real data
zmap = squeeze((diff_h0 - mean(diff_perm_t_values, 3)) ./ std(diff_perm_t_values, [], 3));
zmap(abs(zmap) < norminv(1 - preclust_pval)) = 0;

clust_info = bwconncomp(zmap);
clust_size = cellfun(@numel, clust_info.PixelIdxList);
clust_th = prctile(clust_max, 100 - clust_pval * 100);
clust_rem = find(clust_size < clust_th);
for i = 1:length(clust_rem)
    zmap(clust_info.PixelIdxList{clust_rem(i)}) = 0;
end
zmap(isnan(zmap)) = 0;
zmap = logical(zmap);

% Save file
save('dlpth_vs_norm_results.mat', 'diff_h0', 'diff_perm_t_values', 'clust_max', 'zmap');

% Data creation
erp_dlpth = cell(num_subjects_dlpth, 1);
for s = 1:num_subjects_dlpth
    erp_dlpth{s} = erp_mean_dlpth{s}(11, :);
end
erp_dlpth = cell2mat(erp_dlpth);

erp_norm = cell(num_subjects_norm, 1);
for s = 1:num_subjects_norm
    erp_norm{s} = erp_mean_norm{s}(11, :);
end
erp_norm = cell2mat(erp_norm);

% Number of permutations
num_perms = 1000;
% Alpha value
alpha = 0.05;

% Calcola i valori t per ogni timepoint tra data_ct e data_pth
[~, ~, ~, stats] = ttest2(erp_dlpth, erp_norm);
t_values_erp_dlpthvsnorm = stats.tstat;

% Permuta i dati e calcola i valori t per ogni permutazione
perm_t_values_erp_dlpthvsnorm = zeros(num_perms, Nt);
for p = 1:num_perms
    perm_data = [erp_dlpth; erp_norm];
    perm_data = perm_data(randperm(size(perm_data, 1)), :);
    perm_data_dlpth = perm_data(1:size(erp_dlpth, 1), :);
    perm_data_norm = perm_data(size(erp_dlpth, 1)+1:end, :);
    
    for sample = 1:Nt
        [~, ~, ~, stats] = ttest2(perm_data_dlpth(:, sample), perm_data_norm(:, sample));
        perm_t_values_erp_dlpthvsnorm(p, sample) = stats.tstat;
    end
end

% Calcola i p-value
p_values_erp_dlpthvsnorm = sum(perm_t_values_erp_dlpthvsnorm > repmat(t_values_erp_dlpthvsnorm, num_perms, 1)) / num_perms;

% Confronta i p-value con il livello alpha
significantly_different_erp_dlpthvsnorm = p_values_erp_dlpthvsnorm < alpha;

% Stampa i risultati
for sample = 1:Nt
    if significantly_different_erp_dlpthvsnorm(sample)
        fprintf('Il campione %d è significativamente diverso tra le due condizioni (p = %.4f).\n', sample, p_values_erp_dlpthvsnorm(sample));
    else
        fprintf('Il campione %d non è significativamente diverso tra le due condizioni (p = %.4f).\n', sample, p_values_erp_dlpthvsnorm(sample));
    end
end

% Soglia per l'identificazione dei cluster - qui sto usando la soglia che hai calcolato precedentemente
threshold_erp_dlpthvsnorm = quantile(perm_t_values_erp_dlpthvsnorm(:), 0.975);

% Identifica i cluster di campioni adiacenti per i quali il valore t supera la soglia
clusters_pos_erp_dlpthvsnorm = bwconncomp(t_values_erp_dlpthvsnorm > threshold_erp_dlpthvsnorm);
clusters_neg_erp_dlpthvsnorm = bwconncomp(t_values_erp_dlpthvsnorm < -threshold_erp_dlpthvsnorm);

% Calcola la statistica del cluster per ogni cluster come la somma dei valori t nel cluster
masses_pos_erp_dlpthvsnorm = cellfun(@(I) sum(t_values_erp_dlpthvsnorm(I)), clusters_pos_erp_dlpthvsnorm.PixelIdxList);
masses_neg_erp_dlpthvsnorm = cellfun(@(I) sum(abs(t_values_erp_dlpthvsnorm(I))), clusters_neg_erp_dlpthvsnorm.PixelIdxList);

% Trova il cluster positivo/negativo più grande
[max_mass_pos_erp_dlpthvsnorm, ~] = max(masses_pos_erp_dlpthvsnorm);
[max_mass_neg_erp_dlpthvsnorm, ~] = max(masses_neg_erp_dlpthvsnorm);

% Assicurati che time_axis_corrected_ct abbia la stessa lunghezza di t_values
assert(length(time_axis_corrected) == length(t_values_erp_dlpthvsnorm), 'Time axis and t values must have the same length');

% Trova l'indice del cluster positivo/negativo più grande
[~, idx_max_mass_pos_erp_dlpthvsnorm] = max(masses_pos_erp_dlpthvsnorm);
[~, idx_max_mass_neg_erp_dlpthvsnorm] = max(masses_neg_erp_dlpthvsnorm);

% Save variables for plot
save('ERP_cluster_dlpth_vs_norm.mat', 'time_axis_corrected', 't_values_erp_dlpthvsnorm', 'threshold_erp_dlpthvsnorm', 'clusters_pos_erp_dlpthvsnorm', 'clusters_neg_erp_dlpthvsnorm', 'idx_max_mass_pos_erp_dlpthvsnorm', 'idx_max_mass_neg_erp_dlpthvsnorm');

% Data creation
MI_dlpth = cell(num_subjects_dlpth, 1);
for s = 1:num_subjects_dlpth
    MI_dlpth{s} = I_values_dlpth{s};
end
MI_dlpth = cell2mat(MI_dlpth);

MI_norm = cell(num_subjects_norm, 1);
for s = 1:num_subjects_norm
    MI_norm{s} = I_values_norm{s};
end
MI_norm = cell2mat(MI_norm);

% Number of permutations
num_perms = 1000;
% Alpha value
alpha = 0.05;

% Calcola i valori t per ogni timepoint tra data_ct e data_pth
[~, ~, ~, stats] = ttest2(MI_dlpth, MI_norm);
t_values_MI_dlpthvsnorm = stats.tstat;

% Permuta i dati e calcola i valori t per ogni permutazione
perm_t_values_MI_dlpthvsnorm = zeros(num_perms, Nt);
for p = 1:num_perms
    perm_data = [MI_dlpth; MI_norm];
    perm_data = perm_data(randperm(size(perm_data, 1)), :);
    perm_data_dlpth = perm_data(1:size(MI_dlpth, 1), :);
    perm_data_norm = perm_data(size(MI_dlpth, 1)+1:end, :);
    
    for sample = 1:Nt
        [~, ~, ~, stats] = ttest2(perm_data_dlpth(:, sample), perm_data_norm(:, sample));
        perm_t_values_MI_dlpthvsnorm(p, sample) = stats.tstat;
    end
end

% Calcola i p-value
p_values_MI_dlpthvsnorm = sum(perm_t_values_MI_dlpthvsnorm > repmat(t_values_MI_dlpthvsnorm, num_perms, 1)) / num_perms;

% Confronta i p-value con il livello alpha
significantly_different_MI_dlpthvsnorm = p_values_MI_dlpthvsnorm < alpha;

% Stampa i risultati
for sample = 1:Nt
    if significantly_different_MI_dlpthvsnorm(sample)
        fprintf('Il campione %d è significativamente diverso tra le due condizioni (p = %.4f).\n', sample, p_values_MI_dlpthvsnorm(sample));
    else
        fprintf('Il campione %d non è significativamente diverso tra le due condizioni (p = %.4f).\n', sample, p_values_MI_dlpthvsnorm(sample));
    end
end

% Soglia per l'identificazione dei cluster - qui sto usando la soglia che hai calcolato precedentemente
threshold_MI_dlpthvsnorm = quantile(perm_t_values_MI_dlpthvsnorm(:), 0.975);

% Identifica i cluster di campioni adiacenti per i quali il valore t supera la soglia
clusters_pos_MI_dlpthvsnorm = bwconncomp(t_values_MI_dlpthvsnorm > threshold_MI_dlpthvsnorm);
clusters_neg_MI_dlpthvsnorm = bwconncomp(t_values_MI_dlpthvsnorm < -threshold_MI_dlpthvsnorm);

% Calcola la statistica del cluster per ogni cluster come la somma dei valori t nel cluster
masses_pos_MI_dlpthvsnorm = cellfun(@(I) sum(t_values_MI_dlpthvsnorm(I)), clusters_pos_MI_dlpthvsnorm.PixelIdxList);
masses_neg_MI_dlpthvsnorm = cellfun(@(I) sum(abs(t_values_MI_dlpthvsnorm(I))), clusters_neg_MI_dlpthvsnorm.PixelIdxList);

% Trova il cluster positivo/negativo più grande
[max_mass_pos_MI_dlpthvsnorm, ~] = max(masses_pos_MI_dlpthvsnorm);
[max_mass_neg_MI_dlpthvsnorm, ~] = max(masses_neg_MI_dlpthvsnorm);

% Assicurati che time_axis_corrected_ct abbia la stessa lunghezza di t_values
assert(length(time_axis_corrected) == length(t_values_MI_dlpthvsnorm), 'Time axis and t values must have the same length');

% Trova l'indice del cluster positivo/negativo più grande
[~, idx_max_mass_pos_MI_dlpthvsnorm] = max(masses_pos_MI_dlpthvsnorm);
[~, idx_max_mass_neg_MI_dlpthvsnorm] = max(masses_neg_MI_dlpthvsnorm);

% Save variables for plot
save('MI_cluster_dlpth_vs_norm.mat', 'time_axis_corrected', 't_values_MI_dlpthvsnorm', 'threshold_MI_dlpthvsnorm', 'clusters_pos_MI_dlpthvsnorm', 'clusters_neg_MI_dlpthvsnorm', 'idx_max_mass_pos_MI_dlpthvsnorm', 'idx_max_mass_neg_MI_dlpthvsnorm');

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
plot(time_plot, t_values_erp_dlpthvsnorm, 'k');
hold on;
title('Valori t per ogni punto temporale');
xlabel('Tempo');
ylabel('Valore t');

line([min(time_axis_corrected), max(time_axis_corrected)], [0, 0], 'Color', 'k');
line([min(time_axis_corrected), max(time_axis_corrected)], [threshold_erp_dlpthvsnorm, threshold_erp_dlpthvsnorm], 'Color', 'k', 'LineStyle', '--');
line([min(time_axis_corrected), max(time_axis_corrected)], [-threshold_erp_dlpthvsnorm, -threshold_erp_dlpthvsnorm], 'Color', 'k', 'LineStyle', '--');

try
    [~, idx_max_mass_pos_erp_dlpthvsnorm] = max(masses_pos_erp_dlpthvsnorm);
    for idx = clusters_pos_erp_dlpthvsnorm.PixelIdxList{idx_max_mass_pos_erp_dlpthvsnorm}
        area(time_axis_corrected(idx), t_values_erp_dlpthvsnorm(idx), 'FaceColor', 'r', 'EdgeColor', 'none');
    end
catch
    disp('No positive clusters found');
end

try
    [~, idx_max_mass_neg_erp_dlpthvsnorm] = max(masses_neg_erp_dlpthvsnorm);
    for idx = clusters_neg_erp_dlpthvsnorm.PixelIdxList{idx_max_mass_neg_erp_dlpthvsnorm}
        area(time_axis_corrected(idx), t_values_erp_dlpthvsnorm(idx), 'FaceColor', 'b', 'EdgeColor', 'none');
    end
catch
    disp('No negative clusters found');
end

% Imposta i limiti dell'asse x
xlim([-100, 1000]);

hold off;

% Third subplot: Mutual Information as a function of time
subplot(5, 5, [1 6 11 16]);
plot(time_axis_corrected, t_values_MI_dlpthvsnorm, 'k');
hold on;
xlabel('Tempo');
ylabel('Valore t');

line([min(time_axis_corrected), max(time_axis_corrected)], [0, 0], 'Color', 'k');
line([min(time_axis_corrected), max(time_axis_corrected)], [threshold_MI_dlpthvsnorm, threshold_MI_dlpthvsnorm], 'Color', 'k', 'LineStyle', '--');
line([min(time_axis_corrected), max(time_axis_corrected)], [-threshold_MI_dlpthvsnorm, -threshold_MI_dlpthvsnorm], 'Color', 'k', 'LineStyle', '--');

try
    [~, idx_max_mass_pos_MI_dlpthvsnorm] = max(masses_pos_MI_dlpthvsnorm);
    for idx = clusters_pos_MI_dlpthvsnorm.PixelIdxList{idx_max_mass_pos_MI_dlpthvsnorm}
        area(time_axis_corrected(idx), t_values_MI_dlpthvsnorm(idx), 'FaceColor', 'r', 'EdgeColor', 'none');
    end
catch
    disp('No positive clusters found');
end

try
    [~, idx_max_mass_neg_MI_dlpthvsnorm] = max(masses_neg_MI_dlpthvsnorm);
    for idx = clusters_neg_MI_dlpthvsnorm.PixelIdxList{idx_max_mass_neg_MI_dlpthvsnorm}
        area(time_axis_corrected(idx), t_values_MI_dlpthvsnorm(idx), 'FaceColor', 'b', 'EdgeColor', 'none');
    end
catch
    disp('No negative clusters found');
end

% Imposta i limiti dell'asse x
xlim([-100, 1000]);
set(gca, 'CameraUpVector', [-1 0 0]);
hold off;
