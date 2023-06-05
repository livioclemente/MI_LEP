% Initialize parameters
preclust_pval = 0.1;  % Set this as needed
clust_pval = 0.1;  % Set this as needed
nperm = 1000; % Adjust this to your needs

% Inizializzazione delle variabili
num_subjects_ct = size(Iint_values_ct, 1);
num_subjects_pth = size(Iint_values_pth, 1);
num_subjects_dlpth = size(Iint_values_dlpth, 1);
num_subjects_norm = size(Iint_values_norm, 1);
Iint_values_tot = cell(num_subjects_ct + num_subjects_pth + num_subjects_dlpth + num_subjects_norm, 1);

% Combinazione dei dati dei quattro gruppi in Iint_values_tot
Iint_values_tot(1:num_subjects_ct) = Iint_values_ct;
Iint_values_tot(num_subjects_ct+1:num_subjects_ct+num_subjects_pth) = Iint_values_pth;
Iint_values_tot(num_subjects_ct+num_subjects_pth+1:num_subjects_ct+num_subjects_pth+num_subjects_dlpth) = Iint_values_dlpth;
Iint_values_tot(num_subjects_ct+num_subjects_pth+num_subjects_dlpth+1:end) = Iint_values_norm;


% Inizializzazione delle matrici e variabili
h0_diff = zeros(Nt, Nt);
hp_diff = zeros(Nt, Nt, nperm);

% Calcolo dell'ipotesi nulla (differenza media reale tra i gruppi)
for t1 = 1:Nt
    for t2 = (t1 + 1):Nt
        h0_diff(t1, t2) = mean(group_Iint_values_pth(t1, t2), 'all') - mean(group_Iint_values_ct(t1, t2), 'all') + mean(group_Iint_values_dlpth(t1, t2), 'all') - mean(group_Iint_values_norm(t1, t2), 'all');
    end
end

% Inizializzazione delle matrici per le permutazioni
num_subjects_total = num_subjects_ct + num_subjects_pth + num_subjects_dlpth + num_subjects_norm;
Iint_perm_ct_3D = zeros(Nt, Nt, num_subjects_ct);
Iint_perm_pth_3D = zeros(Nt, Nt, num_subjects_pth);
Iint_perm_dlpth_3D = zeros(Nt, Nt, num_subjects_dlpth);
Iint_perm_norm_3D = zeros(Nt, Nt, num_subjects_norm);

% Permutazioni
for iperm = 1:nperm
    % Permuta l'ordine dei soggetti in tutti e quattro i gruppi
    labels_perm = randperm(num_subjects_total);

    % Divide i dati permutati nei quattro gruppi
    Iint_perm_ct = Iint_values_tot(labels_perm(1:num_subjects_ct));
    Iint_perm_pth = Iint_values_tot(labels_perm(num_subjects_ct+1:num_subjects_ct+num_subjects_pth));
    Iint_perm_dlpth = Iint_values_tot(labels_perm(num_subjects_ct+num_subjects_pth+1:num_subjects_ct+num_subjects_pth+num_subjects_dlpth));
    Iint_perm_norm = Iint_values_tot(labels_perm(num_subjects_ct+num_subjects_pth+num_subjects_dlpth+1:end));

    % Assegna i dati permutati alle matrici 3D
    for i = 1:num_subjects_ct
        Iint_perm_ct_3D(:, :, i) = Iint_perm_ct{i};
    end
    for i = 1:num_subjects_pth
        Iint_perm_pth_3D(:, :, i) = Iint_perm_pth{i};
    end
    for i = 1:num_subjects_dlpth
        Iint_perm_dlpth_3D(:, :, i) = Iint_perm_dlpth{i};
    end
    for i = 1:num_subjects_norm
        Iint_perm_norm_3D(:, :, i) = Iint_perm_norm{i};
    end

    % Calcola la differenza media per ciascuna permutazione
    for t1 = 1:Nt
        for t2 = t1 + 1:Nt
            hp_diff(t1, t2, iperm) = mean(Iint_perm_pth_3D(t1, t2, :), 'all') - mean(Iint_perm_ct_3D(t1, t2, :), 'all') + mean(Iint_perm_dlpth_3D(t1, t2, :), 'all') - mean(Iint_perm_norm_3D(t1, t2, :), 'all');
        end
    end
end

% Initialize parameters
preclust_pval = 0.05;  % Set this as needed
clust_pval = 0.05;  % Set this as needed

% Collect the largest suprathreshold clusters
clust_max = zeros(nperm, 1);
for iperm = 1:nperm
    perms = true(1,nperm);
    perms(iperm) = 0;
    zvals = squeeze((hp_diff(:,:,iperm) - mean(hp_diff(:,:,perms),3)) ./ std(hp_diff(:,:,perms),[],3));
    zvals(abs(zvals) < norminv(1 - preclust_pval)) = 0; 
    
    clust_info = bwconncomp(zvals);
    clust_max(iperm) = max([0 cellfun(@numel, clust_info.PixelIdxList)]);
end


h0_diff = h0_diff + h0_diff';
for iperm = 1:nperm
    hp_diff(:,:,iperm) = hp_diff(:,:,iperm) + hp_diff(:,:,iperm)';
end

% Identify significant clusters in real data
zmap = squeeze((h0_diff - mean(hp_diff, 3)) ./ std(hp_diff, [], 3));
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

% Convert Nt to time in milliseconds assuming your sampling rate
time_plot = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, Nt);

figure; 
imagesc(time_plot, time_plot, h0_diff);
hold on;

% Plot significant clusters
contour(time_plot, time_plot, zmap, 1, 'linecolor', 'k', 'LineWidth', 1);

colormap(brewermap([], '*RdBu'));
colorbar;
axis square;

% Adjust color limits to center around zero
lim = max(abs(min(min(h0_diff))), max(max(h0_diff)));
caxis([-lim, lim]);

xlabel('Time (ms)');
ylabel('Time (ms)');
title('Interaction Information with Cluster-based Correction');
