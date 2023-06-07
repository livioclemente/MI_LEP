% Initialization
Nt1 = size(x_ct, 2);
Nt2 = size(x_pth, 2);
Nt_cluster = min(Nt1, Nt2);  % If x_pth1 and x_pth2 have different numbers of time points, use the smaller number

nperm = 100;  % Adjust this to your needs

% Get data
x1 = x_ct;
x2 = x_pth;
stim1 = stim_ct;
stim2 = stim_pth;
n1 = size(x1, 1);
n2 = size(x2, 1);

% Calculate interaction information for each group
h0_ct = zeros(Nt_cluster, Nt_cluster);
h0_pth = zeros(Nt_cluster, Nt_cluster);

for t1 = 1:Nt_cluster
    for t2 = (t1+1):Nt_cluster
        JMI_tot1 = mi_gg([x1(:, t1) x1(:, t2)], stim1, false);
        II_tot1 = JMI_tot1 - mi_gg(x1(:, t1), stim1, false) - mi_gg(x1(:, t2), stim1, false);
        h0_ct(t1, t2) = II_tot1;

        JMI_tot2 = mi_gg([x2(:, t1) x2(:, t2)], stim2, false);
        II_tot2 = JMI_tot2 - mi_gg(x2(:, t1), stim2, false) - mi_gg(x2(:, t2), stim2, false);
        h0_pth(t1, t2) = II_tot2;
    end
end

% Calculate the difference between the groups
diff_h0 = h0_ct - h0_pth;

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

% Calculate the p-values for the difference
diff_p_values = sum(abs(diff_perm_t_values) >= abs(diff_h0), 3) / nperm;

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

% Convert Nt_cluster to time in milliseconds assuming your sampling rate
time_plot = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, Nt_cluster);

% Plot the difference in interaction information with significance clusters
figure;
imagesc(time_plot, time_plot, diff_h0);
hold on;
contour(time_plot, time_plot, zmap, 1, 'linecolor', 'k', 'LineWidth', 1);
colormap(brewermap([], '*RdBu'));
colorbar;
axis square;
lim = max(abs(min(min(diff_h0))), max(max(diff_h0)));
caxis([-lim, lim]);
xlabel('Time (ms)');
ylabel('Time (ms)');
title('Difference in Interaction Information between Groups with Cluster-based Correction');