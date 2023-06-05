% Initialization
Nt = size(x_pth, 2);
nperm = 100; % Adjust this to your needs
h0_pth = zeros(Nt, Nt);
hp_pth = zeros(Nt, Nt, nperm);

% Get data
x = x_pth;
stim = stim_pth;
n = size(x, 1);

for t1 = 1:Nt
    for t2 = (t1+1):Nt
        % Calculate total joint mutual information and interaction information
        JMI_tot = mi_gg([x(:, t1) x(:, t2)], stim, false);
        II_tot = JMI_tot - mi_gg(x(:, t1), stim, false) - mi_gg(x(:, t2), stim, false);

        % Store interaction information
        h0_pth(t1, t2) = II_tot;

        % Permutation test
        for i = 1:nperm
            stim_perm = stim(randperm(n));
            JMI_perm = mi_gg([x(:, t1) x(:, t2)], stim_perm, false);
            II_perm = JMI_perm - mi_gg(x(:, t1), stim_perm, false) - mi_gg(x(:, t2), stim_perm, false);

            % Store permutation results
            hp_pth(t1, t2, i) = II_perm;
        end
    end
end

% Store results
h0_values_pth = h0_pth;
hp_values_pth = hp_pth;

% Initialize parameters
preclust_pval = 0.05;  % Set this as needed
clust_pval = 0.05;  % Set this as needed

% Fill in the other half of the matrices
h0_pth = h0_pth + h0_pth';
for iperm = 1:nperm 
    hp_pth(:,:,iperm) = hp_pth(:,:,iperm) + hp_pth(:,:,iperm)';
end 

% Collect the largest suprathreshold clusters
clust_max_pth = zeros(nperm, 1);
for iperm = 1:nperm
    perms = true(1,nperm);
    perms(iperm) = 0;
    zvals = squeeze((hp_pth(:,:,iperm) - mean(hp_pth(:,:,perms),3)) ./ std(hp_pth(:,:,perms),[],3));
    zvals(abs(zvals) < norminv(1 - preclust_pval)) = 0; 
    
    clust_info = bwconncomp(zvals);
    clust_max_pth(iperm) = max([0 cellfun(@numel, clust_info.PixelIdxList)]);
end

% Identify significant clusters in real data
zmap_pth = squeeze((h0_pth - mean(hp_pth, 3)) ./ std(hp_pth, [], 3));
zmap_pth(abs(zmap_pth) < norminv(1 - preclust_pval)) = 0;

clust_info = bwconncomp(zmap_pth);
clust_size = cellfun(@numel, clust_info.PixelIdxList);
clust_th = prctile(clust_max_pth, 100 - clust_pval * 100);
clust_rem = find(clust_size < clust_th);
for i = 1:length(clust_rem)
    zmap_pth(clust_info.PixelIdxList{clust_rem(i)}) = 0;
end
zmap_pth(isnan(zmap_pth)) = 0;
zmap_pth = logical(zmap_pth);

% Convert Nt to time in milliseconds assuming your sampling rate
time_plot = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, Nt);

figure; 
imagesc(time_plot, time_plot, h0_pth);
hold on;

% Plot significant clusters
contour(time_plot, time_plot, zmap_pth, 1, 'linecolor', 'k', 'LineWidth', 1);

colormap(brewermap([], '*RdBu'));
colorbar;
axis square;

% Adjust color limits to center around zero
lim = max(abs(min(min(h0_pth))), max(max(h0_pth)));
caxis([-lim, lim]);

xlabel('Time (ms)');
ylabel('Time (ms)');
title('Interaction Information with Cluster-based Correction');

