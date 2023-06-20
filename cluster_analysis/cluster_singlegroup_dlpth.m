clear
clc

load("group_dlpth.mat");

% Initialization
Nt = size(x_dlpth, 2);
nperm = 1000; % Adjust this to your needs
h0_dlpth = zeros(Nt, Nt);
hp_dlpth = zeros(Nt, Nt, nperm);

% Get data
x = x_dlpth;
stim = stim_dlpth;
n = size(x, 1);

for t1 = 1:Nt
    for t2 = (t1+1):Nt
        % Calculate total joint mutual information and interaction information
        JMI_tot = mi_gg([x(:, t1) x(:, t2)], stim, false);
        II_tot = JMI_tot - mi_gg(x(:, t1), stim, false) - mi_gg(x(:, t2), stim, false);

        % Store interaction information
        h0_dlpth(t1, t2) = II_tot;

        % Permutation test
        for i = 1:nperm
            stim_perm = stim(randperm(n));
            JMI_perm = mi_gg([x(:, t1) x(:, t2)], stim_perm, false);
            II_perm = JMI_perm - mi_gg(x(:, t1), stim_perm, false) - mi_gg(x(:, t2), stim_perm, false);

            % Store permutation results
            hp_dlpth(t1, t2, i) = II_perm;
        end
    end
end

% Store results
h0_values_dlpth = h0_dlpth;
hp_values_dlpth = hp_dlpth;

% Initialize parameters
preclust_pval = 0.05;  % Set this as needed
clust_pval = 0.05;  % Set this as needed

% Fill in the other half of the matrices
h0_dlpth = h0_dlpth + h0_dlpth';
for iperm = 1:nperm 
    hp_dlpth(:,:,iperm) = hp_dlpth(:,:,iperm) + hp_dlpth(:,:,iperm)';
end 

% Collect the largest suprathreshold clusters
clust_max_dlpth = zeros(nperm, 1);
for iperm = 1:nperm
    perms = true(1,nperm);
    perms(iperm) = 0;
    zvals = squeeze((hp_dlpth(:,:,iperm) - mean(hp_dlpth(:,:,perms),3)) ./ std(hp_dlpth(:,:,perms),[],3));
    zvals(abs(zvals) < norminv(1 - preclust_pval)) = 0; 
    
    clust_info = bwconncomp(zvals);
    clust_max_dlpth(iperm) = max([0 cellfun(@numel, clust_info.PixelIdxList)]);
end

% Identify significant clusters in real data
zmap_dlpth = squeeze((h0_dlpth - mean(hp_dlpth, 3)) ./ std(hp_dlpth, [], 3));
zmap_dlpth(abs(zmap_dlpth) < norminv(1 - preclust_pval)) = 0;

clust_info = bwconncomp(zmap_dlpth);
clust_size = cellfun(@numel, clust_info.PixelIdxList);
clust_th = prctile(clust_max_dlpth, 100 - clust_pval * 100);
clust_rem = find(clust_size < clust_th);
for i = 1:length(clust_rem)
    zmap_dlpth(clust_info.PixelIdxList{clust_rem(i)}) = 0;
end
zmap_dlpth(isnan(zmap_dlpth)) = 0;
zmap_dlpth = logical(zmap_dlpth);

% General info
t_start = -0.1;
t_end = 1;
Fs = 256;
time_plot = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, Nt);
lw = 1.5;
xl = [min(time_plot), max(time_plot)];

% Create a figure with three subplots
figure('Position', [0, 0, 1200, 800])


% Main subplot: Mutual Information
axm = subplot(5, 5, [2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
imagesc(time_plot, time_plot, h0_dlpth);
hold on;

% Plot significant clusters
contour(time_plot, time_plot, zmap_dlpth, 1, 'linecolor', 'k', 'LineWidth', 1);

colormap(brewermap([], '*RdBu'));
colorbar;
axis square;

% Adjust color limits to center around zero
lim = max(abs(min(min(h0_dlpth))), max(max(h0_dlpth)));
clim([-lim, lim]);

xlabel('Time (ms)');
ylabel('Time (ms)');
title('Interaction Information with Cluster-based Correction');

% Second subplot: Average ERP of all conditions
subplot(5, 5, [22.45 23 24.2 24.34]);
plot(time_axis_corrected, group_erp_mean_dlpth(11, :), 'k','LineWidth', lw);
axis tight;
xlim([-100, 1000]);
xlabel('Time (ms)');
ylabel('ERP');
box off;

% Third subplot: Mutual Information as a function of time
subplot(5, 5, [1 6 11 16]);
plot(time_axis_corrected, group_I_values_dlpth, 'k', 'LineWidth', lw);
axis tight;
box off;
xlabel('Time (ms)');
ylabel('II');
set(gca, 'CameraUpVector', [-1 0 0]);