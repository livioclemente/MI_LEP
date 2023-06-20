load("group_ct.mat");
load("group_dlpth.mat");
load("group_norm.mat");
load("group_pth.mat");

% General info
t_start = -0.1;
t_end = 1;
Fs = 256;
time_plot = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, Nt);
lw = 1.5;
xl = [min(time_plot), max(time_plot)];

%% figure ct

% Create a figure with three subplots
figure('Position', [0, 0, 1200, 800])

% Main subplot: Mutual Information
axm = subplot(5, 5, [2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
imagesc(time_axis_corrected, time_axis_corrected, group_Iint_values_ct);
% Imposta i limiti della colorbar in modo simmetrico
clim_max = max(abs(group_Iint_values_ct(:)));
clim([-clim_max, clim_max]);
colormap(bluewhitered);
xlabel('Time (ms)');
ylabel('Time (ms)');
colorbar;
title('Interaction Information');

% Second subplot: Average ERP of all conditions
subplot(5, 5, [22 23 24.5 24.8]);
plot(time_axis_corrected, group_erp_mean_ct(11, :), 'k','LineWidth', lw);
axis tight;
xlim([-100, 1000]);
xlabel('Time (ms)');
ylabel('ERP');
box off;

% Third subplot: Mutual Information as a function of time
subplot(5, 5, [1 6 11 16]);
plot(time_axis_corrected, group_I_values_ct, 'k', 'LineWidth', lw);
axis tight;
box off;
xlabel('Time (ms)');
ylabel('II');
set(gca, 'CameraUpVector', [-1 0 0]);

%% figure norm

% Create a figure with three subplots
figure('Position', [0, 0, 1200, 800])

% Main subplot: Mutual Information
axm = subplot(5, 5, [2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
imagesc(time_axis_corrected, time_axis_corrected, group_Iint_values_norm);
% Imposta i limiti della colorbar in modo simmetrico
clim_max = max(abs(group_Iint_values_norm(:)));
clim([-clim_max, clim_max]);
colormap(bluewhitered);
xlabel('Time (ms)');
ylabel('Time (ms)');
colorbar;
title('Interaction Information');

% Second subplot: Average ERP of all conditions
subplot(5, 5, [22 23 24.5 24.8]);
plot(time_axis_corrected, group_erp_mean_norm(11, :), 'k','LineWidth', lw);
axis tight;
xlim([-100, 1000]);
xlabel('Time (ms)');
ylabel('ERP');
box off;

% Third subplot: Mutual Information as a function of time
subplot(5, 5, [1 6 11 16]);
plot(time_axis_corrected, group_I_values_norm, 'k', 'LineWidth', lw);
axis tight;
box off;
xlabel('Time (ms)');
ylabel('II');
set(gca, 'CameraUpVector', [-1 0 0]);

%% figure dlpth

% Create a figure with three subplots
figure('Position', [0, 0, 1200, 800])

% Main subplot: Mutual Information
axm = subplot(5, 5, [2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
imagesc(time_axis_corrected, time_axis_corrected, group_Iint_values_dlpth);
% Imposta i limiti della colorbar in modo simmetrico
clim_max = max(abs(group_Iint_values_dlpth(:)));
clim([-clim_max, clim_max]);
colormap(bluewhitered);
xlabel('Time (ms)');
ylabel('Time (ms)');
colorbar;
title('Interaction Information');

% Second subplot: Average ERP of all conditions
subplot(5, 5, [22 23 24.5 24.8]);
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

%% figure pth

% Create a figure with three subplots
figure('Position', [0, 0, 1200, 800])

% Main subplot: Mutual Information
axm = subplot(5, 5, [2 3 4 5 7 8 9 10 12 13 14 15 17 18 19 20]);
imagesc(time_axis_corrected, time_axis_corrected, group_Iint_values_pth);
% Imposta i limiti della colorbar in modo simmetrico
clim_max = max(abs(group_Iint_values_pth(:)));
clim([-clim_max, clim_max]);
colormap(bluewhitered);
xlabel('Time (ms)');
ylabel('Time (ms)');
colorbar;
title('Interaction Information');

% Second subplot: Average ERP of all conditions
subplot(5, 5, [22 23 24.5 24.8]);
plot(time_axis_corrected, group_erp_mean_pth(11, :), 'k','LineWidth', lw);
axis tight;
xlim([-100, 1000]);
xlabel('Time (ms)');
ylabel('ERP');
box off;

% Third subplot: Mutual Information as a function of time
subplot(5, 5, [1 6 11 16]);
plot(time_axis_corrected, group_I_values_pth, 'k', 'LineWidth', lw);
axis tight;
box off;
xlabel('Time (ms)');
ylabel('II');
set(gca, 'CameraUpVector', [-1 0 0]);
