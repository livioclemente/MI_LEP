% Load data
addpath('/Users/livioclemente/Documents/MATLAB/eeglab2023.0/');
addpath('/Users/livioclemente/Documents/MATLAB/EKG_EEG-master/');
addpath('/Users/livioclemente/Documents/MATLAB//bluewhitered/');

% Calculate baseline indices
Fs = 256;
baseline_duration = 100;
baseline_samples = round(baseline_duration / 1000 * Fs);
t_start = -0.1; % start time in seconds
t_end = 1; % end time in seconds
time_vector = t_start:1/Fs:t_end-1/Fs;

% Path tracciati
data_path1 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/pth/ginocchio/';
data_path2 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/pth/piede/';
data_path3 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/pth/mano/';
data_path4 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/ct/ginocchio/';
data_path5 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/ct/piede/';
data_path6 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/ct/mano/';
data_path7 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/dl_pth/ginocchio/';
data_path8 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/dl_pth/piede/';
data_path9 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/dl_pth/mano/';
data_path10 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/normali/ginocchio/';
data_path11 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/normali/piede/';
data_path12 = '/Users/livioclemente/Documents/prova_connettivita/64_ch_3_site/MI/conditions/normali/mano/';

% List of subject files
subject_files1 = dir([data_path1 '*.set']);
subject_files2 = dir([data_path2 '*.set']);
subject_files3 = dir([data_path3 '*.set']);
subject_files4 = dir([data_path4 '*.set']);
subject_files5 = dir([data_path5 '*.set']);
subject_files6 = dir([data_path6 '*.set']);
subject_files7 = dir([data_path7 '*.set']);
subject_files8 = dir([data_path8 '*.set']);
subject_files9 = dir([data_path9 '*.set']);
subject_files_10 = dir([data_path10 '*.set']);
subject_files11 = dir([data_path11 '*.set']);
subject_files12 = dir([data_path12 '*.set']);

% numero soggetti
num_subjects_pth = length(subject_files1);
num_subjects_ct = length(subject_files4);
num_subjects_dlpth = length(subject_files7);
num_subjects_norm = length(subject_files_10);

% Inizializzazione variabili
% pth
subj_data_pth = cell(num_subjects_pth, 1);
subj_data_baseline_corrected_pth = cell(num_subjects_pth, 1);
erp_mean_pth = cell(num_subjects_pth, 1);
erp_mean_baseline_corrected_pth = cell(num_subjects_pth, 1);
I_values_pth = cell(num_subjects_pth, 1);
Iint_values_pth = cell(num_subjects_pth, 1);
% ct
subj_data_ct = cell(num_subjects_ct, 1);
subj_data_baseline_corrected_ct = cell(num_subjects_ct, 1);
erp_mean_ct = cell(num_subjects_ct, 1);
erp_mean_baseline_corrected_ct = cell(num_subjects_ct, 1);
I_values_ct = cell(num_subjects_ct, 1);
Iint_values_ct = cell(num_subjects_ct, 1);
% dlpth
subj_data_dlpth = cell(num_subjects_dlpth, 1);
subj_data_baseline_corrected_dlpth = cell(num_subjects_dlpth, 1);
erp_mean_dlpth = cell(num_subjects_dlpth, 1);
erp_mean_baseline_corrected_dlpth = cell(num_subjects_dlpth, 1);
I_values_dlpth = cell(num_subjects_dlpth, 1);
Iint_values_dlpth = cell(num_subjects_dlpth, 1);
% norm
subj_data_norm = cell(num_subjects_norm, 1);
subj_data_baseline_corrected_norm = cell(num_subjects_norm, 1);
erp_mean_norm = cell(num_subjects_norm, 1);
erp_mean_baseline_corrected_norm = cell(num_subjects_norm, 1);
I_values_norm = cell(num_subjects_norm, 1);
Iint_values_norm = cell(num_subjects_norm, 1);

%%

% Pazienti PTH

for s = 1:num_subjects_pth

    EEG1 = pop_loadset(subject_files1(s).name, data_path1);
    EEG2 = pop_loadset(subject_files2(s).name, data_path2);
    EEG3 = pop_loadset(subject_files3(s).name, data_path3);
    subj_data_pth{s} = cat(3, EEG1.data, EEG2.data, EEG3.data);

    
    % Subtract baseline mean from each trial for each condition
    subj_data_baseline_corrected_pth{s} = subj_data_pth{s} - mean(subj_data_pth{s}(:, 1:baseline_samples, :), 2);

    [electrodes, Nt, trials1] = size(subj_data_baseline_corrected_pth{s});
    x_pth = squeeze(subj_data_baseline_corrected_pth{s}(11, :, :)); % Select one electrode

    % stim vector for all subjects
    trials1 = size(EEG1.data, 3);
    trials2 = size(EEG2.data, 3);
    trials3 = size(EEG3.data, 3);

    stim_pth = [ones(1, trials1), 2 * ones(1, trials2), 3 * ones(1, trials3)];

    x_pth = x_pth';
    stim_pth = stim_pth';
    x_pth = copnorm(x_pth);
    stim_pth = copnorm(stim_pth);

    for ti = 1:Nt
        I_pth(ti) = mi_gg(x_pth(:, ti), stim_pth(:, 1), true, true);
    end

    Nintt_pth = Nt;
    Iint_pth = zeros(Nintt_pth, Nintt_pth);
    noise=.00000005*randn(size(x_pth,1),1);
    for t1_pth = 1:Nintt_pth
        for t2_pth = (t1_pth+1):Nintt_pth
            Ijoint_pth = mi_gg([x_pth(:, t1_pth) x_pth(:, t2_pth)+noise], stim_pth(:, 1), true, true);
            Iint_pth(t1_pth, t2_pth) = Ijoint_pth - I_pth(t1_pth) - I_pth(t2_pth);
        end
    end
    Iint_pth = Iint_pth + Iint_pth';

    % Store the I and Iint values for the current subject
    I_values_pth{s} = I_pth;
    Iint_values_pth{s} = Iint_pth;

    % average ERP for each condition
    erp_mean_pth{s} = mean(subj_data_baseline_corrected_pth{s}, 3);
    n_samples_pth = size(erp_mean_pth{s}, 2);
    time_axis_corrected_pth = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, n_samples_pth);

    % baseline correction
    baseline_range_pth = find(time_vector >= -100 & time_vector <= 0);
    erp_mean_baseline_pth{s} = mean(erp_mean_pth{s}(:, baseline_range_pth), 2);

    erp_mean_baseline_corrected_pth{s} = erp_mean_pth{s} - repmat(erp_mean_baseline_pth{s}, 1, size(erp_mean_pth{s}, 2));

end

% group average for Mutual Information, Iint, and ERP
group_I_values_pth = mean(cat(3, I_values_pth{:}), 3);
group_Iint_values_pth = mean(cat(3, Iint_values_pth{:}), 3);
group_erp_mean_baseline_corrected_pth = mean(cat(3, erp_mean_baseline_corrected_pth{:}), 3);

%%

%Controlli CT

for s = 1:num_subjects_ct

    EEG4 = pop_loadset(subject_files4(s).name, data_path4);
    EEG5 = pop_loadset(subject_files5(s).name, data_path5);
    EEG6 = pop_loadset(subject_files6(s).name, data_path6);
    subj_data_ct{s} = cat(3, EEG4.data, EEG5.data, EEG6.data);

    % Subtract baseline mean from each trial for each condition
    subj_data_baseline_corrected_ct{s} = subj_data_ct{s} - mean(subj_data_ct{s}(:, 1:baseline_samples, :), 2);

    [electrodes, Nt_ct, trials4] = size(subj_data_baseline_corrected_ct{s});
    x_ct = squeeze(subj_data_baseline_corrected_ct{s}(11, :, :)); % Select one electrode

    % stim vector for all subjects
    trials4 = size(EEG4.data, 3);
    trials5 = size(EEG5.data, 3);
    trials6 = size(EEG6.data, 3);

    stim_ct = [ones(1, trials4), 2 * ones(1, trials5), 3 * ones(1, trials6)];

    x_ct = x_ct';
    stim_ct = stim_ct';
    x_ct = copnorm(x_ct);
    stim_ct = copnorm(stim_ct);

    for ti_ct = 1:Nt_ct
        I_ct(ti_ct) = mi_gg(x_ct(:, ti_ct), stim_ct(:, 1), true, true);
    end

    Nintt_ct = Nt_ct;
    Iint_ct = zeros(Nintt_ct, Nintt_ct);
    noise=.00000005*randn(size(x_ct,1),1);
    for t1_ct = 1:Nintt_ct
        for t2_ct = (t1_ct+1):Nintt_ct
            Ijoint_ct = mi_gg([x_ct(:, t1_ct) x_ct(:, t2_ct)+noise], stim_ct(:, 1), true, true);
            Iint_ct(t1_ct, t2_ct) = Ijoint_ct - I_ct(t1_ct) - I_ct(t2_ct);
        end
    end
    Iint_ct = Iint_ct + Iint_ct';

    % Store the I and Iint values for the current subject
    I_values_ct{s} = I_ct;
    Iint_values_ct{s} = Iint_ct;

    % average ERP for each condition
    erp_mean_ct{s} = mean(subj_data_baseline_corrected_ct{s}, 3);
    n_samples_ct = size(erp_mean_ct{s}, 2);
    time_axis_corrected_ct = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, n_samples_ct);

    % baseline correction
    baseline_range_ct = find(time_vector >= -100 & time_vector <= 0);
    erp_mean_baseline_ct{s} = mean(erp_mean_ct{s}(:, baseline_range_ct), 2);

    erp_mean_baseline_corrected_ct{s} = erp_mean_ct{s} - repmat(erp_mean_baseline_ct{s}, 1, size(erp_mean_ct{s}, 2));

end

% group average for Mutual Information, Iint, and ERP
group_I_values_ct = mean(cat(3, I_values_ct{:}), 3);
group_Iint_values_ct = mean(cat(3, Iint_values_ct{:}), 3);
group_erp_mean_baseline_corrected_ct = mean(cat(3, erp_mean_baseline_corrected_ct{:}), 3);

%%

%Pazienti DL_PTH

for s = 1:num_subjects_dlpth
    % Load data for the current subject
    EEG7 = pop_loadset(subject_files7(s).name, data_path7);
    EEG8 = pop_loadset(subject_files8(s).name, data_path8);
    EEG9 = pop_loadset(subject_files9(s).name, data_path9);
    subj_data_dlpth{s} = cat(3, EEG7.data, EEG8.data, EEG9.data);

    % Subtract baseline mean from each trial for each condition
    subj_data_baseline_corrected_dlpth{s} = subj_data_dlpth{s} - mean(subj_data_dlpth{s}(:, 1:baseline_samples, :), 2);

    % Proceed with the rest of the script using the baseline-corrected data
    [electrodes, Nt_dlpth, trials7] = size(subj_data_baseline_corrected_dlpth{s});
    x_dlpth = squeeze(subj_data_baseline_corrected_dlpth{s}(11, :, :)); % Select one electrode

    % Create the stim vector with condition labels for all subjects
    trials7 = size(EEG7.data, 3);
    trials8 = size(EEG8.data, 3);
    trials9 = size(EEG9.data, 3);

    stim_dlpth = [ones(1, trials7), 2 * ones(1, trials8), 3 * ones(1, trials9)];

    x_dlpth = x_dlpth';
    stim_dlpth = stim_dlpth';
    x_dlpth = copnorm(x_dlpth);
    stim_dlpth = copnorm(stim_dlpth);

    for ti_dlpth = 1:Nt_dlpth
        I_dlpth(ti_dlpth) = mi_gg(x_dlpth(:, ti_dlpth), stim_dlpth(:, 1), true, true);
    end

    Nintt_dlpth = Nt_dlpth;
    Iint_dlpth = zeros(Nintt_dlpth, Nintt_dlpth);
    noise=.00000005*randn(size(x_dlpth,1),1);
    for t1_dlpth = 1:Nintt_dlpth
        for t2_dlpth = (t1_dlpth+1):Nintt_dlpth
            Ijoint_dlpth = mi_gg([x_dlpth(:, t1_dlpth) x_dlpth(:, t2_dlpth)+noise], stim_dlpth(:, 1), true, true);
            Iint_dlpth(t1_dlpth, t2_dlpth) = Ijoint_dlpth - I_dlpth(t1_dlpth) - I_dlpth(t2_dlpth);
        end
    end
    Iint_dlpth = Iint_dlpth + Iint_dlpth';

    % Store the I and Iint values for the current subject
    I_values_dlpth{s} = I_dlpth;
    Iint_values_dlpth{s} = Iint_dlpth;

    % Calculate the average ERP for each condition
    erp_mean_dlpth{s} = mean(subj_data_baseline_corrected_dlpth{s}, 3);
    n_samples_dlpth = size(erp_mean_dlpth{s}, 2);
    time_axis_corrected_dlpth = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, n_samples_dlpth);

    % Apply baseline correction
    baseline_range_dlpth = find(time_vector >= -100 & time_vector <= 0);
    erp_mean_baseline_dlpth{s} = mean(erp_mean_dlpth{s}(:, baseline_range_dlpth), 2);

    erp_mean_baseline_corrected_dlpth{s} = erp_mean_dlpth{s} - repmat(erp_mean_baseline_dlpth{s}, 1, size(erp_mean_dlpth{s}, 2));

end

% Calculate the group average for Mutual Information, Iint, and ERP
group_I_values_dlpth = mean(cat(3, I_values_dlpth{:}), 3);
group_Iint_values_dlpth = mean(cat(3, Iint_values_dlpth{:}), 3);
group_erp_mean_baseline_corrected_dlpth = mean(cat(3, erp_mean_baseline_corrected_dlpth{:}), 3);

%%

%Pazienti normali NORM

for s = 1:num_subjects_norm
    % Load data for the current subject
    EEG10 = pop_loadset(subject_files_10(s).name, data_path10);
    EEG11 = pop_loadset(subject_files11(s).name, data_path11);
    EEG12 = pop_loadset(subject_files12(s).name, data_path12);
    subj_data_norm{s} = cat(3, EEG10.data, EEG11.data, EEG12.data);

    % Subtract baseline mean from each trial for each condition
    subj_data_baseline_corrected_norm{s} = subj_data_norm{s} - mean(subj_data_norm{s}(:, 1:baseline_samples, :), 2);

    % Proceed with the rest of the script using the baseline-corrected data
    [electrodes, Nt_norm, trials10] = size(subj_data_baseline_corrected_norm{s});
    x_norm = squeeze(subj_data_baseline_corrected_norm{s}(11, :, :)); % Select one electrode

    % Create the stim vector with condition labels for all subjects
    trials10 = size(EEG10.data, 3);
    trials11 = size(EEG11.data, 3);
    trials12 = size(EEG12.data, 3);

    stim_norm = [ones(1, trials10), 2 * ones(1, trials11), 3 * ones(1, trials12)];

    x_norm = x_norm';
    stim_norm = stim_norm';
    x_norm = copnorm(x_norm);
    stim_norm = copnorm(stim_norm);

    for ti_norm = 1:Nt_norm
        I_norm(ti_norm) = mi_gg(x_norm(:, ti_norm), stim_norm(:, 1), true, true);
    end

    Nintt_norm = Nt_norm;
    Iint_norm = zeros(Nintt_norm, Nintt_norm);
    noise=.00000005*randn(size(x_norm,1),1);
    for t1_norm = 1:Nintt_norm
        for t2_norm = (t1_norm+1):Nintt_norm
            Ijoint_norm = mi_gg([x_norm(:, t1_norm) x_norm(:, t2_norm)+noise], stim_norm(:, 1), true, true);
            Iint_norm(t1_norm, t2_norm) = Ijoint_norm - I_norm(t1_norm) - I_norm(t2_norm);
        end
    end
    Iint_norm = Iint_norm + Iint_norm';

    % Store the I and Iint values for the current subject
    I_values_norm{s} = I_norm;
    Iint_values_norm{s} = Iint_norm;

    % Calculate the average ERP for each condition
    erp_mean_norm{s} = mean(subj_data_baseline_corrected_norm{s}, 3);
    n_samples_norm = size(erp_mean_norm{s}, 2);
    time_axis_corrected_norm = linspace(t_start * 1000, t_end * 1000 - 1 / Fs * 1000, n_samples_norm);

    % Apply baseline correction
    baseline_range_norm = find(time_vector >= -100 & time_vector <= 0);
    erp_mean_baseline_norm{s} = mean(erp_mean_norm{s}(:, baseline_range_norm), 2);

    erp_mean_baseline_corrected_norm{s} = erp_mean_norm{s} - repmat(erp_mean_baseline_norm{s}, 1, size(erp_mean_norm{s}, 2));

end

% Calculate the group average for Mutual Information, Iint, and ERP
group_I_values_norm = mean(cat(3, I_values_norm{:}), 3);
group_Iint_values_norm = mean(cat(3, Iint_values_norm{:}), 3);
group_erp_mean_baseline_corrected_norm = mean(cat(3, erp_mean_baseline_corrected_norm{:}), 3);

%%

%Salvataggio dei dati

% pth
all_subj_data_pth = cell(num_subjects_pth, 1);
for s = 1:num_subjects_pth
    % ...

    % Store all subject data for each condition
    all_subj_data_pth{s} = subj_data_baseline_corrected_pth{s};

    % ...
end
save('data_pth.mat', 'group_I_values_pth', 'group_Iint_values_pth', 'group_erp_mean_baseline_corrected_pth', 'I_values_pth', 'Iint_values_pth', 'erp_mean_baseline_corrected_pth', 'all_subj_data_pth', 'hdr');

% ct
all_subj_data_ct = cell(num_subjects_ct, 1);
for s = 1:num_subjects_ct
    % ...

    % Store all subject data for each condition
    all_subj_data_ct{s} = subj_data_baseline_corrected_ct{s};

    % ...
end
save('data_ct.mat', 'group_I_values_ct', 'group_Iint_values_ct', 'group_erp_mean_baseline_corrected_ct', 'I_values_ct', 'Iint_values_ct', 'erp_mean_baseline_corrected_ct', 'all_subj_data_ct', 'hdr');

% dlpth
all_subj_data_dlpth = cell(num_subjects_dlpth, 1);
for s = 1:num_subjects_dlpth
    % ...

    % Store all subject data for each condition
    all_subj_data_dlpth{s} = subj_data_baseline_corrected_dlpth{s};

    % ...
end
save('data_dlpth.mat', 'group_I_values_dlpth', 'group_Iint_values_dlpth', 'group_erp_mean_baseline_corrected_dlpth', 'I_values_dlpth', 'Iint_values_dlpth', 'erp_mean_baseline_corrected_dlpth', 'all_subj_data_dlpth', 'hdr');

% norm
all_subj_data_norm = cell(num_subjects_norm, 1);
for s = 1:num_subjects_norm
    % ...

    % Store all subject data for each condition
    all_subj_data_norm{s} = subj_data_baseline_corrected_norm{s};

    % ...
end
save('data_norm.mat', 'group_I_values_norm', 'group_Iint_values_norm', 'group_erp_mean_baseline_corrected_norm', 'I_values_norm', 'Iint_values_norm', 'erp_mean_baseline_corrected_norm', 'all_subj_data_norm', 'hdr');
