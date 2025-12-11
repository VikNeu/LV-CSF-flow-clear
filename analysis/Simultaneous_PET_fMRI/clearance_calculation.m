function  [ventricular_slope,brain_slope,ventricle_volume, mean_striatum_value_suv_normalized, brain_volume_integer, etiv_integer, mean_seed_suv_normalized_cut] = clearance_calculation(options, first_frame_for_derivative_calculation,last_frame)

%% Script to calculate ventricular clearance from masks in PET images %%

% Add Cosmo-software
addpath(genpath(['/data_august/pro_brain_clearance_scz/software/CoSMoMVPA-master']));
    % info: all of the COSMO commands run with "fmri...", they also work for PET as we use NIFTI anyways!

%% SETTINGS %%
    % load data and get indices above set threshold

% Get the masks we are using: seed, brain, wm, gm 
seed_msk = cosmo_fmri_dataset(options.seed_msk_dir);
seed_idx = find(seed_msk.samples > options.seed_msk_thrsh);
brain_msk = cosmo_fmri_dataset(options.brain_msk_dir);
brain_idx = find(brain_msk.samples > options.brain_msk_thrsh);
striatum_msk = cosmo_fmri_dataset(options.striatum_mask_dir);
striatum_idx = find(striatum_msk.samples > 0);


% Get the pet data we are using
pet_cut = cosmo_fmri_dataset(options.pet_cut_dir);


%% DO THE SUVR NORMALISATION NEEDED FOR CALCULATION %%
    % Divide PET image through whole brain SUV of the first 5 minutes after
    % the arterial and ventricular values get constant
mean_tc_brain_cut = mean(pet_cut.samples(:,brain_idx'),2);
brain_suv_5_min = mean_tc_brain_cut(1,1); % Mean brain value at 20-25min
pet_cut_suv_normalized.samples = pet_cut.samples/brain_suv_5_min;




% time = 1:(timepoints-frames_to_delete);
% full_time = 1:timepoints;

%% CALCULATE THE CSF CLEARANCE OUT OF ELVM (ERODED LATERAL VENTRICLE MASK) %%
    % Transform the timepoints into actual time in seconds using the time
    % scale of the frames
    % Overlay a linear model fit onto the data and use the slope as a proxy
    % for clearance (Interpret in percent)
time_seconds = [0;30;45;60;75;90;105;120;135;150;165;180;200;220;240;300;360;480;600;900;1200;1500;1800;2100;2400;2700;3000;3300;3600;3900];
timepoints_cut = time_seconds(first_frame_for_derivative_calculation:last_frame);

% Do the linear model
mean_seed_suv_normalized_cut = mean(pet_cut_suv_normalized.samples(:,seed_idx'),2);
matrix_seed_suv_normalized_cut = [timepoints_cut, mean_seed_suv_normalized_cut];
mdl_seed_suv_normalized_cut = fitlm(matrix_seed_suv_normalized_cut(:, 1), matrix_seed_suv_normalized_cut(:, 2));
ventricular_slope = mdl_seed_suv_normalized_cut.Coefficients.Estimate(2); %get slope


%% CALCULATE THE PET CLEARANCE OUT OF WHOLE BRAIN
mean_brain_suv_normalized_cut = mean(pet_cut_suv_normalized.samples(:,brain_idx'),2);
matrix_brain_suv_normalized_cut = [timepoints_cut, mean_brain_suv_normalized_cut];
mdl_brain_suv_normalized_cut = fitlm(matrix_brain_suv_normalized_cut(:, 1), matrix_brain_suv_normalized_cut(:, 2));
brain_slope = mdl_brain_suv_normalized_cut.Coefficients.Estimate(2);

%% Calculate ventricular volume
fileID_freesurfer = fopen(options.freesurfer_stats);
freesurfer_stats_data = textscan(fileID_freesurfer, '%s', 'Delimiter', '\n');
freesurfer_stats = freesurfer_stats_data{1,1};
brain_volume = freesurfer_stats{15,1};
parts_brain_volume = split(brain_volume);
brain_volume_string = parts_brain_volume{8,1};
brain_volume_integer = str2double(brain_volume_string(1:end-1));
etiv = freesurfer_stats{34,1};
parts_etiv = split(etiv);
etiv_string = parts_etiv{9,1};
etiv_integer = str2double(etiv_string);
left_ventricle_volume = freesurfer_stats{80,1};
right_ventricle_volume = freesurfer_stats{98,1};
parts_left_ventricle_volume = split(left_ventricle_volume);
parts_right_ventricle_volume = split(right_ventricle_volume);
left_ventricle_volume_string = parts_left_ventricle_volume{4,1};
left_ventricle_volume_integer = str2double(left_ventricle_volume_string);
right_ventricle_volume_string = parts_right_ventricle_volume{4,1};
right_ventricle_volume_integer = str2double(right_ventricle_volume_string);
ventricle_volume = left_ventricle_volume_integer + right_ventricle_volume_integer;

%% Calculate Mean striatum value (Striatum SUVr calculation)
mean_striatum_value_suv_normalized = mean(pet_cut_suv_normalized.samples(:,striatum_idx),"all");