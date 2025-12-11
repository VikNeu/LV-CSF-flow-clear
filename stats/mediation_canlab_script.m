%% Import toolboxes: TODO: EDIT THESE TWO LINES WITH YOUR SPECIFIC DIRECTORY
addpath(genpath('.../software/CanlabCore-master'));
addpath(genpath('.../software/canlab-MediationToolbox-19917e7'));

%% Set directories and filenames
% Path setup
thisFile = mfilename('fullpath');
thisFolder = fileparts(thisFile);
projectRoot = fileparts(thisFolder);

Simultaneous_PET_fMRI_dir = fullfile(projectRoot, 'data', 'Simultaneous_PET_FMRI');
flicker_dir = fullfile(projectRoot, 'data', 'Visual_stimulation');
sleep_dir = fullfile(projectRoot, 'data', 'Sleep');

Input_table_dopa = readtable(fullfile(DOPA_dir,'mediation_table_LONG.csv'));
Input_table_dopa = Input_table_dopa(2:22,1:end); % First participant has no PET therefore don't use it
Input_table_trimage = readtable(fullfile(TRIMAGE_dir,'mediation_table_TRIMAGE.csv'));
Input_table_trimage = Input_table_trimage(1:24,1:end);
Input_table_flicker = readtable(fullfile(flicker_dir,'mediation_table.csv'));
Input_table_sleep_sleep = readtable(fullfile(sleep_dir,'mediation_table_only_sleep.csv'));
Input_table_sleep_awake = readtable(fullfile(sleep_dir,'mediation_table_only_wake.csv'));


%% Extract timecourses for mediation
Input_table_rest = vertcat(Input_table_trimage,Input_table_dopa);
subjectlist = Input_table_rest.subID;
subjectlist_double = [];
subject_covar = subjectlist';
Input_table_rest = removevars(Input_table_rest, 'subID');
Input_matrix_rest = table2array(Input_table_rest);

gm_matrix_rest = Input_matrix_rest(:,1:(size(Input_matrix_rest,2)/3))';
pfi_matrix_rest = Input_matrix_rest(:,((size(Input_matrix_rest,2)/3)+1):((size(Input_matrix_rest,2)/3)*2))';
csf_matrix_rest = Input_matrix_rest(:,((size(Input_matrix_rest,2)/3)*2+1):((size(Input_matrix_rest,2)/3)*3))';


[path_rest,stats_rest] = mediation(gm_matrix_rest(2:end-1,:),csf_matrix_rest(3:end,:),pfi_matrix_rest(1:end-2,:));

%% Same for visual stimulation dataset
Input_table_flicker = removevars(Input_table_flicker, 'subID');
Input_matrix_flicker = table2array(Input_table_flicker);
gm_matrix_flicker = Input_matrix_flicker(:,1:(size(Input_matrix_flicker,2)/3))';
pfi_matrix_flicker = Input_matrix_flicker(:,((size(Input_matrix_flicker,2)/3)+1):((size(Input_matrix_flicker,2)/3)*2))';
csf_matrix_flicker = Input_matrix_flicker(:,((size(Input_matrix_flicker,2)/3)*2+1):((size(Input_matrix_flicker,2)/3)*3))';

[path_flicker,stats_flicker] = mediation(gm_matrix_flicker(2:end-1,:),csf_matrix_flicker(3:end,:),pfi_matrix_flicker(1:end-2,:));

%% Same for the sleep dataset (only the timepoints where the subjects where actually sleeping)
Input_table_sleep_sleep = removevars(Input_table_sleep_sleep, 'subID');
Input_matrix_sleep_sleep = table2array(Input_table_sleep_sleep);
gm_matrix_sleep_sleep = Input_matrix_sleep_sleep(:,1:(size(Input_matrix_sleep_sleep,2)/3))';
pfi_matrix_sleep_sleep = Input_matrix_sleep_sleep(:,((size(Input_matrix_sleep_sleep,2)/3)+1):((size(Input_matrix_sleep_sleep,2)/3)*2))';
csf_matrix_sleep_sleep = Input_matrix_sleep_sleep(:,((size(Input_matrix_sleep_sleep,2)/3)*2+1):((size(Input_matrix_sleep_sleep,2)/3)*3))';

[path_sleep_sleep,stats_sleep_sleep] = mediation(gm_matrix_sleep_sleep(2:end,:),csf_matrix_sleep_sleep(2:end,:),pfi_matrix_sleep_sleep(1:end-1,:));

%% Same for the sleep dataset (only the timepoints where the subjects where actually awake)
Input_table_sleep_awake = removevars(Input_table_sleep_awake, 'subID');
Input_matrix_sleep_awake = table2array(Input_table_sleep_awake);
gm_matrix_sleep_awake = Input_matrix_sleep_awake(:,1:(size(Input_matrix_sleep_awake,2)/3))';
pfi_matrix_sleep_awake = Input_matrix_sleep_awake(:,((size(Input_matrix_sleep_awake,2)/3)+1):((size(Input_matrix_sleep_awake,2)/3)*2))';
csf_matrix_sleep_awake = Input_matrix_sleep_awake(:,((size(Input_matrix_sleep_awake,2)/3)*2+1):((size(Input_matrix_sleep_awake,2)/3)*3))';

[path_sleep_awake,stats_sleep_awake] = mediation(gm_matrix_sleep_awake(2:end,:),csf_matrix_sleep_awake(2:end,:),pfi_matrix_sleep_awake(1:end-1,:));