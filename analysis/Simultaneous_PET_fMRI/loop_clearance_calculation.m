
%% Loop to calculate ventricular tracer clearance in pet %%

% Add Cosmo-software
addpath(genpath(['/data_august/pro_brain_clearance_scz/software/CoSMoMVPA-master']));

clear;clc;

%% SETTINGS %%
% call set pathnahmes to get necessary directories
dataset = 'SCZ_DOPA_LONG';
projectDir = fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/data',dataset);
seed = 'lateral_ventricle';
brain = 'brain';
number_of_sessions = 1;
tracer = '18FDOPA';
first_frame_for_derivative_calculation = 21;
first_frame_for_realignment = 5;
last_frame = 30;

Cfg = set_pathnames_calculation(projectDir);
% define the thresholds
options.seed_msk_thrsh = 0.9; 
options.brain_msk_thrsh = 0.9;
Output_table = table;
num_iterations = 1;
%% SUBJECT AND SESSION LOOP FOR CLEARANCE CALCULATION %%
for session = 1:number_of_sessions
    session = num2str(session);
    for sub = 1:length(Cfg.subNames)
        subject = Cfg.subNames{sub};
        fprintf('Processing %s\n', subject);
        
        
        options.seed_msk_dir = fullfile(Cfg.masks,subject,strcat('ses-',session),'lateral_ventricles_mask_reduced.nii.gz');
        options.brain_msk_dir = fullfile(Cfg.masks,subject,strcat('ses-',session),'brain_pet_space.nii.gz');
        options.full_ventricle_msk_dir = fullfile (Cfg.masks,subject,strcat('ses-',session),'lateral_ventricles_mask_pet_space.nii.gz');
        options.striatum_mask_dir = fullfile(Cfg.masks,subject,strcat('ses-',session),'striatum_mask.nii');


        % Define the pet images
        options.pet_cut_dir = fullfile(Cfg.preprocDir,subject,['ses-',session],'/pet/',[subject,'_ses-',session,'_trc-',tracer,'_rec-acdyn_pet_derivative_calculation_SUVcorr.nii']);
        % Define Freesurfer Output for potential covariates
        options.freesurfer_stats = fullfile(Cfg.preprocDir,subject,['ses-',session],'/anat/freesurfer',subject,'stats/aseg.stats');
        % Define the output for results table
        options.Output_table_filepath = fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results',dataset,"pet_results_27_05_25_with_timecourse.xlsx");
        options.Output_table_filepath_mat = fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results',dataset,"pet_results_27_05_25_with_timecourse.mat");
        
        % Only start the calculation script when pet and lateral ventricle mask exist
        if isfile(options.seed_msk_dir) && isfile(options.pet_cut_dir)
            % Call function
            [ventricular_slope,brain_slope, ventricle_volume, mean_striatum_value_suv_normalized, brain_volume, etiv, mean_seed_suv_normalized_cut] = clearance_calculation(options, first_frame_for_derivative_calculation,last_frame); 


            Output_table.subID(num_iterations) = string(strcat(subject, '-',num2str(session)));
            Output_table.ventricular_slope(num_iterations) = ventricular_slope;
            Output_table.brain_slope(num_iterations) = brain_slope;
            Output_table.ventricle_volume(num_iterations) = ventricle_volume;
            Output_table.brain_volume(num_iterations) = brain_volume;
            Output_table.etiv(num_iterations) = etiv;
            Output_table.mean_striatum_value_suv_normalized(num_iterations) = mean_striatum_value_suv_normalized;
            Output_table.ventricular_tc(num_iterations) = {mean_seed_suv_normalized_cut};
            num_iterations = num_iterations+1;
        else
            disp(strcat(subject,'-',session,':Is not fully preprocessed yet'))
        end
    end
end
writetable(Output_table,options.Output_table_filepath);
save(options.Output_table_filepath_mat, "Output_table");
