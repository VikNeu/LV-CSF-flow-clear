%% Add toolboxes and set directories
addpath(genpath('/data_august/pro_brain_clearance_scz/software/CoSMoMVPA-master'))
file_path_subjectlist = fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/data/flicker/progs/subjectlist_bids.txt'); 
fileID_subjectlist = fopen(file_path_subjectlist, 'r'); 
projectDir = '/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/data/flicker';

Output_table = table();
mediation_table = table();
plot_timecourse_table = table();
plot_timecourse_table_rest = table();
num_iterations = 1;

%% Loop through sessions and subjects
for session = 1:12
    session_id = num2str(session);
    frewind(fileID_subjectlist)
    while ~feof(fileID_subjectlist)
         sub_id = fgetl(fileID_subjectlist); 

         task_fmri_filename = fullfile(projectDir,'preprocData/',sub_id,['ses-',session_id],'func/frs001_adv_search.feat/filtered_func_data.nii.gz');
         task_aroma_fmri_filename = fullfile(projectDir,'preprocData/',sub_id,['ses-',session_id],'func/frs001_adv_search.feat/ICA_AROMA/denoised_func_data_nonaggr.nii.gz');
         json_filename = fullfile(projectDir,'rawData/',sub_id,['ses-',session_id],['func/',sub_id,'_ses-',session_id,'_task-rest_bold.json']); 
         mean_task_fmri_filename = fullfile(projectDir,'preprocData/',sub_id,['ses-',session_id],'func/frs001_adv_search.feat/mean_func.nii.gz');
         ventricular_border_mask_filename = fullfile(projectDir,'masks/',sub_id,['ses-',session_id],'lateral_ventricles_mask_fmri_space_border.nii.gz');
         eroded_ventricle_mask_filename = fullfile(projectDir,'masks/',sub_id,['ses-',session_id],'lateral_ventricles_mask_fmri_space_eroded.nii.gz');
         occipital_cortex_mask_filename = fullfile(projectDir,'masks/',sub_id,['ses-',session_id],'occipital_cortex_mask_fmri_space.nii.gz');
         csf_mask_filename = fullfile(projectDir,'masks/',sub_id,'func_csf_mask_universal.nii');
         frontal_cortex_mask_filename = fullfile(projectDir,'masks/',sub_id,['ses-',session_id],'frontal_cortex_mask_fmri_space.nii.gz');
         anterior_border_mask_filename = fullfile(projectDir,'masks/',sub_id,['ses-',session_id],'anterior_ventricular_border_fmri_space.nii.gz');
         posterior_border_mask_filename = fullfile(projectDir,'masks/',sub_id,['ses-',session_id],'posterior_ventricular_border_fmri_space.nii.gz');
         extraction_timeframe_filename = fullfile(projectDir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/timeframe_borders.txt');

         if exist(anterior_border_mask_filename) && exist(ventricular_border_mask_filename) && exist(task_fmri_filename) && exist(task_aroma_fmri_filename)
             
             %% Load fMRIs and masks
             task_fmri = cosmo_fmri_dataset(task_fmri_filename);
             task_aroma_fmri = cosmo_fmri_dataset(task_aroma_fmri_filename);
             mean_task_fmri = cosmo_fmri_dataset(mean_task_fmri_filename);
             ventricular_border_mask = cosmo_fmri_dataset(ventricular_border_mask_filename);
             ventricular_border_mask_idx = find(ventricular_border_mask.samples > 0.9);
             eroded_ventricle_mask = cosmo_fmri_dataset(eroded_ventricle_mask_filename);
             eroded_ventricle_mask_idx = find(eroded_ventricle_mask.samples > 0.9);
             occipital_cortex_mask = cosmo_fmri_dataset(occipital_cortex_mask_filename);
             occipital_cortex_mask_idx = find(occipital_cortex_mask.samples > 0.9);
             frontal_cortex_mask = cosmo_fmri_dataset(frontal_cortex_mask_filename);
             frontal_cortex_mask_idx = find(frontal_cortex_mask.samples > 0.9);
             csf_mask = cosmo_fmri_dataset(csf_mask_filename);
             csf_mask_idx = find(csf_mask.samples > 0.9);
             anterior_border_mask = cosmo_fmri_dataset(anterior_border_mask_filename);
             posterior_border_mask = cosmo_fmri_dataset(posterior_border_mask_filename);
             anterior_border_mask_idx = find(anterior_border_mask.samples > 0.1);
             posterior_border_mask_idx = find(posterior_border_mask.samples > 0.1);
             
             % Anterior/Posterior ventricular border mask are the voxels
             % that are in the ventricular border mask and in the anterior
             % border mask from the MNI space
             anterior_ventricular_border_mask_idx = intersect(ventricular_border_mask_idx,anterior_border_mask_idx);
             posterior_ventricular_border_mask_idx = intersect(ventricular_border_mask_idx,posterior_border_mask_idx);
             
             anterior_border_mask_fmri = mean_task_fmri;
             anterior_border_mask_fmri.samples(:,:) = 0;
             anterior_border_mask_fmri.samples(:,anterior_ventricular_border_mask_idx) = 1;

             posterior_border_mask_fmri = mean_task_fmri;
             posterior_border_mask_fmri.samples(:,:) = 0;
             posterior_border_mask_fmri.samples(:,posterior_ventricular_border_mask_idx) = 1;


             % Import TR from json file
             json_file = jsondecode(fileread(json_filename));
             TR = json_file.RepetitionTime;
             

            


            %Create timecourses for ventricular border (different regions) and gray matter (different regions) and bandpass them
            anterior_ventricular_border_tc = mean(task_fmri.samples(:,anterior_ventricular_border_mask_idx),2);
            anterior_ventricular_border_tc_detrend = detrend(anterior_ventricular_border_tc);
            anterior_ventricular_border_tc_bandpass = bandpass(anterior_ventricular_border_tc_detrend,[0.01,0.1],1/TR);
            posterior_ventricular_border_tc = mean(task_fmri.samples(:,posterior_ventricular_border_mask_idx),2);
            posterior_ventricular_border_tc_detrend = detrend(posterior_ventricular_border_tc);
            posterior_ventricular_border_tc_bandpass = bandpass(posterior_ventricular_border_tc_detrend,[0.01,0.1],1/TR);
            posterior_ventricular_border_derivative = gradient(posterior_ventricular_border_tc_bandpass);
            eroded_ventricle_tc = mean(task_fmri.samples(:,eroded_ventricle_mask_idx),2);
            eroded_ventricle_tc_detrend = detrend(eroded_ventricle_tc);
            eroded_ventricle_tc_bandpass = bandpass(eroded_ventricle_tc_detrend,[0.01,0.1],1/TR);
            csf_tc = mean(task_fmri.samples(:,csf_mask_idx),2);
            csf_tc_detrend = detrend(csf_tc);
            csf_tc_bandpass = bandpass(csf_tc_detrend,[0.01,0.1],1/TR);

            occipital_gm_tc = mean(task_aroma_fmri.samples(:,occipital_cortex_mask_idx),2);
            occipital_gm_tc_detrend = detrend(occipital_gm_tc);
            occipital_gm_tc_bandpass = bandpass(occipital_gm_tc_detrend,[0.01,0.1],1/TR);
            occipital_gm_derivative = gradient(occipital_gm_tc_bandpass);
            frontal_gm_tc = mean(task_aroma_fmri.samples(:,frontal_cortex_mask_idx),2);
            frontal_gm_tc_detrend = detrend(frontal_gm_tc);
            frontal_gm_tc_bandpass = bandpass(frontal_gm_tc_detrend,[0.01,0.1],1/TR);


            [r_anterior_occipital_gm,p] = corr(anterior_ventricular_border_tc_bandpass(1:end-1),occipital_gm_tc_bandpass(2:end),type='Spearman')
           
            [r_posterior_occipital_gm,p] = corr(posterior_ventricular_border_tc_bandpass(1:end-1),occipital_gm_tc_bandpass(2:end),type='Spearman')
          
            [r_eroded_ventricle_occipital_gm,p] = corr(eroded_ventricle_tc_bandpass(1:end-1),occipital_gm_tc_bandpass(2:end),type='Spearman')
           
          
            [r_posterior_frontal_gm,p] = corr(posterior_ventricular_border_tc_bandpass(1:end-1),frontal_gm_tc_bandpass(2:end),type='Spearman')

            [r_anterior_frontal_gm,p] = corr(anterior_ventricular_border_tc_bandpass(1:end-1),frontal_gm_tc_bandpass(2:end),type='Spearman')
        
            

            %Cross correlations of only posterior ventricle

            [x_corr_posterior_ventricle_csf,lags_posterior_ventricle_csf] = xcorr(posterior_ventricular_border_derivative,csf_tc_bandpass,15,'normalized');
            xcorr_idx_posterior_ventricle_csf = find(lags_posterior_ventricle_csf == 0);
            [max_corr_value,max_corr_idx] = max(x_corr_posterior_ventricle_csf);
            posterior_partial_volume_csf_lag = lags_posterior_ventricle_csf(max_corr_idx);

            [x_corr_posterior_ventricle_occipital_gm,lags_posterior_ventricle_occipital_gm] = xcorr(posterior_ventricular_border_tc_bandpass,occipital_gm_tc_bandpass,15,'normalized');
            xcorr_idx_posterior_ventricle_occipital_gm = find(lags_posterior_ventricle_occipital_gm == 0);
            [max_corr_value,max_corr_idx] = max(x_corr_posterior_ventricle_occipital_gm);
            posterior_partial_volume_occipital_gm_lag = lags_posterior_ventricle_occipital_gm(max_corr_idx);

            [x_corr_csf_occipital_gm,lags_csf_occipital_gm] = xcorr(occipital_gm_derivative,csf_tc_bandpass,15,'normalized');
            xcorr_idx_csf_occipital_gm = find(lags_csf_occipital_gm == 0);
            [max_corr_value,max_corr_idx] = max(x_corr_csf_occipital_gm);
            csf_occipital_gm_lag = lags_csf_occipital_gm(max_corr_idx);
           
            % Calculate minmax amplitude (normalized by mean timecourse)
            % for Occipital gm and posterior PFI
            occipital_gm_amplitude = max(occipital_gm_tc) - min (occipital_gm_tc);
            occipital_gm_amplitude_normalized = occipital_gm_amplitude/mean(occipital_gm_tc);
            partial_volume_effect_amplitude_normalized = (max(posterior_ventricular_border_tc) - min(posterior_ventricular_border_tc))/mean(posterior_ventricular_border_tc);

            %Correlations are all adjusted to the crosscorrelation results
            [r_posterior_ventricle_csf,p] = corr(posterior_ventricular_border_derivative(1:end-2),csf_tc_bandpass(3:end), type = 'Spearman');

            Output_table.r_anterior_occipital_gm(num_iterations) = r_anterior_occipital_gm;
            Output_table.r_anterior_frontal_gm(num_iterations) = r_anterior_frontal_gm;
            Output_table.r_posterior_occipital_gm(num_iterations) = r_posterior_occipital_gm;
            Output_table.r_eroded_ventricle_occipital_gm(num_iterations) = r_eroded_ventricle_occipital_gm;
            Output_table.r_posterior_frontal_gm(num_iterations) = r_posterior_frontal_gm;
            Output_table.occipital_gm_amplitude_normalized(num_iterations) = occipital_gm_amplitude_normalized;
            Output_table.partial_volume_effect_amplitude_normalized(num_iterations) = partial_volume_effect_amplitude_normalized;
            Output_table.r_posterior_ventricle_csf(num_iterations) = r_posterior_ventricle_csf;
            Output_table.xcorr_posterior_ventricle_occipital_gm(num_iterations) = {x_corr_posterior_ventricle_occipital_gm(xcorr_idx_posterior_ventricle_occipital_gm-10:xcorr_idx_posterior_ventricle_occipital_gm+10)};
            Output_table.xcorr_posterior_ventricle_csf(num_iterations) = {x_corr_posterior_ventricle_csf(xcorr_idx_posterior_ventricle_csf-10:xcorr_idx_posterior_ventricle_csf+10)};
            Output_table.sub_id(num_iterations) = string(sub_id);
            Output_table.session(num_iterations) = session;

            %Save the timecourses in a table for mediation analysis later
            mediation_table.subID(num_iterations) = string(strcat(sub_id,'-',session_id));
            mediation_table.occipital_gm_derivative_tc(num_iterations) = {occipital_gm_derivative};
            mediation_table.posterior_tfi_derivative_tc(num_iterations) = {posterior_ventricular_border_derivative};
            mediation_table.csf_tc(num_iterations) = {csf_tc_bandpass};
            
            %Save the timecourses for single subject plotting
            plot_timecourse_table.subID(num_iterations) = string(strcat(sub_id,'-',session_id));
            plot_timecourse_table.occipital_gm_derivative_tc(num_iterations) = {occipital_gm_derivative};
            plot_timecourse_table.occipital_gm_tc(num_iterations) = {occipital_gm_tc_bandpass};
            plot_timecourse_table.posterior_pfi_derivative_tc(num_iterations) = {posterior_ventricular_border_derivative};
            plot_timecourse_table.posterior_pfi_tc(num_iterations) = {posterior_ventricular_border_tc_bandpass};
            %Add the remaining data, as the csv gets too big otherwise
            plot_timecourse_table_rest.subID(num_iterations) = string(strcat(sub_id,'-',session_id));
            plot_timecourse_table_rest.anterior_pfi_tc(num_iterations) = {anterior_ventricular_border_tc_bandpass};
            plot_timecourse_table_rest.csf_tc(num_iterations) = {csf_tc_bandpass};
            
            num_iterations = num_iterations + 1;


         end
    end
end

writetable(Output_table,'/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker/flicker_subtraction_approach_25_05_25.xlsx')
writetable(mediation_table,'/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker/mediation_table_25_05_25.csv')
writetable(plot_timecourse_table,'/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker/plot_timecourse_table_25_05_25.csv')
writetable(plot_timecourse_table_rest,'/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker/plot_timecourse_table_rest_25_05_25.csv')