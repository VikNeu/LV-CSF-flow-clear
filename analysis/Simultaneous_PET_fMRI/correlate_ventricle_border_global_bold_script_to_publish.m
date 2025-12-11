%% Add toolboxes and set directories/files
addpath(genpath('/data_august/pro_brain_clearance_scz/software/CoSMoMVPA-master'));
addpath(genpath('/data_august/pro_brain_clearance_scz/software/CanlabCore-master'));
addpath(genpath('/data_august/pro_brain_clearance_scz/software/canlab-MediationToolbox-19917e7'));
dataset = 'SCZ_TRIMAGE';
project_dir = fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/data',dataset);

file_path_subject_list = fullfile(project_dir,'progs/subjectlist_bids.txt'); 
file_path_results = fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/',dataset,'/ventricular_border_approach_24_07_25.xls');
fileID_good_subjects = fopen(file_path_subject_list, 'r');


% Check if the file was opened successfully
if fileID_good_subjects == -1
    error('Failed to open the file.')
end

% Create tables/arrays and a variable(num_iterations) to later save the
% Output into the correct row

num_iterations = 1;
Output_table = table ();
mediation_table = table();
plot_timecourse_table = table();
mediation_values = [];

%% Start looping over the sessions and subject

for session = 1:2
    frewind(fileID_good_subjects)
    while ~feof(fileID_good_subjects)
        % Read the next line from the file
        sub_id = fgetl(fileID_good_subjects);
        session_id = num2str(session);
        % Check if the line is empty (end of file)
        if isempty(sub_id)
            break; % Exit the loop when the end of the file is reached
        end
        
        %Get filepaths to the subjects fMRIs, masks and supplementary files
        TR = 2;
        fmri_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/filtered_func_data.nii.gz');
        fmri_aroma_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/ICA_AROMA/denoised_func_data_nonaggr.nii.gz');
        fmri_mean_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/mean_func.nii.gz');
        gm_mask_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/gm_mask_fast.nii.gz');
        csf_mask_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/func_csf_mask.nii'); % This is where your bottom slice CSF mask should be put
        ventricular_border_mask_filename = fullfile(project_dir,'masks', sub_id,strcat('ses-',session_id),'lateral_ventricles_mask_fmri_space_border.nii.gz'); % This is where your ventricular edge mask should be put
        fileID_fmri_file = fopen(fmri_filename, 'r');

        if fileID_fmri_file == -1
            disp('Subject is not completely preprocessed')
            continue
        end
        
        if isfile(ventricular_border_mask_filename)
            
            %Load fMRIs and masks
            fmri = cosmo_fmri_dataset(fmri_filename,'nifti_form','sform');
            fmri_aroma = cosmo_fmri_dataset(fmri_aroma_filename,'nifti_form', 'sform');
            fmri_mean = cosmo_fmri_dataset(fmri_mean_filename,'nifti_form','sform');
            gm_mask = cosmo_fmri_dataset(gm_mask_filename,'nifti_form','sform');
            csf_mask = cosmo_fmri_dataset(csf_mask_filename,'nifti_form','sform');
            ventricular_border_mask = cosmo_fmri_dataset(ventricular_border_mask_filename,'nifti_form','sform');
            gm_mask_idx = find(gm_mask.samples > 0.9);
            csf_mask_idx = find(csf_mask.samples > 0.9);
            ventricular_border_idx = find(ventricular_border_mask.samples > 0.9); 


            %% Timecourse extraction out of pipeline results - not used
            timeseries_table = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/timeseries.csv');
            %Read motion file
            data_table = readtable(timeseries_table);

            % Extract columns by column name directly out of pipeline file
%             gm_column = data_table.GM_AROMA; 
%             csf_column = data_table.CSF;


            %% Timecourse extraction by FMRI image - used
            gm_column = mean(fmri_aroma.samples(:,gm_mask_idx),2); 
            csf_column = mean(fmri.samples(:,csf_mask_idx),2);


            %% over 
            csf_column_z_transform = zscore(csf_column);
            
            %Get indices where there is a big CSF inflow (based on the
            %inflow effect of the bs_CSF mask timecourse
            inflow_indices = csf_column_z_transform > 1;
            
            %Get indices where ther is a CSF outflow (no CSF inflow)
            outflow_indices = csf_column_z_transform < 0;

            %Bandpass filter timeseries in the frequency spectrum of the
            %neuronal activity and take gradient of the gm_timeseries
            gm_column_detrend = detrend(gm_column);
            gm_column_band_timeseries = bandpass(gm_column_detrend, [0.01, 0.1], (1/TR)); 
            gm_derivative = gradient(gm_column_band_timeseries);
    
            csf_column_detrend = detrend(csf_column);
            csf_column_band_timeseries = bandpass(csf_column_detrend, [0.01, 0.1], (1/TR)); 
            
            % Create fMRI where the value of each voxels is the (Pearson) Correlation
            % value with the gm timecourse
            fmri_correlation = fmri_mean;
            for i = 1:size(fmri.samples,2)
                voxel = fmri.samples(:,i);
                correlation_value = corr(voxel,gm_column); 
                fmri_correlation.samples(1,i) = correlation_value;
            end
            
            
            %Preallocate subtraction_fmri
            subtraction_fmri = fmri_mean;
            
            %Subtract all images where there is a csf inflow from all
            %images where no CSF inflow happens
            subtraction_fmri.samples = mean(fmri.samples(outflow_indices,:)) - mean(fmri.samples(inflow_indices,:));
            
            %Take all voxels that don't correlate with GM (<0.1) and have a
            %negative value in the subtraction image for the further
            %analysis of the PFI timecourse
            subtraction_partial_volume_effect_idx = find(fmri_correlation.samples < 0.1 & subtraction_fmri.samples < 0);
            partial_volume_voxel_indices = intersect(ventricular_border_idx,subtraction_partial_volume_effect_idx);
            
            
            if ~isempty(partial_volume_voxel_indices)
                
                %Get mean timecourse of partial_volume_voxels bandpass and
                %take the gradient of the timecourse
                mean_partial_volume_voxel_tc = mean(fmri.samples(:,partial_volume_voxel_indices),2);
                partial_volume_voxel_tc_detrend = detrend(mean_partial_volume_voxel_tc);
                partial_volume_voxel_tc_band_timeseries = bandpass(partial_volume_voxel_tc_detrend, [0.01, 0.1], (1/TR));
                partial_volume_derivative = gradient(partial_volume_voxel_tc_band_timeseries);


                % Correlation of the GM, PFI and CSF timecourse with each
                % other: If GM or PFI are correlated with CSF we take the
                % gradient of the GM or PFI timecourse because the CSF
                % timecourse is inherently a signal for velocity (first
                % derivative)

                %GM and CSF timecourse are shifted by one TR due to
                %cross-correlation results

                [r_gm_csf,p] = corr(gm_derivative,csf_column_band_timeseries, type = 'Spearman');
                disp('Correlation gm and csf:')
                disp(r_gm_csf)
                disp(p)
        
                [r_ventricle_gm,p] = corr(partial_volume_voxel_tc_band_timeseries(1:end-1),gm_column_band_timeseries(2:end), type = 'Spearman');
                disp('Correlation gm and ventricular_border:')
                disp(r_ventricle_gm)
                disp(p)
        
                [r_ventricle_csf,p] = corr(partial_volume_derivative(1:end-2),csf_column(3:end), type = 'Spearman');
                disp('Correlation ventricular border and csf:')
                disp(r_ventricle_csf)
                disp(p)
                
           
                
                
                %% Compute cross correlation of PFI to CSF and GM timecourse
                [x_corr_ventricle_csf,lags_ventricle_csf] = xcorr(partial_volume_derivative,csf_column_band_timeseries,16,"normalized");
                xcorr_idx_ventricle_csf = find(lags_ventricle_csf == 0);
                [max_corr_value,max_corr_idx] = max(x_corr_ventricle_csf);
                partial_volume_csf_lag = lags_ventricle_csf(max_corr_idx);

                
    
                [x_corr_ventricle_gm,lags_ventricle_gm] = xcorr(partial_volume_voxel_tc_band_timeseries,gm_column_band_timeseries,16,"normalized");
                xcorr_idx_ventricle_gm = find(lags_ventricle_gm == 0);
                [min_corr_value,min_corr_idx] = min(x_corr_ventricle_gm);
                partial_volume_gm_lag = lags_ventricle_gm(min_corr_idx);
                
                [x_corr_gm_csf,lags_gm_csf] = xcorr(gm_derivative,csf_column_band_timeseries,16,"normalized");
                xcorr_idx_gm_csf = find(lags_gm_csf == 0);
                [min_corr_value,min_corr_idx] = min(x_corr_gm_csf);
                gm_csf_lag = lags_gm_csf(min_corr_idx);

%                 plot(mean_partial_volume_voxel_tc)
%                 plot(partial_volume_voxel_tc_detrend)
%                 plot(partial_volume_voxel_tc_band_timeseries)
%                 close all
                %% Compute standard deviation and normalized amplitude of the PFI and GM signal
                partial_volume_effect_amplitude_normalized = (max(mean_partial_volume_voxel_tc) - min(mean_partial_volume_voxel_tc))/mean(mean_partial_volume_voxel_tc);
                T_power = length(partial_volume_voxel_tc_detrend) * TR;
                partial_volume_effect_mean_power = (1/T_power) * sum(partial_volume_voxel_tc_detrend.^2)*TR;
                partial_volume_effect_confidence_interval_normalized = (prctile(partial_volume_voxel_tc_detrend,97.5) - prctile(partial_volume_voxel_tc_detrend,2.5))/mean(mean_partial_volume_voxel_tc);

                gm_amplitude_normalized = (max(gm_column) - min(gm_column))/mean(gm_column);

                %% Create the output table
                Output_table.subID(num_iterations) = string(strcat(sub_id,'-',session_id));
                Output_table.partial_volume_effect_amplitude_normalized(num_iterations) = partial_volume_effect_amplitude_normalized;
                Output_table.gm_amplitude_normalized(num_iterations) = gm_amplitude_normalized;
                Output_table.partial_volume_effect_confidence_interval_normalized(num_iterations) = partial_volume_effect_confidence_interval_normalized;
                Output_table.partial_volume_effect_mean_power(num_iterations) = partial_volume_effect_mean_power;
                Output_table.r_ventricle_gm(num_iterations) = r_ventricle_gm;
                Output_table.r_ventricle_csf(num_iterations) = r_ventricle_csf;
                Output_table.r_gm_csf(num_iterations) = r_gm_csf;
                Output_table.partial_volume_gm_lag(num_iterations) = partial_volume_gm_lag;
                Output_table.xcorr_ventricle_gm(num_iterations) = {x_corr_ventricle_gm(xcorr_idx_ventricle_gm-10:xcorr_idx_ventricle_gm+10)};
                Output_table.xcorr_ventricle_csf(num_iterations) = {x_corr_ventricle_csf(xcorr_idx_ventricle_csf-10:xcorr_idx_ventricle_csf+10)};
                Output_table.xcorr_gm_csf(num_iterations) = {x_corr_gm_csf(xcorr_idx_gm_csf-10:xcorr_idx_gm_csf+10)};

                %% Create a table with all the timecourse later used for the mediation of the timeseries
                mediation_table.subID(num_iterations) = string(strcat(sub_id,'-',session_id));
                mediation_table.gm_derivative_tc(num_iterations) = {gm_derivative};
                mediation_table.tfi_derivative_tc(num_iterations) = {partial_volume_derivative};
                mediation_table.csf_tc(num_iterations) = {csf_column_band_timeseries};
                
                plot_timecourse_table.subID(num_iterations) = string(strcat(sub_id,'-',session_id));
                plot_timecourse_table.gm_derivative_tc(num_iterations) = {gm_derivative};
                plot_timecourse_table.gm_tc(num_iterations) = {gm_column_band_timeseries};
                plot_timecourse_table.pfi_derivative_tc(num_iterations) = {partial_volume_derivative};
                % plot_timecourse_table.pfi_tc(num_iterations) = {partial_volume_voxel_tc_detrend}; % alternatively: partial_volume_voxel_tc_band_timeseries
                plot_timecourse_table.pfi_tc(num_iterations) = {mean_partial_volume_voxel_tc};
                plot_timecourse_table.csf_tc(num_iterations) = {csf_column_band_timeseries};
                num_iterations = num_iterations + 1;

            end
            
        end
    end
end
writetable(mediation_table,fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/',dataset,'/mediation_table_24_07_25.csv'))
writetable(plot_timecourse_table,fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/',dataset,'/plot_timecourse_table_24_07_25.csv'))
writetable(Output_table,file_path_results);
