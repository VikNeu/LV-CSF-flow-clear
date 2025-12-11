%% Set toolboxes and directories
addpath(genpath('/data_august/pro_brain_clearance_scz/software/CoSMoMVPA-master'));
addpath(genpath('/data_august/pro_brain_clearance_scz/software/CanlabCore-master'));
addpath(genpath('/data_august/pro_brain_clearance_scz/software/canlab-MediationToolbox-19917e7'));
project_dir = '/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/data/sleep';
file_path_subjectlist = fullfile(project_dir,'progs/subjectlist_bids.txt'); 

% Open the file for reading
fileID_good_subjects = fopen(file_path_subjectlist, 'r');
Output_table_file_path = fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/sleep/sleep_amplitudes_90_sec_chunks_06_06_2025.xls');
% Check if the file was opened successfully
if fileID_good_subjects == -1
    error('Failed to open the file.')
end
num_iterations = 1;
num_iterations_mediation = 1;
num_iterations_mediation_only_sleep = 1;
num_iterations_mediation_only_wake = 1;
mediation_table = table();
mediation_table_only_sleep = table();
mediation_table_only_wake = table();
plot_timecourse_table_wake = table();
plot_timecourse_table_sleep = table();
plot_timecourse_table_wake_unfil = table();
plot_timecourse_table_sleep_unfil = table();
Output_table = table ();
for session = 1:12   
    frewind(fileID_good_subjects)
    while ~feof(fileID_good_subjects)
        % Read the next line from the file
        sub_id = fgetl(fileID_good_subjects);
        session_id = num2str(session);
        % Check if the line is empty (end of file)
        if isempty(sub_id)
            break; % Exit the loop when the end of the file is reached
        end
        
        %% Get filenames of fMRIs, masks and supplementary files for every subject
        fmri_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/filtered_func_data.nii.gz');
        fmri_aroma_filename = fullfile(project_dir,'preprocData/',sub_id,['ses-',session_id],'func/frs001_adv_search.feat/ICA_AROMA/denoised_func_data_nonaggr.nii.gz');
        fmri_json_filename = fullfile(project_dir,'/rawData',sub_id,strcat('ses-',session_id),'/func',[sub_id,'_ses-',session_id,'_task-rest_bold.json']);
        fmri_mean_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/mean_func.nii.gz');
        gm_mask_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/gm_mask_fast.nii.gz');
        csf_mask_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/func_csf_mask.nii');
        ventricular_border_mask_filename = fullfile(project_dir,'masks', sub_id,strcat('ses-',session_id),'lateral_ventricles_mask_fmri_space_border.nii.gz');
        extraction_timeframe_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/timeframe_borders.txt');
        fileID_fmri_file = fopen(fmri_filename, 'r');

        sleep_stadium_table_filename = fullfile(project_dir,'sourcedata',[sub_id,'-sleep-stage.txt']);
        %% Import sleep data 
        opts = detectImportOptions(sleep_stadium_table_filename);
        if size(opts.VariableTypes,2) == 4  & (string(sub_id) ~= 'sub-27') %Formatting reasons
            opts.VariableTypes{1,4} = 'char';
        else 
            opts.VariableTypes{1,3} = 'char';
        end
        sleep_stadium_table = readtable(sleep_stadium_table_filename,opts);
        sleep_stadium_table.x30_sec_epoch_sleep_stage(strcmp(sleep_stadium_table.x30_sec_epoch_sleep_stage,'W')) = {'0'};
        sleep_stadium_table.x30_sec_epoch_sleep_stage(strcmp(sleep_stadium_table.x30_sec_epoch_sleep_stage,'Unscorable')) = {'99'};
        sleep_stadium_table.x30_sec_epoch_sleep_stage(strcmp(sleep_stadium_table.x30_sec_epoch_sleep_stage,'Uncertain')) = {'99'};
        sleep_stadium_table.session_number = NaN(height(sleep_stadium_table),1);
        for i = 1:height(sleep_stadium_table)
            lastChar = sleep_stadium_table.session{i}(end);
            lastInt = str2num(lastChar);
            if sleep_stadium_table.session{i}(1:10) == 'task-sleep'
                lastInt = lastInt + 2;
            end
            sleep_stadium_table.session_number(i) = lastInt;
        end
        
        %% Skip if subject is not preprocessed
        if fileID_fmri_file == -1
            disp('Subject is not completely preprocessed')
            continue
        end
        
        

        if isfile(ventricular_border_mask_filename) && isfile(gm_mask_filename) && isfile(csf_mask_filename) && isfile(fmri_aroma_filename)
    
            
            %% Import data of fMRIs and masks
            fmri = cosmo_fmri_dataset(fmri_filename,'nifti_form','sform');
            fmri_aroma = cosmo_fmri_dataset(fmri_aroma_filename,'nifti_form', 'sform');
            fmri_mean = cosmo_fmri_dataset(fmri_mean_filename,'nifti_form','sform');
            fmri_json = jsondecode(fileread(fmri_json_filename));
            TR = fmri_json.RepetitionTime;
            gm_mask = cosmo_fmri_dataset(gm_mask_filename,'nifti_form','sform');
            csf_mask = cosmo_fmri_dataset(csf_mask_filename,'nifti_form','sform');
            ventricle_mask = cosmo_fmri_dataset(ventricular_border_mask_filename,'nifti_form','sform');
            gm_mask_idx = find(gm_mask.samples > 0.9);
            csf_mask_idx = find(csf_mask.samples > 0.9);
            ventricle_idx = find(ventricle_mask.samples > 0.9);
            timeframe_borders = importdata(extraction_timeframe_filename);

            %% Extract timeseries from the timeseries table 
            gm_column = mean(fmri_aroma.samples(:,gm_mask_idx),2); 
            gm_derivative = gradient(gm_column);
            csf_column = mean(fmri.samples(:,csf_mask_idx),2);
            csf_column_z_transform = zscore(csf_column);
            
            %% Get timepoints with big inflow effects or outflow (no inflow) from the bsCSF timecourse
            inflow_indices = csf_column_z_transform > 1;
            outflow_indices = csf_column_z_transform < 0; %prctile(csf_column_z_transform,15);

            %% Detrend and bandpassfilter different timecourses
            gm_column_detrend = detrend(gm_column);
            gm_column_band_timeseries = bandpass(gm_column_detrend, [0.01, 0.1], (1/TR));
    
            csf_column_detrend = detrend(csf_column);
            csf_column_band_timeseries = bandpass(csf_column_detrend, [0.01, 0.1], (1/TR));
    
            
            
            %% Filter supplementary sleep stage file for the session
            sleep_stadium_table_session_wise = sleep_stadium_table(sleep_stadium_table.session_number == session,:);
            
            %% Create an fMRI where every voxel has the value of its correlation to the gray matter
            fmri_correlation = fmri_mean;  
            for i = 1:size(fmri.samples,2)
                voxel = fmri.samples(:,i);
                correlation_value = corr(voxel,gm_column); %include time_without_movement again
                fmri_correlation.samples(1,i) = correlation_value;
            end


            %% Create an fMRI where each voxel has the value of the mean of every outflow timepoint subtracted by the mean of every inflow timepoint
            subtraction_fmri = fmri_mean;
            subtraction_fmri.samples = mean(fmri.samples(outflow_indices,:)) - mean(fmri.samples(inflow_indices,:));
            subtraction_partial_volume_effect_idx = find(fmri_correlation.samples < 0.1 & subtraction_fmri.samples < 0);
            
            
            %% The PFI mask is defined by voxels that are part of the ventricular border mask do not correlate with the GM timcourse more than r = 0.1 and have a negative value in the subtraction fMRI
            partial_volume_voxel_indices = intersect(ventricle_idx,subtraction_partial_volume_effect_idx);


            if ~isempty(partial_volume_voxel_indices) 
                timeslot = 2;
                while timeslot < (height(sleep_stadium_table_session_wise)-3)
                    if strcmp(sleep_stadium_table_session_wise.x30_sec_epoch_sleep_stage{timeslot-1},sleep_stadium_table_session_wise.x30_sec_epoch_sleep_stage{timeslot}) && strcmp(sleep_stadium_table_session_wise.x30_sec_epoch_sleep_stage{timeslot},sleep_stadium_table_session_wise.x30_sec_epoch_sleep_stage{timeslot+1}) && strcmp(sleep_stadium_table_session_wise.x30_sec_epoch_sleep_stage{timeslot},sleep_stadium_table_session_wise.x30_sec_epoch_sleep_stage{timeslot+2}) && strcmp(sleep_stadium_table_session_wise.x30_sec_epoch_sleep_stage{timeslot},sleep_stadium_table_session_wise.x30_sec_epoch_sleep_stage{timeslot+3})
                        %Adjust fmri TRs and Sleep scoring to the same
                        %timescale
                        time_displacement = ((timeframe_borders(1,1)-1)*TR); %If the fmri got cropped because of movement
                        Start_TR = ceil((sleep_stadium_table_session_wise.epoch_start_time_sec(timeslot)-(time_displacement))/TR); %because we cut first 5 TR in preprocessing
                        End_TR = Start_TR + floor(90/TR);
                        
                        %% If the extracted Start_TR and End_TR are not part of the fMRI, because they were cropped due to movement then continue
                        if End_TR > (timeframe_borders(1,2) - timeframe_borders(1,1))
                            timeslot = timeslot + 1;
                            continue
                            
                        end
                        if Start_TR < 1
                            timeslot = timeslot + 1;
                            continue
                        end
                        
                        %Calculate GM amplitude and nonfiltered data
                        gm_column_sessionwise = gm_column(Start_TR:End_TR);
                        csf_column_sessionwise = csf_column(Start_TR:End_TR);
                        gm_amplitude = max(gm_column_sessionwise) - min(gm_column_sessionwise);
                        gm_amplitude_normalized = gm_amplitude/mean(gm_column_sessionwise);
      
                        %% Get the timecourses of the specific sleep session
                        gm_column_band_timeseries_sessionwise = gm_column_band_timeseries(Start_TR:End_TR);
                        gm_derivative_sessionwise = gradient(gm_column_band_timeseries_sessionwise);
                        
                        csf_column_band_timeseries_sessionwise = csf_column_band_timeseries(Start_TR:End_TR);
                        
                        mean_partial_volume_voxel_tc = mean(fmri.samples(:,partial_volume_voxel_indices),2);
                        partial_volume_voxel_tc_detrend = detrend(mean_partial_volume_voxel_tc);
                        partial_volume_voxel_tc_band_timeseries = bandpass(partial_volume_voxel_tc_detrend, [0.01, 0.1], (1/TR));
                        partial_volume_derivative = gradient(partial_volume_voxel_tc_band_timeseries);
                        partial_volume_voxel_tc_band_timeseries_sessionwise = partial_volume_voxel_tc_band_timeseries(Start_TR:End_TR);
                        partial_volume_derivative_sessionwise = gradient(partial_volume_voxel_tc_band_timeseries_sessionwise);
            
                        % Add in for nonfiltered data
                        pfi_column_sessionwise = mean_partial_volume_voxel_tc(Start_TR:End_TR);
                        %% Get correlations of the timecourses

                        [r_ventricle_gm,p] = corr(partial_volume_voxel_tc_band_timeseries_sessionwise(1:end-1),gm_column_band_timeseries_sessionwise(2:end), type='Spearman');
  
                
                        [r_ventricle_csf,p] = corr(partial_volume_derivative_sessionwise(1:end-1),csf_column_band_timeseries_sessionwise(2:end), type='Spearman');

                        %% Get the cross-correlations of the timecourses
                        [x_corr_ventricle_csf,lags_ventricle_csf] = xcorr(partial_volume_derivative_sessionwise,csf_column_band_timeseries_sessionwise, 'normalized');
                        xcorr_idx_ventricle_csf = find(lags_ventricle_csf == 0);
            
                        [x_corr_ventricle_gm,lags_ventricle_gm] = xcorr(partial_volume_derivative_sessionwise,gm_derivative_sessionwise, 'normalized');
                        xcorr_idx_ventricle_gm = find(lags_ventricle_gm == 0);
                        
                     
                        %% Get the PFI Amplitude of the specific sleep session
                        partial_volume_effect_amplitude = max(mean_partial_volume_voxel_tc(Start_TR:End_TR)) - min(mean_partial_volume_voxel_tc(Start_TR:End_TR));
                        partial_volume_effect_amplitude_normalized = partial_volume_effect_amplitude/mean(mean_partial_volume_voxel_tc);
                        
    
                        %% Add data to Output table
                        Output_table.session(num_iterations) = sleep_stadium_table_session_wise.session(timeslot);
                        Output_table.start_time(num_iterations) = sleep_stadium_table_session_wise.epoch_start_time_sec(timeslot);
                        Output_table.sleep_stage(num_iterations) = sleep_stadium_table_session_wise.x30_sec_epoch_sleep_stage(timeslot);
                        Output_table.session_number(num_iterations) = sleep_stadium_table_session_wise.session_number(timeslot);
                        Output_table.partial_volume_amplitude(num_iterations) = partial_volume_effect_amplitude;
                        Output_table.partial_volume_amplitude_normalized(num_iterations) = partial_volume_effect_amplitude_normalized;
                        Output_table.gm_amplitude_normalized(num_iterations) = gm_amplitude_normalized;
                        Output_table.Start_TR(num_iterations) = Start_TR;
                        Output_table.End_TR(num_iterations) = End_TR;
                        Output_table.r_ventricle_gm(num_iterations) = r_ventricle_gm;
                        Output_table.r_ventricle_csf(num_iterations) = r_ventricle_csf;
                        Output_table.xcorr_ventricle_gm(num_iterations) = {x_corr_ventricle_gm(xcorr_idx_ventricle_gm-10:xcorr_idx_ventricle_gm+10)};
                        Output_table.xcorr_ventricle_csf(num_iterations) = {x_corr_ventricle_csf(xcorr_idx_ventricle_csf-10:xcorr_idx_ventricle_csf+10)};
                        Output_table.sub_id(num_iterations) = string(sub_id);
                        
                        %% Add timecourses to a table which is later used for mediation (seperated for awake and sleeping timecourses) and plotting single subject timecourses
                        if Output_table.sleep_stage{num_iterations} == '1' | Output_table.sleep_stage{num_iterations} == '2' 
                            mediation_table_only_sleep.subID(num_iterations_mediation_only_sleep) = string(sub_id);
                            mediation_table_only_sleep.gm_derivative_tc(num_iterations_mediation_only_sleep) = {gm_derivative_sessionwise};
                            mediation_table_only_sleep.pfi_derivative_tc(num_iterations_mediation_only_sleep) = {partial_volume_derivative_sessionwise};
                            mediation_table_only_sleep.csf_tc(num_iterations_mediation_only_sleep) = {csf_column_band_timeseries_sessionwise};

                            plot_timecourse_table_sleep.subID(num_iterations_mediation_only_sleep) = string(sub_id);
                            plot_timecourse_table_sleep.gm_derivative_tc(num_iterations_mediation_only_sleep) = {gm_derivative_sessionwise/mean(gm_derivative_sessionwise)};
                            plot_timecourse_table_sleep.gm_tc(num_iterations_mediation_only_sleep) = {gm_column_band_timeseries_sessionwise/mean(gm_column_band_timeseries_sessionwise)};
                            plot_timecourse_table_sleep.pfi_derivative_tc(num_iterations_mediation_only_sleep) = {partial_volume_derivative_sessionwise/mean(partial_volume_derivative_sessionwise)};
                            plot_timecourse_table_sleep.pfi_tc(num_iterations_mediation_only_sleep) = {partial_volume_voxel_tc_band_timeseries_sessionwise/mean(partial_volume_voxel_tc_band_timeseries_sessionwise)};
                            plot_timecourse_table_sleep.csf_tc(num_iterations_mediation_only_sleep) = {csf_column_band_timeseries_sessionwise/mean(csf_column_band_timeseries_sessionwise)};
                            
                            plot_timecourse_table_sleep_unfil.subID(num_iterations_mediation_only_sleep) = string(sub_id);
                            plot_timecourse_table_sleep_unfil.gm_derivative_tc_unfil(num_iterations_mediation_only_sleep) = {gradient(gm_column_sessionwise)/mean(gradient(gm_column_sessionwise))};
                            plot_timecourse_table_sleep_unfil.gm_tc_unfil(num_iterations_mediation_only_sleep) = {gm_column_sessionwise/mean(gm_column_sessionwise)};
                            plot_timecourse_table_sleep_unfil.pfi_derivative_tc_unfil(num_iterations_mediation_only_sleep) = {gradient(pfi_column_sessionwise)/mean(gradient(pfi_column_sessionwise))};
                            plot_timecourse_table_sleep_unfil.pfi_tc_unfil(num_iterations_mediation_only_sleep) = {pfi_column_sessionwise/mean(pfi_column_sessionwise)};
                            plot_timecourse_table_sleep_unfil.csf_tc_unfil(num_iterations_mediation_only_sleep) = {csf_column_sessionwise/mean(csf_column_sessionwise)};

                            num_iterations_mediation_only_sleep = num_iterations_mediation_only_sleep + 1;
                        end
                        if Output_table.sleep_stage{num_iterations} == '0'
                            mediation_table_only_wake.subID(num_iterations_mediation_only_wake) = string(sub_id);
                            mediation_table_only_wake.gm_derivative_tc(num_iterations_mediation_only_wake) = {gm_derivative_sessionwise};
                            mediation_table_only_wake.pfi_derivative_tc(num_iterations_mediation_only_wake) = {partial_volume_derivative_sessionwise};
                            mediation_table_only_wake.csf_tc(num_iterations_mediation_only_wake) = {csf_column_band_timeseries_sessionwise};

                            plot_timecourse_table_wake.subID(num_iterations_mediation_only_wake) = string(sub_id);
                            plot_timecourse_table_wake.gm_derivative_tc(num_iterations_mediation_only_wake) = {gm_derivative_sessionwise};
                            plot_timecourse_table_wake.gm_tc(num_iterations_mediation_only_wake) = {gm_column_band_timeseries_sessionwise};
                            plot_timecourse_table_wake.pfi_derivative_tc(num_iterations_mediation_only_wake) = {partial_volume_derivative_sessionwise};
                            plot_timecourse_table_wake.pfi_tc(num_iterations_mediation_only_wake) = {partial_volume_voxel_tc_band_timeseries_sessionwise};
                            plot_timecourse_table_wake.csf_tc(num_iterations_mediation_only_wake) = {csf_column_band_timeseries_sessionwise};

                            plot_timecourse_table_wake_unfil.subID(num_iterations_mediation_only_wake) = string(sub_id);
                            plot_timecourse_table_wake_unfil.gm_derivative_tc_unfil(num_iterations_mediation_only_wake) = {gradient(gm_column_sessionwise)/mean(gradient(gm_column_sessionwise))};
                            plot_timecourse_table_wake_unfil.gm_tc_unfil(num_iterations_mediation_only_wake) = {gm_column_sessionwise/mean(gm_column_sessionwise)};
                            plot_timecourse_table_wake_unfil.pfi_derivative_tc_unfil(num_iterations_mediation_only_wake) = {gradient(pfi_column_sessionwise)/mean(gradient(pfi_column_sessionwise))};
                            plot_timecourse_table_wake_unfil.pfi_tc_unfil(num_iterations_mediation_only_wake) = {pfi_column_sessionwise/mean(pfi_column_sessionwise)};
                            plot_timecourse_table_wake_unfil.csf_tc_unfil(num_iterations_mediation_only_wake) = {csf_column_sessionwise/mean(csf_column_sessionwise)};



                            num_iterations_mediation_only_wake = num_iterations_mediation_only_wake + 1;
                        end
                        num_iterations = num_iterations + 1;
                        timeslot = timeslot + 3;
                    else
                        timeslot = timeslot + 1;
                    end
                end
                mediation_table.subID(num_iterations_mediation) = string(sub_id);
                mediation_table.session_ID(num_iterations_mediation) = string(session_id);
                mediation_table.gm_derivative_tc(num_iterations_mediation) = {gm_derivative};
                mediation_table.tfi_derivative_tc(num_iterations_mediation) = {partial_volume_derivative};
                mediation_table.csf_tc(num_iterations_mediation) = {csf_column_band_timeseries};
                num_iterations_mediation = num_iterations_mediation + 1;
            end
            
        end


    end
end

%% Export data
writetable(Output_table,Output_table_file_path)
writetable(mediation_table, fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/sleep/mediation_table_06_06_25.csv'))
writetable(mediation_table_only_sleep,fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/sleep/mediation_table_only_sleep_06_06_25.csv'))
writetable(mediation_table_only_wake,fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/sleep/mediation_table_only_wake_06_06_25.csv'))
writetable(plot_timecourse_table_sleep,fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/sleep/plot_timecourse_table_only_sleep_06_06_25.csv'))
writetable(plot_timecourse_table_wake,fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/sleep/plot_timecourse_table_only_wake_06_06_25.csv'))
writetable(plot_timecourse_table_sleep_unfil,fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/sleep/plot_timecourse_table_only_sleep_unfil_12_06_25.csv'))
writetable(plot_timecourse_table_wake_unfil,fullfile('/media/vneumaier/pfi_paper_25/Data_and_scripts_for_publishing/results/sleep/plot_timecourse_table_only_wake_unfil_12_06_25.csv'))