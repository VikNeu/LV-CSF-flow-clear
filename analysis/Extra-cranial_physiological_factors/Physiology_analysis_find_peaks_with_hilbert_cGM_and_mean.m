%% Set the toolboxes and directories

addpath(genpath('/data_august/pro_brain_clearance_scz/software/CoSMoMVPA-master'));
project_dir = '/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/data/flicker_physiology';

file_path_subjectlist = fullfile(project_dir,'progs/subjectlist_bids.txt'); 

% Open the file for reading
fileID_good_subjects = fopen(file_path_subjectlist, 'r');
% Check if the file was opened successfully
if fileID_good_subjects == -1
    error('Failed to open the file.')
end

%Initialize Output table and a variable that always tells in which row of
%the table the results are written
num_iterations = 1;
Output_table = table ();
Output_table_cycle_length = table();



%% Loop over subjects and sessions    
for session = 11:12
    frewind(fileID_good_subjects)
    while ~feof(fileID_good_subjects)
        matrix_breathing_ventricle = [];
        matrix_cardiac = [];
        % Read the next line from the file
        sub_id = fgetl(fileID_good_subjects);
        session_id = num2str(session);
        % Check if the line is empty (end of file)
        if isempty(sub_id)
            break; % Exit the loop when the end of the file is reached
        end
        
        %% Set files of fMRIs, masks and supplementary information of every specific subject
        fmri_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/filtered_func_data.nii.gz');
        fmri_mean_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/mean_func.nii.gz');
        ventricular_border_mask_filename = fullfile(project_dir,'masks', sub_id,strcat('ses-',session_id),'lateral_ventricles_mask_fmri_space_border.nii.gz');
        csf_mask_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/slice_zero_mask_manual_segmentation.nii.gz');
        gm_mask_filename = fullfile(project_dir,'masks', sub_id,strcat('ses-',session_id),'GM_mask_func_space_bin_no_lvs_or_subcort.nii.gz');
        extraction_timeframe_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/timeframe_borders.txt');
        fileID_fmri_file = fopen(fmri_filename, 'r');
        physiology_table_filename_rest_session = fullfile(project_dir,'/rawData',sub_id,strcat('ses-',session_id),'/func',strcat(sub_id,'_task-rest_run-',sprintf('%02d',session),'_events.physio.tsv'));

        
        if exist(fmri_filename) && exist(physiology_table_filename_rest_session ) && exist(ventricular_border_mask_filename) && exist(csf_mask_filename) && exist(gm_mask_filename)
            
             % Your c:3ode here using sub_id
             
            physiology_table = readtable(physiology_table_filename_rest_session, 'FileType','text','Delimiter','\t');
            
            
            %% Load fMRIs and mask
            fmri = cosmo_fmri_dataset(fmri_filename,'nifti_form','sform');
            ventricular_border_mask = cosmo_fmri_dataset(ventricular_border_mask_filename,'nifti_form','sform');
            csf_mask = cosmo_fmri_dataset(csf_mask_filename,'nifti_form','sform');
            gm_mask = cosmo_fmri_dataset(gm_mask_filename,'nifti_form','sform');
            fmri_mean = cosmo_fmri_dataset(fmri_mean_filename,'nifti_form','sform');
            ventricular_border_mask_idx = find(ventricular_border_mask.samples > 0.9);
            csf_mask_idx = find(csf_mask.samples > 0.9);
            gm_mask_idx = find(gm_mask.samples > 0.9);


            
            %% Specify TR and samling frequencies of physiological signals (specified in Openneuro description) 
            TR = 0.227;
            fs_resp = 1000;
            fs_cardiac = 1000;
            fs_fmri = 1/TR;

            %% Bring the physiology signal and the fMRI signal to the same timeline
            timeframe_borders = importdata(extraction_timeframe_filename);
            
            time_displacement_start = ((timeframe_borders(1,1)-1)*TR); %If the fmri got cropped because of movement
            time_displacement_end = ((timeframe_borders(1,2)-1)*TR);

            time_displacement_start_physiology_TR = (time_displacement_start/(1/fs_cardiac));
            time_displacement_end_physiology_TR = (time_displacement_end/(1/fs_cardiac));
           
            if time_displacement_end_physiology_TR > size(physiology_table,1)
                physiology_table = physiology_table(time_displacement_start_physiology_TR:end,:); 
            else
                physiology_table = physiology_table(time_displacement_start_physiology_TR:round(time_displacement_end_physiology_TR),:);
            end
    
    
            %% Detrend, bandpass and calculate the hilbert function of the respiratory signal
            respiration = table2array(physiology_table(:,2));
            respiration = detrend (respiration);
            filtered_respiration = bandpass(respiration,[0.1,0.4],fs_resp); %was 0.1 and 0.4
            analyticSignal_respiration = hilbert(filtered_respiration);
            respPhase = angle(analyticSignal_respiration);
            
            t_resp = (0:length(filtered_respiration)-1)/fs_resp;
            
            %% Detrend, bandpass and calculate the hilbert function of the cardiac signal
            cardiac = table2array(physiology_table(:,1));
            cardiac = detrend (cardiac);
            filtered_cardiac = bandpass(cardiac,[0.6,1.3],fs_cardiac);
            analyticSignal_cardiac = hilbert(filtered_cardiac);
            cardiacPhase = angle(analyticSignal_cardiac);
            
            t_cardiac = (0:length(filtered_cardiac)-1)/fs_cardiac;
            
            
            %% Get CSF and GM Timecourse
            csf_tc = ((mean(fmri.samples(:,csf_mask_idx),2)/mean(fmri.samples(:,csf_mask_idx),"all")-1)*100); % The minus one is for getting % Signal Change
            %csf_tc = ((mean(fmri.samples(:,csf_mask_idx),2)/median(fmri.samples(:,csf_mask_idx),"all")-1)*100); % Testing MEDIAN!
            gm_tc = ((mean(fmri.samples(:,gm_mask_idx),2)/mean(fmri.samples(:,gm_mask_idx),"all")-1)*100); % The minus one is for getting % Signal Change

            %% PFI mask is the complete ventricular border mask
            partial_volume_voxel_indices = ventricular_border_mask_idx;
            
            %Extract the PFI timecourse and normalize it
            ventricular_border_tc = ((mean(fmri.samples(:,partial_volume_voxel_indices),2)/mean(fmri.samples(:,partial_volume_voxel_indices),"all")-1)*100); % The minus one is for getting % Signal Change

            %% Detrend timecourses and get gradient of PFI timecourse

            csf_tc = detrend(csf_tc);
            gm_tc = detrend(gm_tc);
            gm_tc_gradient = gradient(gm_tc);
            ventricular_border_tc = detrend(ventricular_border_tc);
            ventricular_border_tc_gradient = gradient(ventricular_border_tc);
               
    
    
            %% Resample respiratory signal to match the fMRI TR and find the peaks of the hilbert function to calculate the average breathing cycle length
            t_fmri = (0:length(ventricular_border_tc)-1)/fs_fmri;
            
            respPhase_resampled = interp1(t_resp,respPhase,t_fmri,'linear','extrap');
            [pks_breathing,locs] = findpeaks(respPhase_resampled,'MinPeakDistance',8);
            average_breathing_cycle_length_TR = length(respPhase_resampled)/length(pks_breathing);
            average_breathing_cycle_length = average_breathing_cycle_length_TR * TR;
                

            %% Bin the hilbert output in 20 degree bins and get the mean PFI signal of those bins
    
                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_ventricle_respiration] = histc(respPhase_resampled, phaseEdges);
                binned_fmri_mean_ventricle_respiration = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_ventricle_respiration = ventricular_border_tc(binIndices_ventricle_respiration == bin);
                    if ~isempty(dataInBin_ventricle_respiration)
                        binned_fmri_mean_ventricle_respiration{bin} = dataInBin_ventricle_respiration;
                    else
                        binned_fmri_mean_ventricle_respiration(bin) = NaN;
                    end
                end
                
                %same for CSF
                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_csf_respiration] = histc(respPhase_resampled, phaseEdges);
                binned_fmri_mean_csf_respiration = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_csf_respiration = csf_tc(binIndices_csf_respiration == bin);
                    if ~isempty(dataInBin_csf_respiration)
                        binned_fmri_mean_csf_respiration{bin} = dataInBin_csf_respiration;
                    else
                        binned_fmri_mean_csf_respiration(bin) = NaN;
                    end
                end

                %same for PFI gradient
                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_ventricle_gradient_respiration] = histc(respPhase_resampled, phaseEdges);
                binned_fmri_mean_ventricle_gradient_respiration = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_ventricle_gradient_respiration = ventricular_border_tc_gradient(binIndices_ventricle_gradient_respiration == bin);
                    if ~isempty(dataInBin_ventricle_gradient_respiration)
                        binned_fmri_mean_ventricle_gradient_respiration{bin} = dataInBin_ventricle_gradient_respiration;
                    else
                        binned_fmri_mean_ventricle_gradient_respiration(bin) = NaN;
                    end
                end
    

                %same for GM
                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_gm_respiration] = histc(respPhase_resampled, phaseEdges);
                binned_fmri_mean_gm_respiration = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_gm_respiration = gm_tc(binIndices_gm_respiration == bin);
                    if ~isempty(dataInBin_gm_respiration)
                        binned_fmri_mean_gm_respiration{bin} = dataInBin_gm_respiration;
                    else
                        binned_fmri_mean_gm_respiration(bin) = NaN;
                    end
                end
                
                %same for GM gradient
                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_gm_gradient_respiration] = histc(respPhase_resampled, phaseEdges);
                binned_fmri_mean_gm_gradient_respiration = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_gm_gradient_respiration = gm_tc_gradient(binIndices_gm_gradient_respiration == bin);
                    if ~isempty(dataInBin_gm_gradient_respiration)
                        binned_fmri_mean_gm_gradient_respiration{bin} = dataInBin_gm_gradient_respiration;
                    else
                        binned_fmri_mean_gm_gradient_respiration(bin) = NaN;
                    end
                end
    
    
            %% Now the same for the cardiac signal
                cardiacPhase_resampled = interp1(t_cardiac,cardiacPhase,t_fmri,'linear','extrap');
                
                [pks_cardiac,locs] = findpeaks(cardiacPhase_resampled,'MinPeakDistance',2);
                average_cardiac_cycle_length_TR = length(cardiacPhase_resampled)/length(pks_cardiac);
                average_cardiac_cycle_length = average_cardiac_cycle_length_TR * TR;


                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_ventricle_cardiac] = histc(cardiacPhase_resampled, phaseEdges);
                binned_fmri_mean_ventricle_cardiac = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_ventricle_cardiac = ventricular_border_tc(binIndices_ventricle_cardiac == bin);
                    if ~isempty(dataInBin_ventricle_cardiac)
                        binned_fmri_mean_ventricle_cardiac{bin} = dataInBin_ventricle_cardiac;
                    else
                        binned_fmri_mean_ventricle_cardiac(bin) = NaN;
                    end
                end
                
               
                
                %same for CSF
                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_csf_cardiac] = histc(cardiacPhase_resampled, phaseEdges);
                binned_fmri_mean_csf_cardiac = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_csf_cardiac = csf_tc(binIndices_csf_cardiac == bin);
                    if ~isempty(dataInBin_csf_cardiac)
                        binned_fmri_mean_csf_cardiac{bin} = dataInBin_csf_cardiac;
                    else
                        binned_fmri_mean_csf_cardiac(bin) = NaN;
                    end
                end
                
                 %same for PFI gradient
                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_ventricle_gradient_cardiac] = histc(cardiacPhase_resampled, phaseEdges);
                binned_fmri_mean_ventricle_gradient_cardiac = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_ventricle_gradient_cardiac = ventricular_border_tc_gradient(binIndices_ventricle_gradient_cardiac == bin);
                    if ~isempty(dataInBin_ventricle_gradient_cardiac)
                        binned_fmri_mean_ventricle_gradient_cardiac{bin} = dataInBin_ventricle_gradient_cardiac;
                    else
                        binned_fmri_mean_csf_cardiac(bin) = NaN;
                    end
                end
                
                
                %same for GM
                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_gm_cardiac] = histc(cardiacPhase_resampled, phaseEdges);
                binned_fmri_mean_gm_cardiac = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_gm_cardiac = gm_tc(binIndices_gm_cardiac == bin);
                    if ~isempty(dataInBin_gm_cardiac)
                        binned_fmri_mean_gm_cardiac{bin} = dataInBin_gm_cardiac;
                    else
                        binned_fmri_mean_gm_cardiac(bin) = NaN;
                    end
                end

                %same for GM gradient
                numBins = 18;
                phaseEdges = linspace(-pi,pi,numBins+1);
                [~, binIndices_gm_gradient_cardiac] = histc(cardiacPhase_resampled, phaseEdges);
                binned_fmri_mean_gm_gradient_cardiac = {};
                alpha = 0.05;
                for bin = 1:numBins
                    dataInBin_gm_gradient_cardiac = gm_tc_gradient(binIndices_gm_gradient_cardiac == bin);
                    if ~isempty(dataInBin_gm_gradient_cardiac)
                        binned_fmri_mean_gm_gradient_cardiac{bin} = dataInBin_gm_gradient_cardiac;
                    else
                        binned_fmri_mean_gm_gradient_cardiac(bin) = NaN;
                    end
                end



                Output_table.sub_id(num_iterations) = string(sub_id);
                Output_table.session(num_iterations) = string(session_id);
                Output_table.matrix_breathing_ventricle(num_iterations) = {binned_fmri_mean_ventricle_respiration};
                Output_table.matrix_cardiac_ventricle(num_iterations) = {binned_fmri_mean_ventricle_cardiac};
                Output_table.matrix_breathing_csf(num_iterations) = {binned_fmri_mean_csf_respiration};
                Output_table.matrix_cardiac_csf(num_iterations) = {binned_fmri_mean_csf_cardiac};
                Output_table.matrix_breathing_ventricle_gradient(num_iterations) = {binned_fmri_mean_ventricle_gradient_respiration};
                Output_table.matrix_cardiac_ventricle_gradient(num_iterations) = {binned_fmri_mean_ventricle_gradient_cardiac};
                Output_table.matrix_breathing_gm(num_iterations) = {binned_fmri_mean_gm_respiration};
                Output_table.matrix_cardiac_gm(num_iterations) = {binned_fmri_mean_gm_cardiac};
                Output_table.matrix_breathing_gm_gradient(num_iterations) = {binned_fmri_mean_gm_gradient_respiration};
                Output_table.matrix_cardiac_gm_gradient(num_iterations) = {binned_fmri_mean_gm_gradient_cardiac};
                
                Output_table_cycle_length.average_breathing_cycle_length(num_iterations) = average_breathing_cycle_length;
                Output_table_cycle_length.average_cardiac_cycle_length(num_iterations) = average_cardiac_cycle_length;
                num_iterations = num_iterations + 1;
                disp(sub_id)
                disp(session_id)

                % Get sample data from the first subject
                if strcmp(sub_id, 'sub-01') && strcmp(session_id,'11')

                    Output_matrix_timecourse = [t_fmri',cardiacPhase_resampled',respPhase_resampled',ventricular_border_tc_gradient, csf_tc,gm_tc_gradient];
                    filename_exemplary_tc = strcat('exemplary_tc_sub-01_ses-01.xls');
                    writematrix(Output_matrix_timecourse,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_exemplary_tc))
                elseif strcmp(sub_id, 'sub-05') && strcmp(session_id,'11')
                    Output_matrix_timecourse = [t_fmri',cardiacPhase_resampled',respPhase_resampled',ventricular_border_tc_gradient, csf_tc,gm_tc_gradient];
                    filename_exemplary_tc = strcat('exemplary_tc_sub-05_ses-01.xls');
                    writematrix(Output_matrix_timecourse,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_exemplary_tc))
                end
        end
        
    end

end


% Shuffle the data to create null distribution for every condition
n_subjects = size(Output_table,1);
shuffles = 10000;

% For Cardiac CSF
null_dist_cardiac_csf = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_cardiac_csf;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_cardiac_csf(shuffle_iter, :) = shuffled_mean_signal;
end

% For Cardiac PFI
null_dist_cardiac_pfi = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_cardiac_ventricle;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_cardiac_pfi(shuffle_iter, :) = shuffled_mean_signal;
end

% For Cardiac d/dt PFI
null_dist_cardiac_d_dt_pfi = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_cardiac_ventricle_gradient;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_cardiac_d_dt_pfi(shuffle_iter, :) = shuffled_mean_signal;
end

% For Cardiac GM
null_dist_cardiac_gm = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_cardiac_gm;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_cardiac_gm(shuffle_iter, :) = shuffled_mean_signal;
end

% For Cardiac d/dt GM
null_dist_cardiac_d_dt_gm = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_cardiac_gm_gradient;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_cardiac_d_dt_gm(shuffle_iter, :) = shuffled_mean_signal;
end



%----
% For Respiration CSF
null_dist_respiration_csf = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_breathing_csf;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_respiration_csf(shuffle_iter, :) = shuffled_mean_signal;
end

% For Respiration PFI
null_dist_respiration_pfi = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_breathing_ventricle;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_respiration_pfi(shuffle_iter, :) = shuffled_mean_signal;
end

% For Respiration d/dt PFI
null_dist_respiration_d_dt_pfi = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_breathing_ventricle_gradient;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_respiration_d_dt_pfi(shuffle_iter, :) = shuffled_mean_signal;
end

% For Respiration GM
null_dist_respiration_gm = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_breathing_gm;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_respiration_gm(shuffle_iter, :) = shuffled_mean_signal;
end

% For Respiration d/dt GM
null_dist_respiration_d_dt_gm = zeros(shuffles,numBins);
for shuffle_iter = 1:shuffles

    data_to_shuffle = Output_table.matrix_breathing_gm_gradient;
    shuffled_data = cell(n_subjects,1);

    for subject = 1:n_subjects
        shuffled_bins = randperm(numBins); % Randomly shuffle bin indices
        shuffled_data{subject} = data_to_shuffle{subject}(shuffled_bins);
    end

    shuffled_means = cellfun(@(x) cellfun(@mean, x), shuffled_data, 'UniformOutput', false);
    shuffled_mean_signal = mean(cell2mat(shuffled_means), 1);
    null_dist_respiration_d_dt_gm(shuffle_iter, :) = shuffled_mean_signal;
end

% Save the 6 null distributions
filename_null_cardiac_csf = strcat('null_cardiac_csf_31_07.xls');
filename_null_cardiac_pfi = strcat('null_cardiac_pfi_31_07.xls');
filename_null_cardiac_d_dt_pfi = strcat('null_cardiac_d_dt_pfi_31_07.xls');
filename_null_cardiac_gm = strcat('null_cardiac_gm_31_07.xls');
filename_null_cardiac_d_dt_gm = strcat('null_cardiac_d_dt_gm_31_07.xls');
filename_null_respiration_csf = strcat('null_respiration_csf_31_07.xls');
filename_null_respiration_pfi = strcat('null_respiration_pfi_31_07.xls');
filename_null_respiration_d_dt_pfi = strcat('null_respiration_d_dt_pfi_31_07.xls');
filename_null_respiration_gm = strcat('null_respiration_gm_31_07.xls');
filename_null_respiration_d_dt_gm = strcat('null_respiration_d_dt_gm_31_07.xls');
%--
writematrix(null_dist_cardiac_csf,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_cardiac_csf))
writematrix(null_dist_cardiac_pfi,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_cardiac_pfi))
writematrix(null_dist_cardiac_d_dt_pfi,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_cardiac_d_dt_pfi))
writematrix(null_dist_cardiac_gm,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_cardiac_gm))
writematrix(null_dist_cardiac_d_dt_gm,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_cardiac_d_dt_gm))
writematrix(null_dist_respiration_csf,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_respiration_csf))
writematrix(null_dist_respiration_pfi,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_respiration_pfi))
writematrix(null_dist_respiration_d_dt_pfi,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_respiration_d_dt_pfi))
writematrix(null_dist_respiration_gm,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_respiration_gm))
writematrix(null_dist_respiration_d_dt_gm,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_null_respiration_d_dt_gm))

% Save the output table
filename_output_table_cycle_length = strcat('Output_table_cycle_length_31_07_25.xls');
writetable(Output_table_cycle_length,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/',filename_output_table_cycle_length))




Final_table_respiration = table();
frewind(fileID_good_subjects)
%% Loop over subjects to extract the means
while ~feof(fileID_good_subjects)
    % Read the next line from the file
    sub_id = fgetl(fileID_good_subjects);
    sub_id_table_ventricle = strcat(strrep(sub_id,'-','_'),'_ventricle');
    sub_id_table_csf = strcat(strrep(sub_id,'-','_'),'_csf');
    sub_id_table_ventricle_gradient = strcat(strrep(sub_id,'-','_'),'_ventricle_gradient');
    sub_id_table_gm = strcat(strrep(sub_id,'-','_'),'_gm');
    sub_id_table_gm_gradient = strcat(strrep(sub_id,'-','_'),'_gm_gradient');
    session_id = num2str(session);
    % Check if the line is empty (end of file)
    if isempty(sub_id)
        break; % Exit the loop when the end of the file is reached
    end
    
    % Create new table with only the timecourses of one specific
    % subject
    Output_table_subject = Output_table(Output_table.sub_id == string(sub_id),:);

    if isempty(Output_table_subject)
        continue
    end

    Output_subject_ventricle = [];
    Output_subject_csf = [];
    Output_subject_ventricle_gradient = [];
    Output_subject_gm = [];
    Output_subject_gm_gradient = [];

    for i = 1:size(Output_table_subject,1)
        Output_subject_ventricle = [Output_subject_ventricle; Output_table_subject.matrix_breathing_ventricle{i}];
    end

    for i = 1:size(Output_table_subject,1)
        Output_subject_csf = [Output_subject_csf; Output_table_subject.matrix_breathing_csf{i}];
    end

    for i = 1:size(Output_table_subject,1)
        Output_subject_ventricle_gradient = [Output_subject_ventricle_gradient; Output_table_subject.matrix_breathing_ventricle_gradient{i}];
    end
    
    for i = 1:size(Output_table_subject,1)
        Output_subject_gm = [Output_subject_gm; Output_table_subject.matrix_breathing_gm{i}];
    end

    for i = 1:size(Output_table_subject,1)
        Output_subject_gm_gradient = [Output_subject_gm_gradient; Output_table_subject.matrix_breathing_gm_gradient{i}];
    end

    %% Create cells to be able to mean over the subject
    Output_subject_ventricle_final = cell(1,numBins);
    Output_subject_csf_final = cell(1,numBins);
    Output_subject_ventricle_gradient_final = cell(1,numBins);
    Output_subject_gm_final = cell(1,numBins);
    Output_subject_gm_gradient_final = cell(1,numBins);
    median_tc_ventricular_border_respiration_subject = zeros(1,numBins);
    median_tc_csf_respiration_subject = zeros(1,numBins);
    median_tc_ventricular_border_gradient_respiration_subject = zeros(1,numBins);
    median_tc_gm_respiration_subject = zeros(1,numBins);
    median_tc_gm_gradient_respiration_subject = zeros(1,numBins);
    for i = 1:size(Output_subject_ventricle,2)
        for j = 1:size(Output_subject_ventricle,1)
            Output_subject_ventricle_final{1,i} = [Output_subject_ventricle_final{i};Output_subject_ventricle{j,i}];
            Output_subject_csf_final{1,i} = [Output_subject_csf_final{i};Output_subject_csf{j,i}];
            Output_subject_ventricle_gradient_final{1,i} = [Output_subject_ventricle_gradient_final{i};Output_subject_ventricle_gradient{j,i}];
            Output_subject_gm_final{1,i} = [Output_subject_gm_final{i};Output_subject_gm{j,i}];
            Output_subject_gm_gradient_final{1,i} = [Output_subject_gm_gradient_final{i};Output_subject_gm_gradient{j,i}];

        end
        median_tc_ventricular_border_respiration_subject(1,i) = median(Output_subject_ventricle_final{1,i});
        median_tc_csf_respiration_subject(1,i) = median(Output_subject_csf_final{1,i});
        median_tc_ventricular_border_gradient_respiration_subject(1,i) = median(Output_subject_ventricle_gradient_final{1,i});
        median_tc_gm_respiration_subject(1,i) = median(Output_subject_gm_final{1,i});
        median_tc_gm_gradient_respiration_subject(1,i) = median(Output_subject_gm_gradient_final{1,i});

    end

% because of 18 bins
    for j = 1:18 
        Final_table_respiration.(sub_id_table_ventricle) (j) = median_tc_ventricular_border_respiration_subject(j);
        Final_table_respiration.(sub_id_table_csf) (j) = median_tc_csf_respiration_subject(j);
        Final_table_respiration.(sub_id_table_ventricle_gradient) (j) = median_tc_ventricular_border_gradient_respiration_subject(j);
        Final_table_respiration.(sub_id_table_gm) (j) = median_tc_gm_respiration_subject(j);
        Final_table_respiration.(sub_id_table_gm_gradient) (j) = median_tc_gm_gradient_respiration_subject(j);
    end
    
end
%% Export results
writetable(Final_table_respiration,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/Final_table_respiration_31_07_25.xls'))

%% Do the same for cardiac
frewind(fileID_good_subjects)
Final_table_cardiac = table();
    while ~feof(fileID_good_subjects)
        % Read the next line from the file
        sub_id = fgetl(fileID_good_subjects);
        sub_id_table_ventricle = strcat(strrep(sub_id,'-','_'),'_ventricle');
        sub_id_table_csf = strcat(strrep(sub_id,'-','_'),'_csf');
        sub_id_table_ventricle_gradient = strcat(strrep(sub_id,'-','_'),'_ventricle_gradient');
        sub_id_table_gm = strcat(strrep(sub_id,'-','_'),'_gm');
        sub_id_table_gm_gradient = strcat(strrep(sub_id,'-','_'),'_gm_gradient');
        session_id = num2str(session);
        if isempty(sub_id)
            break; % Exit the loop when the end of the file is reached
        end
        Output_table_subject = Output_table(Output_table.sub_id == string(sub_id),:);

        if isempty(Output_table_subject)
            continue
        end

        Output_subject_ventricle = [];
        Output_subject_csf = [];
        Output_subject_ventricle_gradient = [];
        Output_subject_gm = [];
        Output_subject_gm_gradient = [];

        for i = 1:size(Output_table_subject,1)
            Output_subject_ventricle = [Output_subject_ventricle; Output_table_subject.matrix_cardiac_ventricle{i}];
        end

        for i = 1:size(Output_table_subject,1)
            Output_subject_csf = [Output_subject_csf; Output_table_subject.matrix_cardiac_csf{i}];
        end

        for i = 1:size(Output_table_subject,1)
            Output_subject_ventricle_gradient = [Output_subject_ventricle_gradient; Output_table_subject.matrix_cardiac_ventricle_gradient{i}];
        end

        for i = 1:size(Output_table_subject,1)
            Output_subject_gm = [Output_subject_gm; Output_table_subject.matrix_cardiac_gm{i}];
        end

        for i = 1:size(Output_table_subject,1)
            Output_subject_gm_gradient = [Output_subject_gm_gradient; Output_table_subject.matrix_cardiac_gm_gradient{i}];
        end


        Output_subject_ventricle_final = cell(1,numBins);
        Output_subject_csf_final = cell(1,numBins);
        Output_subject_ventricle_gradient_final = cell(1,numBins);
        Output_subject_gm_final = cell(1,numBins);
        Output_subject_gm_gradient_final = cell(1,numBins);
        median_tc_ventricular_border_cardiac_subject = zeros(1,numBins);
        median_tc_csf_cardiac_subject = zeros(1,numBins);
        median_tc_ventricular_border_gradient_cardiac_subject = zeros(1,numBins);
        median_tc_gm_cardiac_subject = zeros(1,numBins);
        median_tc_gm_gradient_cardiac_subject = zeros(1,numBins);
        for i = 1:size(Output_subject_ventricle,2)
            for j = 1:size(Output_subject_ventricle,1)
                Output_subject_ventricle_final{1,i} = [Output_subject_ventricle_final{i};Output_subject_ventricle{j,i}];
                Output_subject_csf_final{1,i} = [Output_subject_csf_final{i};Output_subject_csf{j,i}];
                Output_subject_ventricle_gradient_final{1,i} = [Output_subject_ventricle_gradient_final{i};Output_subject_ventricle_gradient{j,i}];
                Output_subject_gm_final{1,i} = [Output_subject_gm_final{i};Output_subject_gm{j,i}];
                Output_subject_gm_gradient_final{1,i} = [Output_subject_gm_gradient_final{i};Output_subject_gm_gradient{j,i}];
            end
            median_tc_ventricular_border_cardiac_subject(1,i) = median(Output_subject_ventricle_final{1,i});
            median_tc_csf_cardiac_subject(1,i) = median(Output_subject_csf_final{1,i});
            median_tc_ventricular_border_gradient_cardiac_subject(1,i) = median(Output_subject_ventricle_gradient_final{1,i});
            median_tc_gm_cardiac_subject(1,i) = median(Output_subject_gm_final{1,i});
            median_tc_gm_gradient_cardiac_subject(1,i) = median(Output_subject_gm_gradient_final{1,i});
        end       

        for j = 1:18
            Final_table_cardiac.(sub_id_table_ventricle) (j) = median_tc_ventricular_border_cardiac_subject(j);
            Final_table_cardiac.(sub_id_table_csf) (j) = median_tc_csf_cardiac_subject(j);
            Final_table_cardiac.(sub_id_table_ventricle_gradient) (j) = median_tc_ventricular_border_gradient_cardiac_subject(j);
            Final_table_cardiac.(sub_id_table_gm) (j) = median_tc_gm_cardiac_subject(j);
            Final_table_cardiac.(sub_id_table_gm_gradient) (j) = median_tc_gm_gradient_cardiac_subject(j);

        end
        
    end
writetable(Final_table_cardiac,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/flicker_physiology/Final_table_cardiac_31_07_25.xls'))