addpath(genpath('/data_august/pro_brain_clearance_scz/software/CoSMoMVPA-master'));
project_dir = '/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/data/Hypercapnia';
subjectlist = 'subjectlist.txt';

file_path_subjectlist = fullfile(project_dir,'progs',subjectlist); 

% Open the file for reading
fileID_good_subjects = fopen(file_path_subjectlist, 'r');
Output_table = table();
% Check if the file was opened successfully
if fileID_good_subjects == -1
    error('Failed to open the file.')
end
num_iterations = 1;
for session = 1   
    frewind(fileID_good_subjects)
    while ~feof(fileID_good_subjects)
        % Read the next line from the file to loop through the subjectlist
        sub_id = fgetl(fileID_good_subjects);
        session_id = num2str(session);
        % Check if the line is empty (end of file)
        if isempty(sub_id)
            break; % Exit the loop when the end of the file is reached
        end
        
        fmri_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/filtered_func_data_without_inflow_slices.nii.gz');
        mean_fmri_filename = fullfile(project_dir,'/preprocData',sub_id,strcat('ses-',session_id),'/func/frs001_adv_search.feat/mean_func_without_inflow_slices.nii.gz');
        ventricle_mask_filename = fullfile(project_dir,'/masks',sub_id,strcat('ses-',session_id),'/lateral_ventricle_border_mask_fmri_space_with_pve.nii.gz');
        fileID_fmri_file = fopen(fmri_filename, 'r');
        if fileID_fmri_file == -1 || (exist(ventricle_mask_filename) == 0)
            disp('Subject is not completely preprocessed')
            continue
        end
        
        Normocapnia_Start = 25;
        Normocapnia_End = 45;

        Hypercapnia_Start = 65;
        Hypercapnia_End = 85;

        fmri = cosmo_fmri_dataset(fmri_filename,'nifti_form','sform');
        mean_fmri = cosmo_fmri_dataset(mean_fmri_filename,'nifti_form','sform');
        ventricle_border_mask = cosmo_fmri_dataset(ventricle_mask_filename,'nifti_form','sform');
        ventricle_border_idx = find(ventricle_border_mask.samples > 0.9); 
        ventricle_border_tc = mean(fmri.samples(:,ventricle_border_idx),2);
        ventricle_border_tc_normalized = ventricle_border_tc/mean(ventricle_border_tc);

        PFI_median_subtraction_nc_hc_normalized = median(ventricle_border_tc_normalized(Normocapnia_Start:Normocapnia_End))-median(ventricle_border_tc_normalized(Hypercapnia_Start:Hypercapnia_End));
        


        Output_table.sub_id(num_iterations) = string(sub_id);
        Output_table.session_id(num_iterations) = string(session_id);
        Output_table.PFI_median_subtraction_nc_hc_normalized(num_iterations) = PFI_median_subtraction_nc_hc_normalized;
        Output_table.ventricle_border_tc_normalized(num_iterations) = {ventricle_border_tc_normalized};
        num_iterations = num_iterations + 1;
    end
end

writetable(Output_table,fullfile('/media/mbonhoeffer/pfi_paper_25/Data_and_scripts_for_publishing/results/Hypercapnia/Nonramp_PFI_table_23_05_25.csv'))