project_dir = fullfile('...'); % Your project dir here
subjectlist = 'subjectlist.txt';

file_path_good_subjects_txt = fullfile(project_dir,'progs',subjectlist); 
segmentation_table = table();
% Open the file for reading
fileID_good_subjects = fopen(file_path_good_subjects_txt, 'r');

% Check if the file was opened successfully
if fileID_good_subjects == -1
    error('Failed to open the file.')
end
num_iterations = 1;
Output_table = table ();
for session = 1:2   
    frewind(fileID_good_subjects)
    while ~feof(fileID_good_subjects)
        % Read the next line from the file to loop through the subjectlist
        sub_id = fgetl(fileID_good_subjects);
        session_id = num2str(session);
        % Check if the line is empty (end of file)
        if isempty(sub_id)
            break; % Exit the loop when the end of the file is reached
        end
    segmentation_csv_HC_filename = fullfile(project_dir,'preprocData',sub_id,strcat('ses-',session_id),'anat/HC/synthseg.vol.csv');
    segmentation_csv_NC_filename = fullfile(project_dir,'preprocData',sub_id,strcat('ses-',session_id),'anat/NC/synthseg.vol.csv');
    if exist(segmentation_csv_NC_filename)
        segmentation_csv_NC = readtable(segmentation_csv_NC_filename);
        segmentation_csv_NC.subject = strcat(string(sub_id),'-NC');
        segmentation_table = [segmentation_table;segmentation_csv_NC];
    end
    
    if exist(segmentation_csv_HC_filename)
        segmentation_csv_HC = readtable(segmentation_csv_HC_filename);
        segmentation_csv_HC.subject = strcat(string(sub_id),'-HC');
        segmentation_table = [segmentation_table;segmentation_csv_HC];
    end
    end
end


writetable(segmentation_table,'...') % Your output dir and filename here