function [single_trial_data] = F_2023_11_29_ExtractNameDateFromFolder(single_trial_folder_subdir_name)
    
% 1. Extract z_scores and file name properties of z-scores from excel directory
%     cd(master_folder);
    full_trial_folder_name = single_trial_folder_subdir_name;
    
    % A. Get CG and ms starting and stopping indices:
    cg_start_ind = strfind(full_trial_folder_name, 'CG'); 
        if isempty(cg_start_ind)
            cg_start_ind = strfind(full_trial_folder_name, 'cg');    % Get all caps versions of cg, Cg or CG
            if isempty(cg_start_ind)
                cg_start_ind = strfind(full_trial_folder_name, 'Cg');
            end
        end
    ms_start_ind = strfind(full_trial_folder_name, 'ms');    
    ms_stop_ind = strfind(full_trial_folder_name(ms_start_ind:end), '_');
    ms_stop_ind = ms_stop_ind(1) + ms_start_ind - 1;
    cg_ms_base_name = full_trial_folder_name(cg_start_ind:ms_stop_ind);

    % B. Get DateTime information:
    datetime_start_ind = strfind(full_trial_folder_name, 'DA');
    datetime_stop_ind = datetime_start_ind + 12;
    datetime_string = full_trial_folder_name(datetime_start_ind:datetime_stop_ind);
    datetime_string(7) = 'm';
    datetime_string(10) = 'd';

    % C. Get CO - condition information
    cond_start_ind = strfind(full_trial_folder_name, 'CO');
    cond_substring = full_trial_folder_name(cond_start_ind:end);
    cond_stop_ind = strfind(cond_substring, '_');
    cond_stop_ind = cond_start_ind  + cond_stop_ind(1) - 1;
    cond_string = full_trial_folder_name(cond_start_ind:cond_stop_ind);

% 2. Put name together into base-string to read file-names with:
    full_basename_string = strcat(datetime_string, cg_ms_base_name, cond_string);

% 3. Save trial names/properties in trial_list struct of STRINGS:
    single_trial_data.date = datetime_string;
    single_trial_data.mouse = cg_ms_base_name;
    single_trial_data.cond = cond_string;
    single_trial_data.trial_name = full_basename_string;
    single_trial_data.folder = full_trial_folder_name;






    