function [single_trial_data] = F_2023_11_29_ImportSlitTagMat(single_trial_folder, single_trial_data, slittag_bhv_score_present, sec_gap_iso_licks)
% [slittag_cell_matlab, slittag_matrix_matlab, slit_frame_diff_at_end] = F_2023_11_29_ImportSlitTagMat(zscr_excel_folder, full_basename_string, total_movie_frames, len_z_scr_trace, slittag_bhv_score_present, sec_gap_iso_licks)
%F_2023_04_18_ImportSlitTagMat Summary of this function goes here

% 0. Properties:
NUM_FRAMES_BHV_DELAY = 0;
NUM_OF_INVIS = 0;
LICK_NUM = 3;
ISOLATED_LICK_NUM = 5;
FRAME_RATE = 15;
% 0.1 Properties:
total_movie_frames = single_trial_data.props_mov.num_frames;
len_z_scr_trace = single_trial_data.props_xls.len_z_scr_trace;

%   Detailed explanation goes here
data_folder = single_trial_folder;
cd(data_folder);
% Get file name of pkl file:
        slit_tag_dir = dir('*2023Jun1*SlitTags.mat');
        if ~isempty(slit_tag_dir) && slittag_bhv_score_present == 1
            
            % A. Load pkl file and convert to matlab array:
            slit_tag_mat = load(slit_tag_dir(1).name);
            slit_tag_timeline = slit_tag_mat.slit_tag_timeline;

            % B. Categories of behavior that you have imported are: 
                % 1. Head-Entry. %2. Low-threshold Lick. %3. High threshold-lick.
                % 4. Define 4th category of behavior that's lick for the first
                % time in sec_gap_iso_licks seconds (first translate this into index):
                inds_of_licks = find(slit_tag_timeline >= LICK_NUM);
                if size(inds_of_licks, 1) > 1
                    prev_lick_ind_diffs = inds_of_licks(2:end) - inds_of_licks(1:(end-1));
                    % Vertcat with position #1 as the FIRST-LICK recorded, it
                    % should have at least sec_gap_iso_licks distance after previous!
                    prev_lick_sec_diffs = vertcat(sec_gap_iso_licks, prev_lick_ind_diffs/FRAME_RATE);
                    % Now, decide if the previous lick sec differences are less
                    % than the FRAME_RATE:
                    isolated_lick_log = prev_lick_sec_diffs > sec_gap_iso_licks;
                    isolated_lick_inds = inds_of_licks(isolated_lick_log);
                    slit_tag_timeline(isolated_lick_inds) = ISOLATED_LICK_NUM;
                end
 
            % C. STEREO frame mismatch with movie frames:
                % i. Delay at beginning of movie
                    % Now, need to take into account the delay of STEREO tags after
                    % the beginning of the movie! So add n frames of tag
                    % 'invisible' or 'back to camera':
                    matrix_delay_to_add = ones(NUM_FRAMES_BHV_DELAY, 1) * NUM_OF_INVIS;
                    slit_tag_timeline_raw_one = vertcat(matrix_delay_to_add, slit_tag_timeline);
                % ii. Truncation of BHV-TAG/STEREO file at end of movie:
                    slit_frame_diff_at_end = total_movie_frames - length(slit_tag_timeline_raw_one);
                    matrix_end_frames_to_add = ones(slit_frame_diff_at_end, 1) * NUM_OF_INVIS;
                    slittag_matrix_matlab = vertcat(slit_tag_timeline_raw_one, matrix_end_frames_to_add);
        else
%             % D. Create arbitrary pkl file based on length of movie!
              slittag_matrix_matlab = ones(total_movie_frames, 1) * NUM_OF_INVIS;
%             predict_pkl_name = strcat(full_basename_string, 'movs');
              slit_frame_diff_at_end = total_movie_frames - length(slittag_matrix_matlab);
%             % 2: Don't need step 2 since this predict_mat is ALREADY aligned to movie
        end
        
%         % 3. Now, save predictions according to named file:
%         save_name_mat = strcat(predict_pkl_name(1:end-4), '.mat');
%         save(save_name_mat, 'predict_matrix_matlab');
%         save_name_excel = strcat(predict_pkl_name(1:end-4), '.xls');
        
            % 4. Convert matrix into cell array and prepare to save:
                % bhv_letter_tags are the following: 1. 'o' other, 2. 'e' entry
                % 3. lick 'l' lick 4. 'L' - Lick 5. 'I' - Isolated-Lick.
            bhv_letter_tags = {'o', 'e', 'p', 'l', 'L', 'I'};
            slittag_cell_matlab = cell(size(slittag_matrix_matlab));
            for i_bhv = 1:length(bhv_letter_tags)
                curr_tag = bhv_letter_tags{i_bhv};
                ind_value = i_bhv - 1;
                log_correct = (slittag_matrix_matlab == ind_value);
                slittag_cell_matlab(log_correct) = {curr_tag}; 
            end      

            % 5. Assign properties to store in trial-list:
            single_trial_data.slittag_cell_matlab = slittag_cell_matlab;
            single_trial_data.slittag_matrix_matlab = slittag_matrix_matlab;

end

