function [single_trial_data] = F_2023_11_29_ReadMovie(single_trial_folder, single_trial_data, exp_props, movie_read_yes)
% 2021_12_23 Extract traces from excel files:

% 0.1 Enter movie directory, and then get movie files:
cd(single_trial_folder);

% 0.2 Set movie parameters:
movie_frame_total = [];
movie_duration_total = [];
curr_mov_read = [];

    if movie_read_yes == 0
        mov_stats_dir = dir('*mov_stats*');          % Now, try and find a Matlab file storing movie stats that you can load
        if length(mov_stats_dir) >= 1
            mov_to_load_name = mov_stats_dir(1).name;
            mov_stats = load(mov_to_load_name);
            movie_frame_total = mov_stats.movie_frame_total;
            movie_duration_total = mov_stats.movie_duration_total;
        end
    end

    if movie_read_yes == 1 || isempty(movie_frame_total)
        mov_dir = dir('*Cam1*');                     % All Cam1 files in directory to get frame length:
        if length(mov_dir) >= 1
    
            % 1. Establish movie reader:
            file_sizes = mov_dir(:).bytes;
             [~, ind_largest_file] = max(file_sizes);
            mov_to_read_name = mov_dir(ind_largest_file).name;
            curr_mov_read = VideoReader(mov_to_read_name);
            
            % 2. Get total number of frames:
            movie_frame_total = curr_mov_read.NumFrames;   % 
            movie_duration_total = curr_mov_read.Duration;  

            % 3. Now, save these stats in file:
            mov_to_save_name = mov_to_read_name;
            mov_to_save_name(mov_to_save_name == '.') = '_';
            save(strcat(mov_to_save_name, '_mov_stats'), 'movie_frame_total', 'movie_duration_total');

        else
            msgbox('There is no Movie file for this trial', 'Movie-file missing');
        end
    end
    
    % Store trial properties:
    single_trial_data.props_mov.num_frames = movie_frame_total;
    single_trial_data.props_mov.duration_mov = movie_duration_total;
    a = 5;
end