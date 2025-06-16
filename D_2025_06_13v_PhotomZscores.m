%% D_2025_06_13_MultiSite Photometry Signal Analysis for Axonal Correlates of Fiber Photometry paper

% This script will compare correlation values and analyze multi-site
% signals and correlations around peaks
% 
% Input: 
% A) Raw Neurophotometrics Photometry recordings
% B) Optionally: Behavior scores (e.g. ButterSlit-tagged behavior scores).
%  
% Output: Figures for each trial showing: 
%   i) Raw traces of photometry recordings with different pre-processing steps
%   ii) Summary of correlation values
%   iii) All-site signal/correlations at peaks in one channel
%   iv) All-site signal/correlation at peaks/behavioral eventss in one
%   channel

    % Master and function paths
    master_fxn_folder = 'Z:\shared\Manuscripts\AxonalCorrelatesFiberPhotometry\PhotomAnalysisCode\';
    fxn_folder = strcat(master_fxn_folder, '\FUNCTIONS');
    master_folder = 'Z:\shared\Manuscripts\AxonalCorrelatesFiberPhotometry\RawData\iSPN_GC6s';            % Photom012b - iSPN-GC6
    % master_folder = 'Z:\shared\Manuscripts\AxonalCorrelatesFiberPhotometry\RawData\dSPN_GC6s';            % Photom012b - dSPN-GC6
    % master_folder = 'Z:\shared\Manuscripts\AxonalCorrelatesFiberPhotometry\RawData\dSPN_GC8m';             % Photom030 - dSPN-GC8
    addpath(fxn_folder);

    % Set specific folder of data to process:
    input_folder_string = 'InputFiles1';            % Photom012b - iSPN-GC6
    % input_folder_string = 'InputFiles2';          % Photom012b - dSPN-GC6
    % input_folder_string = 'InputFiles3';          % Photom030 - dSPN-GC8

    % Define properties for whole experiment/data-processing:
    exp_props = struct;
    exp_props.side_filt_only = 0;       % Filter all mice so that you only plot 'Right' side!
        % Side filtering:
        exp_props.side_str = '';
        if exp_props.side_filt_only == 1
            exp_props.side_str = 'Right';    % Should be either 'Right' or 'Left'!
        end
    exp_props.combine_hemispheres = 1;       % This should be either 1 or 2, 1 = combine_hemispheres, 2 = separate (default).
    exp_props.stereo_bhv_score_present = 0;   % Should stereo-tags be imported? 0 = no; 1 = yes.
    exp_props.slittag_bhv_score_present = 1;  % Should slittag-tags be imported? 0 = no; 1 = yes.
    
    % Filtering!
    exp_props.event_filt_names = {'Pk', 'STEREO', 'Slit'};
    exp_props.event_filt_cond = 1;            % What to use for event filtering. 1 = photometry z-score, 2 = STEREO/Bhv.
    exp_props.filt_params = {2.0, 20};               % IF empty, then will be defined in either DefinePeakProps or defineBehavioralParms according to: exp_props.event_filt_cond
    %     exp_props.filt_params = [];               % IF empty, then will be defined in either DefinePeakProps or defineBehavioralParms according to: exp_props.event_filt_cond

    % Set initial 'unchanging' parameters for trial_props:
    [exp_props, exp_props.pho] = F_2023_12_04_DefineExpProps(exp_props, input_folder_string);    % Here, set properties containing to z-score/excels/signal-raw-data + behavior
    [exp_props.pk] = F_2023_12_04_DefinePeakProps(exp_props);                     % Here, define peak-table properties
    bhv_choice = 2;                                                     % Set bhv_choice property:
    [exp_props.bhv] = F_2024_04_07_defineBehavioralParams(exp_props, bhv_choice);     % Here, define bhv_props: 1 = STEREO, 2 = Slit-tag, 3 = Combo of STEREO/SLIT-TAG!
    exp_props.output_dir = strcat(master_folder, exp_props.output_date_input, num2str(exp_props.filt_param_str{1})); 

%% 1. Establish data storage arrays for use in program 

    % 1A) File folder paths:
    cd(master_folder);
    folder_dirs.input_master_dir = strcat(master_folder, '\', input_folder_string);
    datetime_ans = datetime;
    datetime_ans.Format = 'yyyy_MM_dd_';
    datetime_str = char(datetime_ans);
    date_str = datetime_str(1:11);
    output_folder_strings = strcat('FigsFrom', input_folder_string,...
        '_', num2str(exp_props.pk.z_prom_params{1}), 'to', num2str(exp_props.pk.z_prom_params{end}), '_', num2str(size(exp_props.pk.z_prom_params, 1)), 'X_');
    
    folder_dirs.output_fig_dir = strcat(master_folder, '\', date_str, '_', output_folder_strings);
    folder_dirs.import_raw_output_fig_dir = strcat(folder_dirs.output_fig_dir, '\', 'RawPlots');

    % 1B) Calculate number of trials to process in designated folder:
    cd(folder_dirs.input_master_dir);
    folder_dirs.name_dir = dir('*DA*');    % "DA" is used to identify folders! Must be in name!
    num_ms_trials = size(folder_dirs.name_dir, 1);
    if num_ms_trials < 1
        fprintf('Not enough correctly labeled files');
        waitforbuttonpress;
    end

%% 2. Import: A) File Name from Folder, B) ZScore(s) C) Movie Reader D) Behavior-Scores % E) Align and Store... 
% all imported z-score traces, bhv_scores, and movie_reader variables into z_trace_bhv_table:

    % A. Set up structs, lists, and tables to store raw data:
        all_trial_data = struct;

        table_var_types = {'double', 'double', 'double', 'double', 'double', 'double', 'double'};
        table_var_names = {'sec_photom_delay_prefilt', 'total_frames_behav', 'total_frames_behav_trunc', 'total_frames_zscore', 'frame_diff_pre', 'sec_photom_delay_postfilt', 'bhv_frame_diff_at_end'};
        time_sync_table = table('Size', [num_ms_trials, length(table_var_names)], 'VariableTypes', table_var_types, 'VariableNames', table_var_names);

    % B. Run through each trial and import and align data for each trial:

        % Make a struct to store trial and corr values:
        filt_opts = [3, 2, 5];  % Which filt-opts to plot (if you want to just use one filtering option, then input a single number):
        filt_opts = 5
        filt_opt_strings = {'PMAT'; 'CurveFit'; 'NoLowPass'; 'MultiStep100sec'; 'MultiStep500sec'};

        % Define structs to contain all_filt_data and all_filt_corrs:
        all_filt_data = struct;
        all_filt_corrs = struct;

        for i_filt_opt = filt_opts  % Which filt-opts to plot (if you want to just use one filtering option, then input a single number):

            % For this particular round of filtering, get the type of filtering:
%             filt_opt = 2;                  filt_opt_strings = {'PMAT'; 'Curve-fit'; 'NoLowPass'; 'MultiStep_100sec'; 'MultiStep_500sec'};
            filt_opt = i_filt_opt;  % Set the filtering option
            type_of_filt = filt_opt_strings{filt_opt};

            % Define trial_corr_struct to contain correlation values:
            corr_field_names = {'corr_vals_L_Str', 'corr_siteorder_L_Str'; 'corr_vals_R_Str', 'corr_siteorder_R_Str';...
                'corr_vals_L_Axon', 'corr_siteorder_L_Axon'; 'corr_vals_R_Axon', 'corr_siteorder_R_Axon'};
            trial_corr_struct_fields = {'trialname'; 'mouse'; 'cond'; 'day'};
            trial_corr_struct_fields = cat(1, trial_corr_struct_fields, corr_field_names(:));
            trial_corr_struct = struct;
            trial_corr_struct.(corr_field_names{1,2}) = {'striatal-left', 'striatal-right', 'axon-left', 'axon-right'};
            trial_corr_struct.(corr_field_names{2,2}) = {'striatal-right', 'striatal-left', 'axon-right', 'axon-left'};
            trial_corr_struct.(corr_field_names{3,2}) = {'axon-left', 'axon-right', 'striatal-left', 'striatal-right'};
            trial_corr_struct.(corr_field_names{4,2}) = {'axon-right', 'axon-left', 'striatal-right', 'striatal-left'};
            
            for i_list = 1:num_ms_trials
    
                % Import variables for each trial:
                single_trial_folder_subdir_name = folder_dirs.name_dir(i_list).name;
                single_trial_folder = strcat(folder_dirs.input_master_dir, '\', single_trial_folder_subdir_name);
    
                %*** FUNCTION - A) Extract Name from folder 
                [single_trial_data] = F_2023_11_29_ExtractNameDateFromFolder(single_trial_folder_subdir_name);
    
                %*** FUNCTION - Bi) Extract_Multisite_ZScores_Excel for Neurophotometrics data!!!
%                     filt_opt = 2;               % 1 = 'PMAT', 2 = 'Curve-fit', 3 = 'NoLowPass', 4 = 'MultiStep_100sec', 5 = 'MultiStep_500sec'
%                     filt_opt_strings = filt_opt_strings = {'PMAT'; 'Curve-fit'; 'NoLowPass'; 'MultiStep_100sec'; 'MultiStep_500sec'};
%                     type_of_filt = filt_opt_strings{filt_opt};

                % Added the following code 2024/12/15:
                folder_dirs.output_fig_dir = strcat(master_folder, '\', output_folder_strings, type_of_filt, '_', date_str);
                folder_dirs.import_raw_output_fig_dir = strcat(folder_dirs.output_fig_dir, type_of_filt, '\', 'RawPlots');
                exp_props.output_dir = folder_dirs.output_fig_dir;

                raw_pks_plot = 0;           % If raw_input_plot == 3, then plot as SVG and Matlab fig, If == 2, then FIG + PNG, If == 1 --> PNG only. If == 0 don't save!
                plt_timeline_struct = struct('num_x_ticks', 9, 'plt_minute_start', 0, 'plt_minute_len', 36);
                pks_or_sigisos = 1;         % Signal+Isos == 1 or Peaks Plot == 2. Only relevant if raw_pks_plot > 0 (which == 'plot on').
                plot_dff_zscr = 'zscr';      % Can either be 'dff' or 'zcr'.  'zcr' is the default.
                fig_output_folder = folder_dirs.import_raw_output_fig_dir;
                site_col_names = {'Region2G', 'Region3G', 'Region4G', 'Region5G'}; 
                color = '470';
                sites_to_plot = [2, 4];     % Which of sites to include on your raw-plots figure
                corr_baseline_site = 2;     % Which of the sites to use as comparison for calculation of correlation values.
                corr_cross = 'no';
                tiled_overlay = 'overlay';     % Should be EITHER 'overlay' OR 'tile'. If not written exactly as 'overlay', default code execution is for 'tile'. If tile, then will also plot baselines
                seg_time = 1;                   % How many minutes to divide whole trace into (basically how many correlation points to have on graph)
                [single_trial_data, trial_corr_struct] = F_2025_06_13b_Extract_Multisite_ZScores_Paper(single_trial_folder, single_trial_data, exp_props, filt_opt, raw_pks_plot, pks_or_sigisos, plot_dff_zscr, site_col_names, sites_to_plot, corr_baseline_site, corr_cross, color, fig_output_folder, trial_corr_struct, corr_field_names, plt_timeline_struct, tiled_overlay, i_list);    
                
                %*** FUNCTION - C) ReadMovie!!!
                movie_read_yes = 0;     % Try to read a 'mov_stats' file, or open movie-reader of 'Cam1' movie: (0 == read stats; 1 == open movie reader).
                [single_trial_data] = F_2023_11_29_ReadMovie(single_trial_folder, single_trial_data, exp_props, movie_read_yes);    % Get new trial_props, with frame and duration
    
                %*** FUNCTION - D) Import STEREO data and TransformPklToBhvMat!!!
                stereo_bhv_score_present = 0;   % bhv_score_present = 1 means to process bhv_tags, if 0 then create according to total_movie_frames
                [single_trial_data] = F_2023_11_29_TransformPklToBhvMat(single_trial_folder, single_trial_data, stereo_bhv_score_present);
    
                %*** FUNCTION - E) Import SLIT-TAG data and TransformPklToBhvMat!!!
                slittag_bhv_score_present = 1;  % bhv_score_present = 1 means to process bhv_tags, if 0 then create according to total_movie_frames
                sec_gap_iso_licks = 5;          % The gap to define (in seconds) between isolated licks
                [single_trial_data] = F_2023_11_29_ImportSlitTagMat(single_trial_folder, single_trial_data, slittag_bhv_score_present, sec_gap_iso_licks);
    
                % Add single trial to all-trials list of trials:
                if i_list == 1
                    all_trial_data = single_trial_data;
                else
                    all_trial_data(i_list) = single_trial_data;
                end
            end

            % Now store the trial-data and the trial-corrs in proper
            % structs:
            all_filt_data.(type_of_filt) = all_trial_data;
            all_filt_corrs.(type_of_filt) = trial_corr_struct;

        end
        
        % Make new struct called trial_list_sync which will have all input
        % sources' timings synchronized according to total frame #'s:
        all_trial_sync = all_trial_data;

        %*** FUNCTION - F) Make plot of correlation values:
        plot_corr_values = 1;   % 
        if plot_corr_values == 1
            % corr_vals_fields_to_proc = {'corr_vals_L_Str'; 'VLS-Left'};
            corr_vals_fields_to_proc = {'corr_vals_R_Str'; 'VLS-Right'};
            % The choice of 'baseline' in corr_vals_fields_to_proc determines what the site of site_to_corr_comp is.
            % So, if corr_vals_fields_to_proc = corr_vals_R_Str, then, the values of site_to_corr_comp are: 1 = R-Str, 2 = L-Str, 3 = R-Axons, 4 = L-Axons
            site_to_corr_comp = 3;  % 1 = Same-site (e.g. striatum) ipsi (self), % 2 = Same-site (e.g. striatum) contra % 3 = Diff site (e.g. axons-SNr) ipsi % 4 =  Diff site (e.g. axons-SNr) contra
            plot_title = strcat('All-Trial-Correlation-Plot', type_of_filt);
            F_2025_06_04_Plot_Correlation_Values(fig_output_folder, all_filt_corrs, filt_opts, filt_opt_strings, exp_props, corr_vals_fields_to_proc, site_to_corr_comp, plot_title)
        end

        run_time_sync = 1;  % Run_time_sync == 0 --> Don't run; run_time_sync == 1 --> Run!!!
        if run_time_sync == 1
            for i_list = 1:num_ms_trials
                %*** FUNCTION - G) Align Behavior Tags and Z-Score Traces to Movie and STORE z-scores and behavior matrices!!!
                % 2023/03/17 photom012 note: len_z_scr_trace is length of z-score trace  starting from beginning of movie if all aligned properly (actually 1sec/1 frames AFTER movie start)
                % 2023/03/17 photom012 note: len(z_scores_allsites) is length of z-score trace starting from beginning of movie if all aligned properly (actually 1sec/1 frames AFTER movie start), minus 450 frames/30sec
                % 2023/03/17 photom012 note: predict_cell_matlab +  predict_matrix_matlab start should already be aligned to length of movie from function D!
                % 2023/03/17 photom012 note: --> So align z_score, behav_score timelines to total_movie_frames from sideCam1 movie!
                [all_trial_sync, time_sync_table] = F_2023_11_29_AlignAndStoreZScrBhv(all_trial_sync, exp_props, time_sync_table, i_list);
            end
        end

        % Create lists according to all trials:
        mice_list = unique(vertcat({all_trial_sync.mouse}),'sorted');
        cond_list = unique(vertcat({all_trial_sync.cond}),'sorted');
        date_list = unique(vertcat({all_trial_sync.date}),'sorted');
        
        if exp_props.combine_hemispheres == 1
            %*** FUNCTION - G) Combine Hemispheres:
            [exp_props, trial_list_sync_hemi, num_ms_trials, time_sync_table] = F_2024_01_29_CombineHemispheres(exp_props, all_trial_sync, time_sync_table);
            all_trial_data_to_event = trial_list_sync_hemi;
        else
            all_trial_data_to_event = all_trial_sync;
        end

        % Finally, delete all variables from Sections 1 & 2 connected to input functions which will
        % clutter workspace for analysis functions:
        clear('all_trial_data', 'filt_opt', 'raw_pks_plot', 'pks_or_sigisos', 'fig_output_folder', 'site_col_names',...
            'color', 'movie_read_yes', 'sec_gap_iso_licks', 'bhv_choice', 'date_str', 'table_var_types');

        % Create hemi_list which is a list of all mouse hemispheres:
        hemi_list = unique(vertcat({all_trial_data_to_event.mouse}),'sorted');
        
%% 3A. Fiber photometry activity at a single site - make peak-control table for target site and make peak-control-table with separate params for other site(s):

% Make tables to store data for peaks and controls:
  cd(master_folder);

for i = 1:1
    % Make a separate peak-table for each of the defined z_prom boundaries:
    % Use separately defined z-prominence boundaries, one for 'target-site' and one 
    % for 'other-site' so you can calculate Jaccard/Success-Rates for separately defined data sets: 
    for i_target_comp = 1:2
        if i_target_comp == 1       % First, create peak-table for target site at defined z-intervals
            z_proms_cell = exp_props.pk.z_prom_params;
        elseif i_target_comp == 2   % Second, create peak-table for comparison site at defined z-intervals
            z_proms_cell = exp_props.pk.z_proms_comp_site;
        end
        filt_min_max = [z_proms_cell{1}, z_proms_cell{2}];
   
        % Make a separate peak-table for each of the defined z_prom boundaries:
        % Use separately defined z-prominence boundaries, one for 'target-site' and one filt_string_names
        % for 'other-site' so you can calculate Jaccard/Success-Rates for separately defined data sets: 
    
        % Loop over all trials
        for i_list = 1:num_ms_trials
    
            % Get a single trial to filter:
            single_trial_traces = all_trial_data_to_event(i_list).z_scores_allsites{:};
            trial_name = all_trial_data_to_event(i_list).trial_name;
            mouse_name = all_trial_data_to_event(i_list).mouse;
    
            % Create a table with peaks defined for each recorded site
            for i_site = 1:exp_props.num_sites
    
                % Run on a for-loop to store each peak table that takes into account
                % all data to create a peak-table defined by defined characteristics:
                target_site = i_site;   % Create a different peak table for each site!
                single_timeline_trace = single_trial_traces(i_site, :);
                
                % *** FUNCTION - A. Now, run makePeakEventTables to filter:
                event_filter_props.event_filt_cond = exp_props.event_filt_cond;
                event_filter_props.filt_string = exp_props.pk.filt_string_names{1};
                event_filter_props.filt_min_max = vertcat(filt_min_max);
                event_filter_props.min_interval_btwn_events = 1;
                event_filter_props.sec_trace_around_event = 6;
                % FUNCTION: F_2023_11_30_makePeakEventTables
                [event_table_single, num_events] = F_2023_11_30_makePeakEventTables(single_timeline_trace, trial_name, mouse_name, exp_props, time_sync_table, event_filter_props, i_list);
                
                %*** FUNCTION - B. makeBhvFiltersForTbl
                % Now add behavioral scores to the table from bhv_predict_mat
                for i_bhv_type = 1:2 
                    if i_bhv_type == 1              % STEREO Data:  
                        behavior_matrix = all_trial_data_to_event(i_list).bhv_predict_mat_allsites{:};
                        bhv_scr_type = 'Stereo';
                    elseif i_bhv_type == 2          % Slit-tag Data:
                        behavior_matrix = all_trial_data_to_event(i_list).bhv_slittag_mat_allsites{:};
                        bhv_scr_type = 'Slittag';
                    end
    %                   2023/08/07: The following line of code is extremely expensive time-wise 
                        % FUNCTION: F_2023_11_30_MakeBhvFiltersForTbl
                        [event_table_single_with_bhv] = F_2023_11_30_MakeBhvFiltersForTbl(event_table_single, behavior_matrix, bhv_scr_type, exp_props);
                        peak_table_bhv_site{i_site} = event_table_single_with_bhv;
                end
                
                if i_target_comp == 1
                    master_event_table_for_target{i_list} = peak_table_bhv_site;
                elseif i_target_comp == 2
%                     master_event_table_for_cmpr{i_list} = peak_table_bhv_site;
                end
            end
        end
    end
end

%% 3B. Compare relationships between fiber-photometry peaks across multiple sites 
% a. Threshold and plot relationships between peaks to prepare for Jaccard plotting:
% b. Compare global correlations:
for i = 1:1

    % Make a correlations_storage vector for global correlations:
    correlations_matrix = zeros(num_ms_trials, exp_props.num_sites);
    corr_lags_matrix = correlations_matrix;

    if exp_props.event_filt_cond == 1
        % Get z_prom_minimum and z_prom maximum:
        z_prom_minimum = exp_props.pk.z_prom_params{1}; 
        z_prom_maximum = exp_props.pk.z_prom_params{2};
        filt_min_max = [z_prom_minimum, z_prom_maximum];
    elseif exp_props.event_filt_cond == 2 || exp_props.event_filt_cond == 3
        filt_min_max = [bhv_scr_filt{1}, bhv_scr_filt{2}];
    end

    for i_list = 1:num_ms_trials

        % Re-load peak-tables for each trial:
        trial_table_site_for_target = master_event_table_for_target{i_list};

        for i_target = 1:exp_props.num_sites

            % A. Import peak-table for both target-site and for all
            % other sites!
            target_site = i_target;                           % This is the number for selecting the target site!
%                 min_prom_of_other_site = z_prom_params{1};
            trial_table_target_site = trial_table_site_for_target(target_site);

            % B. %%% F_2023_04_CompareTwoSiteSignals:
            % This function will both compute correlations:
            corr_reg_cross = 1;         % corr_reg_cross = 1 means --> regular correlation, and corr_reg_cross = 2 --> cross correlation
            sec_lag = 1;                % For cross-correlation - the maximum lag.
            [peak_table_wCorr, event_table_name]...
                = F_2025_06_13_CompareTwoSiteSignals(all_trial_data_to_event, trial_table_target_site, i_target, corr_reg_cross, sec_lag, exp_props, i_list);
            all_trial_data_to_event(i_list).(event_table_name) = {peak_table_wCorr};
        end
    end
end


%% 4. Filter struct according to single mouse or multiple mice and plot!

    % IMPORTANT NOTE: In order to plot Jaccard Index/Success Rate/Failure
    % Rate, MUST plot 'all-mice', which is i_hemi = length(hemi_list) + 1!
    % Otherwise ms_string identification won't work!
    mkdir(exp_props.output_dir)

    for i_hemi = length(hemi_list) + 1
    % for i_hemi = 1:length(hemi_list) + 1
        % Get back into main directory so that subsequent directories created will be created in this directory!:
        cd(exp_props.output_dir);

        % Make mouse-specific directories for saving
        if i_hemi <= length(hemi_list)
            ms_specific_dir = strcat(exp_props.output_dir, '_Ind_Mice\', hemi_list{i_hemi}, exp_props.side_str, exp_props.filt_param_str{1});
    %         bhv_sum_check = cell(3, 3); 

            % Ai. Make logical log_filt, which filters all of
            % all_trial_data_to_event based on 2 conditions: #1) Mouse-name #2) Hemisphere
            ms_string = hemi_list{i_hemi}(1:end);
%             ms_string(ms_string == '_') = '-';
            log_filt = zeros(num_ms_trials, 1);
            for i_list = 1:num_ms_trials
                % Condition #1, include in log_filt if ms_name == ms_string (current mouse in i_hemi loop) 
                    log_filt(i_list) = ~isempty(strfind(all_trial_data_to_event(i_list).mouse, ms_string));
                % Condition #2, set second, check if right hemisphere or not...
                    if exp_props.side_filt_only == 1 && log_filt(i_list) == 1
                        log_filt(i_list) = ~isempty(strfind(all_trial_data_to_event(i_list).mouse, exp_props.side_str)); 
                    end
            end
            all_mice_filtering = 0;     % Send to plotting program the information that you need to plot multiple graphs for each mouse
        end
        % Aii. In this condition DO NOT filter by mouse, but instead run code over ALL mice!
        if i_hemi == (length(hemi_list) + 1)
            ms_specific_dir = strcat(exp_props.output_dir, 'all_mice', exp_props.side_str, exp_props.filt_param_str{1});
            ms_string = 'Allmice';
            log_filt = zeros(num_ms_trials, 1);
            for i_list = 1:num_ms_trials
                log_filt(i_list) = num_ms_trials; 
            end
            % Optionally filter only right-side trials
            if exp_props.side_filt_only == 1
                for i_list = 1:num_ms_trials
                    log_filt(i_list) =~isempty(strfind(all_trial_data_to_event(i_list).mouse, exp_props.side_str)); 
                end
            end
            all_mice_filtering = 1;
        end

        % Now use logical to select only the trials corresponding to the
        % correct mouse:
        data_struct_per_ms = all_trial_data_to_event(logical(log_filt));
        
        % Only execute plotting code if data_struct is populated:
        if ~isempty(data_struct_per_ms)
    
            % Define min and max z-score thresholds for defining the z-score
            % parameters:
            z_prom_minimum = exp_props.pk.z_prom_params{1};
            z_prom_maximum = exp_props.pk.z_prom_params{2};
            
            % Now plot z-peak_window_for_maxscores according to functions, below:
            filt_string_ms = strcat(ms_string, ':Z=', num2str(z_prom_minimum), 'to', num2str(z_prom_maximum), '-');
            colors = {'red', 'magenta', 'blue', 'green'};
            % Choose between the following behav_type(s): 1. 'bhv_at_pk', 2. 'majority_bhv_around_pk', 3. 'initiated_bhv_close_to_pk'
            for i_bhv_scr = 2:2

                % Choose i_bhv_type:
                bhv_scr_str = exp_props.bhv.bhv_scr_type;
                bhv_key = exp_props.bhv.bhv_key;
                bhv_colormap = exp_props.bhv.colormap;

                for i_type = 4
                    if i_type == 1
                        behav_type = strcat(bhv_scr_str, '_bhv_at_pk');
                    elseif i_type == 2
                        behav_type = strcat(bhv_scr_str, '_majority_bhv_around_pk');
                    elseif i_type == 3
                        behav_type = strcat(bhv_scr_str, '_initiated_bhv_close_to_pk');
                    elseif i_type == 4
                        behav_type = strcat(bhv_scr_str, '_max_bhv_close_to_pk');
                    end
    
                    for i_site = 1
                        % Define target site and both temporal alignment
                        % for graph and z-score alignment:
                        target_site = i_site;
                        time_align_string = strcat(exp_props.site_names{target_site}, 'Pk');
                        z_score_denom_string = exp_props.site_names{target_site};
    
                        % Now get the target-site peak_table name to load to send to plotting program:
                        event_table_name = strcat('peak_table_', exp_props.site_names{i_site}, '_', exp_props.filt_param_str{1});
                        event_table_name(event_table_name == '-') = '_';
                        event_table_all_trial_cell = vertcat(data_struct_per_ms(:).(event_table_name));
                        event_table_all_trial_tbl = vertcat(event_table_all_trial_cell{:});
                        if isempty(event_table_all_trial_tbl)
                            a = 5;
                        end
                                
                        % Make site specific directory:
                        site_specific_dir = strcat(ms_specific_dir, '_', exp_props.event_filt_names{exp_props.event_filt_cond},...
                            '_', exp_props.site_names{i_site}, '\');
                        mkdir(site_specific_dir);
                        curr_output_dir = site_specific_dir;
 
                        % Get filt_min_max
                        peak_prom_min_max = [z_prom_minimum, z_prom_maximum];
                        filt_min_max = peak_prom_min_max;
                
                        % B. Now, plot data of this struct:
                        F_2025_06_13_PlotMultiSiteOverlay_Paper(exp_props, curr_output_dir, all_trial_data_to_event, event_table_all_trial_tbl,...
                            colors, target_site, i_bhv_scr, behav_type, bhv_key, bhv_colormap, time_align_string, sec_lag);
                     end
                end
            end
        end
    end
