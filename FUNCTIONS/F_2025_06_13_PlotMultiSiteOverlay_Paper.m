function [] = F_2025_06_13_PlotMultiSiteOverlay_Paper(exp_props, curr_output_dir, trial_data_to_event, event_table_all_trial_tbl, colors, target_site, i_bhv_scr, behav_type, bhv_key, bhv_colormap, time_align_string, sec_lag)

 % This function will plot all z-scores at each peak according 

    % Get experimental properties variables:
    for i_param = 1:1
        z_prom_params = exp_props.pk.z_prom_params;
        z_proms_comp_site = exp_props.pk.z_proms_comp_site;
        num_sites = exp_props.num_sites;
        zscr_strings_plot = exp_props.site_names;
        min_prom_Z_other_site = z_proms_comp_site{1};
        num_bins_for_jacc = exp_props.pk.num_bins_for_jacc;
        zscr_inds_around_peak = exp_props.pk.SECS_PEAK_MOV * exp_props.pk.IMG_FRAME_GAP;
        filt_string = exp_props.event_filt_names(exp_props.event_filt_cond);
        PHOTOM_FR = exp_props.pho.PHOTOM_FR;
    end

    % Get variables from data_struct:
%     z_proms_table = data_struct_per_ms();
    z_scores_allmice = horzcat(trial_data_to_event.z_scores_allsites);
    z_controls_allmice = horzcat(trial_data_to_event.z_controls_allsites);
    signal_raw_allmice = horzcat(trial_data_to_event.signal_raw_allsites);
    isos_raw_allmice = horzcat(trial_data_to_event.isos_raw_allsites);

%     % Make Directory:
%     curr_output_dir = strcat(corr_fig_folder, '\', zscr_strings{target_site}, behav_type, '\');
%     mkdir(curr_output_dir);

    % 0.1 Plot correlation plots:

    means_heat_plot = 1;    % Section 3: Plot means (3A) and heat-maps (3B) - toggle 0/1
    ert_by_bhv_plot = 1;    % Section 4: Plot ERT's by behavior - toggle 0/1
 
    sig_control_graph_nums = 1;  % Section 3b: Decide types of graphs to plot - defined by #'s:
%     sig_control_graph_nums = [1, 2, 3, 4];
                                            % 1: z_scores_allmice
                                            % 2: z_controls_allmice
                                            % 3: signal_raw_allmice
                                            % 4: isos_raw_allmice

    % Define set of behaviors to plot: define all behaviors except
    % back-to-camera and jump
    if i_bhv_scr == 1
        iso_bhvs_to_plot = [1, 2, 3, 4, 5, 7, 8, 10];
    elseif i_bhv_scr == 2
        iso_bhvs_to_plot = 4;
    end

    % Define median z-prominence and y-lim:
    pk_prominences_tbl = event_table_all_trial_tbl(:, 'pk_prom');
    pk_prominences = pk_prominences_tbl{:, :};
    median_prom = median(pk_prominences);
    max_prom = max(pk_prominences);
    peak_size_filt = [z_prom_params{1}, z_prom_params{end}];
    num_filt_conds = size(peak_size_filt, 1);

    % Replace '_' in plotting strings with '-':
    filt_string_plot = filt_string;
    filt_string_plot = strrep(filt_string_plot, '_', '-');
    filt_string_plot = filt_string_plot{:};

    % 0.1 Parameters for making calculations:
    corr_reg_cross = 1;         % corr_reg_cross = 1 means --> regular correlation, and corr_reg_cross = 2 --> cross correlation
    corr_strings_all = {'Corr', 'XCorr'};
    corr_string = corr_strings_all{corr_reg_cross};

    % 0.2 Parameters for plotting the number of peaks:
    SEC_BASELINE_SLIT = 2;
    BASE_MARKER_SIZE = 12;
    MAX_SEC_DIFF = 90;
    MIN_NUM_PKS_BHV = 10;
    y_lim_nums = [-1, 5];
    num_peaks = size(event_table_all_trial_tbl, 1);

    % 0.3 Plotting parameters: 
    FONT_SIZE = 8;
    BHV_NM_COL = 4;
    strings_to_rep = {':',': ', ' ', '-', '='};

    % Section 1 default behavior parameters
    default_bhv_tag = 0;

    % Section 2 Heatmap plotting parameters:
    heat_sort_opt = 1;                      % Sort all heat-map trials by: 1 = prominence of target_site 2 = Sort by max-z-score of i_site, 3 = sort by index/time of max-z-score, 4 = max z-scores of target-site

    % Section 3 mean/heatmap plotting parameters:
    SEC_AROUND_PK = 4;
    MAX_LAG_SEC = 1;
    MAX_VAL_THRESH = 2.5;
    max_lag_inds = ceil(MAX_LAG_SEC * PHOTOM_FR);
    clims_for_heatmap = [-2, 8];
    num_trace_segs = 3;

    % What metric should be used to calculate the baseline? 'median' or 'average'?
    baseline_avg_med = 'avg';   % Should be 'med' or 'avg'.
    
    % Define an array to contain all z-score traces for all peaks:
    z_score_array = NaN(num_sites, (zscr_inds_around_peak*2 +1), num_peaks);
    z_controls_array = NaN(num_sites, (zscr_inds_around_peak*2 +1), num_peaks);

    % Define arrays to store numerical data that you can use for filtering
    % all peaks for plotting (all taken directly from peaks_table)

    % 1. Behavior tag of peak. 2. Z-scr prominence of peak 3. Time
    % isolation of peak:
    behav_tag_pk = zeros(num_peaks, 1);
    z_prom_of_pk = zeros(num_peaks, 1);
    sec_after_prev_peak = zeros(num_peaks, 1);

    %% Section 1: Get all z-score and behavior information in the form you need for subsequent plotting:
    for i_pk = 1:num_peaks

        % 1. Get a z-score across all sites
        i_list = event_table_all_trial_tbl{i_pk, 'i_list_num'};
        z_scores_to_plot = [];
        z_controls_to_plot = [];
        z_control_raw_sig_to_plot = [];
        z_control_raw_isos_to_plot = [];
        z_scores_to_plot(:, :) = z_scores_allmice{i_list};
        z_controls_to_plot(:, :) = z_controls_allmice{i_list};
        z_control_raw_sig_to_plot(:, :) = signal_raw_allmice{i_list};
        z_control_raw_isos_to_plot(:, :) = isos_raw_allmice{i_list};
        
        % 2. Get start, center, and ending indices:
        start_ind = event_table_all_trial_tbl{i_pk, 'ind_zscr_start'};
        stop_ind = event_table_all_trial_tbl{i_pk, 'ind_zscr_stop'};
        if start_ind > 0 && stop_ind < size(z_scores_to_plot, 2)
            z_scores_around_peak = z_scores_to_plot(:, start_ind:stop_ind);
            z_score_array(:, :, i_pk) = z_scores_around_peak;
            z_controls_around_peak = z_controls_to_plot(:, start_ind:stop_ind);
            z_controls_array(:, :, i_pk) = z_controls_around_peak;
            raw_sig_around_peak = z_control_raw_sig_to_plot(:, start_ind:stop_ind);
            raw_sig_array(:, :, i_pk) = raw_sig_around_peak;
            raw_isos_around_peak = z_control_raw_isos_to_plot(:, start_ind:stop_ind);
            raw_isos_array(:, :, i_pk) = raw_isos_around_peak;
            sec_after_prev_peak(i_pk) = event_table_all_trial_tbl{i_pk, 'peak_time_diff'};
        end

        % 3. Get behavior-score across all sites:
        tags_present_q = ismember (behav_type, event_table_all_trial_tbl.Properties.VariableNames);
        if tags_present_q == 1
            behav_tag_pk(i_pk) = event_table_all_trial_tbl{i_pk, behav_type};
        else
            behav_tag_pk(i_pk) = default_bhv_tag;
        end

        % 4. Get peak-prominence:
        z_prom_of_pk(i_pk) = event_table_all_trial_tbl{i_pk, 'pk_prom'};
    end

    % Get overall z_score min and max for subsequent plotting:
    max_z_scores_allsites = max(z_score_array, [], 'all');
    min_z_scores_allsites = min(z_score_array, [], 'all');

    %% Section 2: Now filter all of this raw information to get arrays of filtered/selected info!

    % Now make a cell array which contains the indices that the behavior
    % is positive for a given: A) Behavior B) Peak-size.
    bhvs_to_plot = (0:1:size(bhv_key, 1)-1)';
    bhv_prom_inds = cell(size(bhv_key, 1), num_filt_conds);
    bhv_prom_text = cell(size(bhv_key, 1), num_filt_conds);
    bhv_prom_color = cell(size(bhv_key, 1), num_filt_conds);
    bhv_prom_nums = cell(size(bhv_key, 1), num_filt_conds);

    for i_bhv_loop = 1:size(bhv_key, 1) + 1
        for i_prom = 1:num_filt_conds
            % Use a separate looping # to go over all behaviors in index since
            % a possible index for behavior is "0" - according to STEREO, and
            % you want to loop in indices:
            if i_bhv_loop <= size(bhv_key, 1)
                i_bhv = i_bhv_loop - 1;             % The i_bhv # encoded is '1 less' than the loop #.
            elseif i_bhv_loop > size(bhv_key, 1)
                i_bhv = 0:(size(bhv_key, 1) - 1);
            end
            
            % Now, store prominences and behavior values:
            % 1. Behavior condition met?
            bhv_check = (behav_tag_pk == i_bhv);
            bhv_check_max = max(bhv_check, [], 2);      % Need bhv_check_max to be a sum of all dimensions
            % 2. Peak-prominence condition met? Now for this behavior, group the peaks by their prominences:
                % Is peak prominence between lower-bound and upper-bound?
                % (If upper-bound not-set, upper bound = infinity):
                peak_prom_low_bound = peak_size_filt(i_prom, 1);
                peak_prom_upper_bound = peak_size_filt(i_prom, 2);
                prom_check_low = peak_prom_low_bound <= z_prom_of_pk;
                prom_check_high = peak_prom_upper_bound >= z_prom_of_pk;
                prom_check = prom_check_low.* prom_check_high;
            % 3. Now combine these two metrics to produce one set of logicals
            % with the indices where the peaks have met the both
            % conditions: 1) Behavior and 2) Peak prominence.
            bhv_prom_check = logical(prom_check).*bhv_check_max;
            % 4. Now, store logicals of indices indicating where these two
            % conditions (1. Behavior 2. Prominence) have been met;
            bhv_prom_inds{i_bhv_loop, i_prom} = logical(bhv_prom_check);
            bhv_prom_nums{i_bhv_loop, i_prom} = nnz(bhv_prom_check);
            % Now, store text strings and colors for plotting:
            if i_bhv_loop <= size(bhv_key, 1)
                bhv_descr = bhv_key{i_bhv_loop, 4};
                bhv_prom_color{i_bhv_loop, i_prom} = bhv_colormap(i_bhv_loop, :);
            elseif i_bhv_loop > size(bhv_key, 1)
                bhv_descr = 'All behaviors';
                bhv_prom_color{i_bhv_loop, i_prom} = [0, 0, 0];
            end
            text_string1 = strcat(bhv_descr, '. n=', num2str(bhv_prom_nums{i_bhv_loop, i_prom}));
            text_string2 = strcat('Z:', num2str(peak_prom_low_bound), 'to',...
                num2str(peak_prom_upper_bound));
            bhv_prom_text{i_bhv_loop, i_prom} = [text_string1 newline text_string2];
        end
    end

    % Secton 2b: shared graph parameters:
    xlims_on_ax = [1, length(z_controls_array(1,:,1))];

    %% Section 3: Calculate means along the entire matrix and plot:

    if means_heat_plot == 1

         for plot_mean_or_ind = 2
            
            if plot_mean_or_ind == 1
                % Create axes
                h_fig_overlay = figure();
                overlay_ax = axes(h_fig_overlay);
                for i_site = 1:num_sites
                    z_scores_one_site = z_score_array(i_site, :, :);
                    z_score_to_mean = permute(z_scores_one_site, [3, 2, 1]);
                    z_score_mean = nanmean(z_score_to_mean, 1);
                    h_plot = plot(overlay_ax, z_score_mean, colors{i_site});
                    hold on;
                    z_controls_one_site = z_controls_array(i_site, :, :);
                    z_controls_to_mean = permute(z_controls_one_site, [3, 2, 1]);
                    z_controls_mean = nanmean(z_controls_to_mean, 1);
                    h_plot = plot(overlay_ax, z_controls_mean, colors{i_site}, 'LineStyle', '--');
                end

                % Now Label and Save Figure:
                main_title_string = {strcat(filt_string_plot, ': Time0=', time_align_string), strcat(' ZscrsOf-AllSites')};
                save_title_string = strcat(main_title_string{1}, main_title_string{2}, '_', corr_string);
                for i_str = 1:length(strings_to_rep)
                    save_title_string = strrep(save_title_string, strings_to_rep{i_str}, '_');
                end
                save_title_string = strrep(save_title_string, '.', 'p');
                plot_type_string = 'Z-Scores-Overlay';
                xlabel_string = 'Time (sec)';
                ylabel_string = 'ZScore';
                xtick_create = 1;
                xlims_on_ax = [1, length(z_controls_mean)];
                fig_position = [50 101 400 850];
                F_2025_04_03_LabelSaveFigure(h_fig_overlay, curr_output_dir, main_title_string, save_title_string, plot_type_string, xlabel_string, ylabel_string, xtick_create, xlims_on_ax, zscr_inds_around_peak, PHOTOM_FR, [], []);

            elseif plot_mean_or_ind == 2
                for i_site = 1:num_sites

                    %  Now - make a figure and plot on heatmap:
                    h_fig_heat = figure();
                    h_ax_heat = axes(h_fig_heat);
                    
                    % Set variables and plot on heat-map:
                    z_score_mat_all_sites = z_score_array;

                    % Get global correlation to print:
                    corr_global_col_name = strcat(zscr_strings_plot{i_site}, 'global_corrs');
                    corrs_global = event_table_all_trial_tbl.(corr_global_col_name);
                    corr_global_mean = mean(corrs_global);

                    % Calculate correlation values for both standard and cross-correlation
                    [ind_threshold, parc_indices, proms_sorted_to_corr, z_scores_i_sorted_to_heat, z_scores_i_sorted_to_corr, z_scores_target_sorted_to_heat, z_scores_target_sorted_to_corr, corr_standard_all_trials, corr_standard_mean, corr_cross_all_trials, corr_cross_mean, corr_cross_lag_all_trials, corr_cross_lag_mean] =...
                        F_2024_12_22_Create_Corr_Dataset(pk_prominences, z_score_mat_all_sites, target_site, i_site, heat_sort_opt, num_trace_segs, PHOTOM_FR, min_prom_Z_other_site, sec_lag);
                    % Now use correlation values - and plot
                    [title_text] = F_2024_12_18_Create_Heat_Map(h_ax_heat, z_score_mat_all_sites, target_site, i_site, zscr_strings_plot, corr_reg_cross, heat_sort_opt, clims_for_heatmap, PHOTOM_FR, SEC_AROUND_PK, min_prom_Z_other_site,...
                        num_trace_segs, ind_threshold, parc_indices, z_scores_i_sorted_to_heat, z_scores_target_sorted_to_heat, corr_standard_all_trials, corr_standard_mean, corr_cross_all_trials, corr_cross_mean, corr_cross_lag_all_trials, corr_cross_lag_mean, corr_global_mean);
                    % Define threshold text for plotting/saving:
                    threshold_text = {strcat('Red-Line  Threshold = ', num2str(min_prom_Z_other_site), ' Zscr', zscr_strings_plot{i_site})};

                    % Now Label and Save Figure:
                    main_title_string = {strcat(filt_string_plot, ': Time0=', time_align_string), zscr_strings_plot{i_site}};
                    save_title_string = strcat(main_title_string{1}, main_title_string{2}, '_', corr_string, num2str(heat_sort_opt));
                    for i_str = 1:length(strings_to_rep)
                        save_title_string = strrep(save_title_string, strings_to_rep{i_str}, '_');
                    end
                    save_title_string = strrep(save_title_string, '.', 'p');
                    plot_type_string = 'HeatMap';
                    xlabel_string = 'Time (sec)';
                    ylabel_string = 'Trial';
%                     clr_bar_leg = 'Z-Score';
                    xtick_create = 1;
                    xlims_on_ax = [1, size(z_score_mat_all_sites, 2)];
                    main_title_string_w_thr = horzcat(main_title_string, threshold_text, title_text);
                    clr_bar_leg = 'ZScores';
                    fig_position = [50 101 400 850];
                    F_2025_04_03_LabelSaveFigure(h_fig_heat, curr_output_dir, main_title_string_w_thr, save_title_string, plot_type_string, xlabel_string, ylabel_string, xtick_create, xlims_on_ax, zscr_inds_around_peak, PHOTOM_FR, fig_position, clr_bar_leg);
                end
            end
        
        end
    end


    %% Section 4: Now, plot z-scores from all 4 sites according to filtered peaks.
    
    if ert_by_bhv_plot == 1

        a = 5;

        % Create two graphs
        for plot_mean_or_ind = 1:2

            if plot_mean_or_ind == 1
            
            % Define graph_type:
            graph_type = 'ZScoresOverlay';
            
            % Peaks are filtered by: 1) Behavior. 2) Z-Score prominence.
            h_fig_ert_by_bhv = figure();
            curr_fig = h_fig_ert_by_bhv;
        
            % Define dimensions of subplots:
            sub_rows = size(bhv_prom_inds, 2);
            sub_cols = length(iso_bhvs_to_plot);
            
            for i_bhv_loop = 1:length(iso_bhvs_to_plot)
                if i_bhv_loop <= length(iso_bhvs_to_plot)
                    i_bhv = iso_bhvs_to_plot(i_bhv_loop);

                    for i_prom = 1:size(bhv_prom_inds, 2)
            
                        % Access the right subplot:
                        overlay_ax = axes();
                        subplot(sub_rows, sub_cols, ((i_prom-1)*sub_cols + i_bhv_loop), overlay_ax);
        
                        % Get indexing array:
                        index_array = bhv_prom_inds{i_bhv, i_prom};
            
                        for i_hide_plot = 1:1
                            for i_site = 1:num_sites
                                z_scores_one_site = z_score_array(i_site, :, index_array);
                                z_score_to_mean = permute(z_scores_one_site, [3, 2, 1]);
                                z_score_mean = nanmean(z_score_to_mean, 1);
                                    % Need to normalize each z_scores baseline by subtracting the first
                                    % two seconds of the z_scores values:
                                    z_scores_bhv_avg = z_score_mean;
                                    baseline_avg = nanmean(z_score_mean(1:PHOTOM_FR * SEC_BASELINE_SLIT));
                                    z_scores_bhv_avg_norm = z_scores_bhv_avg - baseline_avg;
                                    z_scores_bhv_avg_norm_allsites(i_site, :) = z_scores_bhv_avg_norm;
                                h_plot = plot(overlay_ax, z_scores_bhv_avg_norm, colors{i_site}, 'LineWidth', 3);
                                hold on;
                            end
                
                            % Create time-axes around peak:
                            x_tick_vals = (-zscr_inds_around_peak:PHOTOM_FR:zscr_inds_around_peak);
                            x_tick_secs = x_tick_vals/PHOTOM_FR;
                            x_tick_str = cell(size(x_tick_secs, 2), 1);
                            for i_pt = 1:size(x_tick_secs, 2)
                                x_tick_str{i_pt} = num2str(x_tick_secs(i_pt));
                            end
                            x_tick_inds = x_tick_vals + zscr_inds_around_peak;
                            xlim([x_tick_inds(1), x_tick_inds(end)]);
                            xticks(x_tick_inds);
                            xticklabels(x_tick_str);
                            xlabel('Time - Seconds');    
                            ylabel('ZScores');
                            ylim(y_lim_nums);
                            title(bhv_prom_text{i_bhv, i_prom});
                        end
                    end
                end
            end

            % Add Title to whole figure:
            title_string = strcat(filt_string_plot, zscr_strings_plot{target_site}, ':');
            sgtitle(strcat(title_string, newline, graph_type));
            curr_fig.Position = [50 101 400 850];
        
            % Save figure:
            main_title_string = {strcat(filt_string{:}, ': Time0=', time_align_string, ': ZscrsOf', zscr_strings_plot{i_site})};
            save_name_root = char(main_title_string);
            save_name_root(save_name_root == ':') = '_';
            save_name_root(save_name_root == '=') = '';
            save_name_root(save_name_root == '-') = '_';
            save_name_root(save_name_root == ' ') = '_';
            save_name_root(save_name_root == '.') = 'p';
            save_name = strcat(curr_output_dir, save_name_root);
            saveas(curr_fig, strcat(save_name, 'ByBehavior'));
            saveas(curr_fig, strcat(save_name, 'ByBehavior.png'));

        elseif plot_mean_or_ind == 2

            % Define descriptors for plots:
            graph_types = {'HeatMap', 'Control-IsosProc', 'Control-Raw Signal', 'Control-Raw Isos'};

            % To see controls, change sig_control_graph_nums
            for i_sig_control = sig_control_graph_nums

                for i_site = 1:num_sites
                    % Define graph_type:
                    graph_type = graph_types{i_sig_control};
        
                    % Peaks are filtered by: 1) Behavior. 2) Z-Score prominence.
                    h_fig_heat_by_bhv = figure();
                    curr_fig = h_fig_heat_by_bhv;
                
                    % Define dimensions of subplots:
                    sub_rows = size(bhv_prom_inds, 2);
                    sub_cols = length(iso_bhvs_to_plot);
                    h_fig_tiles = tiledlayout(sub_rows, sub_cols, 'TileSpacing','compact');
                    
                    for i_bhv_loop = 1:length(iso_bhvs_to_plot)
                        if i_bhv_loop <= length(iso_bhvs_to_plot)
                            i_bhv = iso_bhvs_to_plot(i_bhv_loop);
                            for i_prom = 1:size(bhv_prom_inds, 2)
                    
                                % Access the right subplot:
    %                             overlay_ax = axes();
                                h_ax_heat = nexttile(((i_prom-1)*sub_cols + i_bhv_loop));
                
                                % Get indexing array and filter z_score array:
                                index_array = bhv_prom_inds{i_bhv, i_prom};
                                pk_prominences_target_site = pk_prominences(index_array);
                                if i_sig_control == 1
                                    z_score_mat_all_sites = z_score_array(:, :, index_array);
                                elseif  i_sig_control == 2
                                    controls_one_site = z_controls_array(i_site, :, index_array);
                                elseif i_sig_control == 3
                                    controls_one_site = raw_sig_array(i_site, :, index_array);
                                elseif i_sig_control == 4
                                    controls_one_site = raw_isos_array(i_site, :, index_array);
                                end

                                % Run function to create Heat_Map:
                                % Calculate correlation values for both standard and cross-correlation
                                heat_sort_opt = 4;                         % Sort all heat-map trials by: 1 = prominence of target_site 2 = Sort by max-z-score of i_site, 3 = sort by index/time of max-z-score                
                                [ind_threshold, parc_indices, proms_sorted_to_corr, z_scores_i_sorted_to_heat, z_scores_i_sorted_to_corr, z_scores_target_sorted_to_heat, z_scores_target_sorted_to_corr, corr_standard_all_trials, corr_standard_mean, corr_cross_all_trials, corr_cross_mean, corr_cross_lag_all_trials, corr_cross_lag_mean] =...
                                    F_2024_04_07_Create_Corr_Dataset(pk_prominences_target_site, z_score_mat_all_sites, target_site, i_site, heat_sort_opt, num_trace_segs, PHOTOM_FR, min_prom_Z_other_site, sec_lag);
                                % Now use correlation values - and plot
                                [title_text] = F_2023_05_10_Create_Heat_Map(h_ax_heat, z_score_mat_all_sites, target_site, i_site, zscr_strings_plot, corr_reg_cross, heat_sort_opt, clims_for_heatmap, PHOTOM_FR, SEC_AROUND_PK, min_prom_Z_other_site,...
                                    num_trace_segs, ind_threshold, parc_indices, z_scores_i_sorted_to_heat, z_scores_target_sorted_to_heat, corr_standard_all_trials, corr_standard_mean, corr_cross_all_trials, corr_cross_mean, corr_cross_lag_all_trials, corr_cross_lag_mean);
                                % Define threshold text for plotting/saving:
                                threshold_text = {strcat('Red-Line Threshold = ', num2str(min_prom_Z_other_site), ' Zscr', zscr_strings_plot{i_site})};

                                % Now Label Axes:
                                main_title_string = {bhv_prom_text{i_bhv, i_prom}, title_text};
                                xlabel_string = 'Time (sec)';
                                ylabel_string = 'Trial';
                                clr_bar_leg = [];
                                x_tick_inds = [];
                                F_2023_01_08_LabelAxes(overlay_ax, main_title_string, xlabel_string, ylabel_string, zscr_inds_around_peak, PHOTOM_FR, x_tick_inds, clr_bar_leg);
                            end
                        end
                    end
            
                    % Add Title to whole figure:
                    clr_bar_leg = 'ZScores';
                    fig_position = [50 101 400 850];
                    plot_type_string = graph_type;
                    main_title_string = {strcat(filt_string{:}, ': Time0=', time_align_string, ': ZscrsOf', zscr_strings_plot{i_site}),...
                        strcat('Correlation-Of-', zscr_strings_plot{i_site}, '-Zscr-To-', zscr_strings_plot{target_site}, 'DefinedZscr'), newline};
                    save_title_string = strcat(filt_string_plot, zscr_strings_plot{target_site}, 'FiltPeaks-vs-', zscr_strings_plot{i_site}, graph_type, '_', corr_string);
                    save_title_string = strrep(save_title_string,'.','p');
                    save_title_string = strrep(save_title_string,':','');
                    save_title_string = strrep(save_title_string,'=','eq');
                    xtick_create = 0;
                    main_title_string_w_thr = horzcat(main_title_string, threshold_text);
                    [h_fig_tiles_2] = F_2025_04_03_LabelSaveFigure(curr_fig, curr_output_dir, main_title_string_w_thr, save_title_string, plot_type_string, xlabel_string, ylabel_string, xtick_create, xlims_on_ax, zscr_inds_around_peak, PHOTOM_FR, fig_position, clr_bar_leg);
                    close('all');
                
                    % Figure out the number of each behavior's tagged peaks:
                    for i_bhv = 0:max(behav_tag_pk)
                        sum_pks_per_bhv((i_bhv + 1), 1) = i_bhv;
                        sum_pks_per_bhv((i_bhv + 1), 2) = nnz(behav_tag_pk == i_bhv);
                    end
        
                    end
                end
            end
        end
    end

end

