function [single_trial_data, trial_corr_struct] = F_2025_06_13b_Extract_Multisite_ZScores_Paper(single_trial_folder, single_trial_data, exp_props, filt_opt, raw_pks_plot, pks_or_sigisos, plot_dff_zscr, site_col_names, sites_to_plot, corr_baseline_site, corr_cross, color, fig_output_folder, trial_corr_struct, corr_field_names, plt_timeline_struct, tiled_overlay, i_list)

%% 0. Properties Definition section
    % 0.1 Properties for extracting information about raw signal from each
    % recorded site in neurophotometrics raw data file:
    COL_TIMESTAMP = 2;      % Time-stamp column in file
    LED_STATE_SIG = 2;      % LED-State signal channel
    LED_STATE_ISOS = 1;     % LED-State isos channel
    SEC_TO_SUBTRACT = 30;       % Sec to subtract from the beginning of recording

    % 0.2 Get/extract properties from exp_prop about the names of recording
    % sites and peak-definition:
    site_names = exp_props.site_names;
    num_sites = exp_props.num_sites;
    min_prom_pk = exp_props.pk.z_prom_params{1};
    % Define full_basename_string:
    full_basename_string = single_trial_data.trial_name;

    % 0.3 Pre-processing parameters
    POLYFIT_N_VAL = 1;

    % 0.4 Plotting properties
    % tiled_overlay = 'overlay';          % Should be either 'tile' or 'overlay'
    line_widths = [0.2, 2, 0.2, 1];         % If overlay; if 
    line_colors = {[.5 .5 .5]; 'black'; [.75 .75 .75]; 'magenta'};
    num_tm_pts = plt_timeline_struct.num_x_ticks;
    plt_minute_start = plt_timeline_struct.plt_minute_start;
    plt_minute_len = plt_timeline_struct.plt_minute_len;
    plt_signal_or_pks = 'signal';               % Should be a string. If == 'pks', will plot using 'findpeaks'. If any other value (default) plot signal.
    num_tiles = length(sites_to_plot);          % Get num_tiles from sites_to_plot 

%% 1. Importing raw data from neurophotometrics photometry data files
    % 1A. Enter zscore directory, and then get and define paths to interleaved-channel raw-photometry data:
        zscr_excel_folder = single_trial_folder;
        cd(zscr_excel_folder);
        interleaved_file_dir = dir('*Interleaved*');                % Get Interleaved file directory
    
    % 1B. If one or more interleaved file was found in the folder, read the data from the longest recording file in the folder:
    if length(interleaved_file_dir) >= 1
        file_sizes = interleaved_file_dir(:).bytes;
         [~, ind_largest_file] = max(file_sizes);
        filename_to_read = interleaved_file_dir(ind_largest_file).name;
        interleaved_matrix = readmatrix(filename_to_read);
        opts = detectImportOptions(filename_to_read);
        var_names = vertcat({opts.VariableOptions.Name})';
        cols_to_process = zeros(length(site_col_names), 1);
        for i_col_name = 1:length(site_col_names)
            str_to_find = site_col_names{i_col_name};
            cols_to_process(i_col_name) = find(strcmp(var_names, str_to_find));
        end
        signal_raw_only_allcols = interleaved_matrix(interleaved_matrix(:, 3) == LED_STATE_SIG, :);    % Get matrix of all 470_only values
        isos_raw_only_allcols = interleaved_matrix(interleaved_matrix(:, 3) == LED_STATE_ISOS, :);     % Get matrix of all 415_only values
    else
        msgbox('There is no interleaved file for this trial', 'Interleaved-file missing');
    end

    % 1C. Get photometry frame rate and other properties from interleaved file:
        % Load the Excel file
        data = readtable(filename_to_read);
        % Extract the Timestamp column
        timestamps = data.Timestamp;
        % Extract the LED column:
        ledstates = data.LedState;
        % The number of channels imaged will be equal to unique_ledstates - 1,
        % because ledstates always starts with '0', and then alternates between channels 1, 2, and 4.
        unique_ledstates = unique(ledstates);
        channels = unique_ledstates(unique_ledstates > 0);

    % 1D. Compute the difference between the last and first timestamp
        timestamp_difference = timestamps(end) - timestamps(1);
        approx_frame_rate = length(timestamps)/timestamp_difference;
        PHOTOM_FR = round(approx_frame_rate)./length(channels);
        % Calculate frame rate stats:
        frames_to_subtract = PHOTOM_FR * SEC_TO_SUBTRACT;

%% Section 2: Data pre-processing for each site:

    % 2A. Define arrays to store z-score information for one trial
    for i = 1:1
        z_scr_length_min = 0;
        z_scores_allsites = [];
        dff_allsites = [];
        z_timelines_allsites = [];
        z_controls_allsites = [];
        signal_raw_allsites = [];
        isos_raw_allsites = [];
        z_scr_lengths = [];
    end
 
    % 2B. Loop over all sites and import and pre-process signal for the photometry recording at that site::
    for i_site = 1:num_sites

        % i. Now, read the raw signal and isosbestic channel data at all recording sites:
            % Isolate the raw signal at recording_site_i:
            signal_raw_site = signal_raw_only_allcols(:, cols_to_process(i_site))';
            % Isolate the isosbestic at recording_site_i: 
            isos_raw_site = isos_raw_only_allcols(:, cols_to_process(i_site))';
            % Get the recording timeline
            z_timeline = signal_raw_only_allcols(:, COL_TIMESTAMP)' - signal_raw_only_allcols(1, COL_TIMESTAMP);
            % How many indices is the recording trace
            trace_length = size(z_timeline, 2);
            % Frame rate calculation:
            photom_FR = trace_length/( z_timeline(end) - z_timeline(1));            % This is the actual calculate frame rate for photometry!
            targeted_FR = (round(photom_FR * 100)/100);                             % This is the targeted integer frame rate (e.g. 10, 15, 30).    

        
        % ii. Smooth signal and process and plot:
            signal_cut = signal_raw_site(frames_to_subtract:end);
            isos_cut = isos_raw_site(frames_to_subtract:end);

        if filt_opt == 1
            % From PMAT - process data:
            span = 10;                   % Span = 5 indices, which at PHOTOM_FR of 15 is 1/3 of a second, which is ~timecourse of GCaMP6 signal change
                                         % Span = 15 indices, which at PHOTOMR_FR is 1 second!
            Fsignal=smooth(signal_cut, span, 'lowess'); 
            Fisos=smooth(isos_cut, span, 'lowess');
            if length(Fsignal) > length(Fisos)
                Fsignal = Fsignal(1:end-1);
            end
            bls=polyfit(Fisos, Fsignal, POLYFIT_N_VAL);
            % "p = polyfit(x,y,n) returns the coefficients for a polynomial p(x)...
            % of degree n that is a best fit (in a least-squares sense) for the data in y...
            % The coefficients in p are in descending powers, and the length of p
            % is n+1"
            Y_Fit=bls(1).*Fisos+bls(2);

            % ii. Now, make and plot the z-scores of deltaF/F with flattened signal and control.
                deltaF = ((Fsignal - Y_Fit)./Y_Fit)';
                control_signal = Y_Fit';
                z_scores = normalize(deltaF);

        elseif filt_opt == 2
            % Curve-fit/GUPPY - process data:

            % i. Smooth the data using filtfilt:
            filter_window = 10;
            a = 1;
            b = (ones(filter_window, 1)/filter_window);
            Fsignal=filtfilt(b, a, signal_cut); 
            Fisos=filtfilt(b, a, isos_cut); 
            if length(Fsignal) > length(Fisos)
                Fsignal = Fsignal(1:end-1);
            elseif length(Fisos) > length(Fsignal)
                Fisos = Fisos(1:end-1);
            end

            % ii. Now, fit the control to the signal: 
            bls=polyfit(Fisos, Fsignal, POLYFIT_N_VAL);
            % "p = polyfit(x,y,n) returns the coefficients for a polynomial p(x)...
            % of degree n that is a best fit (in a least-squares sense) for the data in y...
            % The coefficients in p are in descending powers, and the length of p
            % is n+1"
            Y_Fit=bls(1).*Fisos+bls(2);
            deltaF = ((Fsignal - Y_Fit)./Y_Fit);
            control_signal = Y_Fit;
            z_scores = normalize(deltaF);

        elseif filt_opt == 3
            % B. Curve-fit WITHOUT low-pass filtering!
            Fsignal=signal_cut; 
            Fisos=isos_cut; 
            if length(Fsignal) > length(Fisos)
                Fsignal = Fsignal(1:end-1);
            elseif length(Fisos) > length(Fsignal)
                Fisos = Fisos(1:end-1);
            end

            % ii. Now, fit the control to the signal: 
            bls=polyfit(Fisos, Fsignal, POLYFIT_N_VAL);
            % "p = polyfit(x,y,n) returns the coefficients for a polynomial p(x)...
            % of degree n that is a best fit (in a least-squares sense) for the data in y...
            % The coefficients in p are in descending powers, and the length of p
            % is n+1"
            Y_Fit=bls(1).*Fisos+bls(2);
            deltaF = ((Fsignal - Y_Fit)./Y_Fit);
            control_signal = Y_Fit;
            z_scores = normalize(deltaF);

        elseif filt_opt == 4 || filt_opt == 5
            % Multi-step filtering based on Simpson et al (adapted
            % photom_pp_Nov24 to Matlab code:)
           
            % i. Smooth the data using filtfilt:
            filter_window = 10;
            a = 1;
            b = (ones(filter_window, 1)/filter_window);
            % [b, a] = butter(2, 10/(sr/2), 'low');
            Fsignal=filtfilt(b, a, signal_cut); 
            Fisos=filtfilt(b, a, isos_cut); 
            if length(Fsignal) > length(Fisos)
                Fsignal = Fsignal(1:end-1);
            elseif length(Fisos) > length(Fsignal)
                Fisos = Fisos(1:end-1);
            end
            
            % ii. Fitting a curve for the photobleaching decay
            if filt_opt == 4
                bleaching_cutoff_freq = 0.01 / 15;
            elseif filt_opt == 5
                bleaching_cutoff_freq = 0.002 / 15;
            end
            [b, a] = butter(2, bleaching_cutoff_freq, 'low');
            sig_decfit = filtfilt(b, a, Fsignal);
            iso_decfit = filtfilt(b, a, Fisos);

            % iii. Removing the decay
            Fsignal_noBleach = Fsignal - sig_decfit;
            Fisos_noBleach = Fisos - iso_decfit;

            % iv-a. Calculation of motion artifact by regression of
            % isosbestic channel against signal channel:
            % Create a table with the data
            tbl = table(Fisos_noBleach(:), Fsignal_noBleach(:), 'VariableNames', {'X', 'Y'});
            % Perform linear regression
            mdl = fitlm(tbl, 'Y ~ X');
            % Extract the coefficients and statistics
            slope = mdl.Coefficients.Estimate(2);
            intercept = mdl.Coefficients.Estimate(1);
            r_value = mdl.Rsquared.Ordinary;
            p_value = mdl.Coefficients.pValue(2);
            std_err = mdl.Coefficients.SE(2);

            % iv-b. Now calculate a motion-corrected Fsignal
            sig_est_motion = intercept + slope*Fisos_noBleach;
            Fsignal_motion_corrected = Fsignal_noBleach-sig_est_motion;

            % v-other-scores-to-store:
            control_signal = sig_decfit;

            % v. Calculating df/f and z-score
            deltaF = Fsignal_motion_corrected./sig_decfit;
            z_scores = normalize(deltaF);

        end
    
        % Now, calculate timing of photometry start-recording, duration, and stop recording:
            photom_zscr_sec_delay = SEC_TO_SUBTRACT;         % This will be the time-value of the first index of the processed, Z-score file!
            photom_zscr_sec_end = ceil(z_timeline(end));        % This is the duration of the whole recording according to the photometry file
            photom_zscr_sec_duration = photom_zscr_sec_end - photom_zscr_sec_delay;            
            num_inds_start_missing = round(photom_zscr_sec_delay * PHOTOM_FR);

        % Make sure all traces are the same length:
        for i = 1:1
            z_scr_lengths = [size(z_scores, 2), size(z_timeline, 2), size(control_signal, 2)];
            z_scr_length_min = min(z_scr_lengths);
        end

        % Set values of all z_scores and traces for each site:
        dff_allsites(i_site, :) =  deltaF(1, 1:z_scr_length_min);
        z_scores_allsites(i_site, :) =  z_scores(1, 1:z_scr_length_min);
        z_timelines_allsites(i_site, :) = z_timeline(1, 1:z_scr_length_min);
        z_controls_allsites(i_site, :) = control_signal(1, 1:z_scr_length_min);
        
        % This is the default length of a z-score, which will be used by
        % subsequent read-pkl-bhv function to create a trace. This should
        % be the mininmum length of the z_scores plus the number of inds
        % missing from start:
        len_z_scr_trace = round(min(z_scr_length_min) + num_inds_start_missing);
        len_photom_raw = round(length(signal_raw_only_allcols));

        % Here, there is not the missing part of the signal at the
        % beginning of the movie, so subtract from the beginning so that
        % all z-score filtered and raw traces are aligned in time:
        signal_raw_allsites(i_site, :) = signal_raw_site(1, num_inds_start_missing:(z_scr_length_min+num_inds_start_missing-1));
        isos_raw_allsites(i_site, :) = isos_raw_site(1, num_inds_start_missing:(z_scr_length_min+num_inds_start_missing-1));
        
    end

    %% 3. Get correlation values:
    
    % Loop over different 'orders' of sites, with each site taking turns as
    % the main-site: Str_L, Str_R, Axon_L, Axon_R:
    trial_corr_struct(i_list).trialname = single_trial_data.trial_name(1:end-1);
    trial_corr_struct(i_list).mouse = single_trial_data.mouse(1:end-1);
    trial_corr_struct(i_list).cond = single_trial_data.cond(1:end-1);
    trial_corr_struct(i_list).day = single_trial_data.date(1:end-1);

    for i_site_perm = 1:num_sites
        if i_site_perm == 1
            site_order = [1, 2, 3, 4];
        elseif i_site_perm == 2
            site_order = [2, 1, 4, 3];
        elseif i_site_perm == 3
            site_order = [3, 4, 1, 2];
        elseif i_site_perm == 4
            site_order = [4, 3, 2, 1];
        end
        corr_val_name_string = corr_field_names{i_site_perm, 1};
        corr_vals_array = zeros(num_sites, 1);
        for i_site = 1:num_sites
            baseline_site = site_order(1);
            if strcmpi(corr_cross, 'yes')
                [corr_vals, lags] = xcorr(z_scores_allsites(baseline_site, :), z_scores_allsites(site_order(i_site), :), exp_props.pho.PHOTOM_FR, 'unbiased');
                [corr_max, corr_max_ind] = max(corr_vals);
                corr_max_lag = lags(corr_max_ind)/exp_props.pho.PHOTOM_FR;
                corr_add_txt = strcat('Xcorr=', num2str(corr_max), 'Xcorr-Time-Shift-to-Target=', num2str(corr_max_lag));
            else
                corr_vals = corr(z_scores_allsites(corr_baseline_site, :)', z_scores_allsites(site_order(i_site), :)');
                corr_max = corr_vals;
                corr_max_lag = 0;
                corr_add_txt = strcat('Corr=', num2str(corr_max), 'Standard-Corr-Time-Shift=', num2str(corr_max_lag));
            end
            corr_vals_array(i_site) = corr_max;
        end
        trial_corr_struct(i_list).(corr_val_name_string) = corr_vals_array; 
    end


    %% 4. Plotting signals:

    if raw_pks_plot > 0
        % Now, enter the raw_data_plot figure directory and plot all z-scores
        % with find-peaks:
        mkdir(fig_output_folder);
        cd(fig_output_folder);
        raw_plot_string = {};

            % Figure creation depending on plotting options:
            h_fig_raw_plt = figure();
            h_fig_raw_plt.Position = [50 101 1600 850];
            if strcmp(tiled_overlay, 'overlay') ~= 1
%                 h_fig_raw_plt = tiledlayout("vertical", "TileSpacing", "compact");
                h_fig_raw_plt = tiledlayout(num_tiles,1,"TileSpacing","compact");
                h_fig_raw_plt.Title.String = full_basename_string(1:end-1);
            else
                title(full_basename_string(1:end-1));
            end

            % h_fig_raw_plt.Position = [50 101 1600 850];
            % h_ax_site = axes(h_fig_raw_plt);
%             h_fig_raw_plt.Title.String = full_basename_string(1:end-1);

            % Set y-lim parameters:
            y_lim_max = 1;
            y_lim_min = -.5;
            y_lim_diff = y_lim_max - y_lim_min;

            % Set x-lim parameters:
            timeline_trace = z_timelines_allsites(1, :);
            starting_time_ind = plt_minute_start * 60 * PHOTOM_FR + 1;
            length_time_inds = plt_minute_len * 60 * PHOTOM_FR;
            timeline_trace_trim_inds = starting_time_ind:(starting_time_ind+length_time_inds);
            if (length(timeline_trace) > timeline_trace_trim_inds(end))
                timeline_trace_trim_inds = starting_time_ind:(starting_time_ind+length_time_inds);
                % if length(timeline_trace_trim_inds) > length(timeline_trace)
                %     timeline_trace_trim_inds = timeline_trace_trim_inds(1):length(timeline_trace);
                % end
            else 
                timeline_trace_trim_inds = starting_time_ind:length(timeline_trace);
            end
            % if timeline_trace_trim_inds(end) > length(timeline_trace)
            %     timeline_trace_trim_inds = timeline_trace_trim_inds(1):length(timeline_trace);
            % end
            if timeline_trace(end) == 18185
                a = 5;
            end
            timeline_trace_trim = timeline_trace(1, timeline_trace_trim_inds);
            xticks_plot_inds = floor(1:length(timeline_trace_trim)/num_tm_pts:length(timeline_trace_trim));
            cell(num_tm_pts, 1);
            for i_tm = 1:num_tm_pts
                xticklabels_plot{i_tm} = round(timeline_trace_trim(xticks_plot_inds(i_tm))/(60),2);
            end
            xtickangle_plot = 30;

            for i_site = sites_to_plot

                if strcmp(tiled_overlay, 'overlay') == 1
                    h_ax_site = gca;
                    % Line-width and line-color:
                    line_width = line_widths(i_site);
                    line_color = line_colors{i_site};
                else
                    % Default condition is to 'tile' plot:
                    h_ax_tile(i_site) = nexttile;
                    h_ax_site = h_ax_tile(i_site);
                    line_width = 0.4;
                    line_color = 'blue';
                end

                % Get z-score for this site:
                % z_score_to_plot = z_scores_allsites(i_site, :)/max(z_scores_allsites(i_site, :));
                dff_to_plot = dff_allsites(i_site, timeline_trace_trim_inds);
                z_score_to_plot = z_scores_allsites(i_site, timeline_trace_trim_inds);
                signal_smooth_to_plot = signal_smooth_allsites(i_site, timeline_trace_trim_inds)/max(signal_smooth_allsites(i_site, timeline_trace_trim_inds));
                isos_smooth_to_plot = isos_smooth_allsites(i_site, timeline_trace_trim_inds)/max(isos_smooth_allsites(i_site, timeline_trace_trim_inds));
                z_timeline_to_plot = z_timelines_allsites(i_site, timeline_trace_trim_inds);
                corr_add_txt = '';

                if pks_or_sigisos == 1
                    hold on;
                    if strcmp(plt_signal_or_pks, 'plt_pks') == 1
                        findpeaks(z_score_to_plot, "MinPeakProminence", min_prom_pk);
                        raw_plot_string = 'Pks';
                    else
                        if strcmp(plot_dff_zscr, 'dff') == 1
                            plot(h_ax_site, dff_to_plot, 'Color', line_color, 'LineWidth', line_width);
                            ylabelstring = 'Dff';
                            % ylim_vals = [-0.02, .13];
                            ylim_vals = [-0.06, .13];
                        else
                            plot(h_ax_site, z_score_to_plot, 'Color', line_color, 'LineWidth', line_width);
                            ylabelstring = 'Zscore';
                            ylim_vals = [-2, 8];
                        end
                        raw_plot_string = 'Signal';
                    end
                    if strcmpi(corr_cross, 'yes')
                        [corr_vals, lags] = xcorr(z_scores_allsites(corr_baseline_site, :), z_scores_allsites(i_site, :), exp_props.pho.PHOTOM_FR, 'unbiased');
                        [corr_max, corr_max_ind] = max(corr_vals);
                        corr_max_lag = lags(corr_max_ind)/exp_props.pho.PHOTOM_FR;
                        corr_add_txt = strcat('Xcorr=', num2str(corr_max), 'Xcorr-Time-Shift-to-Target=', num2str(corr_max_lag));
                    else
                        corr_vals = corr(z_scores_allsites(corr_baseline_site, :)', z_scores_allsites(i_site, :)');
                        corr_max = corr_vals;
                        corr_max_lag = 0;
                        corr_add_txt = strcat('Corr=', num2str(corr_max), 'Standard-Corr-Time-Shift=', num2str(corr_max_lag));
                    end
                    hold off;
                elseif pks_or_sigisos > 2
                    hold on;
                    h_leg(1) = plot(z_score_to_plot, 'k', 'LineWidth', 3);
                    h_leg(2) = plot(signal_smooth_to_plot, 'b');
                    h_leg(3) = plot(isos_smooth_to_plot, 'r');
                    raw_plot_string = 'Raw_Isos';
                    hold off;
                end

%                 sublegend;
                % h_ax_site(i_site).YLim = [y_lim_min, y_lim_max];
                h_ax_site.YLim = ylim_vals;
                h_ax_site.XTick = xticks_plot_inds;
                h_ax_site.XTickLabel = xticklabels_plot;
                h_ax_site.XTickLabelRotation = xtickangle_plot;
                y_coors = [y_lim_min, y_lim_max];
                title(strcat(site_col_names{i_site},'-', site_names{i_site}), corr_add_txt);
                ylabel(ylabelstring);
            end

        % Create file_save_name
        file_save_name = strcat(full_basename_string(1:end-1), '_Tsite', num2str(corr_baseline_site), '_', color, raw_plot_string,...
        '_',tiled_overlay, '_',num2str(plt_minute_start),'_',num2str(plt_minute_len), '_', ylabelstring);

        % Save figures:
        saveas(h_fig_raw_plt, file_save_name, 'png');
        if raw_pks_plot > 1
            saveas(h_fig_raw_plt, file_save_name, 'fig');
        end
        if raw_pks_plot > 2
            % Ensure text is stored as text, not outlines
            set(groot, 'defaultAxesFontName', 'none');
            % Save as SVG using 'print' instead of 'saveas'
            set(findall(h_fig_raw_plt, '-property', 'FontName'), 'FontName', 'Arial');
            print(h_fig_raw_plt, file_save_name, '-dsvg', '-vector');
        end
        close('all');
    end

    % Store all z-score values for all sites in trial_list:
    single_trial_data.z_scores_allsites = z_scores_allsites;
    single_trial_data.z_timelines_allsites = z_timelines_allsites;
    single_trial_data.z_controls_allsites = z_controls_allsites;
    single_trial_data.signal_raw_allsites = signal_raw_allsites;
    single_trial_data.isos_raw_allsites = isos_raw_allsites;

    % Store all values for all 
    single_trial_data.props_xls.len_z_scr_trace = len_z_scr_trace;
    single_trial_data.props_xls.len_photom_raw = len_photom_raw;
    single_trial_data.props_xls.photom_zscr_sec_delay = photom_zscr_sec_delay;
    single_trial_data.props_xls.photom_zscr_sec_end = photom_zscr_sec_end;
    single_trial_data.props_xls.photom_zscr_sec_duration = photom_zscr_sec_duration;
    single_trial_data.props_xls.num_inds_start_missing = num_inds_start_missing;
    single_trial_data.props_xls.photom_FR = photom_FR;
    single_trial_data.props_xls.targeted_FR = targeted_FR;

end







