function [] = F_2025_06_04_Plot_Correlation_Values(fig_output_folder, all_filt_corrs, filt_opts, filt_opt_strings, exp_props, corr_vals_fields_to_proc, site_to_corr_comp, plot_title)
% F_2024_12_23_Plot_Correlation_Values will plot correlation values on a
% scatter-plot, according to mouse, condition, and trial as a bar graph,
% and then in a plot with all on the same column:

% 0. Properties:
bar_plot_y = 0;             % Do bar-plot: 0 = No. 1 = yes!
dot_plot_y = 1;             % Make a dot-plot. 0 = No. 1 = yes!
uni_dim_plot_y =1;          % Make unidimensional plot. 0 = No. 1 = yes!
corr_filt_cmpr_plot_y = 1;  % Make a plot comparing values across filtering conditions, like for unidimensional plot. 0 = No. 1 = yes!
marker_size = 104;
ylim_low_lim = 0;

% 0.1 Calculated properties:
sample_type_filt = filt_opt_strings{filt_opts(1)};
trial_corr_struct_sample = all_filt_corrs.(sample_type_filt);
num_trials = size(trial_corr_struct_sample, 2);
num_filt_conds = length(filt_opts);

% 0.2 Define arrays to store corr-data-pt-strings and corr-values:
corr_data_pt_names = cell(num_trials, num_filt_conds);
corr_data_pt_vals = zeros(num_trials, num_filt_conds);

% 0.3 Full trial-names:
trial_hemi_names = cell(num_trials, 1);

% 1. Loop through trials and get trial-name and corr-value
for i_filt = 1:length(filt_opts)
    filt_index = filt_opts(i_filt);
    type_of_filt = filt_opt_strings{filt_index};
    trial_corr_struct_curr_filt = all_filt_corrs.(type_of_filt);

    % Create trial loop-counter;
    loop_count = 1;

    for i_trial = 1:num_trials
        % Get trial:
        trial_name = trial_corr_struct_curr_filt(i_trial).trialname;
        % Get mouse:
        trial_corr_struct_curr_filt(i_trial).mouse;
        % Get hemisphere-name:
        site_hemi_descript = corr_vals_fields_to_proc{2, 1};
        % Get corr-values field in struct for this hemisphere:
        corr_field_name = corr_vals_fields_to_proc{1, 1};
        % Create trial-string
        corr_data_pt_string = strcat(trial_name, '-', site_hemi_descript);
        % Store the trial-ms-hemi descriptors for plotting:
        trial_hemi_names{loop_count, i_filt} = corr_data_pt_string;
        corr_data_pt_string('_') = '-';
        % Get correlation values:
        corr_values = trial_corr_struct_curr_filt(i_trial).(corr_field_name);
        corr_value_cmpr = [];
        corr_value_cmpr = corr_values(site_to_corr_comp);
        % Store data_pt_name:
        corr_data_pt_string = corr_data_pt_string(7:end);
        corr_data_pt_string = strrep(corr_data_pt_string, 'CGA262', '');
        corr_data_pt_string = strrep(corr_data_pt_string, 'CG638', '');
        corr_data_pt_string = strrep(corr_data_pt_string, 'CG650', '');
        corr_data_pt_string = strrep(corr_data_pt_string, 'CG780', '');
        corr_data_pt_string = strrep(corr_data_pt_string, 'CG784', '');
        corr_data_pt_string = strrep(corr_data_pt_string, 'CO', '');
        corr_data_pt_string = strrep(corr_data_pt_string, '_', '-');
        corr_data_pt_names{loop_count, i_filt} = corr_data_pt_string;
        % Store values:
        corr_data_pt_vals(loop_count, i_filt) = corr_value_cmpr;
        % Increase loop counter
        loop_count = loop_count + 1;
    end
end

% 2. Prepare for dot-plotting:
    % Now set up cell_arrays to contain characteristics for plotting:
    % Get condition array - based on trial-name from first filtering
    % condition:
    for i_nm = 1:size(trial_hemi_names, 1)
        [cg_ind] = strfind(trial_hemi_names{i_nm, 1}, 'CG');
        trial_ms_array{i_nm} = trial_hemi_names{i_nm, 1}(cg_ind:cg_ind+8);
        [co_ind] = strfind(trial_hemi_names{i_nm, 1}, 'CO');
        cond_array{i_nm} = trial_hemi_names{i_nm, 1}(co_ind:end);
    end

    % Now, sort the array by 'condition':
    [~, sortedIndices] = sort(cond_array);
    % Now, resort arrays;
    trial_hemi_names_sort = trial_hemi_names(sortedIndices, :);
    corr_data_pt_name_sort = corr_data_pt_names(sortedIndices, :);
    corr_data_pt_vals_sort = corr_data_pt_vals(sortedIndices, :);

    % Convert to a cell array of unique values
    uniqueMsArray = unique(trial_ms_array);

    % Now get shape from the shape array:
    shape_array = cell(size(trial_hemi_names, 1), num_filt_conds);
    markerShapes = {'o', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};   % Add more shapes if needed
    fill_array = cell(size(trial_hemi_names, 1), num_filt_conds);
    color_array = cell(size(trial_hemi_names, 1), num_filt_conds);

    for i_filt = 1:length(filt_opts)
        for i_nm = 1:size(trial_hemi_names_sort)
            for i_ms = 1:length(uniqueMsArray)  % Loop ovr all the mice in the cohort
                if strfind(trial_hemi_names_sort{i_nm}, uniqueMsArray{i_ms}) > 0   % This particular mouse was indeed found
                    shape_array{i_nm, i_filt} = markerShapes{i_ms};
                end
                if strfind(trial_hemi_names_sort{i_nm}, 'Right') > 0
                    fill_array{i_nm, i_filt} = 'fill';
                end
                if strfind(trial_hemi_names_sort{i_nm}, 'CGA262') > 0
                    color_array{i_nm, i_filt} = 'red';
                elseif strfind(trial_hemi_names_sort{i_nm}, 'CG780') > 0
                    color_array{i_nm, i_filt} = 'blue';
                elseif strfind(trial_hemi_names_sort{i_nm}, 'CG784') > 0
                    color_array{i_nm, i_filt} = 'blue';
                else
                    color_array{i_nm, i_filt} = 'green';
                end
            end
        end
    end

    % color_array = cell(size(trial_hemi_names, 1), 1);
    % shape_array = cell(size(trial_hemi_names, 1), 1);
    % fill_array = cell(size(trial_hemi_names, 1), 1);


    % Highlight dot-plots by making a unique, distinguishing characteristic
    % for each mouse, hemisphere, condition.

    % Color for dSPN/iSPN and reporter?

    % Mouse = shape

    % Define marker shapes for each subject

    % Hemisphere: Left = empty; Right = filled.
    % Condition = position on x-axis.


% % 3. Plot bar graph with x-ticks as trial/side = corr-data-pt-names on the
% % x-axis:
% if bar_plot_y == 1
%     h_corr_bar = figure();
%     h_corr_bar.Position = [50 101 1600 850];
%     h_ax = axes();
%     bar(corr_data_pt_vals_sort);
%     
%     % Setting the x-ticks and their labels
%     set(h_ax, 'XTick', 1:length(corr_data_pt_names), 'XTickLabel', corr_data_pt_names);
%     
%     % Rotating the x-tick labels
%     xtickangle(h_ax, 80);
%     
%     % Adding labels and title for clarity
%     xlabel(h_ax, 'Trial Name');
%     ylabel(h_ax, 'Striatal-Axonal Correlation');
%     
%     % Now, adding title;
%     title(plot_title);
% 
%     % Now save plot:
%     datetime_curr = datetime;
%     datetime_curr.Format = 'yyyy-MM-dd';
%     save_name = strcat(fig_output_folder, char(datetime_curr), plot_title, 'bar');
%     saveas(h_corr_bar, save_name, 'fig');
% end
% 
% if dot_plot_y == 1
%     % Make figure:
%     h_corr_dotplot = figure();
%     h_corr_dotplot.Position = [50 101 1600 850];
%     h_ax_dot = axes();
% 
%     % Make scatter-plot:
%     x_pts = 1:length(corr_data_pt_vals_sort);
%     for i_pt = 1:size(corr_data_pt_vals_sort, 1)
%         hold on;
%         if strcmp(fill_array{i_pt}, 'fill')
%             scatter(x_pts(i_pt), corr_data_pt_vals_sort(i_pt), marker_size, color_array{i_pt}, "filled", 'Marker', shape_array{i_pt});
%             % scatter(x_pts(i_pt), corr_data_pt_vals(i_pt), 24, 'black', 'Marker', '+');
%         else
%             scatter(x_pts(i_pt), corr_data_pt_vals_sort(i_pt), marker_size, color_array{i_pt}, 'Marker', shape_array{i_pt});
%         end
%         % scatter(x_pts(i_pt), corr_data_pt_vals(i_pt), 'Color', color_array{i_pt}, 'MarkerShape', shape_array{i_pt});
%     end
% 
%     % 
%     xlim([0.5, (x_pts(end) + 1)]);
%     ylim([ylim_low_lim, 1]);
% 
%     % Setting the x-ticks and their labels
%     set(h_ax_dot, 'XTick', 1:length(corr_data_pt_name_sort), 'XTickLabel', corr_data_pt_name_sort);
%     
%     % Rotating the x-tick labels
%     xtickangle(h_ax_dot, 80);
%     
%     % Adding labels and title for clarity
%     xlabel(h_ax_dot, 'Trial Name');
%     ylabel(h_ax_dot, 'Striatal-Axonal Correlation');
%     
%     % Now, adding title;
%     title(plot_title);
% 
%     % Now save plot:
%     datetime_curr = datetime;
%     datetime_curr.Format = 'yyyy-MM-dd';
%     save_name = strcat(fig_output_folder, char(datetime_curr), plot_title, 'scatt');
%     saveas(h_corr_dotplot, strcat(save_name, '.fig'), 'fig');
% end
% 
if uni_dim_plot_y == 1
    % Make figure:
    h_uni_dim_plot = figure();
    h_uni_dim_plot.Position = [50 101 400 650];
    h_axh_uni_dim_plot = axes();

    % Set x-axis:
    x_pts = ones(length(corr_data_pt_vals_sort), 1);

    % Calculate mean and S.E.M. line:
    % Calculate the mean
    meanData = mean(corr_data_pt_vals_sort);

    % Calculate the standard error of the mean (SEM)
    semData = std(corr_data_pt_vals_sort) / sqrt(length(corr_data_pt_vals_sort)); 

    % Set xlim and ylim:
    xlim_lower = 0.5;
    xlim_upper = 1.5;

    % % Plot mean
    % Plot shaded SEM area
    x_coors_to_fill = [xlim_lower, xlim_lower, xlim_upper, xlim_upper];
    y_coors_to_fill = [meanData + semData, fliplr(meanData - semData), fliplr(meanData - semData), meanData + semData];
    fill(x_coors_to_fill, y_coors_to_fill, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    % errorbar([xlim_lower, xlim_lower], [meanData, meanData], [semData, semData], 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
    hold on;
    ylim([0, 1]);
    plot([xlim_lower, xlim_upper], [meanData, meanData], 'black', 'LineWidth', 1);
    for i_pt = 1:length(x_pts)
        scatter(x_pts, corr_data_pt_vals_sort, color_array{i_pt}, 'd');
    end

    % Labels and title
    % xlabel('Group');
    ylabel('Trial Correlation');
    title('Mean ± SEM of Trial Correlations');
    grid on;
    hold off;

    % Now save plot:
    datetime_curr = datetime;
    datetime_curr.Format = 'yyyy-MM-dd';
    save_name = strcat(fig_output_folder, char(datetime_curr), plot_title, 'scatt_line_sem');
    saveas(h_uni_dim_plot, save_name, 'fig');
    % Ensure text is stored as text, not outlines
    set(groot, 'defaultAxesFontName', 'none');
    % Save as SVG using 'print' instead of 'saveas'
    set(findall(h_uni_dim_plot, '-property', 'FontName'), 'FontName', 'Arial');
    print(h_uni_dim_plot, save_name, '-dsvg', '-vector');

end

if corr_filt_cmpr_plot_y == 1

    % Make figure:
    h_uni_dim_plot = figure();
    h_uni_dim_plot.Position = [50 101 400 650];
    h_axh_uni_dim_plot = axes();

    for i_filt = 1:num_filt_conds

        % Get corr_data_pt_vals_sigfilt
        corr_data_pt_vals_sigfilt = corr_data_pt_vals_sort(:, i_filt);

        % Set x-axis:
        x_pts = ones(size(corr_data_pt_vals_sigfilt, 1), 1);
    
        % Calculate mean and S.E.M. line:
        % Calculate the mean
        meanData = mean(corr_data_pt_vals_sigfilt);
    
        % Calculate the standard error of the mean (SEM)
        semData = std(corr_data_pt_vals_sigfilt) / sqrt(length(corr_data_pt_vals_sigfilt)); 
    
        % Set xlim and ylim:
        xlim_lower = (i_filt-1) + 0.5;
        xlim_upper = (i_filt-1) + 1.5;
    
        % % Plot mean
        % Plot shaded SEM area
        x_coors_to_fill = [xlim_lower, xlim_lower, xlim_upper, xlim_upper];
        y_coors_to_fill = [meanData + semData, fliplr(meanData - semData), fliplr(meanData - semData), meanData + semData];
        fill(x_coors_to_fill, y_coors_to_fill, 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        % errorbar([xlim_lower, xlim_lower], [meanData, meanData], [semData, semData], 'o', 'MarkerSize', 8, 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
        hold on;
        ylim([0, 1]);
        plot([xlim_lower, xlim_upper], [meanData, meanData], 'black', 'LineWidth', 1);
%         for i_pt = 1:length(x_pts)
%             scatter(x_pts, corr_data_pt_vals_sort, color_array{i_pt}, 'd');
%         end
    end
    hold off;

    % Now plot line-plot across different signal filtering conditions:
    x_pts = 1:length(filt_opts);
    hold on;
    for i_trial = 1:num_trials
        y_pts = (corr_data_pt_vals_sort(i_trial, :))';
        plt_trial(i_trial) = plot(x_pts, y_pts, 'Color', color_array{1, 1});
%         scatter(x_pts, corr_data_pt_vals_sort, color_array{i_pt}, 'd');
    end
    % Now set legend for above plot:
    % legend(trial_hemi_names_sort);
    hold off;
    % legend();

    % Labels and title
    % xlabel('Group');
    xticks(x_pts)
    % xlabelstrings = filt_opt_strings(filt_opts);
    % xlabelstrings(xlabelstrings == '_') = '-';
    xticklabels(filt_opt_strings(filt_opts));
    xlabel('Filtering Condition')
    ylabel('Trial Correlation');
    title('Mean ± SEM of Trial Correlations');
    grid on;
    hold off;

    % Now save plot:
    datetime_curr = datetime;
    datetime_curr.Format = 'yyyy-MM-dd';
    % If directory doesn't exist, then create!
    if ~isfolder(fig_output_folder)
        % Create the folder if it doesn't exist
        mkdir(fig_output_folder);
        disp(['Directory "', fig_output_folder, '" created.']);
    else
        disp(['Directory "', fig_output_folder, '" already exists.']);
    end
    save_name = strcat(fig_output_folder, '\', char(datetime_curr), plot_title, 'scatt_line_sem');
    saveas(h_uni_dim_plot, save_name, 'fig');
    % Ensure text is stored as text, not outlines
    set(groot, 'defaultAxesFontName', 'none');
    % Save as SVG using 'print' instead of 'saveas'
    set(findall(h_uni_dim_plot, '-property', 'FontName'), 'FontName', 'Arial');
    print(h_uni_dim_plot, save_name, '-dsvg', '-vector');

end

end