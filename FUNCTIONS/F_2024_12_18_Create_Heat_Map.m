function [title_text] = F_2023_05_21_Create_Heat_Map(h_ax_heat, z_score_mat_all_sites, target_site, i_site, zscr_strings, corr_reg_cross, heat_sort_opt, clims_for_heatmap, PHOTOM_FR, SEC_AROUND_PK, min_prom_Z_other_site,...
    num_trace_segs, ind_threshold, parc_indices, z_scores_i_sorted_to_heat, z_scores_target_sorted_to_heat, corr_standard_all_trials, corr_standard_mean, corr_cross_all_trials, corr_cross_mean, corr_cross_lag_all_trials, corr_cross_lag_mean, corr_global_mean)
% The function will create a Heat-Map given all-sites Z-Scores (that can be
% filtered in a variety of ways):

    % 0. Properties:
    SEC_LAG = 1;        % Number of seconds to define as max window to take cross correlation!
    ind_lag = SEC_LAG * PHOTOM_FR;

    % i. Get z-scores at all sites:
    nan_test = isnan(z_score_mat_all_sites);
    nan_test_mean = mean(mean(nan_test, 2), 1);
    traces_to_include = squeeze(nan_test_mean) == 0;
    z_score_array_chk = z_score_mat_all_sites(:, :, traces_to_include);
        
    % ii. Get z_scores at current site (i_site) and sort
    % according to: 1) 2023/04/27: max-value 2) Can also
    % sort according to center-index.
    z_scores_one_site = z_score_array_chk(i_site, :, :);

    % iii. Now, get pre-sorted z_scores at target site and i_site:
    z_scores_sorted_to_heat = z_scores_i_sorted_to_heat;
    % Calculate and set #'s of inds and rows:
    num_z_inds = size(z_scores_sorted_to_heat, 2);
    num_z_rows = size(z_scores_sorted_to_heat, 1);
    
    % iv. Now, get indices, threshold and correlations for whole time window.
    % Get the z-scores for i_site and target_site:
    corr_section_mean = NaN(num_trace_segs, 1);
    corr_section_lag = NaN(num_trace_segs, 1);
    if corr_reg_cross == 1
        corr_section_mean = corr_standard_mean;
    elseif corr_reg_cross == 2
        corr_section_mean = corr_cross_mean;
        corr_section_lag = corr_cross_lag_mean;
    end
    
    % Now, plot:
    imagesc(h_ax_heat, z_scores_sorted_to_heat, clims_for_heatmap);
    colormap jet;
    hold on;

    % Now draw a boundary on z_scores_sorted_to_heat before
    % peak-area and after peak-area:
    for i_parc = 1:num_trace_segs
        ind_parc_start = parc_indices(i_parc, 1);
        ind_parc_stop = parc_indices(i_parc, 2);
        y_coors_to_plot = [1, size(z_scores_sorted_to_heat, 1)];
        x_coors_to_plot = [ind_parc_start, ind_parc_start];
        plot(x_coors_to_plot, y_coors_to_plot, 'Color', 'black');
        x_coors_to_plot = [ind_parc_stop, ind_parc_stop];
        plot(x_coors_to_plot, y_coors_to_plot, 'Color', 'black');
    end
       
    % If there are more than 1 z_rows, then plot the threshold and
    % correlation values on the graph:
    if num_z_rows > 1
%         y_lim = h_ax_heat.YLim;
%         x_coors = [floor(num_z_inds/6), floor(num_z_inds/2), floor(5*num_z_inds/6)];
%         y_coors_all = 1 + num_z_rows/10;
%         y_coors = [y_coors_all, y_coors_all, y_coors_all];
%         for i_corr = 1:3
%             if corr_reg_cross == 1
%                 corr_text{i_corr} = {'R = ', num2str(corr_section_mean(i_corr), 3)};
%             elseif corr_reg_cross == 2
%                 corr_text{i_corr} = {'X-R = ', num2str(corr_section_mean(i_corr), 3),...
%                 'Lag = ', num2str(corr_section_lag(i_corr), 3)};
%             end
%             text_to_write = strjoin(corr_text{i_corr}, '\n');
% %             text(x_coors(i_corr), y_coors(i_corr), text_to_write, 'Color', 	'magenta', 'FontSize', 11);
%         end

        % Now, plot the threshold line, and description of what
        % that Z-score line is!
        x_coors_to_plot = [1, num_z_inds];
        y_coors_to_plot = [ind_threshold, ind_threshold];
        % if ~isempty(y_coors_to_plot)
        %     plot(x_coors_to_plot, y_coors_to_plot, 'Color', 'red');
        % end
%         text(1, (y_coors_to_plot(1) + num_z_rows/50), text_to_write, 'Color', 'red', 'FontSize', 14);
    end

    % Now, put as title the correlations:
    title_text = [];
    for i_corr = 1:num_trace_segs
        if corr_reg_cross == 1
            r_eq_text = ' R = ';
        elseif corr_reg_cross == 2
            r_eq_text = ' X-R=';
        end
        title_text = strcat(title_text, r_eq_text, num2str(corr_section_mean(i_corr), 3), '||');
        corr_global_str = strcat('Global ', r_eq_text, num2str(corr_global_mean));
    end
%     title_text = title_text(1:(end-2));
    title_text = horzcat(title_text, newline, corr_global_str);
    title(h_ax_heat, title_text);

end