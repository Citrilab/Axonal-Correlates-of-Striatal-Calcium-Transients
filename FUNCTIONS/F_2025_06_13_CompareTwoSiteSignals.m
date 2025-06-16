function [peak_table_wCorr, peak_table_name] = F_2025_06_13_CompareTwoSiteSignals(trial_data_to_event, trial_table_target_site, i_target, corr_reg_cross, sec_lag, exp_props, i_list)

    % A. DEFINE ARRAYS: First get all z_scores for all mice at target-site (site to filter) for global correlation calculations:
    z_scores_all = trial_data_to_event(i_list).z_scores_allsites{:};
    nan_mean_vec = mean(z_scores_all, 1)';
    not_nan_log = ~isnan(nan_mean_vec);
    z_scores_not_nan = z_scores_all(:, not_nan_log);
    zscr_vec_target = z_scores_not_nan(i_target, :)';

    % B. CORRELATIONS: Second, at each reference site, figure out global correlations to
    % each of 3 other sites.
    for i_site = 1:exp_props.num_sites
        % Set names for columns in peak table containing global
        % correlations and global lags:
        global_corr_name = strcat(exp_props.site_names{i_site}, 'global_corrs');
        global_lag_name = strcat(exp_props.site_names{i_site}, 'global_lags');
        % Figure out correlations
        target_i_corr = zeros(exp_props.num_sites, 1);
        target_i_lag = target_i_corr;
        zscr_vec_i = z_scores_not_nan(i_site, :)';
        if corr_reg_cross == 2      % Calculate cross-correlation if corr_reg_cross == 2
            inds_lag = sec_lag * exp_prop.pho.PHOTOM_FR;
            zscr_vec_target_to_corr = zscr_vec_target - mean(zscr_vec_target);
            zscr_vec_i_to_corr = zscr_vec_i - mean(zscr_vec_i);                        
            [corr_vals, lag_vals] = xcorr(zscr_vec_target, zscr_vec_i, inds_lag, 'coeff');
            [corr_max, corr_max_ind] = max(corr_vals);
            target_i_corr(i_site) = corr_max;
            target_i_lag(i_site) = lag_vals(corr_max_ind)/PHOTOM_FR;
            correlations_matrix(i_list, i_target, i_site) = target_i_corr(i_site);
            corr_lags_matrix(i_list, i_target, i_site) = target_i_lag(i_site);
            trial_table_target_site{:}.(global_corr_name) = ones(size(trial_table_target_site{:}, 1), 1) .* target_i_corr(i_site);
            trial_table_target_site{:}.(global_lag_name) = ones(size(trial_table_target_site{:}, 1), 1) .* target_i_lag(i_site);
        else                        % Else calculate standard correlation
            target_i_corr(i_site) = corr(zscr_vec_target, zscr_vec_i);
            correlations_matrix(i_list, i_target, i_site) = target_i_corr(i_site);
            trial_table_target_site{:}.(global_corr_name) = ones(size(trial_table_target_site{:}, 1), 1) .* target_i_corr(i_site);
        end
    end
        
    % Output version of peak-table with peak-compare stats:
    peak_table_wCorr = [trial_table_target_site{:}];    
    
    % C. SAVE! All names to use in saving are below!!!     % Define variable names to save the peak table:
    if exp_props.event_filt_cond == 1                                           % Event to use to filter is peak!
        low_thresh_num = num2str(exp_props.pk.z_prom_params{1});
        low_thresh_num(low_thresh_num == '.') = 'p';
        peak_table_name = strcat('peak_table_', exp_props.site_names{i_target}, '_', low_thresh_num);
        peak_table_name(peak_table_name == '-') = '_';
    elseif exp_props.event_filt_cond == 2 || exp_props.event_filt_cond == 3     % Event to use to filter is behavioral event!
        peak_table_name = strcat('behav_table_', filt_string_names{pk_filt_z_bhv}, '_', num2str(bhv_scr_filt{1}), 'to', num2str(bhv_scr_filt{2}));
        peak_table_name(peak_table_name == '-') = '_';
    end

end