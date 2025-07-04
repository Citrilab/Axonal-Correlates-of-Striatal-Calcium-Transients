function [ind_threshold, parc_indices, proms_sorted_to_corr, z_scores_i_sorted_to_heat, z_scores_i_sorted_to_corr, z_scores_target_sorted_to_heat, z_scores_target_sorted_to_corr, corr_standard_all_trials, corr_standard_mean, corr_cross_all_trials, corr_cross_mean, corr_cross_lag_all_trials, corr_cross_lag_mean] =...
    F_2024_04_07_Create_Corr_Dataset(pk_prominences, z_score_mat_all_sites, target_site, i_site, heat_sort_opt, num_trace_segs, PHOTOM_FR, min_prom_Z_other_site, sec_lag)
% The function will create a Heat-Map given all-sites Z-Scores (that can be
% filtered in a variety of ways): 

    % 0. Properties:
%     SEC_LAG = 1;        % Number of seconds to define as max window to take cross correlation!
    ind_lag = sec_lag * PHOTOM_FR;
 
    % i. Get z-scores at all sites and exclude any z-scores with NaN's in them:
    nan_test = isnan(z_score_mat_all_sites);
    nan_test_mean = mean(mean(nan_test, 2), 1);
    traces_to_include = squeeze(nan_test_mean) == 0;
    z_score_array_chk = z_score_mat_all_sites(:, :, traces_to_include);
    pk_prominences_chk = pk_prominences(traces_to_include);

    % ii. Get z_scores at current site (i_site) and sort them according to
    % specified parameter:
    % Can sort according to: 
    % 1) Prominence of 'independent variable' z-score peak, meaning z-score of target site, which we build peak-table on 
    % 2) Maximum z-score value of 'dependent variable' (that is, i-site) z-score
    % 3) Time-point along time-line of z-score
    z_scores_target_site = z_score_array_chk(target_site, :, :);
    z_scores_i_site = z_score_array_chk(i_site, :, :);
    z_score_target_to_mean = permute(z_scores_target_site, [3, 2, 1]);
    z_score_i_to_mean = permute(z_scores_i_site, [3, 2, 1]);
    [max_values, max_inds] = max(z_score_i_to_mean, [], 2);
    [max_values_target, max_inds_target] = max(z_score_target_to_mean, [], 2);

    % Now sort based on heat_sort_opt condition:
    if heat_sort_opt == 1
        z_scores_i_to_sort = horzcat(pk_prominences_chk, z_score_i_to_mean);
        z_scores_target_to_sort = horzcat(pk_prominences_chk, z_score_target_to_mean);
    elseif heat_sort_opt == 2
        z_scores_i_to_sort = horzcat(max_values, z_score_i_to_mean);
        z_scores_target_to_sort = horzcat(max_values, z_score_target_to_mean);
    elseif heat_sort_opt == 3
        z_scores_i_to_sort = horzcat(max_inds, z_score_i_to_mean);
        z_scores_target_to_sort = horzcat(max_inds, z_score_target_to_mean);
    elseif heat_sort_opt == 4
        z_scores_i_to_sort = horzcat(max_values_target, z_score_i_to_mean);
        z_scores_target_to_sort = horzcat(max_values_target, z_score_target_to_mean);
    else
        % Default is to sort by max-zscore values of i_site:
        z_scores_i_to_sort = horzcat(max_values, z_score_i_to_mean);
        z_scores_target_to_sort = horzcat(max_values, z_score_target_to_mean);
    end

    % Now sort_rows and re-isolate matrix form horzcat sorting vector:
    % i_site
    z_scores_i_sorted = sortrows(z_scores_i_to_sort, 1);
    z_scores_i_sorted_to_heat = z_scores_i_sorted(:, 2:end);
    sort_param_sorted = z_scores_i_sorted_to_heat(:, 1);
    % target-site
    z_scores_target_sorted = sortrows(z_scores_target_to_sort, 1);
    z_scores_target_sorted_to_heat = z_scores_target_sorted(:, 2:end);

    % Now sort z_scores at i_site according to the prominence of target-site z-score
    % peak in order to send for correlation plotting:
    z_scores_i_to_prom_to_sort = horzcat(pk_prominences_chk, z_score_i_to_mean);
    z_scores_i_prom_sorted = sortrows(z_scores_i_to_prom_to_sort, 1);
    z_scores_i_sorted_to_corr = z_scores_i_prom_sorted(:, 2:end);
    proms_sorted_to_corr = z_scores_i_prom_sorted(:, 1);

    % Now sort z_scores at target-site according to the prominence of target-site z-score
    % peak in order to send for correlation plotting:
    z_scores_target_to_prom_to_sort = horzcat(pk_prominences_chk, z_score_target_to_mean);
    z_scores_target_prom_sorted = sortrows(z_scores_target_to_prom_to_sort, 1);
    z_scores_target_sorted_to_corr = z_scores_target_prom_sorted(:, 2:end);

    % iii. Now, calculate the indices above the threshold for values sorted by max-value:
    ind_threshold = [];
    if heat_sort_opt == 1
        max_val_above_thresh = sort_param_sorted > min_prom_Z_other_site;
        inds_above_thresh = find(max_val_above_thresh);
        if ~isempty(inds_above_thresh)
            ind_threshold = inds_above_thresh(1);
        else
            ind_threshold = 0;
        end
    end

    % iv. Now calculate the correlations (either regular correlation or cross-correlation) for peak trace divided into specified segments:
    % Calculate and set #'s of inds and rows:
    num_z_inds = size(z_scores_i_sorted_to_heat, 2);
    num_z_rows = size(z_scores_i_sorted_to_heat, 1);
    
    % Now, cut z-scores traces into x-number of regions to take z-scores
    % over (should ideally be an odd number, and should divide the trace
    % into a integers:
    parcel_ind_num = num_z_inds/num_trace_segs;
    
    % Define variables to store correlations (across each row!):
    corr_standard_all_trials = NaN(num_z_rows, num_trace_segs);
    corr_cross_all_trials = NaN(num_z_rows, num_trace_segs);
    corr_cross_lag_all_trials = NaN(num_z_rows, num_trace_segs);

    % Define variables to store mean correlations:
    corr_standard_mean = NaN(num_trace_segs, 1);
    corr_cross_mean = NaN(num_trace_segs, 1);
    corr_cross_lag_mean = NaN(num_trace_segs, 1);
    parc_indices = NaN(num_trace_segs, 2);

    % Cut scores up into time-segments
    for corr_reg_cross = 1:2
        
        % Define variables to use in calculating correlation:
        corr_answers = NaN(num_z_rows, num_trace_segs);
        corr_lags = NaN(num_z_rows, num_trace_segs);

        for i_parc = 1:num_trace_segs
            % Calculate segment indices:
            ind_parc_start = floor((i_parc-1) * parcel_ind_num) + 1;
            ind_parc_stop = floor((i_parc) * parcel_ind_num);
            parc_indices(i_parc, 1) = ind_parc_start;
            parc_indices(i_parc, 2) = ind_parc_stop;

            % Create target and i-site z-score arrays:
            z_target_to_corr = z_scores_target_sorted_to_corr(:, ind_parc_start:ind_parc_stop);
            z_scores_to_corr = z_scores_i_sorted_to_corr(:, ind_parc_start:ind_parc_stop);
            for i_row = 1:num_z_rows
                z_target_single_row = z_target_to_corr(i_row, :)';
                z_scores_single_row = z_scores_to_corr(i_row, :)';
                if corr_reg_cross == 1
                    corr_standard_all_trials(i_row, i_parc) = corr(z_target_single_row, z_scores_single_row, 'Type', 'Pearson');
                    % corr_standard_mean = mean(corr_standard_all_trials, 1);
                elseif corr_reg_cross == 2
                    z_target_single_row_to_cross = z_target_single_row - mean(z_target_single_row);     % Corr. regular removes the mean as part of fxn; need to do here to compare directly to Pearson's
                    z_scores_single_row_to_cross = z_scores_single_row - mean(z_scores_single_row);     % Corr. regular removes the mean as part of fxn; need to do here to compare directly to Pearson's
                    [cross_corr_val, cross_corr_all_inds] = xcorr(z_target_single_row_to_cross, z_scores_single_row_to_cross, ind_lag, 'coeff');   
                    [max_corr_val, max_corr_ind] = max(cross_corr_val);
                    corr_cross_all_trials(i_row, i_parc) = max_corr_val;
                    corr_cross_lag_all_trials(i_row, i_parc) = cross_corr_all_inds(max_corr_ind)/PHOTOM_FR;
                    % corr_cross_mean = mean(corr_cross_all_trials, 1);
                    % corr_cross_lag_mean = mean(corr_cross_lag_all_trials, 1);
                end 
            end
            % 2024/04/07: Just to double-check corrs were calculcated correctly within
            % loop:
            if corr_reg_cross == 1
                corr_standard_mean = mean(corr_standard_all_trials, 1);
            elseif corr_reg_cross == 2
                corr_cross_mean = mean(corr_cross_all_trials, 1);
                corr_cross_lag_mean = mean(corr_cross_lag_all_trials, 1);
            end
        end

    end

    a = 5;

    end