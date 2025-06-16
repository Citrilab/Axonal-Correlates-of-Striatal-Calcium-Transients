function [exp_props, trial_list_sync_hemi, num_ms_trials, time_sync_table] = F_2024_01_29_CombineHemispheres(exp_props, trial_list_sync, time_sync_table)
%F_2023_11_29_COMBINEHEMISPHERES Summary of this function goes here
%   Detailed explanation goes here

    % 0. Properties:
    left_hemi_order = [1, 3, 2, 4];
    right_hemi_order = [2, 4, 1, 3];

    % 1. Rename site names to be either 'Ipsi' or 'Contra', with either
    % 'left' first or 'right' first.
%     zscr_strings_plot_hemi = exp_props.pho.zscr_strings_plot(left_hemi_order);
%     for i_site = 1:num_sites
%         zscr_strings_plot_hemi{i_site, 1} = strrep(zscr_strings_plot_hemi{i_site, 1}, '-Left', 'Ipsi');
%         zscr_strings_plot_hemi{i_site, 1} = strrep(zscr_strings_plot_hemi{i_site, 1}, '-Right', 'Contra');
%     end
%     exp_props.pho.zscr_strings_plot = zscr_strings_plot_hemi;


    % 2. Code
    for i_hemi = 1:2

        if i_hemi == 1
            temp_data_struct = trial_list_sync;
            hemi_string = 'Left';
            site_order = left_hemi_order;
        elseif i_hemi == 2
            temp_data_struct = trial_list_sync;
            hemi_string = 'Right';
            site_order = right_hemi_order;
        end

        % Now go through entire list and import cg_ms_base_name and
        % all trials vectors:
        for i_list = 1:length(trial_list_sync)
            % Get current mouse name:
            curr_ms_name = temp_data_struct(i_list).mouse;
            curr_ms_name_hemi = strcat(curr_ms_name, hemi_string);
            temp_data_struct(i_list).mouse = curr_ms_name_hemi;
            % Now, if combine_hemispheres == 1, and you should recombine
            % hemispheres, then change the order of all raw data to store as
            % VLS-ipsi, VLS-contra, SNR/downstream-ipsi, SNR/downstream-contra:
            % Now re-order the hemispheres:
            var_names = {'z_scores_allsites', 'z_timelines_allsites', 'z_controls_allsites', 'signal_raw_allsites', 'isos_raw_allsites'};
            for i_var = 1:5
                var_name = var_names{i_var};
                temp_data = trial_list_sync(i_list).(var_names{i_var}){:};
                temp_data_rearrange = temp_data(site_order, :);
                temp_data_struct(i_list).(var_names{i_var}) = {temp_data_rearrange};
            end
        end

        if i_hemi == 1
            data_struct_left = temp_data_struct;
        elseif i_hemi == 2
            data_struct_right = temp_data_struct;
            trial_list_sync_hemi_iso = horzcat(data_struct_left, data_struct_right);
            time_sync_table = repmat(time_sync_table, [2, 1]);
            trial_list_sync_hemi = trial_list_sync_hemi_iso;
            num_ms_trials = length(trial_list_sync_hemi_iso);
        end

    end

    % Rename site strings (do this only once, and with left_hemi_order, since the original order was left-first):
    site_order = left_hemi_order;
    zscr_strings_plot_hemi = exp_props.site_names(site_order);
    for i_site = 1:exp_props.num_sites
        zscr_strings_plot_hemi{i_site, 1} = strrep(zscr_strings_plot_hemi{i_site, 1}, '-Left', 'Ipsi');
        zscr_strings_plot_hemi{i_site, 1} = strrep(zscr_strings_plot_hemi{i_site, 1}, '-Right', 'Contra');
    end
    exp_props.site_names = zscr_strings_plot_hemi;

end

