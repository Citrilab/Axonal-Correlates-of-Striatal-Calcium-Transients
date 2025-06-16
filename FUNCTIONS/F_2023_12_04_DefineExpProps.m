function [exp_props, photom_props] = F_2023_12_04_DefineExpProps(exp_props, input_folder_string)

% Input: exp_props which have user toggled properties
% Output: modified exp_props which have additionally:

    % 0. Define NUMBERS/CONSTANTS to be used in plot-program:
    for i = 1:1
        % A. Photometry properties:
        photom_props.PEAK_Z_MIN = 3;
        photom_props.MOVIE_FR = 15;            % Movie FRAME RATE -  crucial: different from photom_FR
        photom_props.PHOTOM_FR = 15;
        photom_props.PHOTOM_SEC_START_SUBTRACT = 30;
        photom_props.PHOTOM_SEC_END_SUBTRACT = 30;

        % B. Behavior properties/constants:
        exp_props.bhv.BHV_SCR_FR = 15;
        exp_props.bhv.bhv_scr_min_max = {0, 6};
        exp_props.bhv.bhv_scr_filt = {2, 6};        
    end

    % 1. Define STRINGS/NAMES to be used in plotting:
    for i = 1:1
        % A. Number of photometry sites to be processed for each ms_trial for use in program
        num_sites = 4;
        zscr_strings_input = cell(num_sites, 1);
        zscr_strings_input(1) = {'VLS_Left'};
        zscr_strings_input(2) = {'VLS_Right'};
        zscr_strings_input(3) = {'SNR_Left'};
        zscr_strings_input(4) = {'SNR_Right'};
        zscr_strings_plot = zscr_strings_input;
    
        % B. Modify zscr 'strings' and change from 'dSPN' strings to 'iSPN':
        for i_site = 1:num_sites
            zscr_strings_plot{i_site, 1}(zscr_strings_plot{i_site, 1}=='_') = '-';
            if contains(input_folder_string, 'iSPN')
                zscr_strings_plot{i_site, 1} = strrep(zscr_strings_plot{i_site, 1}, 'SNR', 'GPe');
            end
        end
    end

    % 2. Get current datetime
    datetime_str = char(datetime);
    date_str = datetime_str(1:11);

    % 3. Param-strings
    param_one_str = num2str(exp_props.filt_params{1});
    param_one_str(param_one_str == '.') = 'p';
    param_two_str = num2str(exp_props.filt_params{2});
    param_two_str(param_one_str == '.') = 'p';

    % 4. Set all exp_props:
    exp_props.filt_param_str = {param_one_str, param_two_str};
    exp_props.num_sites = num_sites;
    exp_props.site_names = zscr_strings_plot;
    exp_props.output_date_input = strcat(date_str, input_folder_string);