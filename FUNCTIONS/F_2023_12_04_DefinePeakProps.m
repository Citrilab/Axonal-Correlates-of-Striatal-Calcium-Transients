function [peak_props] = F_2023_12_04_DefinePeakProps(exp_props)

% Input: exp_props which have user toggled properties
% Output: modified exp_props which have additionally:

    % 0. Peak-table constants:
    peak_props.SECS_PEAK_MOV = 6;      % This is the # of seconds before/after the peak to isolate! It is the size of the window in either direction. Whole window is twice this amount.
    peak_props.IMG_FRAME_GAP = 15;

    % 1. What is the 'event' that you use to filter peak-table with:
    peak_props.filt_string_names = {'Peak', 'STEREO', 'Slit'};
    peak_props.pk_filt_z_bhv = 1;   % If pk_filt_z_bhv = 1 --> filter according to peak,...
                                      % If pk_filt_z_bhv = 2 --> filter according to stereo-predict behavior
                                      % If pk_filt_z_bhv = 3 --> filter according to slittag behavior
                          
    % 2. Define z-score parameters for program:
    
    % If not empty, then define filtering z_prom_params for BASE-COMPARE-SITE:
    if isempty(exp_props.filt_params) && exp_props.event_filt_cond == 1
        peak_props.z_prom_params = {2.5, 20.0};
    else 
        peak_props.z_prom_params = exp_props.filt_params;
    end

    % Should be a single interval, a cell array of 1 x 2:
    % define filtering z_prom_params for OTHER-COMPARE-SITE:
    peak_props.z_proms_comp_site = {2.0, 20};

    % 3. Define the minimum prominence (or Z Score) of the other site, to be use for i) Jacquard calculations ii) Treshold line on graph:
    peak_props.min_prom_Z_other_site = 2.0; 
    peak_props.import_plot_z_peak_min_prom = peak_props.z_prom_params{1};
    peak_props.jacc_same_z_or_low_thresh = 1;          % Same z- = 0, means that you calculate Jaccard with other site z-params the same as target site, z = 1 means low_thresh params

    % 4. Make strings of above parameters for plotting:
    for i_prom_cmpr = 1:size(peak_props.z_proms_comp_site, 1)
        z_prom_cmpr_minimum = peak_props.z_proms_comp_site{i_prom_cmpr, 1};
        z_prom_cmpr_min_str = num2str(z_prom_cmpr_minimum);
        z_prom_cmpr_min_str(z_prom_cmpr_min_str == '.') = 'p';
        % Save strings of z_prom_cmpr:
        peak_props.z_prom_cmpr_strings{i_prom_cmpr} = z_prom_cmpr_min_str;
    end
    
    % 5. Define plotting parameters related to jaccard:
    num_bins_for_jacc = 4;
    peak_props.num_bins_for_jacc = num_bins_for_jacc;
    jacc_pk_tbl_groups = zeros(2, num_bins_for_jacc);
    jacc_pk_tbl_groups(2, :) = (1:num_bins_for_jacc)/num_bins_for_jacc;
    jacc_pk_tbl_groups(1, :) = jacc_pk_tbl_groups(2, :) - jacc_pk_tbl_groups(2, 1);
    
    % Define Jaccard peak table groups:
    peak_props.jacc_groups = jacc_pk_tbl_groups;

    % Now - save output folder:
    if exp_props.event_filt_cond == 1
        exp_props.output_dir = strcat(exp_props.output_date_input, peak_props.z_prom_cmpr_strings);
    end

end