function [bhv_props] = F_2024_04_07_defineBehavioralParams(exp_props, bhv_choice)
% Choose behavioral properties and enter them 

    % 0. Properties:

    % 1: Bhv score methods:
    for i_bhv = bhv_choice
        if i_bhv == 1
            [stereo_bhv_key, stereo_bhv_colormap] = F_2022_08_21_STEREO_loadBehavKeyCMap();
            bhv_props.bhv_key = stereo_bhv_key;
            bhv_props.colormap = stereo_bhv_colormap; 
            bhv_props.bhv_scr_type = 'STEREO';
        elseif i_bhv == 2
            [slit_tag_bhv_key, slit_tag_bhv_colormap] = F_2023_04_24_SLIT_TAG_loadBehavKeyCMap();
            bhv_props.bhv_key = slit_tag_bhv_key;
            bhv_props.colormap = slit_tag_bhv_colormap;
            bhv_props.bhv_scr_type = 'Slittag';
            bhv_props.bhv_filt_params = [0, 5];
        elseif i_bhv == 3
            % Future: combination of stereo and slit-tag!
        end
    end

    % 2: Define filtering props:
    % If not empty, then define filtering z_prom_params for BASE-COMPARE-SITE:
    if isempty(exp_props.filt_params) && exp_props.event_filt_cond == 2
        bhv_props.bhv_filt_params = {0, size(bhv_props.bhv_key, 1)};
    else 
        bhv_props.bhv_filt_params = exp_props.filt_params;
    end

end