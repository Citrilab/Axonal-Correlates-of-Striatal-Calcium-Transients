function [h_ax] = F_2023_01_08_LabelAxes(h_ax, main_title_string, xlabel_string, ylabel_string, zscr_inds_around_peak, PHOTOM_FR, x_tick_inds, clr_bar_leg)
% F_2023_01_08_LabelSaveFigure: This function will label and save a figure:

    % Create title of graph
    title(main_title_string);
    
    % Create time-axes around peak:
    x_tick_vals = (-zscr_inds_around_peak:PHOTOM_FR:zscr_inds_around_peak);
    x_tick_secs = x_tick_vals/PHOTOM_FR;
    x_tick_str = cell(size(x_tick_secs, 2), 1);
    for i_pt = 1:size(x_tick_secs, 2)
        x_tick_str{i_pt} = num2str(x_tick_secs(i_pt));
    end
    if isempty(x_tick_inds)
        x_tick_inds = x_tick_vals + zscr_inds_around_peak;
    end
    xlim([x_tick_inds(1), x_tick_inds(end)]);
    xticks(x_tick_inds);
    xticklabels(x_tick_str);
    xlabel(xlabel_string);    
    ylabel(ylabel_string);

    % Add optional colorbar
    if ~isempty(clr_bar_leg)
        h_ax = gca;
        h_clr_bar = colorbar(h_ax);
        h_clr_bar.Label.String = clr_bar_leg;
    end
                    
end