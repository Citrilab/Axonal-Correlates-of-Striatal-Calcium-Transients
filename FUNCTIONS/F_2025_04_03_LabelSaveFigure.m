function [curr_fig] = F_2025_04_03_LabelSaveFigure(curr_fig, site_specific_dir, main_title_string, save_title_string, plot_type_string, xlabel_string, ylabel_string, xtick_create, xlims_on_ax, zscr_inds_around_peak, PHOTOM_FR, fig_position, clr_bar_leg)
% F_2023_01_08_LabelSaveFigure: This function will label and save a figure:

    % Input all label strings and titles. Create xtick timeline for
    % displaying as x-axis on graphs.

    % Query if the figure is a tiled figure or not: % 1 if tiled, 0 if not
    % tiled!
    if ~isempty(curr_fig.Children)
        is_graph_tiled = strfind(curr_fig.Children(1).Type, 'tiled');
        tiled_yes = ~isempty(is_graph_tiled);
        if tiled_yes == 0
            curr_ax = gca;
        elseif tiled_yes == 1
            curr_ax = curr_fig.Children;
        end
    else
       curr_ax = gca;
    end

    % Label axis with title, xlabel, and ylabel
        title(curr_ax, main_title_string, 'FontSize', 18);
        xlabel(curr_ax, xlabel_string, 'FontSize', 18);
        ylabel(curr_ax, ylabel_string, 'FontSize', 18);
     
    % Single graph options: optional (if parameter exists), create x-axis
    if xtick_create == 1        % Create timeline on x-axis
        xticks = xlims_on_ax(1):PHOTOM_FR:xlims_on_ax(end);
        xtick_nums = (-zscr_inds_around_peak/PHOTOM_FR):1:(zscr_inds_around_peak/PHOTOM_FR);
        for i_tick = 1:length(xtick_nums)
            xtick_str{i_tick} = num2str(xtick_nums(i_tick));
        end
        curr_ax.XTick = xticks;
        curr_ax.XTickLabel = xtick_str;
    elseif xtick_create == 2    % Create x-axis of prominence values
        xticks = xlims_on_ax(1):1:xlims_on_ax(end);
        xtick_nums = xticks;
        for i_tick = 1:length(xtick_nums)
            xtick_str{i_tick} = num2str(xtick_nums(i_tick));
        end
        curr_ax.XTick = xticks;
        curr_ax.XTickLabel = xtick_str;
    end

    % Tiled graph options: add optional (if parameter exists) colorbar
    if ~isempty(clr_bar_leg)
        h_clr_bar = colorbar;
        h_clr_bar.Label.String = clr_bar_leg;
        h_clr_bar.Label.FontSize = 18;
        if ~isempty(h_clr_bar.Layout)
            h_clr_bar.Layout.Tile = 'east';
        end
    end

    % Optional (if parameter exists), change position on screen of figure:
    if ~isempty(fig_position)
        curr_fig.Position = fig_position;
    end

    % Now Save Figure:
    % Save figure: 
    save_name_root = char(save_title_string);
    save_name_root(save_name_root == ':') = '_';
    save_name_root(save_name_root == '=') = '';
    save_name_root(save_name_root == '-') = '_';
    save_name_root(save_name_root == ' ') = '_';
    save_name_root(save_name_root == '.') = 'p';
    save_name = strcat(site_specific_dir, save_name_root);
    % Save as matlab figure
    saveas(curr_fig, strcat(save_name, plot_type_string));

    % Now save as svg with modifiable text properties:
    saveas(curr_fig, save_name, 'fig');
    % Ensure text is stored as text, not outlines
    set(groot, 'defaultAxesFontName', 'none');
    % Save as SVG using 'print' instead of 'saveas'
    set(findall(curr_fig, '-property', 'FontName'), 'FontName', 'Arial');
    print(curr_fig, save_name, '-dsvg', '-vector');               
end