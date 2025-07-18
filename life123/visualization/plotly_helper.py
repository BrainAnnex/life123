import numpy as np
import pandas as pd
import math
import plotly.express as px
import plotly.graph_objects as pgo
import plotly.subplots as sp
from life123.visualization.colors import Colors
from typing import Union



class PlotlyHelper:
    """
    Static class to assist in the use of the plotly library.
    For color management, see the separate class "Colors"

    TODO: turn the class methods into static methods?
    TODO: improve consistency in argument names; also across UniformCompartment
    TODO: maybe rename to "VisualizationHelper" or "GraphicsHelper"
    """


    #####################################################################################################

    '''                                    ~   PLOTS   ~                                           '''

    def ________PLOTS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    @classmethod
    def plot_curves(cls, x, y, title="", range_x=None, x_label="", y_label="", curve_labels=None, legend_title=None,
                    colors=None, show=False) -> pgo.Figure:
        """
        Plot one or more 2D curves.  If your data is a in a Pandas dataframe, consider using plot_pandas()

        :param x:           A Numpy array, with the (common) x-axis values
        :param y:           Either a Numpy array, or a list/tuple of them, with the y-axis values of the curve(s)
        :param title:       [OPTIONAL] Title to use for the overall plot
        :param range_x:     [OPTIONAL] list of the form [t_start, t_end], to initially only show a part of the timeline.
                                Note: it's still possible to zoom out, and see the excluded portion
        :param x_label:     [OPTIONAL] Caption to use for the x-axis
        :param y_label:     [OPTIONAL] Caption to use for the y-axis
        :param curve_labels:[OPTIONAL] String, or list of strings.
                                Label(s) to use for the various curves in the legend and in the hover boxes.
                                If missing, and there's only 1 curve, the legend box won't be shown
        :param legend_title:[OPTIONAL] String to show at the top of the legend box.
                                Ignored if curve_labels wasn't set.
                                If not provided, and the legend box is shown, it will appear as "variable"
        :param colors:      [OPTIONAL] Either a single color (string with standard plotly CSS name, such as "red"),
                                or list of names to use, in the same order as the y variables; if None, then use the hardwired defaults
        :param show:        If True, the plot will be shown
                                Note: in JupyterLab, simply returning a plot object (without assigning it to a variable)
                                leads to it being automatically shown

        :return:            A plotly "Figure" object
        """
        if type(y) == list or type(y) == tuple:
            number_of_curves = len(y)
        else:
            number_of_curves = 1
            if curve_labels:
                y = [y]     # This will trigger the legend to be shown
                if type(curve_labels) == str:
                    curve_labels = [curve_labels]

        if colors is None:
            colors = Colors.assign_default_colors(number_of_curves)
        elif type(colors) == str:
            colors = [colors]
        # TODO: if any color is missing, assign default ones


        fig = px.line(x=x, y=y, color_discrete_sequence=colors)     # y can be one array or a list

        if not x_label:
            x_label = "x"

        if not y_label:
            y_label = "y"

        if number_of_curves == 1:
            hovertemplate = f"{x_label}=%{{x}}<br>{y_label}=%{{y}}<extra></extra>"
        else:
            hovertemplate = f"{x_label}=%{{x}}<br>value=%{{y}}<extra></extra>"
            # The "<extra></extra>" appears to suppress a redundant name showing in the hoverbox

        for index in range(number_of_curves):
            if curve_labels is not None:
                label = curve_labels[index]
                fig.data[index]["name"] = label
                fig.data[index]["hovertemplate"] = f"{label} :<br>{hovertemplate}"    # variable :
            else:
                fig.data[index]["hovertemplate"] = hovertemplate


        fig.update_layout(title=title,
                          xaxis_title=x_label,
                          yaxis_title=y_label)

        if legend_title:
            fig.update_layout(legend={"title": legend_title})

        if range_x:
            fig.update_layout(xaxis={"range": range_x})


        if show:
            fig.show()  # Actually display the plot

        return fig



    @classmethod
    def plot_pandas(cls, df :pd.DataFrame, x_var="SYSTEM TIME", fields=None,
                    log_y=False,
                    colors=None, title=None, title_prefix=None,
                    range_x=None, range_y=None,
                    x_label=None, y_label="Y", legend_header="Plot",
                    vertical_lines_to_add=None, show_intervals=False,
                    smoothed=False, show=False) -> pgo.Figure:
        """
        Using plotly, draw line plots from the values in the given dataframe.
        One column supplies the values for the x-axis,
        and one or more other columns supply the y-axis values for one or more line plots.

        Note: if this plot is to be later combined with others, use PlotlyHelper.combine_plots()

        :param df:              Pandas dataframe with the data for the plot
        :param x_var:           Name of column with the independent variable for the x-axis
        :param fields:          Name, or list of names, of the dataframe columns whose values are to be plotted;
                                    if a list is passed, also display a figure legend;
                                    if None, then display all columns EXCEPT the one that
                                    was declared as the independent variable thru argument `x_var`
        :param log_y:           [OPTIONAL] If True, a log scale is used for the y-axis
        :param colors:          [OPTIONAL] Either a single color (string with standard plotly name, such as "red"),
                                    or list of names to use, in order; some of the entries may be None.
                                    If None, then the hardwired default colors are used
        :param title:           [OPTIONAL] Title for the plot
        :param title_prefix:    [OPTIONAL] String to prefix (automatically followed by " <br>") to the title
        :param range_x:         [OPTIONAL] list of the form [t_start, t_end], to initially show only a part of the timeline.
                                    Note: it's still possible to zoom out, and see the excluded portion
        :param range_y:         [OPTIONAL] list of the form [y_min, y_max], to initially show only a part of the y values.
                                    Note: it's still possible to zoom out, and see the excluded portion
        :param x_label:         [OPTIONAL] Caption to use for the x-axis
        :param y_label:         [OPTIONAL] Caption to use for the y-axis.  Default: "Y"
        :param legend_header:   [OPTIONAL] Caption to use at the top of the legend box.
                                            Only applicable if more than 1 curve is being shown.
        :param vertical_lines_to_add:  [OPTIONAL] Ignored if the argument `show_intervals` is specified.
                                    Value, or list, or tuple, or Numpy array, or Pandas series,
                                    of x-coordinate(s) at which to draw thin vertical dotted gray lines.
                                    If the number of vertical lines is so large as to overwhelm the plot,
                                    only a sample of them is shown.
                                    Note that vertical lines, if requested, go into the plot's "layout";
                                    as a result they might not appear if this plot is later combined with another one.
        :param show_intervals:  [OPTIONAL] If True, it over-rides any value passed to the `vertical_lines` argument,
                                    and draws thin vertical dotted gray lines at all the x-coords
                                    of the data points in the saved history data;
                                    also, it adds a comment to the title.
        :param smoothed:        [OPTIONAL] If True, a spline is used to smooth the lines;
                                    otherwise (default), line segments are used
        :param show:            If True, the plot will be shown
                                    Note: on JupyterLab, simply returning a plot object (without assigning it to a variable)
                                          gets it automatically shown

        :return:                A plotly "Figure" object
        """
        MAX_NUMBER_VERTICAL_LINES = 150     # Used to avoid extreme clutter in the plot, in case
                                            # a very large number of vertical lines is requested;
                                            # if this value is exceeded, then the vertical lines are sampled
                                            # infrequently enough to bring the total number below this value

        col_list = list(df.columns)
        assert x_var in col_list, \
            f"plot_pandas(): the value of the argument `x_var` ({x_var}) must be the name of one of the columns in the dataframe"

        # Prevent obscure error messages that arise when y_label is not a valid string (due to the fact that this string is used
        # to shows values in the hover boxes when more than 1 curve is being shown)
        assert type(y_label) == str, \
            f"plot_pandas(): the value of the argument `y_label` must be a string; the passed values is of type {type(y_label)}"


        if fields is None:
            number_of_curves = len(col_list) - 1  # All the field but one (since one is the independent variable)
            fields = col_list
            fields.remove(x_var)
            if len(fields) == 1:
                fields = fields[0]      # Take the only element, if there's only one (this has the effect of hiding the legend)
        else:
            number_of_curves = len(fields)

        if colors is None:
            # Entirely use default colors
            colors = Colors.assign_default_colors(number_of_curves)
        elif type(colors) == str:
            # Turn colors into a list, if it was a single entry
            colors = [colors]
        else:
            # If we get here, we were given a list; replace any missing (None) entry with a default color
            number_none = colors.count(None)    # Number of None entries
            if number_none > 0:
                replacement_colors = Colors.assign_default_colors(number_none)   # Get all the replacements in bulk
                colors_adjusted = []
                i = 0
                for c in colors:
                    if c is None:
                        colors_adjusted.append(replacement_colors[i])
                        i += 1
                    else:
                        colors_adjusted.append(c)

                colors = colors_adjusted


        if title_prefix is not None:
            title = f"{title_prefix} <br>{title}"

        if show_intervals:
            vertical_lines_to_add = df[x_var]   # Make use of the simulation times
            title += " (time steps shown in dashed lines)"


        # Create the main plot
        line_shape = "spline" if smoothed else "linear"

        fig = px.line(data_frame=df, x=x_var, y=fields,
                      title=title, range_x=range_x, range_y=range_y,
                      log_y=log_y,
                      color_discrete_sequence = colors,
                      labels={"value": y_label, "variable": legend_header},
                      line_shape=line_shape)

        if type(fields) == str:     # Somehow, the `labels` argument in px.line, above, is ignored when `fields` is just a string
            fig.update_layout(yaxis_title=y_label)   # This line will remedy the above issue

        if x_label is not None:
            fig.update_layout(xaxis_title=x_label)   # Over-ride the default naming of the x-axis


        if vertical_lines_to_add is not None:   # User requested to add vertical lines to the plot
            if isinstance(vertical_lines_to_add, (float, int)):
                vertical_lines_to_add = [vertical_lines_to_add]
            else:
                assert (type(vertical_lines_to_add) == list) or (type(vertical_lines_to_add) == tuple) \
                       or (type(vertical_lines_to_add) == np.ndarray) or (type(vertical_lines_to_add) == pd.core.series.Series), \
                            "plot_pandas(): the argument `vertical_lines`, " \
                            "if not None or a number, must be a list or tuple or Numpy array or Pandas series of numbers (x-axis coords)"

            vline_list = []
            if range_x:
                step = 1    # Always show all vertical lines if a range on the x-axis was specified
            else:
                # Possibly limit the number of vertical lines shown
                step = 1 + len(vertical_lines_to_add) // MAX_NUMBER_VERTICAL_LINES
                if step > 1:
                    print(f"plot_pandas() NOTICE: Excessive number of vertical lines ({len(vertical_lines_to_add)}) - only showing 1 every {step} lines")

            for xi in vertical_lines_to_add[::step] :  # Notice that if step > 1 then we're sampling a subset of the elements
                # The following is the internal data structure used by Plotly Express,
                # for each of the vertical lines
                vline = {  'line': {'color': 'gray', 'dash': 'dot', 'width': 1},
                           'type': 'line',
                           'x0': xi,
                           'x1': xi,
                           'xref': 'x',
                           'y0': 0,
                           'y1': 1,
                           'yref': 'y domain'
                           }
                vline_list.append(vline)
                # Strangely, a direct call to fig.add_vline(), as done below, dramatically slows things down in case
                # of a large number of vertical lines; so, we'll be directly modifying the data structure of the "fig" dictionary
                #fig.add_vline(x=xi, line_width=1, line_dash="dot", line_color="gray")
            # END for

            fig['layout']['shapes'] = vline_list    # The vertical lines are regarded by Plotly Express as "shapes"
                                                    # that are stored in the figure's "layout"
        if show:
            fig.show()  # Display the plot

        return fig



    @classmethod
    def combine_plots(cls, fig_list :Union[list, tuple], layout_index=None,
                      title="", x_label=None, y_label=None,
                      xrange=None, legend_title=None, curve_labels=None, show=False) -> pgo.Figure:
        """
        Combine together several existing Plotly plots into a single one (with combined axes)

        EXAMPLE:
                    from life123 import PlotlyHelper
                    plot_1 = PlotlyHelper.plot_pandas(various args, show=False)
                    plot_2 = PlotlyHelper.plot_pandas(various args, show=False)
                    PlotlyHelper.combine_plots([plot_1, plot_2], other optional args)

        :param fig_list:    List or tuple of plotly "Figure" objects (as returned by several functions)
        :param layout_index:[OPTIONAL] If given, the layout of the "Figure" object with the given index
                                (in `fig_list`) is used as is - and all the layout parameters below are ignored
        :param title:       [OPTIONAL] The title to use for the overall plot
        :param x_label:     [OPTIONAL] Caption to use for the x-axis; if not specified, use that of the 1st plot
        :param y_label:     [OPTIONAL] Caption to use for the y-axis; if not specified, use that of the 1st plot
        :param xrange:      [OPTIONAL] list of the form [t_start, t_end], to initially only show a part of the timeline.
                                Note: it's still possible to zoom out, and see the excluded portion
        :param legend_title:[OPTIONAL] String to show at the top of the legend box
        :param curve_labels:[OPTIONAL] List of labels to use for the various curves in the legend
                                and in the hover boxes; if not specified, use the titles of the individual plots
        :param show:        [OPTIONAL] If True, the plot will be shown
                                Note: on JupyterLab, simply returning a plot object (without assigning it to a variable)
                                      leads to it being automatically shown
        :return:            A plotly "Figure" object for the combined plot
        """
        assert (type(fig_list) == list) or (type(fig_list) == tuple), \
            "combine_plots(): the argument fig_list must be a list or tuple"


        representative_fig = fig_list[0]    # The axes titles of the first plot are used as default values
        if x_label is None:
            x_label = representative_fig.layout.xaxis.title.text
        if y_label is None:
            y_label = representative_fig.layout.yaxis.title.text


        # Put together the data from all the various individual plots
        combined_data = []
        for fig in fig_list:
            combined_data += fig.data      # concatenating lists


        if layout_index is not None:
            assert 0 <= layout_index < len(fig_list), \
                f"combine_plots(): argument {layout_index} must be an integer between 0 and {len(fig_list)-1}, inclusive"

            figure_providing_layout = fig_list[layout_index]
            all_fig = pgo.Figure(data=combined_data, layout = figure_providing_layout.layout)
        else:
            all_fig = pgo.Figure(data=combined_data)    # Note that the + is concatenating lists

            all_fig.update_layout(title=title,
                                  xaxis_title=x_label,
                                  yaxis_title=y_label)

            if legend_title:
                all_fig.update_layout(legend={"title": legend_title})

            if xrange:
                all_fig.update_layout(xaxis={"range": xrange})


            for i, fig in enumerate(fig_list):
                all_fig['data'][i]['showlegend'] = True
                if curve_labels:
                    all_fig['data'][i]['name'] = curve_labels[i]
                    all_fig.data[i]['hovertemplate'] = f"{curve_labels[i]}<br>" + all_fig.data[i]['hovertemplate']
                else:
                    all_fig['data'][i]['name'] = fig.layout.title.text
                    all_fig.data[i]['hovertemplate'] = f"{fig.layout.title.text}<br>" + all_fig.data[i]['hovertemplate']


        if show:
            all_fig.show()  # Actually display the plot

        return all_fig





    #####################################################################################################

    '''                                    ~   SUBPLOTS   ~                                           '''

    def ________SUBPLOTS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    @classmethod
    def combine_in_vertical_grid(cls, fig1, fig2, title1 :str, title2 :str,
                                 title_combined :str, height=900) -> pgo.Figure:
        """
        Combine into a vertical grid 2 plotly graph objects (now treated as subplots)

        :param fig1:
        :param fig2:
        :param title1:
        :param title2:
        :param title_combined:
        :param height           [OPTIONAL] The overall height of the combined plots
        :return:                A combined plot (a plotly "figure" object)
        """

        # Make a visual grid with 2 rows and 1 column; adequate vertical spacing is used to clearly shows axes labels
        combined_fig = sp.make_subplots(rows=2, cols=1,
                                        subplot_titles=(title1, title2),
                                        vertical_spacing=0.15)  # Increase vertical spacing (as a fraction of the plot height)

        # Populate the grid with all the inner parts ("traces") of each individual plot ("figure" object)
        for trace in fig1.data:
            combined_fig.add_trace(trace, row=1, col=1)
        for trace in fig2.data:
            combined_fig.add_trace(trace, row=2, col=1)

        # The axes labels must be copied over separately
        combined_fig.update_xaxes(title_text=fig1.layout.xaxis.title.text, row=1, col=1)
        combined_fig.update_yaxes(title_text=fig1.layout.yaxis.title.text, row=1, col=1)

        combined_fig.update_xaxes(title_text=fig2.layout.xaxis.title.text, row=2, col=1)
        combined_fig.update_yaxes(title_text=fig2.layout.yaxis.title.text, row=2, col=1)

        # Add a title to the combined layout, and set its overall height
        combined_fig.update_layout(title=title_combined, height=height)    # height controls the total height of the entire figure

        return combined_fig



    @classmethod
    def _optimal_subplot_grid_size(cls, ncells :int, max_n_cols = 4) -> (int, int):
        """
        Given the specified number of "cells" (graphic elements) to combine into a subplot grid,
        and given the desired max number of columns to use,
        suggest a visually-pleasing, least crowded, number of rows and columns to use.
        For the needed number of rows, attempt to use fewer columns if possible.
        Note: NO graphic elements are actually produced; this is a helper function

        EXAMPLES, with max_n_cols = 4 :
            3 cells -> 1 row, 3 cols
            4 cells -> 1 row, 4 cols
            5 cells -> 2 rows, 3 cols   (grid of 3 + 2, rather than 4 + 1)
            7 cells -> 2 rows, 4 cols   (grid of 4 + 3)

        :param ncells:      The number of graphic elements intended to being assembled into a subplot grid
        :param max_n_cols:
        :return:            A pair with the recommended number of rows and columns to use
        """

        assert max_n_cols >= 1, \
            "optimal_subplot_grid(): the number of maximum columns must be at least 1" \

        assert ncells >= 1, \
            "optimal_subplot_grid(): the number of cells must be at least 1" \

        ncols = min(ncells, max_n_cols)     # Tentative number of columns to use (no more than ncells nor max_n_cols)
        nrows = math.ceil(ncells / ncols)   # Number of rows needed for the above number of columns
        ncols_revised = math.ceil(ncells / nrows)   # Now that the number of rows has been fixed, attempt to use the least number of columns
        #print(f"{ncols} columns [Revised to {ncols_revised}], over {nrows} rows")

        return nrows, ncols_revised




    #####################################################################################################

    '''                                    ~   HEATMAPS   ~                                           '''

    def ________HEATMAPS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    @classmethod
    def heatmap_stack_1D(cls, data_matrix, labels :[str],
                         title="", data_name="value", entity_name="SAMPLE",
                         colors=None, color_borders=False, monochromatic=False, text_format=None,
                         height=None, grid_vert_spacing=None, colorbar_extend=1.2,
                         xgap=None, ygap=2,
                         barriers=None, barrier_color="brown", show=False) -> pgo.Figure:
        """
        Prepare and return a Plotly Figure object containing a series of vertically-stacked 1-D heatmaps,
        each with its own colorbar, with data from a separate variable.
        All variables share the same number of bins.
        Optional barriers can be places between bins, consistently across all the heatmaps.

        :param data_matrix:     NumPy 2-D array of floats, of dimension: (n_variables) x (n_bins)
                                    Each row represents a different variable, and contains its bin data
        :param labels:          List of names for the variables, in the same order as they appear in argument `data_matrix`.
                                    The number of labels must match the number of rows of `data_matrix`.
                                    Labels will be shown to the left of their respective heatmaps,
                                    and will also appear when the mouse hovers over heatmap cells
                                    EXAMPLE: ['A', 'B']
        :param title:           [OPTIONAL] Caption to place at the top of the overall plot

        :param data_name:       [OPTIONAL] The name of the numerical quantity to visualize as heatmaps.
                                    It will appear in the colorbars, and in the hover boxes.
                                    EXAMPLES: "value" (default), "Concentration", "Temperature", "Pressure"
        :param entity_name:     [OPTIONAL] To specify what data_name is referring to.
                                    This will appear in the hover boxes.
                                    EXAMPLES: "SAMPLE" (default), "Chemical", "Gas"

        :param height:          [OPTIONAL] Height, in pixels, for the overall plot.  If the requested value is considered
                                    too small to fit the plot, it will be automatically increased
                                    (and a message to this effect will be displayed)

        :param grid_vert_spacing: [OPTIONAL] Only applicable if more than 1 heatmap.
                                    Fractional (relative to overall plot height of 1.0) spacing between grid rows;
                                    increase to make the overall plot more spacious vertically.
                                    Note that grid_vert_spacing cannot be greater than (1 / (rows - 1))

        :param colorbar_extend: [OPTIONAL] Multiplicative factor to expand or shrink the colorbar by (Default 1.2)

        :param xgap:            [OPTIONAL] Horizontal gap (in pixels) between adjacent bins

        :param colors:          [OPTIONAL] List of standard color names : one for each of the variables being included, in order.
                                    Colors, if specified, are used to draw a segment to the left of each heatmap,
                                    as well as (if `color_borders` is True) to draw borders around the bins
                                    EXAMPLE: ['turquoise', 'red']
        :param color_borders:   [OPTIONAL] Only applicable if `colors` is set.
                                    Boolean indicating whether to draw colored borders around the bins.  (Default, False)
        :param monochromatic:     [OPTIONAL] Only applicable if `colors` is set.
                                     Boolean (default False) indicating whether all the heatmap values should use a gray scale
                                     (smallest value = black; highest = white) rather than have each row use
                                     a scale based on its corresponding value in argument `colors`

        :param text_format:     [OPTIONAL] String with standard python formatting for numbers,
                                    affecting the values shown on cells.
                                    If not specified, the level of detail is automatically set based on the number of bins.
                                    EXAMPLES: ".3g", ".2f", ".0f"

        :param barriers:        [OPTIONAL] List of integers, in any order,
                                    specifying the positions of vertical barriers between heatmap bins.
                                    Use the bin number of the bin to the RIGHT of the desired barrier
                                    (or 1 + the last bin number if at the far right)
                                    The values must be between 0 and n_bins, both inclusive, where n_bins is the number of bins.
                                    EXAMPLE: [0, 12, 3, 15]
        :param barrier_color:   [OPTIONAL] Only applicable if `barriers` is set.
                                    Standard color name to use for the barriers (default, 'brown')

        :param show:            If True, the plot will be shown
                                    Note: on JupyterLab, simply returning a plot object (without assigning it to a variable)
                                          gets it automatically shown

        :return:                A Plotly "Figure" object containing a stack of Heatmaps
        """
        assert type(labels) == list, \
            f"The argument `labels` must be a list.  The passed argument was of type {type(labels)}"


        if barriers is not None:
            assert type(barriers) == list, "The argument `barriers`, if specified, must be a list of integers"

        nrows, n_bins = data_matrix.shape
        #print("nrows: ", nrows)
        #print("n_bins: ", n_bins)

        assert len(labels) == nrows, \
            f"The length of the list passed to argument `labels` must be {nrows} (same as number of rows in `data_matrix`)"


        # Enforce a minimum plot height based on the number of rows (if too small, besides looking terrible,
        # if can produce an obscure Plotly error)
        requested_height = height

        MIN_ROW_HEIGHT = 250

        min_plot_height = max(MIN_ROW_HEIGHT, (nrows+1)*100)

        if height is None:
            height = 1.5 * min_plot_height       # If not specified by user,
                                                 # use a somewhat bigger value than the minimum height
        else:
            height = max(height, min_plot_height)

        if (requested_height is not None) and (height > requested_height):
            print(f"Height increased to {height} because the requested value of {requested_height} is insufficient to fit the plot")


        if (nrows > 1) and (grid_vert_spacing is None):
            grid_vert_spacing = .135 / (nrows - 1)              # Empirical choice for "compact" but not "squished vertically"
            #print("grid_vert_spacing: ", grid_vert_spacing)

        if nrows == 1:
            row_height = height
        else:
            row_height = height * ( 1 - grid_vert_spacing * (nrows - 1) ) / nrows

        #print("Height of each graphic row: ", row_height)


        delta_border = 0.5          # This value affects the vertical extension
                                    # of the color-coded edges to add to the heatmap cells
        if colors is not None:
            # For smaller plots, it esthetically looks better to make these boxes a tad smaller
            if row_height <= 200:
                delta_border = 0.485
            elif row_height <= 300:
                delta_border = 0.49
            elif row_height <= 750:
                delta_border = 0.495


        delta_barrier = 0.5         # This value affects the vertical extension
                                    # of the color-coded edges to add to the heatmap cells
        if barriers is not None:
            # For smaller plots, it esthetically looks better to make these edges a tad smaller
            if row_height <= 200:
                delta_barrier = 0.48
            elif row_height <= 250:
                delta_barrier = 0.485
            elif row_height <= 750:
                delta_barrier = 0.495


        # Prepare the blank grid for the subplots
        # Note that the first subplot (row=1, col=1) appears at the top of the stack
        # No titles are used for the subplots, in the interest of compactness, and to avoid intrusion into other subplots
        fig = sp.make_subplots(rows=nrows, cols=1,
                               horizontal_spacing=0, vertical_spacing=grid_vert_spacing)

        # Show less cell text if there's a large number of bins
        if text_format is not None:
            texttemplate = f"%{{z:{text_format}}}"
        else:
            if n_bins <= 5:
                texttemplate = "%{z:.4g}"
            elif n_bins <= 10:
                texttemplate = "%{z:.3g}"
            elif n_bins <= 15:
                texttemplate = "%{z:.2g}"
            elif n_bins <= 20:
                texttemplate = "%{z:.1g}"
            else:
                texttemplate = None

        if xgap is None:
            if n_bins > 50:
                xgap = 0    # If there are lots of bins, the heatmap generally looks better without gaps
            else:
                xgap = 2

        for i, hm_label in enumerate(labels):   # For each row of stacked heatmaps (starting at the top)
            row = i + 1   # row values start at 1

            # Get the subplot domain in figure coordinates
            layout_XAxis, layout_YAxis = fig.get_subplot(row=row, col=1)
            #x_domain = layout_XAxis.domain   # EXAMPLE: [0.55, 1.0]
            y_domain = layout_YAxis.domain    # EXAMPLE: [0.0, 0.425]
            domain_height = y_domain[1] - y_domain[0]      # vertical height of subplot row (as a fraction of all plot)
            #print("domain_height: ", domain_height)       # Should be same as (row_height / height)
            #print("domain_height double check: ", row_height / height)


            # Compute colorbar position to the right of the subplot
            cb_y = (y_domain[0] + y_domain[1]) / 2      # Center vertically
            cb_len = domain_height * colorbar_extend    # Empirically, a value of about 1.2 for colorbar_extend
                                                        #   makes the colorbars tall, but unlikely to run into those of other rows


            # Prepare the plotly Heatmap object
            data_row = data_matrix[i]          # A 1-D array
            hm_matrix = np.array([data_row])   # Turn into a 1-row matrix

            hovertemplate = f"{data_name}: %{{z}}<br>Bin #: %{{x}}<br>{entity_name}: %{{y}}<extra>{hm_label}</extra>"

            if colors is not None:
                row_color = colors[i]   # Color assigned to this row of stacked heatmaps
            else:
                row_color = None

            if monochromatic or (row_color is None):
                color_scale = "gray_r"
            else:
                lighter_color = "white"
                color_scale = [
                    [0.0, lighter_color],   # Lighter tint (possibly white)
                    [1.0, row_color],       # Full color
                ]

            hm = pgo.Heatmap(z=hm_matrix,
                             y = [hm_label],
                             xgap=xgap, ygap=ygap,
                             hovertemplate=hovertemplate,
                             texttemplate = texttemplate,
                             colorscale=color_scale,
                             colorbar=dict(
                                            title=data_name,
                                            y=cb_y,
                                            len=cb_len
                                        )
                            )

            # Add the newly-created Heatmap object to the grid of subplots
            fig.add_trace(hm, row=row, col=1)   # row and col values start at 1

            if row == nrows:
                # Add x-axis title for this subplot (identified by row and col numbers)
                fig.update_xaxes(
                                    title_text="Bin number",
                                    row=row,
                                    col=1
                                )


            if row_color is not None:
                # Add color elements to this heatmap row, using plotly "shapes" elements

                # Add annotations for row color-coding, as a vertical bar to the left of the individual heatmap
                fig.add_shape(
                    type="rect",
                    x0 = -0.7,      x1 = -0.6,       # Slightly outside the heatmap
                    #y0 = i - delta_border, y1 = i + delta_border ,
                    y0 = - delta_border, y1 = delta_border ,
                    xref="x", yref="y",
                    line=dict(width=0),
                    fillcolor=row_color,
                    row=row, col=1
                )

                if color_borders:
                    # Add color-coded edges around all heatmap cells
                    for j in range(n_bins):
                        fig.add_shape(
                            type="rect",
                            x0=j - 0.5,   x1=j + 0.5,                  # Reduce x0 to extend to left; increase x1 to extend to right
                            #y0=i - delta_border, y1=i + delta_border,  # Reduce y0 to extend to bottom; increase y1 to extend to top
                            y0 = - delta_border, y1 = delta_border,
                            xref="x", yref="y",
                            line=dict(color=row_color, width=3),        # Edge color and width (note there's NO fill color)
                            row=row, col=1
                        )
            # END "if row_color is not None"


            if barriers is not None:
                # Add barriers (if applicable): red separator bars to the right of the specified cells, using plotly "shapes"
                # Barriers are painted last, to appear at the top
                for separator_bin in barriers:   # separator_bin is the index of the bin to the RIGHT of the barrier
                    assert 0 <= separator_bin <= n_bins, \
                        f"All elements of the list passed to argument `barriers` must be integers between 0 and {n_bins}, inclusive.  " \
                        f"{separator_bin} falls outside of this range"

                    m = separator_bin-1
                    fig.add_shape(
                            type="rect",
                            x0 = m + 0.43,  x1 = m + 0.57,     # Reduce x0 to extend bar's left edge to left;
                                                               #    increase x1 to extend bar's right edge to right
                            y0 = -delta_barrier, y1 = delta_barrier,    # Reduce y0 to extend to bottom; increase y1 to extend to top
                            xref="x", yref="y",
                            line = dict(width=0),
                            fillcolor = barrier_color,
                            row=row, col=1
                        )

            # END of addition of individual heatmaps to the grid of subplots


        fig.update_layout(title=title,
                          height=height)

        if show:
            fig.show()  # Display the plot

        return fig



    @classmethod
    def heatmap_grid(cls, array_list :list, labels :[str], title ="Grid of Heatmaps",
                     height=None, colors=None, z_name="z", max_n_cols=4, cartesian=True) -> pgo.Figure:
        """
        Prepare and return a Plotly Figure object containing a grid of heatmaps (up to a max of 12)

        :param array_list:  List of 2-D Numpy arrays, up to a maximum of 12;
                                each element contains the data for one of the heatmaps
        :param labels:      List of labels for each of the heatmaps; its length must match that of the data
        :param title:       [OPTIONAL] Overall title to show at the top of the grid of heatmaps
        :param height:      [OPTIONAL] Height of the overall grid of heatmaps
        :param colors:      [OPTIONAL] List of CSS color names for each of the heatmaps.
                                If provided, its length must match that of the data;
                                if not provided, default colors are used
        :param z_name:      [OPTIONAL] Name of the quantity being visualized in the heatmaps; e.g. "Conc." or "Temp."
        :param max_n_cols:  [OPTIONAL] The maximum number of columns to use in the grip (at most 4)
        :param cartesian:   If True (default) a Cartesian grid coordinate is used, with y-bin numbers increasing up
        :return:            A Plotly "Figure" object
        """

        assert max_n_cols <= 4, "heatmap_grid(): Can only handle up to 4 columns"

        n_cells=len(array_list)
        assert n_cells <= 12, "Can only handle up to 12 heatmaps"
        assert n_cells == len(labels), "The number of heatmaps and labels must match"

        if colors is None:
            colors = Colors.assign_default_heatmap_colors(n_cells)
        else:
            assert n_cells == len(colors), "heatmap_grid(): The number of labels and colors must match"


        # Create subplots for all the individual heatmaps

        nrows, ncols = cls._optimal_subplot_grid_size(ncells=n_cells, max_n_cols = max_n_cols) # Number of rows/columns in subplot grid

        grid_hor_spacing = 0.15     # TODO: would probably be good to increase this a tad
        grid_vert_spacing = 0.15 if nrows < 3 else 0.12  # A little more spacious vertically if few rows

        # Pick esthetically-pleasing values for x0 and x1, used for horizontal positioning of the subplots
        if ncols == 1:
            x0=1.0
            x1=0    # N/A in this case
        elif ncols == 2:
            x0=0.43
            x1=0.5745
        elif ncols == 3:
            x0=0.235
            x1=0.3834
        elif ncols == 4:
            x0=0.14
            x1=0.2875
        else:
            raise Exception(f"Cannot handle this {ncols} columns")

        # Pick esthetically-pleasing values for ln (Length of the colorbars) and
        # y0, y1, and recommended_height (used for vertical positioning of the subplots)
        if nrows == 1:
            ln=1.2
            y0=0.505
            y1=0    # N/A in this case
            recommended_height=500
        elif nrows == 2:
            ln=0.52
            y0=0.792
            y1=0.573
            recommended_height=850
        elif nrows == 3:
            ln=0.3
            y0=0.895
            y1=0.373   # 0.383
            recommended_height=950
        else:
            raise Exception(f"Cannot handle this {nrows} rows")


        if height is None:
            height = recommended_height

        if nrows == 1:
            assert height >= 250, f"Overall height cannot be less than 250 for {n_cells} heatmaps"
        elif nrows == 2:
            assert height >= 600, f"Overall height cannot be less than 600 for {n_cells} heatmaps"
        elif nrows == 3:
            assert height >= 750, f"Overall height cannot be less than 750 for {n_cells} heatmaps"

        fig = sp.make_subplots(rows=nrows, cols=ncols, subplot_titles=[f'{c}' for c in labels],
                               horizontal_spacing=grid_hor_spacing, vertical_spacing=grid_vert_spacing)

        row = 1
        col = 1
        for i, hm_label in enumerate(labels):
            #print("************* row, col : ", row, col)
            color_name = colors[i]
            if color_name is None:
                color_scale = "gray_r"
            else:
                #lighter_color = Colors.lighten_color(color_name, factor=.96)
                lighter_color = "white"
                color_scale = [
                    [0.0, lighter_color],   # Light tint
                    [1.0, color_name],      # Full color
                ]

            # Create the Heatmap plotly object
            hovertemplate = f"{z_name}: %{{z}}<br>x bin: %{{x}}<br>y bin: %{{y}}<extra>{hm_label}</extra>"
            # EXAMPLE of hovertemplate: "Conc.: %{z}<br>x bin: %{x}<br>y bin: %{y}<extra>Enzyme X</extra>"

            hm = pgo.Heatmap(z=array_list[i],
                             xgap=2, ygap=2,
                             hovertemplate=hovertemplate,
                             texttemplate = '%{z:.2f}',
                             colorscale=color_scale,
                             colorbar=dict(
                                        title=z_name,
                                        len=ln,                     # Length of the colorbar
                                        y = y0 - (row - 1) * y1,    # Adjust vertical position in subplot
                                        x = x0 + x1 *(col-1)        # Adjust horizontal position
                                    )
                             )

            # Add the newly-created Heatmap to the grid of subplots
            fig.add_trace(hm, row=row, col=col)   # row and col values start at 1

            # Add x-axis and y-axis titles for this subplot (identified by row and col numbers)
            fig.update_xaxes(
                                title_text="x bin",
                                row=row,
                                col=col
                            )

            fig.update_yaxes(
                                title_text="y bin",
                                row=row,
                                col=col
                            )

            if not cartesian:
                # Invert the y-axis for this subplot
                fig.update_yaxes(
                                    autorange="reversed",
                                    row=row,
                                    col=col
                                )

            col += 1
            if col > ncols:
                col = 1
                row += 1

        # END of addition of individual heatmaps to the grid of subplots


        # Update layout, with the title, and height adjustment
        fig.update_layout(
            title=title,
            height=height,
            showlegend=False
        )

        return fig
