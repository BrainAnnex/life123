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
    Static class to assist in the use of the plotly library

    TODO: improve consistency in argument names; also across UniformCompartment
    TODO: rename as "VisualizationHelper" or "GraphicsHelper"
    """



    #####################################################################################################

    '''                                    ~   COLORS   ~                                           '''

    def ________COLORS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    @classmethod
    def get_default_colors(cls, n :int) -> [int]:
        """
        Return a list of n colors, specified by their standard plotly names;
        meant for situations when 1 or more default colors are needed.

        The choice of default colors is hardwired in this function.

        :param n:   Desired number of default colors
        :return:    A list of n standard (CSS) color names
        """
        # TODO: provide multiple, user-selectable, harmonious assortments of default colors

        default_colors = ['darkturquoise', 'green', 'brown', 'red', 'gray', 'blue',
                          'orange', 'purple', 'cyan', 'darkorange', 'navy',
                          'darkred', 'black', 'mediumspringgreen']
        '''
        # Available color names:
                aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, saddlebrown, salmon, sandybrown,
                seagreen, seashell, sienna, silver, skyblue,
                slateblue, slategray, slategrey, snow, springgreen,
                steelblue, tan, teal, thistle, tomato, turquoise,
                violet, wheat, white, whitesmoke, yellow,
                yellowgreen
        '''
        colors = default_colors[:n]      # Pick the first n default colors; TODO: rotate if needing more

        return colors



    @classmethod
    def get_default_heatmaps_colors(cls, n :int) -> [int]:
        """
        Return a list of n colors, specified by their standard plotly names;
        meant for situations when 1 or more default colors are needed.

        The choice of default colors is hardwired in this function.

        :param n:   Desired number of default colors
        :return:    A list of n standard (CSS) color names
        """
        default_colors = ["yellow", "green", "blue", "red",
                          "purple", "teal", "black", "brown",
                          "deeppink", "midnightblue", "darkolivegreen", "darkorange"]

        colors = default_colors[:n]      # Pick the first n default colors; TODO: rotate if needing more

        return colors






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
            colors = cls.get_default_colors(number_of_curves)
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
                                    if None, then display all columns except the one that was declared as the independent variable
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
                                    If the number of vertical line is so large as to overwhelm the plot,
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
            colors = PlotlyHelper.get_default_colors(number_of_curves)
        elif type(colors) == str:
            # Turn colors into a list, if it was a single entry
            colors = [colors]
        else:
            # If we get here, we were given a list; replace any missing (None) entry with a default color
            number_none = colors.count(None)    # Number of None entries
            if number_none > 0:
                replacement_colors = PlotlyHelper.get_default_colors(number_none)   # Get all the replacements in bulk
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
            fig.show()  # Actually display the plot

        return fig



    @classmethod
    def combine_plots(cls, fig_list :Union[list, tuple], title="", x_label=None, y_label=None,
                      xrange=None, legend_title=None, curve_labels=None, show=False) -> pgo.Figure:
        """
        Combine together several existing plotly plots into a single one (with combined axes)

        EXAMPLE:
                    from life123 import PlotlyHelper
                    p1 = PlotlyHelper.plot_pandas(various args, show=False)
                    p2 = PlotlyHelper.plot_pandas(various args, show=False)
                    PlotlyHelper.combine_plots([p1, p2], other optional args)

        :param fig_list:    List or tuple of plotly "Figure" objects (as returned by several functions)
        :param title:       [OPTIONAL] The title to use for the overall plot
        :param x_label:     [OPTIONAL] Caption to use for the x-axis; if not specified, use that of the 1st plot
        :param y_label:     [OPTIONAL] Caption to use for the y-axis; if not specified, use that of the 1st plot
        :param xrange:      [OPTIONAL] list of the form [t_start, t_end], to initially only show a part of the timeline.
                                Note: it's still possible to zoom out, and see the excluded portion
        :param legend_title:[OPTIONAL] String to show at the top of the legend box
        :param curve_labels:[OPTIONAL] List of labels to use for the various curves in the legend
                                and in the hover boxes; if not specified, use the titles of the individual plots
        :param show:        If True, the plot will be shown
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
    def heatmap_grid(cls, array_list :list, labels :[str], title ="Grid of Heatmaps",
                     height=None, colors=None, z_name="z", max_n_cols=4, cartesian=True) -> pgo.Figure:
        """
        Prepare and return a Plotly Figure object containing a grid of heatmaps (up to a max of 12)

        :param array_list:  List of 2-D Numpy arrays, up to a maximum of 12,
                                containing the data for each of the heatmaps
        :param labels:      List of labels for each of the heatmaps; its length must match that of the data
        :param title:       [OPTIONAL] Overall title to show at the top of the grid of heatmaps
        :param height:      [OPTIONAL] Height of the overall grid of heatmaps
        :param colors:      [OPTIONAL] List of CSS color names for each of the heatmaps.
                                If provided, its length must match that of the data; otherwise, default colors are used
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
            colors = cls.get_default_heatmaps_colors(n_cells)
        else:
            assert n_cells == len(colors), "The number of labels and colors must match"


        # Create subplots for all the individual heatmaps

        nrows, ncols = cls._optimal_subplot_grid_size(ncells=n_cells, max_n_cols = max_n_cols) # Number of rows/columns in subplot grid

        grid_hor_spacing = 0.15     # TODO: would probably be good to increase this a tad
        grid_vert_spacing = 0.15 if nrows < 3 else 0.12  # A little more spacious vertically if few rows

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

            # Create the Heatmap object
            hovertemplate=f"{z_name}: %{{z}}<br>x bin: %{{x}}<br>y bin: %{{y}}<extra>{hm_label}</extra>"
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

            fig.add_trace(hm, row=row, col=col)   # row and col values start at 1

            # Add x-axis and y-axis titles for this subplot
            fig.update_xaxes(
                                title_text=f"x bin",
                                row=row,
                                col=col
                            )

            fig.update_yaxes(
                                title_text=f"y bin",
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


        # Update layout
        fig.update_layout(
            title=title,
            height=height,
            showlegend=False
        )

        return fig
