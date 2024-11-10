import plotly.express as px
import plotly.graph_objects as pgo
import numpy as np
import pandas as pd
from typing import Union



class PlotlyHelper:
    """
    Static class to assist in the use of the plotly library

    TODO: improve consistency in argument names; also across UniformCompartment
    """

    @classmethod
    def get_default_colors(cls, n :int) -> [int]:
        """
        Return a list of n colors, specified by their standard plotly names;
        meant for situations when 1 or more default colors are needed.

        The choice of default colors is hardwired in this function.

        :param n:   Desired number of default colors
        :return:
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
        colors = default_colors[:n]      # Pick the first default colors; TODO: rotate if needing more

        return colors



    @classmethod
    def plot_curves(cls, x, y, title="", range_x=None, x_label="", y_label="", curve_labels=None, legend_title=None,
                    colors=None, show=False) -> pgo.Figure:
        """
        Plot one or more 2D curves

        :param x:           A Numpy array, with the (common) x-axis values
        :param y:           Either a Numpy array, or a list/tuple of them, with the y-axis values of the curve(s)
        :param title:       [OPTIONAL] Title to use for the overall plot
        :param range_x:      [OPTIONAL] list of the form [t_start, t_end], to initially only show a part of the timeline.
                                Note: it's still possible to zoom out, and see the excluded portion
        :param x_label:      [OPTIONAL] Caption to use for the x-axis
        :param y_label:      [OPTIONAL] Caption to use for the y-axis
        :param curve_labels:[OPTIONAL] String, or list of strings.
                                Label(s) to use for the various curves in the legend and in the hover boxes.
                                If missing, and there's only 1 curve, the legend box won't be shown
        :param legend_title:[OPTIONAL] String to show at the top of the legend box.
                                Ignored if curve_labels wasn't set.
                                If not provided, and the legend box is shown, it will appear as "variable"
        :param colors:      [OPTIONAL] Either a single color (string with standard plotly name, such as "red"),
                                or list of names to use, in order; if None, then use the hardwired defaults
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
                    vertical_lines_to_add=None,
                    show_intervals=False, show=False) -> pgo.Figure:
        """
        Using plotly, draw line plots from the values in the given dataframe.

        Note: if this plot is to be later combined with others, use PlotlyHelper.combine_plots()

        :param df:              Pandas dataframe with the data for the plot
        :param x_var:           Name of column with independent variable for the x-axis
        :param fields:          Name, or list of names, of the dataframe columns whose values are to be plotted;
                                    if a list is passed, also display a figure legend;
                                    if None, then display all columns except the one that is declared as the independent variable
        :param colors:          (OPTIONAL) Either a single color (string with standard plotly name, such as "red"),
                                    or list of names to use, in order; if None, then the hardwired default colors are used
        :param title:           (OPTIONAL) Title for the plot
        :param title_prefix:    (OPTIONAL) Strint to prefixed (automatically followed by " <br>") to the title
        :param range_x:         (OPTIONAL) list of the form [t_start, t_end], to initially show only a part of the timeline.
                                    Note: it's still possible to zoom out, and see the excluded portion
        :param range_y:         (OPTIONAL) list of the form [y_min, y_max], to initially show only a part of the y values.
                                    Note: it's still possible to zoom out, and see the excluded portion
        :param x_label:         (OPTIONAL) Caption to use for the x-axis
        :param y_label:         (OPTIONAL) Caption to use for the y-axis.  Default: "Y"
        :param legend_header:   (OPTIONAL) Caption to use at the top of the legend box.
                                            Only applicable if more than 1 curve is being shown.
        :param vertical_lines_to_add:  (OPTIONAL) Ignored if the argument `show_intervals` is specified.
                                    Value, or list, or tuple, or Numpy array, or Pandas series,
                                    of x-coordinate(s) at which to draw thin vertical dotted gray lines.
                                    If the number of vertical line is so large as to overwhelm the plot,
                                    only a sample of them is shown.
                                    Note that vertical lines, if requested, go into the plot's "layout";
                                    as a result they might not appear if this plot is later combined with another one.
        :param show_intervals:  (OPTIONAL) If True, it over-rides any value passed to the `vertical_lines` argument,
                                    and draws thin vertical dotted gray lines at all the x-coords
                                    of the data points in the saved history data;
                                    also, it adds a comment to the title.
        :param show:            If True, the plot will be shown
                                    Note: on JupyterLab, simply returning a plot object (without assigning it to a variable)
                                          leads to it being automatically shown

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
            colors = PlotlyHelper.get_default_colors(number_of_curves)
        elif type(colors) == str:
            colors = [colors]


        if title_prefix is not None:
            title = f"{title_prefix} <br>{title}"

        if show_intervals:
            vertical_lines_to_add = df[x_var]   # Make use of the simulation times
            title += " (time steps shown in dashed lines)"


        # Create the main plot
        fig = px.line(data_frame=df, x=x_var, y=fields,
                      title=title, range_x=range_x, range_y=range_y,
                      color_discrete_sequence = colors,
                      labels={"value": y_label, "variable": legend_header})

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
        Combine together several existing plotly plots

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
