import plotly.express as px
import plotly.graph_objects as pgo
from typing import Union



class PlotlyHelper:
    """
    To assist in the use of the plotly library
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

        default_colors = ['blue', 'green', 'brown', 'red', 'gray',
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
    def plot_curves(cls, x, y, title="", xrange=None, xlabel="", ylabel="", curve_labels=None, legend_title=None,
                    colors=None, show=False) -> pgo.Figure:
        """
        Plot one or more 2D curves

        :param x:           A Numpy array, with the (common) x-axis values
        :param y:           Either a Numpy array, or a list/tuple of them, with the y-axis values of the curve(s)
        :param title:       (OPTIONAL) Title to use for the overall plot
        :param xrange:      (OPTIONAL) list of the form [t_start, t_end], to initially only show a part of the timeline.
                                Note: it's still possible to zoom out, and see the excluded portion
        :param xlabel:      (OPTIONAL) Caption to use for the x-axis
        :param ylabel:      (OPTIONAL) Caption to use for the y-axis
        :param curve_labels:(OPTIONAL) String, or list of strings.
                                Label(s) to use for the various curves in the legend and in the hover boxes.
                                If missing, and there's only 1 curve, the legend box won't be shown
        :param legend_title:(OPTIONAL) String to show at the top of the legend box.
                                Ignored if curve_labels wasn't set.
                                If not provided, and the legend box is shown, it will appear as "variable"
        :param colors:      (OPTIONAL) Either a single color (string with standard plotly name, such as "red"),
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

        if not xlabel:
            xlabel = "x"

        if not ylabel:
            ylabel = "y"

        if number_of_curves == 1:
            hovertemplate = f"{xlabel}=%{{x}}<br>{ylabel}=%{{y}}<extra></extra>"
        else:
            hovertemplate = f"{xlabel}=%{{x}}<br>value=%{{y}}<extra></extra>"
            # The "<extra></extra>" appears to suppress a redundant name showing in the hoverbox

        for index in range(number_of_curves):
            if curve_labels is not None:
                label = curve_labels[index]
                fig.data[index]["name"] = label
                fig.data[index]["hovertemplate"] = f"{label} :<br>{hovertemplate}"    # variable :
            else:
                fig.data[index]["hovertemplate"] = hovertemplate


        fig.update_layout(title=title,
                          xaxis_title=xlabel,
                          yaxis_title=ylabel)

        if legend_title:
            fig.update_layout(legend={"title": legend_title})

        if xrange:
            fig.update_layout(xaxis={"range": xrange})


        if show:
            fig.show()  # Actually display the plot

        return fig



    @classmethod
    def combine_plots(cls, fig_list :Union[list, tuple], title="", xlabel=None, ylabel=None,
                      xrange=None, curve_labels=None, legend_title=None, show=False) -> pgo.Figure:
        """
        Combine together several existing plotly plots

        :param fig_list:    List or tuple of plotly "Figure" objects (as returned by several functions)
        :param title:       (OPTIONAL) The title to use for the overall plot
        :param xlabel:      (OPTIONAL) Caption to use for the x-axis; if not specified, use that of the 1st plot
        :param ylabel:      (OPTIONAL) Caption to use for the y-axis; if not specified, use that of the 1st plot
        :param xrange:      (OPTIONAL) list of the form [t_start, t_end], to initially only show a part of the timeline.
                                Note: it's still possible to zoom out, and see the excluded portion
        :param curve_labels:(OPTIONAL) List of labels to use for the various curves in the legend
                                and in the hover boxes; if not specified, use the titles of the individual plots
        :param legend_title:(OPTIONAL) String to show at the top of the legend box
        :param show:        If True, the plot will be shown
                                Note: on JupyterLab, simply returning a plot object (without assigning it to a variable)
                                      leads to it being automatically shown
        :return:            A plotly "Figure" object for the combined plot
        """
        assert (type(fig_list) == list) or (type(fig_list) == tuple), \
            "combine_plots(): the argument fig_list must be a list or tuple"


        representative_fig = fig_list[0]    # The axes titles of the first plot are used as default values
        if xlabel is None:
            xlabel = representative_fig.layout.xaxis.title.text
        if ylabel is None:
            ylabel = representative_fig.layout.yaxis.title.text


        # Put together the data from all the various individual plots
        combined_data = []
        for fig in fig_list:
            combined_data += fig.data      # concatenating lists


        all_fig = pgo.Figure(data=combined_data)    # Note that the + is concatenating lists

        all_fig.update_layout(title=title,
                              xaxis_title=xlabel,
                              yaxis_title=ylabel)

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
