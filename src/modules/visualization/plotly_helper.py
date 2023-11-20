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

        :param n:   Desire number of default colors
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
    def plot_curves(cls, x, y, title="", xaxis="", yaxis="", curve_labels=None, legend_title=None, colors=None, suppress=False) \
                            -> Union[None, pgo.Figure]:
        """

        :param x:           A Numpy array
        :param y:           Either a Numpy array, or a list/tuple of them
        :param title:       (OPTIONAL) Title to use for the overall plot
        :param xaxis:       (OPTIONAL) Caption to use for the x-axis
        :param yaxis:       (OPTIONAL) Caption to use for the y-axis
        :param curve_labels:(OPTIONAL) Ignored if just 1 curve.  List of labels to use for the various curves in the legend
                                and in the hover boxes
        :param legend_title:(OPTIONAL) Ignored if just 1 curve.  String to show at the top of the legend.
        :param colors:      (OPTIONAL) Either a single color (string with standard plotly name, such as "red"), or list of names
        :param suppress:    If True, nothing gets shown - and a plotly "Figure" object gets returned instead;
                                this is useful to make additional tweaks, or combine multiple plots, etc
        :return:            None or a plotly "Figure" object, depending on the "suppress" flag
        """
        if type(y) == list or type(y) == tuple:
            number_of_curves = len(y)
        else:
            number_of_curves = 1

        if colors is None:
            colors = cls.get_default_colors(number_of_curves)
        elif type(colors) == str:
            colors = [colors]


        fig = px.line(x=x, y=y, color_discrete_sequence=colors)     # y can be one array or a list

        if not xaxis:
            xaxis = "x"

        if not yaxis:
            yaxis = "y"

        if number_of_curves == 1:
            hovertemplate = f"{xaxis}=%{{x}}<br>{yaxis}=%{{y}}<extra></extra>"
        else:
            hovertemplate = f"{xaxis}=%{{x}}<br>value=%{{y}}<extra></extra>"
            # The "<extra></extra>" appears to suppress a redundant name showing in the hoverbox

        for index in range(number_of_curves):
            if curve_labels is not None:
                label = curve_labels[index]
                fig.data[index]["name"] = label
                fig.data[index]["hovertemplate"] = f"{label}<br>{hovertemplate}"    # variable :
            else:
                fig.data[index]["hovertemplate"] = hovertemplate
        '''
        
        if (number_of_curves > 1) and (curve_labels is not None):
            for index, l in enumerate(curve_labels):
                fig.data[index]["name"] = l
                fig.data[index]["hovertemplate"] = f"variable : {l}<br>{xaxis}=%{{x}}<br>value=%{{y}}<extra></extra>"
                # The "<extra></extra>" appears to suppress a redundant name showing in the hoverbox
        '''

        if legend_title:
            fig.update_layout(title=title,
                              xaxis_title=xaxis,
                              yaxis_title=yaxis,
                              legend={"title": legend_title})
        else:
            fig.update_layout(title=title,
                              xaxis_title=xaxis,
                              yaxis_title=yaxis)

        if suppress:
            return fig      # Return the plot data (without actually displaying the plot)
        else:
            fig.show()      # Actually display the plot
