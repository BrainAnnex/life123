class SVGhelper
/* 	VERSION 1.2

    Helper library to generate either individual SVG tags, or larger SVG constructs.

	All coordinates are screen-coordinate integers (origin at top/left; x-axis point right; y-axis pointing DOWN)

	Screen coordinates are referred to as:  (Sx, Sy)

	Styling is generally relegated to CSS.
 */
{

    /*
        SHAPES
     */

    point(Sx, Sy, r = 2, color = 'black')
    // Create a point, with the specified coordinates, radius, and color
    {
         return `<circle cx='${Sx}' cy='${Sy}' r='${r}' stroke='${color}' stroke-width='1' fill='black'/>`;
    }


    line(Sx1, Sy1, Sx2, Sy2, stroke)
    /* 	Create a line (segment) between the two specified points: (Sx1, Sy1) and (Sx2, Sy2),
        given in screen coordinates.
        The "stroke" argument (with a color name or value) is optional - but, if missing,
        a stroke color must be specified with CSS, or nothing will show up!
        EXAMPLE of CSS:  line {stroke: black;}
     */
    {
         if (stroke === undefined)  // Detect missing argument
            return `<line x1='${Sx1}' y1='${Sy1}' x2='${Sx2}' y2='${Sy2}'/>`;
         else
            return `<line x1='${Sx1}' y1='${Sy1}' x2='${Sx2}' y2='${Sy2}' stroke='${stroke}'/>`;
    }



    /*
        GROUPS
     */

    start_group(class_name)
	// The SVG opening tag for a group of graphic elements, with an optional CSS class name
	{
		if (class_name === undefined)  // Detect missing argument
			return "<g>";
		else
			return `<g class='${class_name}'>`;
	}

	end_group()
	// A SVG end tag for a group of graphic elements
	{
		return "</g>";
	}



    /*
        TEXT
     */

    text(label, Sx, Sy, attributes = "", dx = 0, dy = 0)
    /*  Create a text label, with the given attributes.
        The x-position of the text refers to its left side, unless CSS
        directives such as "text-anchor: middle"  are used.
        The y-position refers to the BOTTOM of the text.

        dx and dy specify, respectively, a rightward/downward shift proportional to the font size,
        to clear the ticks and vertically-align the label regardless of font.
     */
    {
        if (dx != 0)
            attributes += " dx='" + dx + "em'";
        if (dy != 0)
             attributes += " dy='" + dy + "em'";

        if (attributes)
            return `<text x='${Sx}' y='${Sy}' ${attributes}>${label}</text>`;
        else
            return `<text x='${Sx}' y='${Sy}'>${label}</text>`;

        return this;
    }



    /*
        COMMENT
     */

    comment(text)
    {
        return `<!-- ${text} -->`;
    }



    /*
        MISC
     */

    translate(x, y)
    /*  Return a string suitable as an SVG attribute, to indicate a translation by <x,y>
        EXAMPLE:   "translate(10, 50)"
    */
    {
        return `translate(${x}, ${y})`;
    }



    /*
        TICKS and LABELS (for Axes)
     */


	vertical_tick(Sx, Sy, extension_above, extension_below)
	/* 	Create a vertical plot tick mark (i.e. meant for a horizontal axis) at the point (Sx,Sy).
		extension_above and extension_below are, respectively the lengths of the tick portions
		shown above and below the axis
	 */
	{
		const Sy_tick_min = Sy - extension_above;	// Highest point in the vertical segment forming the tack
		const Sy_tick_max = Sy + extension_below;	// Lowest point in the vertical segment forming the tack

        return this.line(Sx, Sy_tick_min, Sx, Sy_tick_max);
	}

	horizontal_tick(Sx, Sy, extension_left, extension_right)
	/* 	Create a horizontal plot tick mark (i.e. meant for a vertical axis) at the point (Sx,Sy).
		extension_left and extension_right are, respectively the lengths of the tick portions
		shown to the left and to the right of the axis
	 */
	{
		const Sx_tick_min = Sx - extension_left;	// Leftmost point in the vertical segment forming the tack
		const Sx_tick_max = Sx + extension_right;	// Rightmost point in the vertical segment forming the tack

        return this.line(Sx_tick_min, Sy, Sx_tick_max, Sy);
	}



    /*
        AXES
     */

    axis_left_scaleLinear( {    y_scale_func,
                                n_intervals,
                                y_min,
                                y_max,
                                Sx_axis,
                                tick_right=0,
                                tick_left=6  } )
    /*  Create a vertical axis line meant to be placed to the left of a plot that used d3.scaleLinear() for the y-axis.

        The y-axis is divided up into n_intervals equal parts,
        resulting in n_intervals+1 tick marks.
        The y-coordinates (in graph coordinates) are used as tick labels.

        Use CSS for styling.  The following classes are created:
                            'axis-line', 'ticks', 'tick-labels-left-axis'

        EXAMPLE of CSS:  g.tick-labels-left-axis { translate(0, 10px); }  to shift down the tick labels
     */
    {
        const bin_width = (y_max - y_min) / n_intervals;
        //console.log("bin_width: ", bin_width);

        // Determine the screen coordinates corresponding to, respectively,
        // the bottom and top point on the y-axis
        // NOTE: in SVG screen coordinates, which start at the top,
        //       the top has a SMALLER value than the bottom
        const Sy_bottom = y_scale_func(y_min);
        const Sy_top = y_scale_func(y_max);
        //console.log("Sy_bottom: ", Sy_bottom);
        //console.log("Sy_top: ", Sy_top);

        let svg = "";

        // Vertical axis line (an SVG group consisting of 1 line element)
        svg += this.start_group("axis-line");
        svg += this.line(Sx_axis, Sy_bottom, Sx_axis, Sy_top);
        svg += this.end_group();


        /* Handle the TICKS (kept together as an SVG group)
         */
        svg += this.start_group("ticks");	            // Pass a class used for all the ticks

        for (let item_index = 0;  item_index <= n_intervals;  ++item_index)  {
            // Determine the vertical position of the tick
            let y_graph = y_min + item_index * bin_width;
            let Sy_tick = y_scale_func(y_graph);
            //console.log("Sy_tick: ", Sy_tick);
            svg += this.horizontal_tick(Sx_axis, Sy_tick, tick_left, tick_right);
        }

        svg += this.end_group();                        // end of the ticks



        /* Handle the TICK LABELS (kept together as an SVG group)
         */
        svg += this.start_group("tick-labels-left-axis");         // Pass a class used for all the labels

        for (let item_index = 0;  item_index <= n_intervals;  ++item_index)  {
            let y_graph = y_min + item_index * bin_width;
            let label = y_graph.toString();    // The y-coordinates (in graph coordinates) are used as tick labels
            // First, determine the vertical position of the tick
            let Sy_tick = y_scale_func(y_graph);
            // Now, position the label
            svg += this.text(label, Sx_axis - tick_left, Sy_tick, "", -0.75 -0.25 * label.length, 0.3);
            // Note: the last 2 arguments specify, respectively, a rightward/downward shift proportional to the font size,
            //       to clear the ticks and vertically-align the label regardless of font.
            //       Extra control can be achieved with CSS
        }

        svg += this.end_group();                                // end of the labels


        return svg;

    } // axis_left_scaleLinear



    /*  TODO: maybe rename axis_left_scaleBand() */
    axis_left( {    y_scale_func,
                    Sx_axis,
                    categorical_labels,
                    tick_right=0,
                    tick_left=6  } )
    /*  Create a vertical axis line meant to be placed to the left of a plot that used d3.scaleBand() for the y-axis.

        Use CSS for styling.  The following classes are created:
                            'axis-line', 'ticks', 'tick-labels-left-axis'

        EXAMPLE of CSS:  g.tick-labels-left-axis { translate(0, 10px); }  to shift down the tick labels

        y_scale_func:       function produced with d3.scaleBand()
        Sx_axis:            x-coord of axis, in screen coordinates
        categorical_labels  List of desired label names (equally-spaced, at the center of their intervals)
        tick_right:         Amount by which ticks stick to the right of the axis, in screen coordinates
        tick_left:          Amount by which ticks stick to the left of the axis, in screen coordinates
     */
    {
        const n_items = categorical_labels.length;

        const bin_width = y_scale_func.bandwidth();     // TODO: move to caller function

        const Sy1 = y_scale_func(categorical_labels[0]);
        const Sy2 = y_scale_func(categorical_labels[n_items-1]);

        const Symin = Math.min(Sy1, Sy2);
        const Symax = Math.max(Sy1, Sy2) + bin_width;

        let svg = "";

        // Vertical axis line
        svg += this.start_group("axis-line");
        svg += this.line(Sx_axis, Symin, Sx_axis, Symax);
        svg += this.end_group();


        /* Handle the TICKS (kept together as an SVG group)
         */
        svg += this.start_group("ticks");	            // Pass a class used for all the ticks

        for (let item_index = 0;  item_index < n_items;  ++item_index)  {
            // Determine the vertical position of the tick
            let Sy_tick = y_scale_func(categorical_labels[item_index]) + 0.5 * bin_width;
            svg += this.horizontal_tick(Sx_axis, Sy_tick, tick_left, tick_right);
        }

        svg += this.end_group();                        // end of the ticks


        /* Handle the TICK LABELS (kept together as an SVG group)
         */
        svg += this.start_group("tick-labels-left-axis");         // Pass a class used for all the labels

        for (let item_index = 0;  item_index < n_items;  ++item_index)  {
            let label = categorical_labels[item_index];
            // First, determine the vertical position of the tick
            let Sy_tick = y_scale_func(categorical_labels[item_index]) + 0.5 * bin_width;
            // Now, position the label
            svg += this.text(label, Sx_axis - tick_left, Sy_tick, "", -1 -0.2 * label.length, 0.3);
            // Note: the last 2 arguments specify, respectively, a rightward/downward shift proportional to the font size,
            //       to clear the ticks and vertically-align the label regardless of font.
            //       Extra control can be achieved with CSS
        }

        svg += this.end_group();                                // end of the labels

        return svg;

    } // axis_left



    /*  TODO: maybe rename axis_bottom_scaleBand() */
    axis_bottom( {  x_scale_func,
                    Sy_axis,
                    categorical_labels,
                    tick_above=0,
                    tick_below=6  } )
    /*  Create a horizontal axis line meant to be placed below a plot that used d3.scaleBand() for the x-axis.

        Use CSS for styling.  The following classes are created:
                            'axis-line', 'ticks', 'tick-labels'

        EXAMPLE:  g.tick-labels { translate(0, 10px); }  to shift down the tick labels

        x_scale_func:       function produced with d3.scaleBand()
        Sy_axis:            y-coord of axis, in screen coordinates
        categorical_labels  List of desired label names (equally-spaced, at the center of their intervals)
        tick_above:         Amount by which ticks stick above axis, in screen coordinates
        tick_below:         Amount by which ticks stick below axis, in screen coordinates
     */
    {
        const n_items = categorical_labels.length;

        const bin_width = x_scale_func.bandwidth();          // TODO: move to caller function

        // Determine the screen coordinates corresponding to, respectively,
        // the first and last point on the x-axis
        const Sx1 = x_scale_func(categorical_labels[0]);
        const Sx2 = x_scale_func(categorical_labels[n_items-1]);


        return this.axis_bottom_common(
            {   x_scale_func: x_scale_func,
                Sx1: Sx1,
                Sx2: Sx2,
                n_items: n_items,
                bin_width: bin_width,
                Sy_axis: Sy_axis,
                categorical_labels: categorical_labels,
                tick_above: tick_above,
                tick_below: tick_below
            });

    } // axis_bottom



    axis_bottom_scaleLinear_bins( {  x_scale_func,
                                n_items,
                                Sy_axis,
                                tick_above=0,
                                tick_below=6  } )
    /*  MEANT FOR PLOTS BASED ON BINS.  THE BINS' X-VALUES ARE MID-BIN.

        Create a horizontal axis line meant to be placed below a plot that used d3.scaleLinear() for the x-axis.

        Use CSS for styling.  The following classes are created:
                            'axis-line', 'ticks', 'tick-labels'

        EXAMPLE:  g.tick-labels { translate(0, 10px); }  to shift down the tick labels

        x_scale_func:       function produced with d3.scaleLinear()
        Sy_axis:            y-coord of axis, in screen coordinates
        tick_above:         Amount by which ticks stick above axis, in screen coordinates
        tick_below:         Amount by which ticks stick below axis, in screen coordinates
     */
    {
        // Determine the screen coordinates corresponding to, respectively,
        // the leftmost point on the x-axis of the first (0-th) bin, and the rightmost point of the last (n-1) bin
        const Sx1 = x_scale_func(-0.5);             // The 0-th bin is centered at 0, and extends: [-0.5, +0.5]
        const Sx2 = x_scale_func(n_items-0.5);      // The last, (n-1)-th, bin is centered at (n-1), and extends: [n-1.5, n-0.5]

        return this.axis_bottom_shared(
            {   x_scale_func: x_scale_func,
                Sx1: Sx1,
                Sx2: Sx2,
                Sy_axis: Sy_axis,
                n_subdivisions: n_items,
                is_bins: true,
                tick_above: tick_above,
                tick_below: tick_below
            });

    } // axis_bottom_scaleLinear_bins



    axis_bottom_shared( {   x_scale_func,
                            Sx1,
                            Sx2,
                            Sy_axis,
                            n_subdivisions,
                            is_bins,
                            categorical_labels=null,
                            tick_above=0,
                            tick_below=6  } )
    /*  Create a horizontal axis line meant to be placed below a plot.      TODO: newest function, to later absorb some older ones

        Use CSS for styling.  The following classes are created:
                            'axis-line', 'ticks', 'tick-labels'

        EXAMPLE:  g.tick-labels { translate(0, 10px); }  to shift down the tick labels

        x_scale_func:       Function produced with D3; for example d3.scaleLinear() or d3.scaleBand(),
                                to map plot coordinates into screen coordinates
        Sx1:                Screen coordinates corresponding to the first point on the x-axis
        Sx2:                Screen coordinates corresponding to the last point on the x-axis

        Sy_axis:            y-coord of axis, in screen coordinates

        n_subdivisions:     The number of parts into which to divide up the interval
        is_bins:            A boolean indicating whether the subdivisions are bins (to label in their midpoints)
        categorical_labels  (OPTIONAL) List of desired label names (equally-spaced, at the center of their respective intervals);
                                typically, only used if is_bins is true

        tick_above:         Amount by which ticks stick above the axis, in screen coordinates
        tick_below:         Amount by which ticks stick below the axis, in screen coordinates
     */
    {
        const Sxmin = Math.min(Sx1, Sx2);
        const Sxmax = Math.max(Sx1, Sx2);

        let svg = "";

        // Horizontal axis line
        svg += this.start_group("axis-line");
        svg += this.line(Sxmin, Sy_axis, Sxmax, Sy_axis);
        svg += this.end_group();


        /* Handle the TICKS (kept together as an SVG group)
         */
        svg += this.start_group("ticks");	            // Pass a class used for all the ticks

        var n_ticks, label, Sx_tick;
        if (is_bins)
            n_ticks = n_subdivisions;
        else
            n_ticks = n_subdivisions + 1;

        let tick_Sx_list = [];
        let tick_label_list = [];

        for (let item_index = 0;  item_index < n_ticks;  ++item_index)  {
            if (is_bins)  {
                if (categorical_labels && categorical_labels.length)    // Use categorical labels if provided; otherwise, fall back to the index
                    label = categorical_labels[item_index];
                else
                    label = item_index;

                Sx_tick = x_scale_func(item_index);               // Bin labels are placed mid-bin
            }
            else  {
                Sx_tick = x_scale_func(item_index);
                label = Sx_tick;
            }

            svg += this.vertical_tick(Sx_tick, Sy_axis, tick_above, tick_below);

            // Save the position and label
            tick_Sx_list.push(Sx_tick);
            tick_label_list.push(label);
        }

        svg += this.end_group();                        // end of the ticks


        /* Handle the TICK LABELS (kept together as an SVG group)
         */
        svg += this.start_group("tick-labels");         // Pass a class used for all the labels

        for (let item_index = 0;  item_index < n_ticks;  ++item_index)  {
            // position the label
            svg += this.text(tick_label_list[item_index], tick_Sx_list[item_index], Sy_axis + tick_below, "", 0., 0.9);
            // Note: the last argument specifies a downward shift proportional to the font size,
            //       to clear the ticks regardless of font.  Extra control can be achieved with CSS
        }

        svg += this.end_group();                        // end of the labels

        return svg;

    } // axis_bottom_shared



    axis_bottom_scaleLinear( {  x_scale_func,
                                n_items,
                                bin_width,
                                Sy_axis,
                                tick_above=0,
                                tick_below=6  } )
    /*  Create a horizontal axis line meant to be placed below a plot that used d3.scaleLinear() for the x-axis.

        Use CSS for styling.  The following classes are created:
                            'axis-line', 'ticks', 'tick-labels'

        EXAMPLE:  g.tick-labels { translate(0, 10px); }  to shift down the tick labels

        x_scale_func:       function produced with d3.scaleLinear()
        Sy_axis:            y-coord of axis, in screen coordinates
        tick_above:         Amount by which ticks stick above axis, in screen coordinates
        tick_below:         Amount by which ticks stick below axis, in screen coordinates
     */
    {
        // Determine the screen coordinates corresponding to, respectively,
        // the first (leftmost) and last (rightmost) point on the x-axis
        const Sx1 = x_scale_func(0);
        const Sx2 = x_scale_func(n_items-1);

        return this.axis_bottom_common(
            {   x_scale_func: x_scale_func,
                Sx1: Sx1,
                Sx2: Sx2,
                n_items: n_items,
                bin_width: bin_width,
                Sy_axis: Sy_axis,
                tick_above: tick_above,
                tick_below: tick_below
            });

    } // axis_bottom_scaleLinear



    axis_bottom_common( {  x_scale_func,
                            Sx1,
                            Sx2,
                            n_items,
                            bin_width,
                            Sy_axis,
                            categorical_labels=null,
                            tick_above=0,
                            tick_below=6  } )
    /*  Create a horizontal axis line meant to be placed below a plot.      TODO: replace by the newer axis_bottom_shared()

        Use CSS for styling.  The following classes are created:
                            'axis-line', 'ticks', 'tick-labels'

        EXAMPLE:  g.tick-labels { translate(0, 10px); }  to shift down the tick labels

        x_scale_func:       Function produced with D3; for example d3.scaleLinear() or d3.scaleBand()
        Sx1:                Screen coordinates corresponding to the first point on the x-axis
        Sx2:                Screen coordinates corresponding to the last point on the x-axis
        n_items:            Number of points of interest on the x-axis
        bin_width:
        Sy_axis:            y-coord of axis, in screen coordinates
        categorical_labels  (OPTIONAL) List of desired label names (equally-spaced, at the center of their intervals)
        tick_above:         Amount by which ticks stick above axis, in screen coordinates
        tick_below:         Amount by which ticks stick below axis, in screen coordinates
     */
    {
        const Sxmin = Math.min(Sx1, Sx2);
        const Sxmax = Math.max(Sx1, Sx2) + bin_width;

        let svg = "";

        // Horizontal axis line
        svg += this.start_group("axis-line");
        svg += this.line(Sxmin, Sy_axis, Sxmax, Sy_axis);
        svg += this.end_group();


        /* Handle the TICKS (kept together as an SVG group)
         */
        svg += this.start_group("ticks");	            // Pass a class used for all the ticks

        var label;
        for (let item_index = 0;  item_index < n_items;  ++item_index)  {
            if (categorical_labels && categorical_labels.length)
                label = categorical_labels[item_index];
            else
                label = item_index;

            let Sx_tick = x_scale_func(label) + 0.5 * bin_width;
            svg += this.vertical_tick(Sx_tick, Sy_axis, tick_above, tick_below);
        }

        svg += this.end_group();                        // end of the ticks


        /* Handle the TICK LABELS (kept together as an SVG group)
         */
        svg += this.start_group("tick-labels");         // Pass a class used for all the labels

        for (let item_index = 0;  item_index < n_items;  ++item_index)  {
            if (categorical_labels && categorical_labels.length)
                label = categorical_labels[item_index];
            else
                label = item_index;

            // First, determine the horizontal position of the tick
            let Sx_tick = x_scale_func(label) + 0.5 * bin_width;
            // Now, position the label
            svg += this.text(label, Sx_tick, Sy_axis + tick_below, "", 0., 0.9);
            // Note: the last argument specifies a downward shift proportional to the font size,
            //       to clear the ticks regardless of font.  Extra control can be achieved with CSS
        }

        svg += this.end_group();                        // end of the labels

        return svg;

    } // axis_bottom_common



    axis_bottom_without_scale( {  Sxmin, Sxmax,
                    Sy_axis,
                    categorical_labels,
                    tick_above=0,
                    tick_below=6} )
    /*  Create a horizontal axis line meant to be placed below a plot.
        ALTERNATE VERSION that doesn't need the xscale function

        Use CSS for styling.  The following classes are created:
                            'axis-line', 'ticks', 'tick-labels'

        EXAMPLE:  g.tick-labels { translate(0, 10px); }  to shift down the tick labels

        Sxmin:              x-coord of left side of axis, in screen coordinates
        Sxmax:              x-coord of right side of axis, in screen coordinates
        Sy_axis:            y-coord of axis, in screen coordinates
        categorical_labels  List of desired label names (equally-spaced, at the center of their intervals)
        tick_above:         Amount by which ticks stick above axis, in screen coordinates
        tick_below:         Amount by which ticks stick below axis, in screen coordinates
     */
    {
        let svg = "";

        // Horizontal axis line
        svg += this.start_group("axis-line");
        svg += this.line(Sxmin, Sy_axis, Sxmax, Sy_axis);
        svg += this.end_group();

        const n_ticks = categorical_labels.length;
        const bin_width = (Sxmax - Sxmin) / n_ticks;

        /* Handle the TICKS (kept together as an SVG group)
         */
        svg += this.start_group("ticks");	            // Pass a class used for all the ticks

        for (let Sx_tick = Sxmin + 0.5 * bin_width;
                                Sx_tick <= Sxmax;
                                Sx_tick += bin_width)
        {
            svg += this.vertical_tick(Sx_tick, Sy_axis, tick_above, tick_below);
        }

        svg += this.end_group();                                // end of the ticks


        /* Handle the TICK LABELS (kept together as an SVG group)
         */
        svg += this.start_group("tick-labels");         // Pass a class used for all the labels

        let label_index = 0;
        for (let Sx_tick = Sxmin + 0.5 * bin_width;
                                Sx_tick <= Sxmax;
                                Sx_tick += bin_width)
        {
            svg += this.text(categorical_labels[label_index], Sx_tick, Sy_axis + tick_below, "dy=0.9em");
            // Note: dy specifies a downward shift proportional to the font size,
            //       to clear the ticks regardless of font.  Extra control can be achieved with CSS

            label_index += 1;
        }

        svg += this.end_group();                                // end of the labels

        return svg;

    } // axis_bottom_without_scale

}