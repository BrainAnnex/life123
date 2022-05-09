class SVGhelper
/* 	Helper library to generate either individual SVG tags, or larger SVG constructs.

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

    text(label, Sx, Sy, attributes = "")
    /*  Create a text label at the specified point.
        The x-position of the text refers to its left side, unless CSS
        directives such as "text-anchor: middle"  are used.
        The y-position refers to the BOTTOM of the text.
     */
    {
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

    axis_left( {    y_scale_func,
                    Sx_axis,
                    categorical_labels,
                    tick_right=0,
                    tick_left=6  } )
    /*  Create a vertical axis line meant to be placed to the left of a plot.

        Use CSS for styling.  The following classes are created:
                            'axis-line', 'ticks', 'tick-labels-left-axis'

        EXAMPLE:  g.tick-labels-left-axis { translate(0, 10px); }  to shift down the tick labels

        y_scale_func:       function produced with d3.scaleBand()
        Sx_axis:            x-coord of axis, in screen coordinates
        categorical_labels  List of desired label names (equally-spaced, at the center of their intervals)
        tick_right:         Amount by which ticks stick to the right of the axis, in screen coordinates
        tick_left:          Amount by which ticks stick to the left of the axis, in screen coordinates
     */
    {
        const n_items = categorical_labels.length;

        const bin_width = y_scale_func.bandwidth();

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
            let Sy_tick = y_scale_func(categorical_labels[item_index]) + 0.5 * bin_width;
            svg += this.horizontal_tick(Sx_axis, Sy_tick, tick_left, tick_right);
        }

        svg += this.end_group();                                // end of the ticks


        /* Handle the TICK LABELS (kept together as an SVG group)
         */
        svg += this.start_group("tick-labels-left-axis");         // Pass a class used for all the labels

        for (let item_index = 0;  item_index < n_items;  ++item_index)  {
            let Sy_tick = y_scale_func(categorical_labels[item_index]) + 0.5 * bin_width;
            svg += this.text(categorical_labels[item_index], Sx_axis - tick_left, Sy_tick, "dx=-0.9em dy=0.2em");
            // Note: dx and dy specify, respectively, a rightward/downward shift proportional to the font size,
            //       to clear the ticks and vertically-align the label regardless of font.
            //       Extra control can be achieved with CSS
        }

        svg += this.end_group();                                // end of the labels

        return svg;

    } // axis_left



    axis_bottom( {  x_scale_func,
                    Sy_axis,
                    categorical_labels,
                    tick_above=0,
                    tick_below=6  } )
    /*  Create a horizontal axis line meant to be placed below a plot.

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

        const bin_width = x_scale_func.bandwidth();

        const Sx1 = x_scale_func(categorical_labels[0]);
        const Sx2 = x_scale_func(categorical_labels[n_items-1]);

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

        for (let item_index = 0;  item_index < n_items;  ++item_index)  {
            let Sx_tick = x_scale_func(categorical_labels[item_index]) + 0.5 * bin_width;
            svg += this.vertical_tick(Sx_tick, Sy_axis, tick_above, tick_below);
        }

        svg += this.end_group();                                // end of the ticks


        /* Handle the TICK LABELS (kept together as an SVG group)
         */
        svg += this.start_group("tick-labels");         // Pass a class used for all the labels

        for (let item_index = 0;  item_index < n_items;  ++item_index)  {
            let Sx_tick = x_scale_func(categorical_labels[item_index]) + 0.5 * bin_width;
            svg += this.text(categorical_labels[item_index], Sx_tick, Sy_axis + tick_below, "dy=0.9em");
            // Note: dy specifies a downward shift proportional to the font size,
            //       to clear the ticks regardless of font.  Extra control can be achieved with CSS
        }

        svg += this.end_group();                                // end of the labels

        return svg;

    } // axis_bottom



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