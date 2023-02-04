Vue.component('vue-heatmap-9',
    /*  A simple heatmap in 2D.  High values are shown as dark (think of ink in water); low values in white.
        An outer box in pale cyan is also shown.
        See: https://julianspolymathexplorations.blogspot.com/2022/01/D3-plus-Vue-visualization-UI.html
        It needs the SVG_helper library for drawing the axes.
     */
    {
        props: ['my_groups', 'my_vars', 'my_data', 'range_min', 'range_max', 'outer_width', 'outer_height', 'margins'],
        /*
            my_groups:      list of x-axis labels.  EXAMPLE: ["A", "B", "C"]
            my_vars:        list of y-axis labels.  EXAMPLE: ["v1", "v2"]
            my_data:        array of objects with 3 keys ('group', 'variable' and 'value')
                            'group' refers to the x-axis, and 'variable' to the y-axis
                            EXAMPLE:  [
                                    { "group": "A", "variable": "v1", "value": "0" },
                                    { "group": "A", "variable": "v2", "value": "82" }
                                   ]
            "range_min":    value of heatmap cell to map to white
            "range_max":    value of heatmap cell to map to black
            "outer_width":  in pixels.  EXAMPLE: 850
            "outer_height": in pixels
            "margins":      Affecting the outside of the container box.
                            EXAMPLE: {"top": 30, "right": 30, "bottom": 30, "left": 30}
         */

        template: `
            <!-- Outer container, serving as Vue-required component root  -->
            <svg v-bind:width="outer_width" v-bind:height="outer_height" class="chart-holder">
                <g v-bind:transform="translate(margins.left, margins.top)">

                    <!-- The main part of the heatmap (a series of rectangles, colored based on the cell value) -->
                    <g class="heatmap">
                        <rect v-for="(item, index) in my_data"
                            v-bind:key="index"
                            v-bind:x="x_scale_func(item.group)"  v-bind:y="y_scale_func(item.variable)"
                            v-bind:width="rect_w"  v-bind:height="rect_h"
                            v-bind:fill="color_scale_func(item.value)"
                            stroke='rgb(200,200,200)' stroke-width='1'
                        >
                        </rect>
                    </g>
                    <!-- End of the main part of the heatmap -->

                    <!--
                        Insert some lines in the plot (the CSS classes are not currently used)
                      -->

                    <!-- A line at the top -->
                    <g class="top-line" v-html="top_line">
                    </g>

                    <!-- A line on the right -->
                    <g class="right-line" v-html="this.svg_helper.line(plot_width, 0, plot_width, plot_height, 'purple')">
                    </g>


                    <!--
                        Add the axes (they make use of CSS)
                      -->
                    <g class="horiz-axis" v-html="X_axis">
                    </g>

                    <g class="vert-axis" v-html="Y_axis">
                    </g>

                </g>
                <!-- End of the translated element -->
            </svg>
            <!-- End of outer container -->
            `,


        data: function() {
            return {
                svg_helper: new SVGhelper(),
                min_val: this.range_min,
                max_val: this.range_max
            }
        }, // data



        // ---------------------------  COMPUTED  ---------------------------
        computed: {     // NOTE: computed methods are only invoked AS NEEDED

            plot_width()
            {
                return this.outer_width - this.margins.left - this.margins.right;
            },

            plot_height()
            {
                return this.outer_height - this.margins.top - this.margins.bottom;
            },


            x_scale_func()
            /*  Create and return a function to build the X scale.
                This function maps a "group" name (in heatmap_data) into an X value in screen coordinates.
                EXAMPLE, if the x_labels are ["A", "B", "C"], and the plot_width is 100:
                            "A" |-> 0  , "B" |-> 33.33 , "C" |-> 66.66

                Example of test in a browser's JS console:
                    >> f = d3.scaleBand().domain(["A", "B", "C"]).range([0, 100])
                    >> f("A") will give 0 ,  f("B") will give 33.33, etc
             */
            {
                const f = d3.scaleBand()
                            .domain(this.my_groups)
                            .range([ 0, this.plot_width ]);     // f is a function
                return f;
            },

            y_scale_func()
            /*  Create and return a function to build the Y scale.
                This function maps a "variable" name (in my_data) into a Y value in screen coordinates
                EXAMPLE:  "v1" |-> 195  , "v2" |-> 0
             */
            {
                const f = d3.scaleBand()
                            .domain(this.my_vars)
                            .range([ this.plot_height, 0 ]);   // f is a function
                return f;
            },

            color_scale_func()
            /*  Create and return a function to build the color scale.
                This returns maps a number ("value" in my_data) into a color code
                EXAMPLES:  0 |-> "rgb(255, 255, 255)"  , 10 |-> "rgb(240, 247, 246)" , 100 |-> "rgb(105, 179, 162)"

                Example of test in JS console:
                    f = d3.scaleLinear().domain([0, 100]).range(["white", "black"])
                    f(50) will give "rgb(128, 128, 128)"
                    Likewise for f("50")
                    It's also acceptable for the domain to contain strings, such as domain(["0", "100"])
             */
            {
                const f = d3.scaleLinear()
                    .domain([this.min_val, this.max_val])
                    .range(["white", "black"]);   // Maps 0 to white, and 100 to black

                return f;   // f is a function
            },


            rect_w()
            /*  Return the width (in pixels) of each rectangle element in the heatmap, based on the previously-set x scale.
                EXAMPLE, if in the earlier call to x_scale_func(),
                    the x_labels is a list with 4 elements, and the plot_width is 600,
                    then rect_w() returns 150
             */
            {
                return this.x_scale_func.bandwidth();
            },

            rect_h()
            /*  Return the height (in pixels) of each rectangle element in the heatmap, based on the previously-set y scale.
                EXAMPLE, in the earlier call to y_scale_func(),
                    if the y_labels is a list with 2 elements, and the plot_height is 400,
                    then rect_h() returns 200
             */
            {
                return this.y_scale_func.bandwidth();
            },

            top_line()
            // Return the SVG code to produce a line at the top of the graph
            {
                return this.svg_helper.line(0, 0, this.plot_width, 0, "blue");
            },

            X_axis()
            // Return the SVG code to produce an x-axis
            {
                return this.svg_helper.axis_bottom(
                            {x_scale_func: this.x_scale_func,
                             Sy_axis: this.plot_height,
                             categorical_labels: this.my_groups
                            }
                        );
                /*
                // Alternative that doesn't use the x scale function
                return this.svg_helper.axis_bottom_without_scale(
                            {Sxmin: 0, Sxmax: this.plot_width,
                             Sy_axis: this.plot_height,
                             categorical_labels: this.my_groups
                            }
                        );
                */
            },

            Y_axis()
            // Return the SVG code to produce a y-axis
            {
                return this.svg_helper.axis_left(
                            {y_scale_func: this.y_scale_func,
                             Sx_axis: 0,
                             categorical_labels: this.my_vars
                            }
                        );
            }

        },  // COMPUTED


        // ---------------------------  METHODS  ---------------------------
        methods: {

            translate(x, y)
            /*  Return a string suitable as an SVG attribute, to indicate a translation by <x,y>
                EXAMPLE:   "translate(10, 50)"
             */
            {
                return `translate(${x}, ${y})`;
            }

        }  // METHODS

    }
); // end component