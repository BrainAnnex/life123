Vue.component('vue-heatmap-9',
    /*  A heatmap.
        It needs the SVG_helper library for drawing the axes.
     */
    {
        props: ['my_groups', 'my_vars', 'my_data', 'range_min', 'range_max', 'outer_width', 'outer_height', 'margins'],
        /*
            my_data:    array of objects with 3 keys ('group', 'variable' and 'value')
         */

        template: `
            <!-- Outer container, serving as Vue-required template root  -->
            <svg v-bind:width="outer_width" v-bind:height="outer_height" class="chart-holder">
                <g v-bind:transform="translate(margins.left, margins.top)">

                    <!-- The main part of the heatmap (a series of rectangles) -->
                    <g class="heatmap">
                        <rect v-for="(item, index) in my_data"
                            v-bind:key="index"
                            v-bind:x="rect_x(item)"  v-bind:y="rect_y(item)"
                            v-bind:width="rect_w"  v-bind:height="rect_h"
                            v-bind:fill="rect_color(item)"
                            stroke='rgb(200,200,200)' stroke-width='1'
                        >
                        </rect>
                    </g>


                    <!--
                        Insert some lines in the plot (no CSS used)
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
                This function maps a "group" name (in my_data) into an X value in screen coordinates.
                EXAMPLE:  "A" |-> 0  , "B" |-> 130 , "C" |-> 260
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
             */
            {
                const f = d3.scaleLog()
                    .domain([this.min_val, this.max_val])
                    .range(["white", "black"]);   // Maps 0 to white, and 100 to black

                return f;   // f is a function
            },


            rect_w()
            {
                return this.x_scale_func.bandwidth();
            },

            rect_h()
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
            },

            rect_x(item)
            // Return a number for the left x-coordinate of the heatmap rectangle
            {
                return this.x_scale_func(item.group);
            },

            rect_y(item)
            // Return a number for the top y-coordinate of the heatmap rectangle
            {
                return this.y_scale_func(item.variable);
            },


            rect_color(item)
            // Return a string such as "rgb(240, 247, 246)", to be used for a heatmap rectangle
            {
                var color_func = this.color_scale_func;   // This is a function
                return color_func(item.value);
            }

        }  // METHODS

    }
); // end component