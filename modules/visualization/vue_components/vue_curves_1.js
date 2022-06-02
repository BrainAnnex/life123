Vue.component('vue_curves_1',
    /*  Line charts and interpolating functions in 1D.
        For now, just 1 dataset at a time.

        Based on the passed dataset, create a group of dots
        connected by segments, plus interpolating curves/segments.
        Clicking on any dot will display its data in the browser console.

        Inspired by https://www.youtube.com/watch?v=CkFktv0p3pw

        DEPENDENCIES:   - the SVG_helper library for drawing the axes
                        - the D3 (v7) libraries
     */
    {
        props: {

            y_labels: {
                // List of y-axis labels.  EXAMPLE: ["Chem 1", "Chem 2"]
                // Note: the x-axis labels are just the bin numbers
            },

            data: {
                /* List of row values.  Each item is a set of bin values from left to right;
                        EXAMPLE:  [20, 85, 100, 50]
                 */
                required: true
            },

            range_min: {
                // NOT SURE WHETHER USING THIS...
                type: Number,
                default: 0
            },

            range_max: {
                // Value corresponding to the top of the graph
                type: Number,
                required: true
            },

            outer_width: {
                // For the container box.  In pixels.  EXAMPLE: 850
                type: Number,
                required: true
            },

            outer_height: {
                // For the container box.  In pixels.  EXAMPLE: 350
                type: Number,
                required: true
            },

            margins: {
                /* Affecting the outside of the container box.
                             EXAMPLE: {"top": 30, "right": 30, "bottom": 30, "left": 30}
                 */
            }
        },


        template: `
            <!-- Outer container, serving as Vue-required component root  -->
            <section>

                <svg v-bind:width="outer_width" v-bind:height="outer_height" class="chart-holder">
                    <g v-bind:transform="translate(margins.left, margins.top)"> <!-- Shift the contained g block below -->

                        <!-- The main part of the plot -->
                        <g class="main-plot">
                            <!-- Show the datapoints as little circles -->
                            <circle r="5"
                                v-for="(val, index) in data"
                                    v-bind:key="index"

                                    v-bind:cx="x_scale_func(index)"
                                    v-bind:cy="y_scale_func(val)"
                                    fill="#111"
                                    @click="show_datapoint_info(index, val)"
                            />

                            <path stroke="yellow" stroke-width="1"
                                fill="none"
                                v-bind:d="path_straight"
                            />

                            <path stroke="gray" stroke-width="1"
                                fill="none"
                                v-bind:d="path_steps"
                            />

                            <path stroke="red" stroke-width="2"
                                fill="none"
                                v-bind:d="path_curve"
                            />

                        </g>
                        <!-- END of the main part of the plot -->


                        <!--
                            Add the axes (they make use of CSS)
                        -->
                        <g class="horiz-axis" v-html="X_axis">
                        </g>
                        <g class="vert-axis" v-html="Y_axis">
                        </g>

                    </g>
                    <!-- END of the translated element -->

                </svg>


                <button @click="cycle_curve_types" style='font-weight:bold; padding:5px'>
                    Cycle between curve types<br>
                    <span style='color:grey; margin-left:10px'>(currently '{{curve_type}}')</span>
                </button>


            </section>
            <!-- End of outer container -->
            `,


        data: function() {
            return {
                svg_helper: new SVGhelper(),

                curve_type: "curveNatural"

                //min_val: this.range_min,
                //max_val: this.range_max,
            }
        }, // data



        // ---------------------------  COMPUTED  ---------------------------
        computed: {     // NOTE: computed methods are only invoked AS NEEDED

            n_bins()
            {
                if (this.data.length == 0)  {
                    console.error("n_bins(): the data is empty");
                    return 0;
                }

                return this.data.length;
            },


            plot_width()
            // For the main part of the graph
            {
                return this.outer_width - this.margins.left - this.margins.right;
            },

            plot_height()
            // For the main part of the graph
            {
                return this.outer_height - this.margins.top - this.margins.bottom;
            },


            x_scale_func()
            /*  Create and return a function to build the X scale.
                This function maps a "bin index" into an X-value in screen coordinates,
                using existing values for the number of bins and the plot width.
                EXAMPLE, if there are 4 bins, and the plot_width is 600:
                            0 |-> 0  , 1 |-> 150 , 2 |-> 300 , 3 |-> 450

                Example of test in a browser's JS console:
                    >> f = d3.scaleLinear().domain([0, 4]).range([0, 600])
                    >> f(0) will give 0 ,  f(1) will give 150, etc
             */
            {
                const f = d3.scaleLinear()
                            .domain([ 0, this.n_bins ])
                            .range([ 0, this.plot_width ]);     // f is a function
                return f;
            },

            y_scale_func()
            /*  Create and return a function to build the Y scale.
                This function maps a "data value" into a Y-value in screen coordinates,
                based on the previously-set maximum y range (the range starts at zero)
                and the plot height.
             */
            {
                const f = d3.scaleLinear()
                            .domain([ 0, this.range_max ])
                            .range([ this.plot_height, 0 ]);     // f is a function
                            // Note the UP/DOWN reversal between data coord and screen coords!
                return f;
            },


            bin_width()
            /*  Return the width (in pixels) of each rectangle element in the heatmap.
                EXAMPLE, if there are 4 bins and a plot_width of 600,
                    then rect_w() returns 150
             */
            {
                return this.plot_width / this.n_bins;
            },


            X_axis()
            // Return the SVG code to produce an x-axis
            {
                return this.svg_helper.axis_bottom_scaleLinear(
                            {x_scale_func: this.x_scale_func,
                             n_items: this.n_bins,
                             bin_width: this.bin_width,
                             Sy_axis: this.plot_height
                            }
                        );
            },

            Y_axis()
            // Return the SVG code to produce a y-axis
            {
                return this.svg_helper.axis_left_scaleLinear(
                            {y_scale_func: this.y_scale_func,
                             n_intervals: 10,
                             y_min: 0,
                             y_max: this.range_max,
                             Sx_axis: 0
                            }
                        );
            },


            path_straight()
            // Connect the data points with segments
            {
                // The x-coord is the array index; the y-coord is the data value
                const line_func = d3.line()
                                    .x((v, i) => this.x_scale_func(i))
                                    .y(v      => this.y_scale_func(v));      // This will be a function

                return line_func(this.data);
            },

            path_steps()
            // Connect the data points with a series of steps
            {
                // The x-coord is the array index; the y-coord is the data value
                const line_func = d3.line()
                                    .curve(d3.curveStepAfter)
                                    .x((v, i) => this.x_scale_func(i))
                                    .y(v      => this.y_scale_func(v));      // This will be a function

                return line_func(this.data);
            },

            path_curve()
            // Connect the data points in interpolating curve
            {
                // The x-coord is the array index; the y-coord is the data value
                const line_func = d3.line()
                                    .curve(d3[this.curve_type])
                                    .x((v, i) => this.x_scale_func(i))
                                    .y(v      => this.y_scale_func(v));      // This will be a function

                return line_func(this.data);
            }

        },  // COMPUTED



        // ---------------------------  METHODS  ---------------------------
        methods: {

            show_datapoint_info(n, val)
            // Print out to the console the point's graph coordinates
            {
                console.log(`Bin ${n}, value ${val}`);
            },


            cycle_curve_types()
            // Cycle among different available types of curve interpolators
            {
                const options = ["curveNatural", "curveBasis", "curveMonotoneX"];
                let pos = options.indexOf(this.curve_type);
                pos += 1;
                if (pos >= options.length)
                    pos = 0;

                this.curve_type = options[pos];
            },


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