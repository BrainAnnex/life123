Vue.component('curves-1',
    /*  Interpolated lots for functions in 1D

        Based on the hardwired dataset, create a group of dots
        connected by segments, plus interpolating curves/segments.
        Clicking on any dot will display its data in the browser console.

        Based on https://www.youtube.com/watch?v=CkFktv0p3pw

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
                // NOT SURE WHETHER USING...
                type: Number,
                default: 0
            },

            range_max: {
                // Value corresponding to the top of the graph
                type: Number,
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
                        </g>
                        <!-- END of the main part of the plot -->


                        <!--
                            Add the axes (they make use of CSS)
                        -->
                        <g class="horiz-axis" v-html="X_axis">
                        </g>

                    </g>
                    <!-- END of the translated element -->

                </svg>



                <svg width="500" height="500">

                    <!-- Show a sample fixed set of datapoints as little circles -->
                    <circle r="10"
                        v-for="(item, index) in dataset"
                            v-bind:key="index + '_old'"
                            v-bind:cx="item[0]"
                            v-bind:cy="item[1]"
                            fill="#555"
                            @click="on_click(item)"
                    />

                    <circle r="5"
                        v-for="(val, index) in data"
                            v-bind:key="index"

                            v-bind:cx="x_scale_func(index)"
                            v-bind:cy="val"
                            fill="#111"
                            @click="show_datapoint_info(index, val)"
                    />

                    <path stroke="green" stroke-width="2"
                        fill="none"
                        v-bind:d="path_data_straight"
                    />

                    <path stroke="red" stroke-width="4"
                        fill="none"
                        v-bind:d="path_data_curve"
                    />

                    <path stroke="yellow" stroke-width="3"
                        fill="none"
                        v-bind:d="path_data_steps"
                    />

                </svg>

                <button @click="switch_curve_type" style='font-weight:bold; padding:5px'>
                    Toggle between Step curves<br>(animated on Chrome)
                </button>
                <span style='color:grey; margin-left:10px'>(currently '{{curve_type}}')</span>

            </section>
            <!-- End of outer container -->
            `,


        data: function() {
            return {
                svg_helper: new SVGhelper(),

                dataset: [ [80, 200], [100, 123], [250, 200] , [250, 300] , [220, 400], [310, 380] ],
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
            {
                return this.outer_width - this.margins.left - this.margins.right;
            },

            plot_height()
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
                using existing values for the number of bins and the plot width.
             */
            {
                const f = d3.scaleLinear()
                            .domain([ 0, this.range_max ])
                            .range([ this.outer_height, 0 ]);     // f is a function
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


            path_data_straight()
            {
                const line_func = d3.line()
                                      .x(v => v[0])
                                      .y(v => v[1]);      // This will be a function

                return line_func(this.dataset);
            },

            path_data_curve()
            {
                const curve_func = d3.line()
                                     .curve(d3[this.curve_type])
                                     .x(v => v[0])
                                     .y(v => v[1]);      // This will be a function

               return curve_func(this.dataset);
            },

            path_data_steps()
            {
                const curve_func = d3.line()
                                     .curve(d3.curveStep)
                                     .x(v => v[0])
                                     .y(v => v[1]);      // This will be a function

               return curve_func(this.dataset);
            }

        },  // COMPUTED


        // ---------------------------  METHODS  ---------------------------
        methods: {

            show_datapoint_info(n, val) {
                console.log(`Bin ${n}, value ${val}`);
            },

            on_click(item) {
                console.log("This item: ", item);
            },

            switch_curve_type()  {
                //this.curve_type = (this.curve_type === "curveBasis" ? "curveNatural" : "curveBasis");
                this.curve_type = (this.curve_type === "curveStepAfter" ? "curveStepBefore" : "curveStepAfter");
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