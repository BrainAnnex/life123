Vue.component('vue_curves_4',
    /*  Line charts and interpolating curves in 1D.
        This version can handle multiple datasets.

        Based on each of the the passed dataset, create a group of dots
        optionally connected by segments, interpolating curves and steps.
        Clicking on any dot will display its data in the browser console and on the UI.

        Inspired by https://www.youtube.com/watch?v=CkFktv0p3pw

        DEPENDENCIES:   - the SVG_helper (v1.2) library for drawing the axes
                        - the D3 (v7) libraries

        CHANGES FROM THE PREVIOUS VERSION:
                - support for multiple datasets, incl. a legend
                - ditched curveNatural from the option for the interpolating curves
                - added a border around the whole component (CSS: section.graphic-component)
                - data from clicking on dots is also shown on the UI
                - checkboxes to let user display segments, interpolating curves and/or steps

        TODO:   - Expand the number of colors (currently, just 4)
                - Allow users to dynamically change y-range
     */
    {
        props: {

            curve_labels: {
                // List of labels for the curves.  EXAMPLE: ["Chem 1", "Chem 2"]
            },

            plot_data: {
                /*  Concentration data for the plots;
                    at the outer level in order of chemical-species index,
                    and at the inner level in bin index order from left to right
                        EXAMPLE:  [[20, 85, 100, 50],     <- 0-th chemical
                                   [14, 99, 5, 65]        <- 1st chemical
                 */
                required: true
            },

            range_min: {
                // NOT USING THIS FOR NOW...
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
            },

            color_mapping: {
                default: {    // TODO: define more colors!
                     0: 'red',
                     1: 'orange',
                     2: 'green',
                     3: 'blue',
                     4: 'pink'
                }
            }
        },


        template: `
            <!-- Outer container, serving as Vue-required component root  -->
            <section class="graphic-component">

                <svg v-bind:width="outer_width" v-bind:height="outer_height" class="chart-holder">
                    <g v-bind:transform="translate(margins.left, margins.top)"> <!-- Shift the contained g blocks below -->

                        <!--
                            Add the axes (they make use of CSS)
                        -->
                        <g class="horiz-axis" v-html="X_axis">
                        </g>
                        <g class="vert-axis" v-html="Y_axis">
                        </g>


                        <!-- The MAIN PART of the plot -->
                        <g class="main-plot">

                            <!-- Connect the data points of each plot with line segments -->
                            <path v-if="show_segments" v-for="(p, index_p) in plot_data"
                                v-bind:key="'seg' + index_p"
                                v-bind:stroke="color_picker(index_p, true)"
                                v-bind:d="path_straight(index_p)"
                                stroke-width="1"
                                fill="none"
                            />

                            <!-- Connect the data points of each plot with a series of steps -->
                            <path v-if="show_steps" v-for="(p, index_p) in plot_data"
                                v-bind:key="'steps' + index_p"
                                v-bind:stroke="color_picker(index_p, true)"
                                v-bind:d="path_steps(index_p)"
                                stroke-width="1"
                                fill="none"
                            />

                            <!-- Connect the data points of each plot with interpolating curves -->
                            <path v-if="show_curves"  v-for="(p, index_p) in plot_data"
                                v-bind:key="'inter' + index_p"
                                v-bind:stroke="color_picker(index_p, false)"
                                v-bind:d="path_curve(index_p)"
                                stroke-width="1"
                                fill="none"
                            />

                            <!-- Show the data points in each plot as little circles (on top of all curves/segments) -->
                            <template v-for="(plot, index_plot) in plot_data">
                                <circle v-for="(val, index) in plot_data[index_plot]"
                                    v-bind:key="index_plot + '-' + index"
                                    v-bind:cx="x_scale_func(index)"
                                    v-bind:cy="y_scale_func(val)"
                                    r="5"
                                    fill="#111"
                                    fill-opacity='0.1'
                                    @click="show_datapoint_info(index, val, index_plot)"
                                />
                            </template>

                        </g>
                        <!-- END of the main part of the plot -->


                    </g>
                    <!-- END of the translated elements -->

                </svg>


                <div style='display: inline-block; vertical-align:top'> <!-- Side box with legend and controls -->

                    <!-- The plot legend-->
                    <svg width="150" v-bind:height="legend_SVG_height">
                        <rect x='5' y='10' width='140' v-bind:height="legend_box_height" stroke="#E4E4E4" stroke-width="1" fill='#F6F6F6'/>

                        <g v-for="(label, index_labels) in curve_labels">
                            <text x="15" v-bind:y="legend_text_Y(index_labels)" v-bind:fill="color_picker(index_labels, false)">{{label}}</text>
                            <line x1="15" v-bind:y1="legend_Y(index_labels)" x2="135" v-bind:y2="legend_Y(index_labels)" v-bind:stroke="color_picker(index_labels, false)" />
                        </g>
                    </svg>

                    <br><span style="color:gray">{{info_line}}</span>

                    <!-- Graph controls -->

                    <!-- Specify what to show -->
                    <br><br>
                    <input type="checkbox" v-model="show_steps"> Show steps<br>
                    <input type="checkbox" v-model="show_segments"> Show segments<br>
                    <input type="checkbox" v-model="show_curves"> Show curves<br>

                    <!-- A button to choose the type of interpolating curves -->
                    <button @click="cycle_curve_types" style='padding:2px'>
                        <span style='font-weight:bold'>Toggle<br>curve type</span><br>
                        <span style='color:grey'>(using '{{curve_type}}')</span>
                    </button>
                </div>


            </section>
            <!-- End of outer container -->
            `,


        data: function() {
            return {
                svg_helper: new SVGhelper(),

                curve_type: "curveMonotoneX",

                show_curves: true,
                show_segments: false,
                show_steps: false,

                info_line: "click points for details"

                //min_val: this.range_min,
                //max_val: this.range_max,
            }
        }, // data of the Vue component



        // ---------------------------  COMPUTED  ---------------------------
        computed: {     // NOTE: computed methods are only invoked AS NEEDED

            n_bins()
            // TODO: should check that all elements have the same length
            {
                if (this.plot_data[0].length == 0)  {
                    console.error("n_bins(): the data is empty");
                    return 0;
                }

                return this.plot_data[0].length;
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
                            .domain([ -0.5, this.n_bins -0.5 ])
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
            /*  Return the width (in pixels) of bin (subdivision) in the system.
                EXAMPLE, if there are 4 bins and a plot_width of 600,
                    then rect_w() returns 150
             */
            {
                return this.plot_width / this.n_bins;
            },


            X_axis()
            // Return the SVG code to produce an x-axis
            {
                return this.svg_helper.axis_bottom_scaleLinear_bins(
                            {x_scale_func: this.x_scale_func,
                             n_items: this.n_bins,
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

            legend_SVG_height()
            {
                return this.legend_box_height + 15;
            },

            legend_box_height()
            {
                return 10 + 25 * this.curve_labels.length;
            },

        },  // COMPUTED



        // ---------------------------  METHODS  ---------------------------
        methods: {

            legend_text_Y(index)
            {
                return 30 + index*25;
            },

            legend_Y(index)
            {
                return 35 + index*25;
            },



            pale_version(color)
            /*  Return the paler version of the given color,
                specified as a string name ("red") or string RBG hex ("#FFAAAA").
                Return the new color as a string RBG hex.
                EXAMPLE: pale_version("red") and pale_version("#FF0000") return "#ffb3b3"
             */
            {
                //console.log(color);
                const c = d3.hsl(color);        // Convert to HSL coordinates
                //console.log(c);
                c.l = 0.8;                      // Set the luminosity relatively close to the max
                const c_new = c.formatHex();    // Convert back to RBG hex
                //console.log(c_new);
                return c_new;
            },

            color_picker(index, pale)
            /*  Given an integer index (0, 1, 2, ...) return a color;
                If the pale flag is true, return a pale version of the color instead.
                A default color (and pale version) are provided if the index isn't recognized.
             */
            {
                var color;

                if (index in this.color_mapping)
                    color = this.color_mapping[index];
                else
                    color = 'black';

                if (pale)
                    return this.pale_version(color);
                else
                    return color;
            },


            path_straight(i)
            // Connect the data points of the i-th graph with line segments
            {
                // The x-coord is the array index; the y-coord is the data value
                const line_func = d3.line()
                                    .x((v, i) => this.x_scale_func(i))
                                    .y(v      => this.y_scale_func(v));      // This will be a function

                return line_func(this.plot_data[i]);
            },

            path_curve(i)
            // Connect the data points of the i-th graph with an interpolating curve
            {
                // The x-coord is the array index; the y-coord is the data value
                const line_func = d3.line()
                                    .curve(d3[this.curve_type])
                                    .x((v, i) => this.x_scale_func(i))
                                    .y(v      => this.y_scale_func(v));      // This will be a function

                return line_func(this.plot_data[i]);
            },

            path_steps(i)
            // Connect the data points of the i-th graph with a series of steps
            {
                // The x-coord is the array index; the y-coord is the data value
                const line_func = d3.line()
                                    .curve(d3.curveStepAfter)
                                    .x((v, i) => this.x_scale_func(i-0.5))
                                    .y(v      => this.y_scale_func(v));      // This will be a function

                return line_func(this.extended_data(i));
            } ,

            extended_data(i)
            // Return an array that is the original data of the i-th graph, with the last entry appended (i.e. repeated twice)
            {
                const last_element = this.plot_data[i][this.n_bins - 1];

                return this.plot_data[i].concat(last_element);
            },


            show_datapoint_info(n, val, plot_no)
            // Print out to the console and to the UI the point's graph coordinates
            {
                var label = this.curve_labels[plot_no];
                var message = `${label}: bin ${n}, value ${val}`
                console.log(message);
                this.info_line = message;
            },


            cycle_curve_types()
            // Cycle among different available types of curve interpolators
            {
                const options = ["curveBasis", "curveMonotoneX"];
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