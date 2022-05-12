Vue.component('vue-heatmap-11',
    /*  A heatmap in 2D, a small change from 'vue-heatmap-10'
        (CHANGED: the x-axis labels are now just the bin numbers; a tooltip added)

        High values are shown as dark (think of ink in water); low values in white.
        An outer box (with border set by CSS "chart-holder") is also shown.
        See: https://julianspolymathexplorations.blogspot.com/2022/01/D3-plus-Vue-visualization-UI.html
        Loosely based on: https://d3-graph-gallery.com/graph/heatmap_tooltip.html

        DEPENDENCIES:   - the SVG_helper library for drawing the axes
                        - the D3 (v7) libraries
     */
    {
        props: {

            y_labels: {
                // List of y-axis labels.  EXAMPLE: ["Chem 1", "Chem 2"]
                // Note: the x-axis labels are just the bin numbers
                required: true
            },

            heatmap_data: {
                /* List of row values.  Each item is a set of bin values from left to right;
                                               consecutive rows are moving UP along the y-axis
                                               EXAMPLE:  [
                                                               [10., 32., 2.6],
                                                               [90., 14.5, 55.1]
                                                         ]
                 */
                required: true
            },

            range_min: {
                // Value of heatmap cell that maps to white
                type: Number,
                default: 0
            },

            range_max: {
                // Value of heatmap cell that maps to black
                type: Number,
                required: true,
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
            <div>  <!-- Outer container, serving as Vue-required component root  -->
                <svg v-bind:width="outer_width" v-bind:height="outer_height" class="chart-holder">
                    <g v-bind:transform="translate(margins.left, margins.top)"> <!-- Shift the contained g block below -->

                        <!-- The main part of the heatmap (a series of rectangles, colored according to the cell value) -->
                        <g class="heatmap" @mouseleave='mouse_leave_heatmap'>
                            <!-- For each row -->
                            <template v-for="(row_data, row_index) in heatmap_data">

                                <!-- For each column, in the row under consideration -->
                                <template v-for="(heatmap_value, col_index) in row_data">
                                    <rect
                                        v-bind:key="row_index + '_' + col_index"

                                        v-bind:x="x_scale_func(col_index)"  v-bind:y="y_scale_func(y_labels[row_index])"
                                        v-bind:width="rect_w"  v-bind:height="rect_h"
                                        v-bind:fill="color_scale_func(heatmap_value)"
                                        stroke="rgb(200,200,200)" stroke-width="1"

                                        @mouseover='mouseover_on_heatmap'
                                        @mousemove='mouse_move(heatmap_value, $event)'
                                    >
                                    </rect>
                                </template>

                            </template>
                        </g>
                        <!-- END of the main part of the heatmap -->


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
                    <!-- END of the translated element -->
                </svg>


                <!-- Movable DIV that shows a tooltip -->
                <div class='tooltip'
                    v-bind:style="{'left': tooltip_left + 'px', 'top': tooltip_top + 'px', 'opacity': tooltip_opacity}"
                    @mouseover='mouseover_on_tooltip' @mouseleave='mouse_leave_tooltip'>
                {{tooltip_value}}
                </div>


                <!--  Slider, to let the user adjust the max value of the heatmap range -->
                <div>
                    <span style="color: #888; margin-right:25px">Heatmap range: </span>
                    <span>{{min_val}}</span> to
                    <input type="range"
                            v-bind:min="slider_min" v-bind:max="original_max_val"
                            v-bind:step="slider_step"
                            v-model:value="max_val">
                    <label>{{max_val}}</label>
                </div>


            </div>
            <!-- End of outer container -->
            `,


        data: function() {
            return {
                svg_helper: new SVGhelper(),
                min_val: this.range_min,
                max_val: this.range_max,
                original_max_val: this.range_max,    // This will become the rightmost value on the slider

                tooltip_value: null,
                tooltip_opacity: 0,
                tooltip_left: 0,
                tooltip_top: 0
            }
        }, // data



        // ---------------------------  COMPUTED  ---------------------------
        computed: {     // NOTE: computed methods are only invoked AS NEEDED

            slider_min()
            /*  To avoid letting the slider go all the way to the minimum value,
                which would causes singularities in the value mappings
             */
            {
                return this.min_val + this.slider_step;
            },

            slider_step()
            // Make slider steps a 20-th of the range
            {
                return (this.original_max_val - this.min_val) / 20.;
            },

            plot_width()
            {
                return this.outer_width - this.margins.left - this.margins.right;
            },

            plot_height()
            {
                return this.outer_height - this.margins.top - this.margins.bottom;
            },

            n_bins()
            {
                if (this.heatmap_data.length == 0)  {
                    console.error("n_bins(): the heatmap data is empty");
                    return 0;
                }

                const first_row = this.heatmap_data[0];
                return first_row.length;
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
                This function maps a "y_labels" entry into a Y-value in screen coordinates.

                EXAMPLE, if the y_labels is ["Chem 1", "Chem 2"], and the plot_height is 400:
                            "Chem 1" |-> 200  , "Chem 2" |-> 0
             */
            {
                const f = d3.scaleBand()
                            .domain(this.y_labels)
                            .range([ this.plot_height, 0 ]);   // f is a function
                return f;
            },


            rect_w()
            /*  Return the width (in pixels) of each rectangle element in the heatmap.
                EXAMPLE, if there are 4 bins and a plot_width of 600,
                    then rect_w() returns 150
             */
            {
                return this.plot_width / this.n_bins;
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


            color_scale_func()
            /*  Create and return a function to build the color scale.
                This returns maps a number ("value" in heatmap_data) into a color code
                EXAMPLES:  0 |-> "rgb(255, 255, 255)"  , 10 |-> "rgb(240, 247, 246)" , 100 |-> "rgb(105, 179, 162)"

                Example of test in a browser's JS console:
                    >> f = d3.scaleLinear().domain([0, 100]).range(["white", "black"])
                    >> f(50) will give "rgb(128, 128, 128)"
                    Likewise for  >>  f("50")
                    It's also acceptable for the domain to contain strings, such as domain(["0", "100"])

                    Can also issue commands such as >> $vm0.color_scale_func(10)
             */
            {
                if ( this.min_val === undefined ) {
                    alert("The range_min parameter is undefined");
                }
                if ( this.max_val === undefined ) {
                    alert("The range_max parameter is undefined");
                }
                const f = d3.scaleLinear()
                    .domain([this.min_val, this.max_val])
                    .range(["white", "black"]);   // Maps the min value to white, and the max value to black

                return f;   // f is a function
            },



            top_line()
            // Return the SVG code to produce a line at the top of the graph
            {
                return this.svg_helper.line(0, 0, this.plot_width, 0, "blue");
            },


            X_axis()
            // Return the SVG code to produce an x-axis
            {
                return this.svg_helper.axis_bottom_scaleLinear(
                            {x_scale_func: this.x_scale_func,
                             n_items: this.n_bins,
                             bin_width: this.rect_w,
                             Sy_axis: this.plot_height
                            }
                        );
            },

            Y_axis()
            // Return the SVG code to produce a y-axis
            {
                return this.svg_helper.axis_left(
                            {y_scale_func: this.y_scale_func,
                             Sx_axis: 0,
                             categorical_labels: this.y_labels
                            }
                        );
            }

        },  // COMPUTED


        // ---------------------------  METHODS  ---------------------------
        methods: {
            mouseover_on_heatmap()
            // Upon detecting the mouse enter the main part of the heatmap
            {
                //console.log("mouseover_on_heatmap");
                // Make tooltip visible
                this.tooltip_opacity = 1;
            },
            mouseover_on_tooltip(ev)
            /* This is needed for times when the mouse has landed on the tooltip DIV;
                a consequence is that it has left the heatmap, and thus been made invisible:
                make it visible again! (though in a dim way that reflects the different modality)
             */
            {
                //console.log("mouseover_on_tooltip");
                // Make tooltip visible, but dimly
                this.tooltip_opacity = 0.3;
            },

            mouse_move(val, ev)
            // Upon detecting the mouse in motion on a heatmap rectangle
            {
                //console.log("In mouse_move", ev);
                //console.log("Value: ", val, " | event.x: ", ev.x, " | event.y: ", ev.y);
                // Show the value of the heatmap rectangle
                this.tooltip_value = val;
                // Position the tooltip box a little to the right and above the mouse
                this.tooltip_left = ev.x + 30;
                this.tooltip_top = ev.y - 30;
            },

            mouse_leave_heatmap()
            // Upon detecting the mouse leave the main part of the heatmap
            {
                //console.log("mouse_leave_heatmap");
                this.tooltip_opacity = 0;
            },
            mouse_leave_tooltip()
            /* This is needed for times when the mouse has landed on the tooltip DIV
               upon leaving the heatmap, and thus is has been made visible again;
               hide it for good when it leaves the tooltip DIV
             */
            {
                //console.log("mouse_leave_tooltip");
                this.tooltip_opacity = 0;
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