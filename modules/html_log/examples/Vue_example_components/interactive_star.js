Vue.component('vue-interactive-star',
    /*  Multipoint star, with a slider to pick how many arms it has.
        Based on https://ultimatecourses.com/blog/create-data-visualization-chart-vuejs-d3
     */
    {
        props: {
                outer_radius: {
                    default: 200
                }
        },

        template: `
            <!-- Outer container, serving as Vue-required template root  -->
            <section>

                <svg v-bind:width="svg_width" v-bind:height="svg_height">
                    <path
                        class="radial"
                        v-bind:d="radial_data"
                        v-bind:transform="this.svg_helper.translate(svg_width/2, svg_height/2)"
                        fill="gold"
                    ></path>
                </svg>

                <div class="range-input">
                    <label for="rays"># of Rays</label>
                    <input name="rays" type="range" min="4" max="60" v-model="number_rays" />
                    <span style="margin-left:10px">{{number_rays}}</span>
                </div>

            </section>
            <!-- End of outer container -->
            `,


        data: function() {
            return {
                number_rays: 8,                     // The number of points in the star
                svg_width: this.outer_radius * 2,
                svg_height: this.outer_radius * 2,
                svg_helper: new SVGhelper()
            }
        }, // data


        // ---------------------------  COMPUTED  ---------------------------
        computed: {

            radial_data()
            /*  Generate and return a string for the "-d" parameter value of the SVG <path> tag.
                EXAMPLE (abridged):   "M0,-200L70.7,-70.7 [...] L0,-200"
                Note: this function could alternatively go under methods - in that case,
                      the HTML would have to be  v-bind:d="radial_data()" , with parentheses
             */
            {
                const radial_line_func = d3.lineRadial();      // This will be a function
                return radial_line_func(this.radial_points());
                // Note: if radial_points had been made as computed rather than a method,
                //       then the call would have been:  radial_line_func(this.radial_points)
            }

        },

        // ---------------------------  METHODS  ---------------------------
        methods: {

            inner_radius() {
                return this.outer_radius * 0.5;
            },

            radial_points()
            /* Generate and return a list of radial coordinates; each one a [Theta, r] pair.
                EXAMPLE: [
                           [0, 300],
                           [0.6283185307179586, 150],
                           [1.2566370614359172, 300],
                            ...
                           [6.283185307179586, 300]
                         ]
             */
            {
                const step = (2 * Math.PI) / (this.number_rays * 2); // The incremental angular step
                const points = [];
                for (let i = 0; i <= this.number_rays * 2; i++) {
                    const current_radius = (i % 2) ? this.inner_radius() : this.outer_radius;
                    points.push([i * step, current_radius]);
                }
                //console.log(points);
                return points;
            }
        }  // METHODS

    }
); // end component