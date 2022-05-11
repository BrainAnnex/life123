# Test 1 of using the heatmap module "vue-heatmap-11" in an HTML log

from modules.html_log.html_log import HtmlLog as log

# Note: paths are from the location of THE LOG FILE
log.config(filename="vue-heatmap-11_1.htm", mode='overwrite',
           use_D3=True,
           Vue_lib = "../../../modules/Vue2_lib/vue2.js",
           js = "../../../modules/SVG_helper/svg_helper.js",
           css="../../../modules/visualization/vue_components/heatmap11.css")


log.write('Example of Heatmap, using the module "vue-heatmap-11":',
          style=log.h2, blanks_after=1)


# EXAMPLE OF DATA VISUALIZATION WITH A HEATMAP (with arbitrary data, hardwired below)

all_data = {
    # Labels of rows
    "y_labels": ["Chem 1", "Chem 2"],

    # Data for the heatmap, by rows, starting with the bottom one
    "heatmap_data": [
        [0, 85, 100],
        [100, 20, 0]
    ],

    # Set the range of values in the heatmap bins
    "range_min": 0,
    "range_max": 100,

    # Set the dimensions and margins of the heatmap
    "outer_width": 900,
    "outer_height": 400,
    "margins": {"top": 30, "right": 30, "bottom": 30, "left": 85}
}

log.export_plot_Vue(data=all_data,
                    component_name="vue-heatmap-11",
                    component_file="../../../modules/visualization/vue_components/heatmap11.js")
