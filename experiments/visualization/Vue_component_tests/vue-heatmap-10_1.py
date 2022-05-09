# Test of the heatmap module "vue-heatmap-10"

from modules.html_log.html_log import HtmlLog as log

# Note: paths are from the location of THE LOG FILE
log.config(filename="vue-heatmap-10_1.htm", mode='overwrite',
           use_D3=True,
           Vue_lib = "../../../modules/Vue2_lib/vue2.js",
           js = "../../../modules/SVG_helper/svg_helper.js",
           css="../../../modules/visualization/vue_components/heatmap10.css")


log.write('Example of Heatmap, using the module "vue-heatmap-10":',
          style=log.h2, blanks_after=1)


# EXAMPLE OF DATA VISUALIZATION WITH A HEATMAP (with arbitrary data, hardwired below)

all_data = {
    "x_labels": ["C0", "C1", "C2"],
    "y_labels": ["Chem 1", "Chem 2"],
    "heatmap_data": [
        [0, 85, 100],
        [100, 20, 0]
    ],
    "range_min": 0,
    "range_max": 100,
    "outer_width": 900,
    "outer_height": 400,
    "margins": {"top": 30, "right": 30, "bottom": 30, "left": 50}
}

log.export_plot_Vue(data=all_data,
                    component_name="vue-heatmap-10",
                    component_file="../../../modules/visualization/vue_components/heatmap10.js")
