"""
This is a test of the heatmap module
"""

from modules.html_log.html_log import HtmlLog as log

log.config(filename="heatmap_test_1.htm", mode='overwrite',
           use_D3=True,
           Vue_lib = "../../../modules/Vue2_lib/vue2.js",
           js = "../../../modules/SVG_helper/svg_helper.js",
           css="../../../modules/visualization/D3_heatmap.css")
# Note: paths are from the location of THE LOG FILE




# EXAMPLE OF DATA VISUALIZATION WITH A HEATMAP (with arbitrary data, hardwired below)

all_data = {
    "my_groups": ["A", "B", "C"],
    "my_vars": ["v1", "v2"],
    "my_data": [
        { "group": "A", "variable": "v1", "value": "0" },
        { "group": "A", "variable": "v2", "value": "82" },
        { "group": "B", "variable": "v1", "value": "37" },
        { "group": "B", "variable": "v2", "value": "50" },
        { "group": "C", "variable": "v1", "value": "100" },
        { "group": "C", "variable": "v2", "value": "13" }
    ],
    "range_min": 0,
    "range_max": 100,
    "outer_width": 850,
    "outer_height": 450,
    "margins": {"top": 30, "right": 30, "bottom": 30, "left": 30}
}

log.write('Example of Heatmap, using the module "vue-heatmap-9"')

log.export_plot_Vue(data=all_data,
                    component_name="vue-heatmap-9",
                    component_file="../../modules/visualization/vue_components/heatmap9.js")
