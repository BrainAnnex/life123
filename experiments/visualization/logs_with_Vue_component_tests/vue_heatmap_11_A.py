# Test 1 of using the heatmap module "vue-heatmap-11" in an HTML log

from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

COMPONENT_NAME = "vue_heatmap_11"       # CHANGE AS NEEDED
FILENAME="vue_heatmap_11_A.htm"         # CHANGE AS NEEDED


# Initialize the HTML logging.
# Note: the relative path is from the location of THE LOG FILE to the project's home
GraphicLog.config(filename=FILENAME,
                  components=COMPONENT_NAME,
                  home_rel_path="../../..")


log.write('Example of Heatmap, using the module "{COMPONENT_NAME}":',
          style=log.h2, blanks_after=1)


# EXAMPLE OF DATA VISUALIZATION WITH A HEATMAP (with arbitrary data, hardwired below)

all_data = {
    # Labels for the rows
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

# Send the plot to the HTML log file
GraphicLog.export_plot(all_data, COMPONENT_NAME)
