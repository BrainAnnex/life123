# Test "A" of using the heatmap module "vue_curves_1" in an HTML log

from src.modules.html_log.html_log import HtmlLog as log
from src.modules.visualization.graphic_log import GraphicLog

COMPONENT_NAME = "vue_curves_1"

# Initialize the HTML logging.
# Note: the relative path is from the location of THE LOG FILE to the project's home
GraphicLog.config(filename="vue_curves_1_A.htm",
                  components=COMPONENT_NAME,
                  home_rel_path="../../../..")


log.write(f'Example of Line plots, using the module "{COMPONENT_NAME}":',
          style=log.h2, blanks_after=1)


# EXAMPLE OF DATA VISUALIZATION WITH A HEATMAP (with arbitrary data, hardwired below)

all_data = {
    # Labels for the rows
    "y_labels": ["v1", "very_long_label"],

    # Concentration data for the plots (for now just 1 chemical species), in index order
    "data": [20, 85, 100, 50],

    # Set the range of values (not sure if being used)
    "range_min": 0,
    "range_max": 100,

    # Set the dimensions and margins of the plot ("outer" means the surrounding box)
    "outer_width": 900,
    "outer_height": 400,
    "margins": {"top": 20, "right": 20, "bottom": 20, "left": 30}
}


# Send the plot to the HTML log file
GraphicLog.export_plot(all_data, COMPONENT_NAME)
