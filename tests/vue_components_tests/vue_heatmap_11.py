# Test of using the graphic module (Vue component) "vue-heatmap-11" in an HTML log

import pathlib
from life123 import HtmlLog as log
from life123 import GraphicLog

COMPONENT_NAME = "vue_heatmap_11"           # CHANGE AS NEEDED


LOG_FILENAME = pathlib.Path(__file__).stem + ".log.htm"
# EXAMPLE:  if this script were named "sample.py", then the log will go to "sample.log.htm" (SAME folder)
#           (Note that this will NOT work in a Jupyter notebook)


# Initialize the HTML logging.
# Note: the relative path is from the location of THE LOG FILE to the project's home
GraphicLog.config(filename=LOG_FILENAME,
                  components=COMPONENT_NAME,
                  local_files=True, home_rel_path="../..")


log.write(f'Example of Heatmap, using the module "{COMPONENT_NAME}":',
          style=log.h2, blanks_after=1)


# EXAMPLE OF DATA VISUALIZATION (with arbitrary data, hardwired below)

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
GraphicLog.export_plot(all_data, COMPONENT_NAME, unpack=True)
