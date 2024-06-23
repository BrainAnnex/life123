# Test of using the graphic module (Vue component)  "vue_curves_3" in an HTML log

import pathlib
from life123 import HtmlLog as log
from life123 import GraphicLog

COMPONENT_NAME = "vue_curves_3"         # CHANGE AS NEEDED


LOG_FILENAME = pathlib.Path(__file__).stem + ".log.htm"
# EXAMPLE:  if this script were named "sample.py", then the log will go to "sample.log.htm" (SAME folder)
#           (Note that this will NOT work in a Jupyter notebook)


# Initialize the HTML logging.
# Note: the relative path is from the location of THE LOG FILE to the project's home
GraphicLog.config(filename=LOG_FILENAME,
                  components=COMPONENT_NAME,
                  local_files=True, home_rel_path="../..")


log.write(f'Example of Line plots, using the module "{COMPONENT_NAME}":',
          style=log.h2, blanks_after=1)


# EXAMPLE OF DATA VISUALIZATION (with arbitrary data, hardwired below)

all_data = {         # CHANGE AS NEEDED
    # Labels for the y-axis
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
GraphicLog.export_plot(all_data, COMPONENT_NAME, unpack=True)
