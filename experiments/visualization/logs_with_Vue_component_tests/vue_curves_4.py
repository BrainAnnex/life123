# Test of using the graphic module (Vue component) "vue_curves_4" in an HTML log

import pathlib
from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

COMPONENT_NAME = "vue_curves_4"         # CHANGE AS NEEDED

LOG_FILENAME = pathlib.Path(__file__).stem + ".log.htm"
# EXAMPLE:  if this script were named "sample.py", then the log will go to "sample.log.htm"
#           (Note that this will NOT work in a Jupyter notebook)


# Initialize the HTML logging.
# Note: the relative path is from the location of THE LOG FILE to the project's home
GraphicLog.config(filename=LOG_FILENAME,
                  components=COMPONENT_NAME,
                  home_rel_path="../../..")


log.write(f'Example of Line plots, using the module "{COMPONENT_NAME}":',
          style=log.h2, blanks_after=1)


# EXAMPLE OF DATA VISUALIZATION (with arbitrary data, hardwired below)

all_data = {         # CHANGE AS NEEDED
    # Labels for the curves
    "curve_labels": ["Chem A", "Hydrogen sulfide", "Element X", "Acetone"],

    # Concentration data for the plots
    #       outer level : order of chemical-species index,
    #       inner level : in bin index order from left to right
    "plot_data": [[20, 85, 100, 50] , [14, 99, 5, 65] , [100, 75, 55, 0] , [0, 15, 45, 90]],

    # Set the range of values
    "range_min": 0,         # (not yet being used)
    "range_max": 100,

    # Set the dimensions and margins of the plot ("outer" means the surrounding box)
    "outer_width": 900,
    "outer_height": 300,
    "margins": {"top": 20, "right": 20, "bottom": 20, "left": 30}
}


# Send the plot to the HTML log file
GraphicLog.export_plot(all_data, COMPONENT_NAME)
