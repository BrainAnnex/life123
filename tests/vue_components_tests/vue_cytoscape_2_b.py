# Same diagram test as "vue_cytoscape_2_a", but using a low-level, hardwired definition of the network diagram

import pathlib
from life123 import HtmlLog as log
from life123 import GraphicLog

COMPONENT_NAME = "vue_cytoscape_2"           # *** CHANGE AS NEEDED ***


LOG_FILENAME = pathlib.Path(__file__).stem + ".log.htm"
# EXAMPLE:  if this script were named "sample.py", then the log will go to "sample.log.htm" (SAME folder)
#           (Note that this will NOT work in a Jupyter notebook)


# Initialize the HTML logging.
# Note: the relative path is from the location of *THE LOG FILE* to the project's home
GraphicLog.config(filename=LOG_FILENAME,
                  components=COMPONENT_NAME,
                  local_files=True, home_rel_path="../..",
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")


log.write(f'Reaction Network Diagram, using the Vue module "{COMPONENT_NAME}"',
          style=log.h3, newline=False)


# Data to define the reaction network
graph_data = {"structure":
                    [
                        {"name": "RXN", "kF": "3", "kR": "2", "delta_G": "-1,005.13", "K": "1.5", "id": "R-0", "labels": ["Reaction"]},
                        {"name": "B", "diff_rate": None, "id": "C-1", "labels": ["Chemical"]},
                        {"name": "produces", "source": "R-0", "target": "C-1", "id": "edge-1", "stoich": 1, "rxn_order": 1},
                        {"name": "A", "diff_rate": None, "id": "C-0", "labels": ["Chemical"]},
                        {"name": "reacts", "source": "C-0", "target": "R-0", "id": "edge-2", "stoich": 1, "rxn_order": 1}
                    ],
                "color_mapping": {"Chemical": "#8DCC92", "Reaction": "#D9C8AD"},
                "caption_mapping": {"Chemical": "name", "Reaction": "name"}
             }

# Send the plot to the HTML log file
GraphicLog.export_plot(graph_data, "vue_cytoscape_2")
