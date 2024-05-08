"""
Test of using the graphic module (Vue component) "vue_cytoscape_2" inside an HTML log
(same filename stem as this script, with ".log.htm" ending)

This test is a simple A <-> B reaction, using high-level functions to define the reaction,
and to send the diagram of the reaction network to the log file.
"""


import pathlib
from src.modules.html_log.html_log import HtmlLog as log
from src.modules.visualization.graphic_log import GraphicLog
from src.modules.reactions.reaction_dynamics import ReactionDynamics


COMPONENT_NAME = "vue_cytoscape_2"           # *** CHANGE AS NEEDED ***

LOG_FILENAME = pathlib.Path(__file__).stem + ".log.htm"
# EXAMPLE:  if this script were named "sample.py", then the log will go to "sample.log.htm"
#           (Note that this will NOT work in a Jupyter notebook)


# Initialize the HTML logging.
# Note: the relative path is from the location of *THE LOG FILE* to the project's home
GraphicLog.config(filename=LOG_FILENAME,
                  components=COMPONENT_NAME,
                  local_files=True, home_rel_path="../..",
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")


log.write(f'Reaction Network Diagram, using the Vue module "{COMPONENT_NAME}"',
          style=log.h3, newline=False)



# Instantiate the simulator and specify the chemicals
dynamics = ReactionDynamics(names=["A", "B"])

# Reaction A <-> B , with 1st-order kinetics in both directions
dynamics.add_reaction(reactants=["A"], products=["B"],
                      forward_rate=3., reverse_rate=2.)


# Send the plot to the HTML log file
dynamics.plot_reaction_network("vue_cytoscape_2")
