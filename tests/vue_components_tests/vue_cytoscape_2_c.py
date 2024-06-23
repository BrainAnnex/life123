# Similar to "vue_cytoscape_2_a", but this time with more reactions, and diffusion rates

import pathlib
from life123 import HtmlLog as log
from life123 import GraphicLog, ChemData, UniformCompartment

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



# Instantiate the simulator and specify the chemicals
chem_data = ChemData(names=["A", "B", "C", "D", "E"],
                     diffusion_rates=[0.1, 0.2, 0.3, 0.4, 0.5])

dynamics = UniformCompartment(chem_data=chem_data)


# Reaction A + 2B <-> C , with 2nd-order kinetics in B
dynamics.add_reaction(reactants=["A", (2, "B", 2)], products=["C"],
                      forward_rate=3., reverse_rate=2.)

# Reaction B + D <-> E , with 1st-order kinetics in all terms
dynamics.add_reaction(reactants=["B", "D"], products=["E"],
                      forward_rate=1., reverse_rate=5.)

# Send the plot to the HTML log file
dynamics.plot_reaction_network("vue_cytoscape_2")
