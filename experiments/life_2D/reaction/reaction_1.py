# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## A <-> B reaction, with 1st-order kinetics in both directions,
# ### taken to equilibrium
#
# Diffusion not done
#
# LAST REVISED: Aug. 19, 2022
# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %%
from experiments.get_notebook_info import get_notebook_basename

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_2D.bio_sim_2d import BioSim2D as bio

import plotly.express as px
from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file
# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js",
                  home_rel_path="../../..")    # relative path is from the location of THE LOG FILE to the project's home

# %%
# Initialize the system
chem_data = chem(names=["A", "B"])     # NOTE: Diffusion not done

rxn = Reactions(chem_data)

# Reaction A <-> B , with 1st-order kinetics in both directions
rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

bio.initialize_system(n_bins=(3,4), chem_data=chem_data, reactions=rxn)

bio.set_bin_conc_all_species(bin_x=0, bin_y=0, conc_list=[10.,50.])
bio.set_bin_conc_all_species(bin_x=0, bin_y=1, conc_list=[20.,35.])
bio.set_bin_conc_all_species(bin_x=2, bin_y=3, conc_list=[5.,100.])

bio.describe_state()

# %%
rxn.describe_reactions()

# %%
# Send the plot to the HTML log file
graph_data = rxn.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown] tags=[]
# ### <a name="sec_2_first_step"></a>First step

# %%
# First step
bio.reaction_step(delta_time=0.1)
print("bio.delta_reactions:\n", bio.delta_reactions)

# %%
bio.system += bio.delta_reactions       # Matrix operation to update all the concentrations
bio.system_time += 0.1

bio.describe_state()
