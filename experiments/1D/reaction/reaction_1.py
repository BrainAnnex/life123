# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# ### A simple A <-> B reaction between 2 species with initial uniform concentrations across 3 bins,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# Diffusion NOT taken into account
#
# LAST REVISED: Aug. 22, 2022
#
# * [First Step](#sec_1)

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

import plotly.express as px
from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% tags=[]
# Initialize the system
chem_data = chem(names=["A", "B"])       # Diffusion NOT taken into account
bio.initialize_system(n_bins=3, chem_data=chem_data)   # We'll specify the reactions later

bio.set_uniform_concentration(species_name="A", conc=10.)
bio.set_uniform_concentration(species_name="B", conc=50.)

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))

# %%
bio.get_history()

# %%
# Specify the reaction
rxn = Reactions(chem_data)

# Reaction A <-> B , with 1st-order kinetics in both directions
rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

bio.set_reactions(rxn)

print("Number of reactions: ", rxn.number_of_reactions())

# %%
rxn.describe_reactions()

# %%
# Send the plot to the HTML log file
graph_data = rxn.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown] tags=[]
# ### <a name="sec_1"></a>First step

# %%
# First step of reaction
bio.react(time_step=0.1, n_steps=1)
bio.describe_state()

# %% [markdown]
# NOTE: the concentration of species A is increasing, while that of species B is decreasing.
# All bins have identical concentrations; so, there's no diffusion (and we're not attempting to compute it):
# [[17. 17. 17.]
#  [43. 43. 43.]]

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# ### <a name="sec_more_steps"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.1, n_steps=10)

bio.describe_state()

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order reactions),
#  the systems settles in the following equilibrium:
#
# [A] = 23.99316406
#  
# [B] = 36.00683594
#

# %%
bio.bin_snapshot(bin_address = 0)

# %%
# Verify that the reaction has reached equilibrium
A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_back_rate(0)}")
print(f"Ratio of equilibrium concentrations: {B_eq / A_eq}")

# %%
# A handy way to do the above
rxn.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['navy', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
# Same plot, but with smooth line
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['navy', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"},
              line_shape="spline")
fig.show()

# %%