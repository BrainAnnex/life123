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
# ## One-bin  2A <-> B reaction, COMPARING 1st-order and 2nd-order kinetics in *forward* direction; reverse direction 1-st order
#
# Diffusion not applicable (just 1 bin)
#
# LAST REVISED: Aug. 22, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %%
from experiments.get_notebook_info import get_notebook_basename

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

import plotly.express as px
from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # INITIALLY, with 1st-order kinetics in both directions

# %%
# Initialize the system
chem_data = chem(names=["A", "B"])     # NOTE: Diffusion not applicable (just 1 bin)

rxn = Reactions(chem_data)

# Reaction  2A <-> B , FOR NOW with 1st-order kinetics in both directions
rxn.add_reaction(reactants=[(2, "A")], products=["B"], forward_rate=5., reverse_rate=2.)

bio.initialize_system(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_all_uniform_concentrations( [3., 5.] )

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %%
rxn.describe_reactions()

# %%
rxn._internal_reactions_data()    # Low-level view of the reactions data

# %%
# Send a header and a plot to the HTML log file
log.write("Reaction 2A <-> B is 1st order in all species:",
          style=log.h2)

graph_data = rxn.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %%
# First step
bio.react(time_step=0.02, n_steps=1)
bio.describe_state()

# %% [markdown]
# Small conc. changes so far:  [A] = 2.8 , [B] = 5.1

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %%
# Numerous more steps, to equilibrium
bio.react(time_step=0.02, n_steps=20)

bio.describe_state()

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the *1st order* reactions),
# the systems settles in the following equilibrium:  
# [A] = 2.16928427 , [B] = 5.41535786

# %%
# Verify that the reaction has reached equilibrium
rxn.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% tags=[]
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="2A <-> B : changes in concentrations",
              color_discrete_sequence = ['navy', 'orange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# A gets depleted, while B gets produced.

# %% [markdown]
# # STARTING OVER, this time with 2nd-order kinetics in the forward reaction

# %%
rxn.clear_reactions()

# %%
# Reaction  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction
rxn.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)

# %%
# RESET the concentrations to their original values
bio.set_all_uniform_concentrations( [3., 5.] )

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0), 
                  caption = "RESET all concentrations to initial values")
bio.get_history()

# %%
rxn.describe_reactions()

# %%
rxn._internal_reactions_data()    # Low-level view of the reactions data

# %%
# Send a header and a plot to the HTML log file
log.write("Reaction 2A <-> B is 2nd order in A, and 1st order in B:",
          style=log.h2)
graph_data = rxn.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %%
# First step
bio.react(time_step=0.02, n_steps=1)
bio.describe_state()

# %% [markdown]
# [A] = 1.6 , [B] = 5.7
# _(Contrast with the counterpart in the 1st order kinetics:  [A] = 2.8 , [B] = 5.1)_

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %%
# Numerous more steps
bio.react(time_step=0.02, n_steps=20)

bio.describe_state()

# %% [markdown]
# The systems settles in the following equilibrium:  [A] = 1.51554944 , [B] = 5.74222528

# %%
# Verify that the reaction has reached equilibrium
rxn.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="2A <-> B : changes in concentrations",
              color_discrete_sequence = ['navy', 'orange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# Compared to first-order kinetics in A, the reaction now takes place much more quickly, and proceeds to almost complete depletion of A

# %%
