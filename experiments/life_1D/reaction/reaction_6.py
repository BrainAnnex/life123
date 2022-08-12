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
# **One-bin 2A + 5B <-> 4C + 3D, with 1st-order kinetics for each species,
# taken to equilibrium**
#
# Diffusion not applicable (just 1 bin)
#
# LAST REVISED: Aug. 11, 2022

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
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js",
                  home_rel_path="../../..")    # relative path is from the location of THE LOG FILE to the project's home

# %%
# Initialize the system
chem_data = chem(names=["A", "B", "C", "D"])     # NOTE: Diffusion not applicable (just 1 bin)

# Specify the reaction
rxn = Reactions(chem_data)

# Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
rxn.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                 forward_rate=5., reverse_rate=2.)

bio.initialize_system(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_all_uniform_concentrations( [4., 7., 5., 2.] )

bio.describe_state()

# %%
rxn.describe_reactions()

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %%
# Send the plot to the HTML log file
graph_data = rxn.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown] tags=[]
# ### <a name="sec_2_first_step"></a>First step

# %%
# First step
bio.react(time_step=0.001, n_steps=1)
bio.describe_state()

# %% [markdown]
# _Early in the reaction :_
# [A] = 3.76 ,  [B] = 6.4 ,  [C] = 5.48 ,  [D] = 2.36

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.001, n_steps=40)

bio.describe_state()

# %% [markdown] tags=[]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:  
# [A] = 2.80284552 , [B] = 4.00711381 , [C] = 7.39430896 , [D] = 3.79573172

# %%
A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
D_eq = bio.bin_concentration(0, 3)
print(f"Ratio of equilibrium concentrations ((C_eq * D_eq) / (A_eq * B_eq)) : {(C_eq * D_eq) / (A_eq * B_eq)}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# A and B get depleted, while C and D get produced.
#
# **2A + 5B <-> 4C + 3D**

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %% tags=[]
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "C", "D"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['navy', 'cyan', 'red', 'orange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
# Same plot, but with smooth line
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "C", "D"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['navy', 'cyan', 'red', 'orange'],
              labels={"value":"concentration", "variable":"Chemical"},
              line_shape="spline")
fig.show()

# %%
