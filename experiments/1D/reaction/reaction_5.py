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
# ## One-bin A <-> 2C + D, with 1st-order kinetics for each species, taken to equilibrium
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

# %%
# Initialize the system
chem_data = chem(names=["A", "C", "D"])     # NOTE: Diffusion not applicable (just 1 bin)

rxn = Reactions(chem_data)

# Reaction A <-> 2C + D , with 1st-order kinetics for each species
rxn.add_reaction(reactants=[("A")], products=[(2, "C") , ("D")],
                 forward_rate=5., reverse_rate=2.)

bio.initialize_system(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_all_uniform_concentrations( [4., 7., 2.] )

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

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
bio.react(time_step=0.2, n_steps=1)
bio.describe_state()

# %% [markdown]
# ---   
#     Note:  the above values are quite inaccurate because of the large time step 0.2
#
#     For example, the value for the concentration of D (0.4) is a wild overshot from the initial 2.0 to the equilibrium value of 1.68941267
#        
#     A more precise calculation with bio.react(time_step=0.1, n_steps=2) gives conc_D(0.2) = 2.304
#        
#     An even more precise calculation with bio.react(time_step=0.05, n_steps=4) gives conc_D(0.2) = 1.69037202
#        
#     I.e. the system is almost at equilibrium already at t=0.2 !
#     
#     TODO: explore the early dynamics of the system in a separate experiment
# ---

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.05, n_steps=30)

bio.describe_state()

# %% [markdown] tags=[]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:  
# [A] = 4.31058733 , [C] = 6.37882534 , [D] = 1.68941267

# %%
# Verify that the reaction has reached equilibrium
rxn.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# C and D get depleted, while A gets produced.
# A wild overshoot is present at t=0.2

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %% tags=[]
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "C", "D"], 
              title="Reaction A <-> 2C + D .  Changes in concentrations",
              color_discrete_sequence = ['navy', 'violet', 'red'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# #### Notice the wild overshoot present at t=0.2

# %%