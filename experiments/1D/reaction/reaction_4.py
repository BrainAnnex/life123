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
# ## One-bin Association/Dissociation reaction A + B <-> C
# ### with 1st-order kinetics for each species, taken to equilibrium
#
# Diffusion not applicable (just 1 bin)
#
# LAST REVISED: Nov. 28, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %%
from experiments.get_notebook_info import get_notebook_basename

from modules.reactions.reaction_data import ReactionData as chem
from modules.reactions.reaction_dynamics import ReactionDynamics
from life_1D.bio_sim_1d import BioSim1D

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
chem_data = chem(names=["A", "B", "C"])     # NOTE: Diffusion not applicable (just 1 bin)



# Reaction A + B <-> C , with 1st-order kinetics for each species
chem_data.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                 forward_rate=5., reverse_rate=2.)

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)
bio.set_uniform_concentration(species_index=2, conc=20.)

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %%
chem_data.describe_reactions()

# %%
# Send the plot to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown] tags=[]
# ### <a name="sec_2_first_step"></a>First step

# %%
# First step
bio.react(time_step=0.002, n_steps=1)
bio.describe_state()

# %% [markdown]
# _Early in the reaction :_
# [A] = 5.08 ,  [B] = 45.08 ,  [C] = [24.92]

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.002, n_steps=40)

bio.describe_state()

# %% [markdown] tags=[]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:  
# [A] = 0.29487741 , [B] = 40.29487741 , [C] = 29.70512259

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# # Note: "A" (now almost completely depleted) is largely the limiting reagent

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "C"], 
              title="Reaction A + B <-> C .  Changes in concentrations",
              color_discrete_sequence = ['navy', 'violet', 'red'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
