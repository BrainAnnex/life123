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
# ## One-bin `2A + 5B <-> 4C + 3D`, with 1st-order kinetics for each species, taken to equilibrium
#
# Diffusion not applicable (just 1 bin)
#
# LAST REVISED: May 28, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.life_1D.bio_sim_1d import BioSim1D

import plotly.express as px
from src.modules.html_log.html_log import HtmlLog as log
from src.modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system
chem_data = chem(names=["A", "B", "C", "D"])     # NOTE: Diffusion not applicable (just 1 bin)

# Specify the reaction

# Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
chem_data.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                       forward_rate=5., reverse_rate=2.)

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_all_uniform_concentrations( [4., 7., 5., 2.] )

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0), caption="Initial state")
bio.get_history()

# %%
chem_data.describe_reactions()

# %%
# Send a header and a plot to the HTML log file
log.write("Reaction 2 A + 5 B <-> 4 C + 3 D",
          style=log.h2)
graph_data = chem_data.prepare_graph_network()
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
bio.add_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.001, n_steps=40, snapshots={"sample_bin": 0})

bio.describe_state()

# %% [markdown] tags=[]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:  
# [A] = 2.80284552 , [B] = 4.00711381 , [C] = 7.39430896 , [D] = 3.79573172

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
df = bio.get_history()
df

# %% [markdown]
# A and B get depleted, while C and D get produced.
#
# **2A + 5B <-> 4C + 3D**

# %% [markdown]
# #### Let's verify that the stoichiometry is being respected

# %%
# We'll check the first two arrays of concentrations, from the run's history
arr0 = bio.reaction_dynamics.get_historical_concentrations(row=0, df=df)
arr1 = bio.reaction_dynamics.get_historical_concentrations(row=1, df=df)
arr0, arr1

# %%
bio.reaction_dynamics.stoichiometry_checker(rxn_index=0, 
                               conc_arr_before = arr0, 
                               conc_arr_after = arr1)

# %% [markdown]
# Indeed, the change in [A] is -2 x 0.12, and the change in [B] is -5 X 0.12,  
#   while the change in [C] is  4 x 0.12, and the change in [D] is  3 X 0.12

# %%
(arr1 - arr0) / 0.12

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %% tags=[]
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "C", "D"], 
              title="Changes in concentrations with time",
              color_discrete_sequence = ['navy', 'cyan', 'red', 'orange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
