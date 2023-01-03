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
# ## Association/Dissociation reaction A + B <-> C
# #### with 1st-order kinetics for each species, taken to equilibrium.
# #### Exploration of debugging and diagnostics options
# (Adaptive variable time resolution is used)
#
# _See also the experiment "1D/reactions/reaction_4"_ 
#
# LAST REVISED: Jan. 1, 2023

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from modules.reactions.reaction_data import ReactionData as chem
from modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.express as px
from modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # Initialize the System
# Specify the chemicals and the reactions

# %% tags=[]
# Specify the chemicals
chem_data = chem(names=["A", "B", "C"])

# Reaction A + B <-> C , with 1st-order kinetics for each species
chem_data.add_reaction(reactants=["A" , "B"], products=["C"],
                       forward_rate=5., reverse_rate=2.)

print("Number of reactions: ", chem_data.number_of_reactions())

# %%
chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# # Start the simulation

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)

# %%
# Initial concentrations of all the chemicals, in index order
dynamics.set_conc([10., 50., 20.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()
#dynamics.verbose_list = [1]      # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

dynamics.single_compartment_react(time_step=0.004, reaction_duration=0.06,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_step=2)      
                                  # Accepting the default:  fast_threshold=5

# %% [markdown]
# ### Note: the argument _dynamic_step=2_ splits the time steps in 2 whenever the reaction is "fast" (as determined using fast_threshold=5)

# %%
df = dynamics.history.get()
df

# %% [markdown]
# ### Notice how the reaction proceeds in smaller steps in the early times, when the concentrations are changing much more rapidly

# %% [markdown]
# ## Note: "A" (now largely depleted) is the limiting reagent

# %% [markdown]
# ### Check the final equilibrium

# %%
dynamics.get_system_conc()

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order reactions),  
#  the systems settles in the following equilibrium:
#
# [A] = 0.29487741 , [B] = 40.29487741 , [C] = 29.70512259
#

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(rxn_index=0, conc=dynamics.get_conc_dict())

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "B", "C"], 
              title="Reaction A + B <-> C .  Changes in concentrations with time",
              color_discrete_sequence = ['red', 'violet', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%

# %% [markdown]
# # Everthing below is just for diagnostic insight 
# # into the adaptive variable time steps

# %%
# This approach, from the run data, is only usable with single-reaction runs
dynamics.examine_run(df=df, time_step=0.004)
# the time step MUST match the value used in call to single_compartment_react()

# %% [markdown]
# # Take a peek at internal diagnostic data from the reactions

# %%
# This approach, from internal diagnostic data, 
# is more generally applicable also to runs with multiple reactions
dynamics.diagnose_variable_time_steps()

# %% [markdown]
# ### The above diagnostics are based on the following diagnostic data

# %%
dynamics.diagnostic_data.get()

# %%
dynamics.diagnostic_data_baselines.get()

# %%
