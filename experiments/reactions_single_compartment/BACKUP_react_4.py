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
# # BACKUP
# ## Association/Dissociation reaction 2A <-> C
# #### with 2nd-order kinetics for A,  
# #### and 1-st order kinetics for C
#
# Taken to equilibrium.  (Adaptive variable time resolution is used)
#
# _See also the experiment "1D/reactions/reaction_7"_ 
#
#
# LAST REVISED: Jan. 20, 2023

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
chem_data = chem(names=["A", "C"])

# Reaction 2A <-> C , with 2nd-order kinetics for A, and 1st-order kinetics for C
chem_data.add_reaction(reactants=[(2, "A", 2)], products=["C"],
                       forward_rate=3., reverse_rate=2.)   
# Note: the first 2 in (2, "A", 2) is the stoichiometry coefficient, while the other one is the order

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
dynamics.set_conc([200., 40.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()
#dynamics.verbose_list = [1]      # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

# The changes of concentrations vary very rapidly early on; so, we'll be using dynamic_steps=4 , i.e. increase time resolution
# by x4 initially, as long as the reaction remains "fast" (based on a threshold of 5% change)
dynamics.single_compartment_react(time_step=0.002, reaction_duration=0.04,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_steps=4)
                                  # Accepting the default:  fast_threshold=5

# %% [markdown]
# ### Note: the argument _dynamic_step=4_ splits the time steps in 4 whenever the reaction is "fast" (as determined using fast_threshold=5)

# %%
df = dynamics.history.get()
df

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Notice how the reaction proceeds in smaller steps in the early times, when the concentrations are changing much more rapidly

# %%
# Let's look at the first two arrays of concentrations, from the run's history
arr0 = dynamics.get_historical_concentrations(0)
arr1 = dynamics.get_historical_concentrations(1)
arr0, arr1

# %%
# Let's verify that the stoichiometry is being respected
dynamics.stoichiometry_checker(rxn_index=0, 
                               conc_arr_before = arr0, 
                               conc_arr_after = arr1)

# %% [markdown]
# #### Indeed, it can be easy checked that the drop in [A] is -119.920000 , twice the 59.96 increase in [C], as dictated by the stoichiometry

# %%
dynamics.diagnostic_data[0].get().loc[0]    # Conveniently seen in the diagnostic data

# %% [markdown]
# ## Note: "A" (now largely depleted) is the limiting reagent

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "C"], 
              title="Reaction 2A <-> C  (2nd order in A).  Changes in concentrations with time",
              color_discrete_sequence = ['red', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# #### For diagnostic insight, uncomment the following lines:

# %%
#dynamics.examine_run(df=df, time_step=0.002)  
# the time step MUST match the value used in call to single_compartment_react()

#dynamics.diagnose_variable_time_steps()

#dynamics.diagnostic_data[0].get()

#dynamics.diagnostic_data_baselines.get()

# %%