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
# ## 2 COUPLED reactions of different speeds:  
# ### A <-> B (fast) and B <-> C (slow)
# All 1st order. Taken to equilibrium. Both reactions are mostly forward.
# The concentration of the intermediate product B manifests 1 oscillation (transient "overshoot")
#
# (Adaptive variable time resolution is used)
#
# LAST REVISED: Jan. 9, 2023

# %% [markdown]
# ## Bathtub analogy:
# A is initially full, while B and C are empty.  
# Tubs are progressively lower (reactions are mostly forward.)  
# A BIG pipe connects A and B: fast kinetics.  A small pipe connects B and C: slow kinetics. 
#
# INTUITION: B, unable to quickly drain into C while at the same time being blasted by a hefty inflow from A,  
# will experience a transient surge, in excess of its final equilibrium level.
#
# * [Compare with the final reaction plot (the red line is B)](#react_5_plot)

# %% [markdown]
# ![2 Coupled Reactions](../../docs/2_coupled_reactions.png)

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

# Reaction A <-> B (fast)
chem_data.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=64., reverse_rate=8.) 

# Reaction B <-> C (slow)
chem_data.add_reaction(reactants=["B"], products=["C"],
                       forward_rate=12., reverse_rate=2.) 

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

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics.set_conc([50., 0, 0.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()
#dynamics.verbose_list = [1, 2, 3]      # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

# The changes of concentrations vary very rapidly early on; so, we'll be using dynamic_step=4 , i.e. increase time resolution
# by x4 initially, as long as the reaction remains "fast" (based on a threshold of 5% change)
dynamics.single_compartment_react(time_step=0.02, reaction_duration=0.4,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_step=10, fast_threshold=15)      

# %% [markdown]
# ### Note: the argument _dynamic_step=10_ splits the time steps in 10 for any reactions that are "fast-changing" (as determined using _fast_threshold=15_ )

# %%
df = dynamics.history.get()
df

# %%
# Let's expand the last part
df.loc[60:]

# %% [markdown]
# ### Notice:
# * the reaction proceeds in smaller steps in the earlier times (until t=0.160, in line 80), when the concentrations are changing much more rapidly 
#
# * between lines 70 and 80, only rection #1 is regarded as fast-changing (based on the fast_threshold we specified in the _simulation run_); previously, both reactions were regarded as fast-changing
#
# * "fast-changing" and "slow-changing" is NOT the same thing as "fast" and "slow" reaction kinetics.  For example, reaction #1, though it has much slower kinetics than reaction #0, involves large relative concentration changes because [C] is small
#
# * after step 80, both reactions are regarded as slow-changing, and no more intermediate steps are used

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium(tolerance=0.2)

# %% [markdown] tags=[]
# ## <a name="react_5_plot"></a>Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "B", "C"], 
              title="Coupled reactions A <-> B and B <-> C",
              color_discrete_sequence = ['blue', 'red', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%

# %%

# %% [markdown]
# # EVERYTHING BELOW IS FOR DIAGNOSTIC INSIGHT
#
# ### Perform some verification

# %%
# Verify that the stoichiometry is respected at every reaction step/substep (NOTE: it requires earlier activation of saving diagnostic data)
dynamics.stoichiometry_checker_entire_run()

# %% [markdown]
# ### Take a peek at the diagnostic data saved during the earlier reaction simulation

# %%
# This table contains a subset of a typical system history
dynamics.diagnostic_data_baselines.get()

# %%
# Concentration increments due to reaction 0 (A <-> B)
# Note that [C] is not affected
dynamics.get_diagnostic_data(rxn_index=0)

# %%
# Expand the last part of the above table
dynamics.get_diagnostic_data(rxn_index=0).loc[60:]

# %%
# Concentration increments due to reaction 1 (B <-> C)
# Note that [A] is not affected
dynamics.get_diagnostic_data(rxn_index=1)

# %%
# Expand the last part of the above table
dynamics.get_diagnostic_data(rxn_index=1).loc[60:]

# %%

# %% [markdown]
# ### Provide a detailed explanation of all the steps/substeps of the reactions, from the saved diagnostic data

# %%
dynamics.explain_reactions()

# %%
