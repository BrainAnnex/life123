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
# ## Association/Dissociation reaction `A + B <-> C`
# #### with 1st-order kinetics for each species, taken to equilibrium.
# #### Exploration of debugging and diagnostics options
# (Adaptive variable time steps are used)
#
# _See also the experiment "1D/reactions/reaction_4"_  
#
# LAST REVISED: May 18, 2023 
#
# # For the 0-th reaction (the only reaction in our case)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.express as px
from src.modules.visualization.graphic_log import GraphicLog

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
dynamics.get_history()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these are the currently the default values, but subject to change
dynamics.set_thresholds(thresholds={"low": 0.5, "high": 0.8, "abort": 1.44, "reduction_factor": 2.})

dynamics.single_compartment_react(initial_step=0.004, reaction_duration=0.06,
                                  variable_steps=True,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"})

# %% [markdown]
# ### Note: the argument _dynamic_step=2_ splits the time steps in 2 whenever the reaction is "fast" (as determined using the given value of _fast_threshold_ )

# %%
df = dynamics.get_history()
df

# %%
dynamics.explain_time_advance()

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
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=['red', 'violet', 'green'], show_intervals=True)

# %%

# %% [markdown]
# # Everthing below is just for diagnostic insight 
# ## into the adaptive variable time steps

# %%
dynamics.get_diagnostic_decisions_data()

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)

# %%
dynamics.get_diagnostic_conc_data()

# %%
dynamics.examine_run(df=df)

# %%
