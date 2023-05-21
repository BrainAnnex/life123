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
# ## Association/Dissociation reaction `2A <-> C`
# #### with 2nd-order kinetics for `A`,  
# #### and 1-st order kinetics for `C`
#
# Taken to equilibrium.  (Adaptive variable time teps are used)
#
# _See also the experiment "1D/reactions/reaction_7"_ 
#
# LAST REVISED: May 21, 2023

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
chem_data = chem(names=["A", "C"])

# Reaction 2A <-> C , with 2nd-order kinetics for A, and 1st-order kinetics for C
chem_data.add_reaction(reactants=[(2, "A", 2)], products=["C"],
                       forward_rate=3., reverse_rate=2.)   
# Note: the first 2 in (2, "A", 2) is the stoichiometry coefficient, while the other one is the order

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
dynamics.get_history()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are the currently the default values... but subject to change
dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
dynamics.set_step_factors(abort=0.5, downshift=0.5, upshift=1.5)

dynamics.single_compartment_react(time_step=0.002, reaction_duration=0.04,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True, explain_variable_steps=True)

# %% [markdown]
# ### Note: the argument _dynamic_step=4_ splits the time steps in 4 whenever the reaction is "fast" (as determined using the given value of _fast_threshold_ )

# %%
df = dynamics.get_history()
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
dynamics.get_diagnostic_rxn_data(rxn_index=0, head=1)    # Conveniently seen in the diagnostic data

# %% [markdown]
# ## Note: "A" (now largely depleted) is the limiting reagent

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=2)

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=['red', 'green'],
                     title="Reaction 2A <-> C  (2nd order in A).  Changes in concentrations with time")

# %% [markdown]
# #### For diagnostic insight, uncomment the following lines:

# %%
#dynamics.get_diagnostic_decisions_data()

#dynamics.get_diagnostic_rxn_data(rxn_index=0)

#dynamics.get_diagnostic_conc_data()

# %%
