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
# ### Adaptive time steps (variable time resolution) for reaction `A <-> B`,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# Same as the experiment _"react_1"_ , but with adaptive variable time steps
#
# LAST REVISED: May 22, 2023

# %%
import set_path      # Importing this module will add the project's home directory +to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import pandas as pd
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
chem_data = chem(names=["A", "B"])

# Reaction A <-> B , with 1st-order kinetics in both directions
chem_data.add_reaction(reactants=["A"], products=["B"], 
                       forward_rate=3., reverse_rate=2.)

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
dynamics.set_conc([10., 50.])

dynamics.describe_state()

# %%
dynamics.get_history()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
dynamics.set_thresholds(norm="norm_B")  # This has the effect of turning off use of "norm_B"
dynamics.set_step_factors(upshift=2.0, downshift=0.5, abort=0.5)
dynamics.set_error_step_factor(0.5)

# %%
dynamics.thresholds

# %% [markdown]
# #### Note how we UNSET the defaults for "norm_B" (i.e. it won't be used)

# %%
dynamics.single_compartment_react(initial_step=0.1, target_end_time=1.2,
                                  variable_steps=True, explain_variable_steps=False,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"}
                                  )

# %% [markdown]
# ## The flag _variable_steps_ automatically adjusts up or down the time step,  whenever the changes of concentrations are, respectively, "slow" or "fast" (as determined using the specified _thresholds_ )

# %%
dynamics.get_history()   # The system's history, saved during the run of single_compartment_react()

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# ### That resulted from passing the flag _variable_steps=True_ to single_compartment_react()

# %% [markdown] tags=[]
# # Scrutinizing some instances of step-size changes:

# %% [markdown] tags=[]
# ### Detailed Example 1: **going from 0.1375 to 0.1875**    

# %%
lookup = dynamics.get_history(t_start=0.1375, t_end=0.1875)
lookup

# %%
delta_concentrations = dynamics.extract_delta_concentrations(lookup, 7, 8, ['A', 'B'])
delta_concentrations

# %% [markdown]
# As expected by the 1:1 stoichiometry, delta_A = - delta_B
#
# The above values coud also be looked up from the diagnostic data, since we only have 1 reaction:

# %%
rxn_data = dynamics.get_diagnostic_rxn_data(rxn_index=0)

# %%
rxn_data[8:12]

# %%
delta_row = dynamics.get_diagnostic_rxn_data(rxn_index=0, t=0.1375) # Locate the row with the interval's start time
delta_row

# %%
delta_row[["Delta A", "Delta B"]].to_numpy()   # Gives same value as delta_concentrations, above

# %%
# Computes a measure of how large delta_concentrations is, and propose a course of action
dynamics.adjust_speed(delta_concentrations)  

# %% [markdown]
# #### The above conclusion is that the time step is on the "high" side, and should be **HALVED** at the next round : that's because the computed norm is > than the "high" value previously given in the argument to _set_thresholds()_ (but doesn't exceed the "abort" threshold)

# %%
dynamics.thresholds   # Consult the previously-set threshold values

# %% [markdown] tags=[]
# ### Detailed Example 2: **going from 0.1875 to 0.2125**   

# %%
lookup = dynamics.get_history(t_start=0.1875, t_end=0.2125)
lookup

# %%
delta_concentrations = dynamics.extract_delta_concentrations(lookup, 8, 9, ['A', 'B'])
delta_concentrations

# %% [markdown]
# Note how substantially smaller _delta_concentrations_ is, compared to the previous example

# %%
dynamics.adjust_speed(delta_concentrations)

# %% [markdown]
# #### The above conclusion is that the time step is on the "low" side, and should be **DOUBLED** at the next round : that's because the computed norm is < than the "low" value previously given in the argument to _set_thresholds()_

# %%
dynamics.thresholds   # Consult the previously-set threshold values

# %% [markdown]
# # Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=['blue', 'orange'])

# %% [markdown]
# ## Note how the left-hand side of this plot is much smoother than it was in experiment `react_1`, where no adaptive time steps were used!

# %%
dynamics.plot_curves(colors=['blue', 'orange'], show_intervals=True)

# %% [markdown]
# #### Compare the above with the fixed step sizes (always 0.1) of experiment `react_1`

# %%
dynamics.plot_step_sizes(show_intervals=True)

# %% [markdown]
# # Diagnostics of the run may be investigated as follows:  
# _(note - this is possible because we make a call to set_diagnostics() prior to running the simulation)_

# %%
dynamics.get_diagnostic_conc_data()   # This will be complete, even if we only saved part of the history during the run

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)      # For the 0-th reaction (the only reaction in our case)

# %% [markdown]
# ### Note that diagnostic data with the DELTA Concentrations - above and below - also record the values that were considered (but not actually used) during ABORTED steps

# %%
dynamics.get_diagnostic_decisions_data()

# %%
