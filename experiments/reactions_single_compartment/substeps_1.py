# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# ## Exploration of errors in the implementation of substeps, for the simulation of the coupled reactions `2 S <-> U` and `S <-> X`   
# Both mostly forward.  1st-order kinetics throughout.   
#
# Same as `variable_steps_1`, but with substeps.
#
# LAST REVISED: Mar. 6, 2023

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.express as px
from src.modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# ### Initialize the system

# %%
# Initialize the system
chem_data = chem(names=["U", "X", "S"])

# Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants=[(2, "S")], products="U",
                       forward_rate=8., reverse_rate=2.)

# Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants="S", products="X",
                       forward_rate=6., reverse_rate=3.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()
dynamics.verbose_list = ["substeps", "variable_steps"]

dynamics.single_compartment_react(time_step=0.01, stop_time=0.3, 
                                  variable_steps=True, thresholds={"low": 0.25, "high": 0.64},
                                  dynamic_substeps=2, abs_fast_threshold=100.)

df = dynamics.get_history()
df

# %%
(transition_times, step_sizes) = dynamics.explain_time_advance(return_times=True)

# %%
np.array(step_sizes)

# %%
np.array(transition_times)    # Note: there will be one more transition time (the end time) than step sizes

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=['green', 'orange', 'blue'])

# %% [markdown]
# # Notice the irregular "bumpiness" on the plot of S.  See comments at very end.

# %%
dynamics.plot_curves(colors=['green', 'orange', 'blue'], show_intervals=True)

# %%
# Show the "critical values", i.e. times when the step size changes
dynamics.plot_curves(colors=['green', 'orange', 'blue'], vertical_lines=transition_times, 
                     title="Critical values of time-step changes for reactions `2 S <-> U` and `S <-> X`")

# %% [markdown]
# ## Note: the dashed lines in the plots immediatly above and below are NOT the steps; they are the "critical values", i.e. times when the step size changes.   
# The time steps were shown in an earlier plots

# %%
dynamics.plot_step_sizes(show_intervals=True)

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=1)

# %%
dynamics.get_diagnostic_conc_data()

# %%
dynamics.get_diagnostic_delta_conc_data()

# %% [markdown]
# #### Notice how the first step got aborted, and re-run, because of the large adjusted L2 value in the concentrations 

# %%
dynamics.get_diagnostic_L2_data()      # TODO: OBSOLETE!

# %% [markdown]
# # At time t=0.185 :   

# %% [markdown]
# ## Data from this simulation (with substeps):   
# [it took 39 computation steps]

# %%
df.iloc[39]

# %%
substep_data = df.iloc[39][['U', 'X', 'S']].to_numpy()
substep_data

# %% [markdown]
# ## Data from experiment `variable_steps_1` (no substeps, and same main steps as in here):   
# [it took 21 computation steps]

# %%
no_substep_data = np.array([54.631805, 71.023480, 19.712910])
no_substep_data

# %% [markdown]
# ## Data from experiment `substeps_2` (with tiny fixed steps, used as a proxy for exact solution):

# %%
exact_data = np.array([54.529071, 71.351769, 19.590088])
exact_data

# %% [markdown]
# ## Comparing the results:

# %%
substep_error = abs(substep_data - exact_data)
substep_error

# %%
no_substep_error = abs(no_substep_data - exact_data)
no_substep_error

# %% [markdown]
# ## Conclusions:  
# In this scenario, the substeps option made the errors WORSE - in spite of a substantial number of extra steps.  
# Also, substeps led to irregularies in the plot of the concentrations of S.   
# At any rates, substeps are far less relevant ever since the introduction of variable (main) step.

# %% [markdown]
# # => The current implementation of substeps is now deprecated - and will be taken out at a later release

# %%