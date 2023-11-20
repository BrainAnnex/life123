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
# ## Exploration of variable time steps in the simulation of the 2 coupled reactions:
# ### `2 S <-> U` and `S <-> X`   
# Both mostly forward.  1st-order kinetics throughout.   
#
# Based on the reactions and initial conditions of the experiment `up_regulate_3`
#
# This experiment gets repeated, with very fine _fixed_ steps (as a proxy for the "exact value"), in `variable_steps_2`
#
# LAST REVISED: July 23, 2023

# %% [markdown]
# ![Adaptive time steps](../../docs/variable_steps.png)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.chemicals.chem_data import ChemData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
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
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.25, high=0.64, abort=1.44)
dynamics.set_thresholds(norm="norm_B")   # We are disabling norm_B (to conform to the original run)
dynamics.set_step_factors(upshift=2.0, downshift=0.5, abort=0.5)    # Note: upshift=2.0 seems to often be excessive.  About 1.4 is currently recommended
dynamics.set_error_step_factor(0.5)

dynamics.single_compartment_react(initial_step=0.01, target_end_time=2.0, 
                                  variable_steps=True, explain_variable_steps=True)

# %%
dynamics.get_history()

# %%
(transition_times, step_sizes) = dynamics.explain_time_advance(return_times=True)

# %%
np.array(step_sizes)

# %%
np.array(transition_times)    # Note: there will be one more transition time (the end time) than step sizes

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_history(colors=['green', 'orange', 'blue'])

# %%
dynamics.curve_intersection("U", "X", t_start=0.3, t_end=0.35)  # Compare with the value from experiment "variable_steps_2"

# %%
dynamics.plot_history(colors=['green', 'orange', 'blue'], show_intervals=True)

# %%
# Show the "critical values", i.e. times when the step size changes
dynamics.plot_history(colors=['green', 'orange', 'blue'], vertical_lines=transition_times,
                      title="Critical values of time-step changes for reactions `2 S <-> U` and `S <-> X`")

# %% [markdown]
# ## Note: the dashed lines in the plots immediatly above are NOT the steps; they are the "critical values", i.e. times when the step size changes.   
# The time steps were shown in an earlier plots

# %%
dynamics.plot_step_sizes(show_intervals=True)

# %%
dynamics.is_in_equilibrium()

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=1)

# %%
dynamics.get_diagnostic_conc_data()

# %%
dynamics.get_diagnostic_decisions_data()

# %% [markdown]
# #### Notice how the first step got aborted, and re-run, because of the large value of `norm_A`

# %%

# %%
dynamics.get_diagnostic_decisions_data_ALT()      # TODO: OBSOLETE!

# %%
