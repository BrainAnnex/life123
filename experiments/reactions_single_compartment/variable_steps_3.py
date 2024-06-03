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
# ## Validation of the independence of the automatically-picked variable time steps from the presence of extraneous chemicals not participating in the dynamical changes.
#
# Just as in experiment `variable_steps_1` : 2 coupled reactions: `2 S <-> U` and `S <-> X`   
#
# However, here:
#
# * in part 1 there are a few extra chemicals in the system that don't participate in any of the reactions  
# * in part 2 there's a hypothetical enzyme (with concentration 1) that catalyzes the first reaction 
#
# In either case, the extra chemicals and the enzyme don't vary in concentration - and thus **get automatically excluded from considerations about the adaptive variable step sizes** , which remain exactly as they were in experiment `variable_steps_1`
#
# LAST REVISED: Dec. 3, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.chemicals.chem_data import ChemData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np

# %%

# %% [markdown]
# # PART 1 : 
# #### Notice the "EXTRA" group of chemicals, that don't participate in any of the reactions 

# %%
# Initialize the system.  
chem_data = chem(names=["EXTRA 1", "U", "EXTRA 2", "X", "S", "EXTRA 3"])

# Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants=[(2, "S", 1)], products="U",
                       forward_rate=8., reverse_rate=2.)

# Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants="S", products="X",
                       forward_rate=6., reverse_rate=3.)

chem_data.describe_reactions()

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0., "EXTRA 1": 0.,  "EXTRA 2": 55.,  "EXTRA 3": 100. })
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

# %% [markdown]
# ## Compare the above printout with its counterpart from the experiment `variable_steps_1`
# Notice the extra lines we have this time, saying  _"Restricting adaptive time step analysis to 3 chemicals only"_

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
dynamics.plot_history(colors=['lightgray', 'green', 'lightgray', 'orange', 'blue', 'lightgray'])

# %%
dynamics.curve_intersect("U", "X", t_start=0.3, t_end=0.35)

# %%

# %% [markdown]
# # PART 2 : 
# #### Notice the fictitious enzyme "E" in the first reaction, with concentration 1 (not affecting any of the kinetic parameters)

# %%
# Initialize the system
chem_data = chem(names=["U", "X", "S", "E"])

# Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants=[(2, "S", 1), "E"], products=["U", "E"],
                       forward_rate=8., reverse_rate=2.)

# Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants="S", products="X",
                       forward_rate=6., reverse_rate=3.)

chem_data.describe_reactions()

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0., "E": 1. })
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

# %% [markdown]
# ## Compare the above printout with its counterpart from the experiment `variable_steps_1`
# Notice the extra lines we have this time, saying  _"Restricting adaptive time step analysis to 3 chemicals only"_

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
dynamics.plot_history(colors=['green', 'orange', 'blue', 'gray'])

# %%
dynamics.curve_intersect("U", "X", t_start=0.3, t_end=0.35)

# %%
