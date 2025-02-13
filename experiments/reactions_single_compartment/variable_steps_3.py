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

# %% [markdown]
# ### TAGS :  "uniform compartment", "under-the-hood"

# %%
LAST_REVISED = "Oct. 11, 2024"
LIFE123_VERSION = "1.0.0.beta.39"   # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData as chem
from life123 import UniformCompartment

import numpy as np

# %%

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
dynamics = UniformCompartment(chem_data=chem_data, preset=None)
dynamics.set_conc(conc={"U": 50., "X": 100., "EXTRA 2": 55.,  "EXTRA 3": 100. })  # The EXTRA's are junk for testing
dynamics.describe_state()

# %%
dynamics.enable_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# These adaptive-time settings (normally specified with a preset) are being set explicitly
dynamics.adaptive_steps.set_thresholds(norm="norm_A", low=0.25, high=0.64, abort=1.44)
dynamics.adaptive_steps.set_step_factors(upshift=2.0, downshift=0.5, abort=0.5, error=0.5)    
                                        # Note: upshift=2.0 seems to often be excessive.  About 1.4 is currently recommended

dynamics.adaptive_steps.show_adaptive_parameters()

# %%
dynamics.single_compartment_react(initial_step=0.01, target_end_time=2.0, 
                                  variable_steps=True, explain_variable_steps=[0, 2])

# %% [markdown]
# ## Compare the above printout with its counterpart from the experiment `variable_steps_1`
# Notice the extra lines we have this time, saying  _"Restricting adaptive time step analysis to 3 chemicals only"_

# %%
dynamics.get_history()

# %%
(transition_times, step_sizes) = dynamics.get_diagnostics().explain_time_advance(return_times=True)

# %%
np.array(step_sizes)

# %%
np.array(transition_times)    # Note: there will be one more transition time (the end time) than step sizes

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_history(colors=['lightgray', 'green', 'lightgray', 'orange', 'darkturquoise', 'lightgray'])

# %%
dynamics.curve_intersect("U", "X", t_start=0.3, t_end=0.35)

# %%

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
dynamics = UniformCompartment(chem_data=chem_data, preset=None)
dynamics.set_conc(conc={"U": 50., "X": 100., "E": 1. })
dynamics.describe_state()

# %%
dynamics.enable_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# These adaptive-time settings (normally specified with a preset) are being set explicitly
dynamics.adaptive_steps.set_thresholds(norm="norm_A", low=0.25, high=0.64, abort=1.44)
dynamics.adaptive_steps.set_step_factors(upshift=2.0, downshift=0.5, abort=0.5, error=0.5)

dynamics.adaptive_steps.show_adaptive_parameters()

# %%
dynamics.single_compartment_react(initial_step=0.01, target_end_time=2.0, 
                                  variable_steps=True, explain_variable_steps=[0,2])

# %% [markdown]
# ## Compare the above printout with its counterpart from the experiment `variable_steps_1`
# Notice the extra lines we have this time, saying  _"Restricting adaptive time step analysis to 3 chemicals only"_

# %%
dynamics.get_history()

# %%
(transition_times, step_sizes) = dynamics.get_diagnostics().explain_time_advance(return_times=True)

# %%
np.array(step_sizes)

# %%
np.array(transition_times)    # Note: there will be one more transition time (the end time) than step sizes

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_history(colors=['green', 'orange', 'darkturquoise', 'gray'])

# %%
dynamics.curve_intersect("U", "X", t_start=0.3, t_end=0.35)

# %%
