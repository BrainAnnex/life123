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
# ## A complex reaction `A <-> C` derived from 2 elementary reactions `A <-> B` and `B <-> C`  
# We are given the time evolution of the complex reaction, and want to estimate the 
# rate constants of the elementary reactions.
# Assume the reaction is known to be 1st order (won't verify that.) 
#
# # IN-PROGRESS
#
# In PART 1, a time evolution of [A] and [B] is obtained by simulation  
# In PART 2, the time functions generated in Part 1 are taken as a _starting point,_ to estimate the rate constants of `A <-> B`  
# In PART 3, we'll repeat what we did in Part 2, but this time showing the full details of how the answer is arrived at
#
# LAST REVISED: Dec. 3, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.modules.visualization.plotly_helper import PlotlyHelper

import numpy as np

# %%

# %% [markdown]
# # PART 1 - We'll generate the time evolution of [A] and [C] by simulating reactions of KNOWN rate constants...
# ## but pretend you don't see this section! (because we later want to estimate those rate constants)

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = ReactionDynamics(names=["A", "B", "C"])

# Reaction A <-> B (slower)
dynamics.add_reaction(reactants="A", products="B",
                      forward_rate=8., reverse_rate=2.) 

# Reaction A <-> B (faster)
dynamics.add_reaction(reactants="B", products="C",
                      forward_rate=12., reverse_rate=1.) 
                                   
dynamics.describe_reactions()

# %% [markdown]
# ### Run the simulation

# %%
dynamics.set_conc([50., 0., 0.], snapshot=True)  # Set the initial concentrations of all the chemicals, in their index order
dynamics.describe_state()

# %%
#dynamics.set_diagnostics()         # To save diagnostic information about the call to single_compartment_react()

# These settings can be tweaked to make the time resolution finer or coarser.  
# Here we use a "mid" heuristic: neither too fast nor too prudent
dynamics.use_adaptive_preset(preset="mid")

dynamics.single_compartment_react(initial_step=0.01, reaction_duration=0.8,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True, explain_variable_steps=False)

# %%
dynamics.plot_history(colors=['blue', 'red', 'green'], show_intervals=True)

# %%
dynamics.is_in_equilibrium(tolerance=15)

# %%

# %% [markdown]
# # PART 2 - This is the starting point of fitting the data from part 1.  
# ### We're given the data of the above curves - i.e. the system history, and we want to estimate the rate constants (forward and reverse) of the reaction `A <-> C`

# %% [markdown]
# Let's start by taking stock of the actual data (saved during the simulation of part 1):

# %%
df = dynamics.get_history() 
df

# %% [markdown]
# ## Column B is NOT given to us.  For example, `B` might be an intermediary we can't measure.  Only [A] and [C] are given to us, on some variable time grid

# %% [markdown]
# #### Let's extract some columns, as Numpy arrays:

# %%
t_arr = df["SYSTEM TIME"].to_numpy()   # The independent variable : Time

# %%
A_conc = df["A"].to_numpy()

# %%
C_conc = df["C"].to_numpy()

# %% [markdown]
# ### Here, we take the easy way out, using a specialized Life123 function!
# (in Part 3, we'll do a step-by-step derivation, to see how it works)

# %%
dynamics.estimate_rate_constants(t=t_arr, reactant_conc=A_conc, product_conc=C_conc, reactant_name="A", product_name="C")

# %% [markdown]
# ### The least-square fit is awful : the complex reaction `A <-> C` doesn't seem to be amenable to being modeled as a simple reaction with some suitable rate constants
# Probably not too surprising given our "secret" knowledge from Part 1 that the complex reaction originates from 2 elementary reactions where one doesn't dominate the other one in terms of reaction kinetics

# %% [markdown]
# ### A glance at the above diagram reveals much-better linear fits, if split into 2 portions, one where A(t) ranges from 0 to about 25, and one from about 25 to 50

# %%
