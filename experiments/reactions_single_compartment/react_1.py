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
# ### A minimalist, "get-stared", demonstration for the reaction `A <-> B`,
# with 1st-order kinetics in both directions, taken to equilibrium.
#
# "No frills!"  For advanced graphics, analysis, diagnostics, fine-tuning, etc, please see other experiments.
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta34.1)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from life123 import UniformCompartment

# %%

# %% [markdown]
# ## Initialize the System

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = UniformCompartment(names=["A", "B"])  

# Reaction A <-> B , with 1st-order kinetics in both directions
dynamics.add_reaction(reactants="A", products="B", 
                       forward_rate=3., reverse_rate=2.)

dynamics.describe_reactions()

# %%
# Set the initial concentrations of all the chemicals
dynamics.set_conc({"A": 80., "B": 10.})

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.single_compartment_react(initial_step=0.1, target_end_time=1.)   # Using defaults for all other parameters

# %%
dynamics.get_history()   # The system's history, saved during the run of single_compartment_react()

# %% [markdown] tags=[]
# ## Plots changes of concentration with time  
# Notice that adaptive variable time steps were automatically taken

# %%
dynamics.plot_history(show_intervals=True)

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%
