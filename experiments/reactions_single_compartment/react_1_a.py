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
# ### A MINIMALIST, "get-started", demonstration for the reaction `A <-> B`,
# with 1st-order kinetics in both directions, taken to equilibrium.
#
# **"No frills!"**  For advanced graphics, analysis, diagnostics, fine-tuning, etc, please see other experiments.

# %% [markdown]
# ### TAGS :   "quick-start", "uniform compartment"

# %%
LAST_REVISED = "Nov. 18, 2024"
LIFE123_VERSION = "1.0.0.rc.0"      # Version this experiment is based on

# %%
#import set_path            # Using MyBinder?  Uncomment this before running the next cell!
                            # Importing this module will add the project's home directory to sys.path

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

import life123

# %%
life123.check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Initialize the System

# %% tags=[]
# Instantiate the simulator and specify the chemicals
uc = life123.UniformCompartment()  

# Reaction A <-> B , with 1st-order kinetics in both directions
uc.add_reaction(reactants="A", products="B", 
                forward_rate=3., reverse_rate=2.)

uc.describe_reactions()

# %%
# Set the initial concentrations of all the chemicals
uc.set_conc({"A": 80., "B": 10.})

# %% [markdown] tags=[]
# ## Run the reaction

# %%
uc.single_compartment_react(initial_step=0.1, target_end_time=1.)   # Using defaults for all other parameters

# %%
uc.get_history()   # The system's history, saved during the run of single_compartment_react()

# %% [markdown] tags=[]
# ## Plots changes of concentration with time  
# Notice that adaptive variable time steps were automatically taken

# %%
uc.plot_history(show_intervals=True)

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium()

# %%