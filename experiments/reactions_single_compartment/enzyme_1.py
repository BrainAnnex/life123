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
# ## Comparing the reaction `A <-> B` , with and without an enzyme
#
# LAST REVISED: Nov. 3, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from src.modules.chemicals.chem_data import ChemData
from src.modules.reactions.reaction_dynamics import ReactionDynamics

# %% [markdown]
# # 1. WITHOUT ENZYME
# ### `A` <-> `B`

# %%
# Initialize the system
chem_data = ChemData(names=["A", "B"])

# Reaction A <-> B , with 1st-order kinetics, favorable thermodynamics in the forward direction, 
# and a forward rate that is much slower than it would be with the enzyme - as seen in part 2, below
chem_data.add_reaction(reactants="A", products="B",
                       forward_rate=1., delta_G=-3989.73)

chem_data.describe_reactions()

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"A": 20., "B": 0.},
                  snapshot=True)
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Take the initial system to equilibrium

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
dynamics.set_step_factors(upshift=1.5, downshift=0.5, abort=0.5)
dynamics.set_error_step_factor(0.5)

dynamics.single_compartment_react(initial_step=0.1, reaction_duration=3.0,
                                  variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['darkorange', 'green'], show_intervals=True, title_prefix="WITHOUT enzyme")

# %% [markdown]
# #### Note how the time steps get automatically adjusted, as needed by the amount of change - including a complete step abort/redo at time=0

# %%
dynamics.curve_intersection("A", "B", t_start=0, t_end=1.0)

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=3)

# %%

# %% [markdown]
# # 2. WITH ENZYME `E`
# ### `A` + `E` <-> `B` + `E`

# %% [markdown]
# ### Note: for the sake of the demo, we'll completely ignore the concomitant (much slower) direct reaction A <-> B
# This in an approximation that we'll drop in later experiments

# %%
# Initialize the system
chem_data = ChemData(names=["A", "B", "E"])

# Reaction A + E <-> B + E , with 1st-order kinetics, and a forward rate that is faster than it was without the enzyme
# Thermodynamically, there's no change from the reaction without the enzyme
chem_data.add_reaction(reactants=["A", "E"], products=["B", "E"],
                       forward_rate=10., delta_G=-3989.73)

chem_data.describe_reactions()  # Notice how the enzyme `E` is noted in the printout below

# %% [markdown]
# ### Notice how, while the ratio kF/kR is the same as it was without the enzyme (since it's dictated by the energy difference), the individual values of kF and kR now are each 10 times bigger than before

# %% [markdown]
# ### Set the initial concentrations of all the chemicals (to what they originally were)

# %%
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"A": 20., "B": 0., "E": 30.},
                  snapshot=True)      # Plenty of enzyme `E`
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Take the initial system to equilibrium

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
dynamics.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4)
dynamics.set_error_step_factor(0.25)

dynamics.single_compartment_react(initial_step=0.1, reaction_duration=0.1,
                                  variable_steps=True, explain_variable_steps=False)

# %% [markdown]
# #### Note how the (proposed) initial step - kept the same as the previous run - is now _extravagantly large_, given the much-faster reaction dynamics.  However, the variable-step engine intercepts and automatically corrects the problem!

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="WITH enzyme")

# %%
dynamics.curve_intersection("A", "B", t_start=0, t_end=0.02)

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown]
# ## Thanks to the (abundant) enzyme, the reaction reaches equilibrium roughtly around t=0.02, far sooner than the roughly t=3.5 without enzyme
# The concentrations of `A` and `B` now become equal (cross-over) at t=0.00246 , rather than t=0.740

# %%
dynamics.get_history()

# %%
