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
# ## Coupled reactions: `S <-> P` , `E + S <-> ES*`, and  `ES* <-> E + P`
# A direct reaction and the same reaction, catalyzed and showing the transient intermediate state
# # Variation of experiment enzyme_3
# ### Re-run from same initial concentrations of S ("Substrate") and P ("Product"), for various concentations of the enzyme `E`: from zero to hugely abundant  
# #### Given the initial [S]=20 and P=[0], the times at which [S] = [P] are computed for various concentrations of `E`
#
# 1. No enzyme `E`  
# 2. [E] = 0.2  (1/100 of the initial [S])  
# 3. [E] = 1  
# 4. [E] = 5  
# 5. [E] = 20  
# 6. [E] = 30  
# 7. [E] = 100  
# 8. [E] = 2000  
#
# LAST REVISED: May 24, 2024

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from src.modules.chemicals.chem_data import ChemData
from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.modules.movies.movies import MovieTabular

import pandas as pd
import plotly.express as px

# %%
# Initialize the system
chem_data = ChemData(names=["S", "P", "E", "ES*"])
                                                     
# Reaction S <-> P , with 1st-order kinetics, favorable thermodynamics in the forward direction, 
# and a forward rate that is much slower than it would be with the enzyme - as seen in the next reaction, below
chem_data.add_reaction(reactants="S", products="P",
                       forward_rate=1., delta_G=-3989.73)

     
# Reaction E + S <-> ES* , with 1st-order kinetics, and a forward rate that is much faster than it was without the enzyme
# Thermodynamically, the forward direction is at a disadvantage (higher energy state) because of the activation barrier in forming the transient state ES*
chem_data.add_reaction(reactants=["E", "S"], products=["ES*"],
                       forward_rate=100., delta_G=2000)                           
                                                      
# Reaction ES* <-> E + P , with 1st-order kinetics, and a forward rate that is much faster than it was without the enzyme
# Thermodynamically, the total energy change of this reaction and the previous one adds up to the same value as the reaction without the enzyme (-3989.73)
chem_data.add_reaction(reactants=["ES*"], products=["E", "P"],
                       forward_rate=200., delta_G=-5989.73)

# NOTE: the forward_rate's of the last 2 reactions (the catalyzed ones) were tweaked, 
#       to lead to a crossover point [S] = [P] at about the same time as in experiment `enzyme3`, in step 2 (when [E] = 0.2)

chem_data.describe_reactions()

# %% [markdown]
# Note that `E` is not labeled as an "enzyme" because it doesn't appear as a catalyst in any of the registered reactions; it only becomes an enzyme in the context of the _compound_ reaction from (2) and (3)

# %%

# %% [markdown]
# # 1. Set the initial concentrations of all the chemicals - starting with no enzyme

# %%
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 0., "ES*": 0.},
                  snapshot=True)      # Initially, no enzyme `E`
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Advance the reactions (for now without enzyme) to equilibrium

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# This setting is a preset that can be adjusted make the time resolution finer or coarser;
# it will stay in effect from now on, unless explicitly changed later
dynamics.use_adaptive_preset(preset="mid")

# Perform the reactions
dynamics.single_compartment_react(duration=4.0,
                                  initial_step=0.1, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=['darkorange', 'green', 'violet', 'red'], show_intervals=True, title_prefix="With ZERO enzyme")

# %% [markdown]
# ### The reactions, lacking enzyme, are proceeding slowly towards equilibrium, just like the reaction that was discussed in part 1 of the experiment "enzyme_1"

# %% [markdown]
# #### Locate the intersection of the curves for [S] and [P]:

# %%
dynamics.curve_intersection("S", "P", t_start=0, t_end=1.0)

# %%
# Verify that the reactions have reached equilibrium
dynamics.is_in_equilibrium(tolerance=2)

# %%

# %% [markdown]
# # 2. Re-start all the reactions from the same initial concentrations - except for now having a tiny amount of enzyme (two orders of magnitude less than the starting [S])

# %%
dynamics = ReactionDynamics(chem_data=chem_data)      # A brand-new simulation, with the same chemicals and reactions as before 

# %%
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 0.2, "ES*": 0.},
                  snapshot=True)      # A tiny bit of enzyme `E`: 1/100 of the initial [S]
dynamics.describe_state()

# %%

# %% [markdown] tags=[]
# ### Re-take the new system (now with a tiny amount of enzyme) to equilibrium

# %%
# These settings can be tweaked to make the time resolution finer or coarser.  
# Here we use a "slower" heuristic: very conservative about taking larger steps
dynamics.use_adaptive_preset(preset="slower")

dynamics.single_compartment_react(duration=1.5, 
                                  initial_step=0.00005, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=['darkorange', 'green', 'violet', 'red'], show_intervals=True, title_prefix="With a tiny amount of enzyme E")

# %% [markdown]
# #### Locate the intersection of the curves for [S] and [P]:

# %%
dynamics.curve_intersection("S", "P", t_start=0, t_end=1.0)

# %% [markdown]
# #### As mentioned earlier, the forward rates of the catalyzed reactions were tweaked, to produce a crossover time (with this initial [E]) close to what we had in the simpler model used in experiment `enzyme_3`

# %% [markdown]
# ## Notice how even a tiny amount of enzyme (1/100 of the initial [S])  makes a very pronounced difference!

# %%
# Verify that the reactions have reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown]
# Interestingly, most of the inital [E] of 0.2 is now, at equilibrium, stored as [ES*]=0.1194; the energy of the "activation barrier" from E + S to ES* might be unrealistically low (2000 Joules).  Zooming in on the very earl part of the plot:

# %%
dynamics.plot_history(chemicals=['E', 'ES*'], colors=['violet', 'red'], show_intervals=True, 
                      title_prefix="With a tiny amount of enzyme E", xrange=[0., 0.002])

# %%

# %%

# %% [markdown] tags=[]
# # 3. Stop when [P] reaches a particular threshold

# %% [markdown]
# The equilibrium concentrations of P is 16.666 and that of S is 3.333  
# We'll stop the simulation when [P] first rises above 70% of its equilibrium concentrations

# %%
P_threshold = 16.666 * 0.7
P_threshold

# %%
dynamics = ReactionDynamics(chem_data=chem_data)   # A brand-new simulation, with the same chemicals and reactions as before 

# %%
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 0.2, "ES*": 0.},
                  snapshot=True)      # A tiny bit of enzyme `E`.  Same as in the previous run
dynamics.describe_state()

# %%
# These settings can be tweaked to make the time resolution finer or coarser.  
# Here we use a "slower" heuristic: very conservative about taking larger steps
dynamics.use_adaptive_preset(preset="slower")

# Instead of specifying a desired duration, like done before, we specify a termination criterion
dynamics.single_compartment_react(stop=("conc_above", ("P", P_threshold)), max_steps=10000, 
                                  initial_step=0.00005, variable_steps=True, explain_variable_steps=False)

# %%
dynamics.plot_history(colors=['darkorange', 'green', 'violet', 'red'], show_intervals=True, title_prefix="With a tiny amount of enzyme")

# %%
dynamics.curve_intersection("S", "P", t_start=0, t_end=1.0)  # This will be the same as in the previous step, as a double-check

# %%
dynamics.get_history()

# %% [markdown]
# #### Notice how the simulation got automatically stopped as soon as [P] rose over P_threshold = 11.6662

# %%

# %% [markdown]
# # WIP BELOW

# %%
dynamics.get_history(t_start=0, t_end=1.0, columns=["SYSTEM TIME", "S", "P"])

# %%

# %%

# %%

# %%

# %%

# %%
