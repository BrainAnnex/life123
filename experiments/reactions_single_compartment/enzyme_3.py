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
# ## Coupled pair of reactions: `S <-> P` , and  `S` + `E` <-> `P` + `E`
# A direct reaction and the same reaction, catalyzed
# ### Re-run from same initial concentrations of S ("Substrate") and P ("Product") for various concentations of the enzyme `E`: from zero to hugely abundant
#
# LAST REVISED: Nov. 7, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from src.modules.chemicals.chem_data import ChemData
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import pandas as pd
import plotly.express as px

# %%
# Initialize the system
chem_data = ChemData(names=["S", "P", "E"])

# Reaction S <-> P , with 1st-order kinetics, favorable thermodynamics in the forward direction, 
# and a forward rate that is much slower than it would be with the enzyme - as seen in the next reaction, below
chem_data.add_reaction(reactants="S", products="P",
                       forward_rate=1., delta_G=-3989.73)

# Reaction S + E <-> P + E , with 1st-order kinetics, and a forward rate that is faster than it was without the enzyme
# Thermodynamically, there's no change from the reaction without the enzyme
chem_data.add_reaction(reactants=["S", "E"], products=["P", "E"],
                       forward_rate=10., delta_G=-3989.73)

chem_data.describe_reactions()     # Notice how the enzyme `E` is noted in the printout below

# %%

# %% [markdown]
# # 1. Set the initial concentrations of all the chemicals - starting with no enzyme

# %%
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 0.},
                  snapshot=True)      # Initially, no enzyme `E`
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Advance the reactions (for now without enzyme) to equilibrium

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
dynamics.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4)
dynamics.set_error_step_factor(0.25)

# Perform the reactions
dynamics.single_compartment_react(reaction_duration=4.0,
                                  initial_step=0.1, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="With ZERO enzyme")

# %% [markdown]
# ### The reactions, lacking enzyme, are proceeding slowly towards equilibrium, just like the reaction that was discussed in part 1 of the experiment "enzyme_1"

# %%
crossover_points = []       # We'll be saving all the crossover points, together with their corresponding enzyme concentration
new_crossover = dynamics.curve_intersection("S", "P", t_start=0, t_end=1.0)
crossover_points.append([new_crossover[0], dynamics.get_chem_conc("E")])
new_crossover

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=2)

# %%

# %% [markdown]
# # 2. Re-start all the reactions from the same initial concentrations - except for now having a tiny amount of enzyme (two orders of magnitude less than the starting [S])

# %%
dynamics = ReactionDynamics(chem_data=chem_data)   # A brand-new simulation  

# %%
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 0.2},
                  snapshot=True)      # A tiny bit of enzyme `E`

# %%
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Re-take the new system (now with a tiny amount of enzyme) to equilibrium

# %%
dynamics.single_compartment_react(reaction_duration=1.5, 
                                  initial_step=0.05, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="With a tiny amount of enzyme")

# %%
new_crossover = dynamics.curve_intersection("S", "P", t_start=0, t_end=0.5)
crossover_points.append([new_crossover[0], dynamics.get_chem_conc("E")])
new_crossover

# %% [markdown]
# ## Notice how even a tiny amount of enzyme (1/100 of the initial [A])  makes a very pronounced difference!

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%

# %% [markdown]
# # 3. Re-start all the reactions from the same initial concentrations - except for now having a more substantial amount of enzyme

# %%
dynamics = ReactionDynamics(chem_data=chem_data)   # A brand-new simulation  

# %%
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 1.},
                  snapshot=True)      # A more substantial amount of enzyme `E`

# %%
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Re-take the new system (now with a more substantial amount of enzyme) to equilibrium

# %%
dynamics.single_compartment_react(reaction_duration=0.5, 
                                  initial_step=0.02, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="With a more substantial amount of enzyme")

# %%
new_crossover = dynamics.curve_intersection("S", "P", t_start=0, t_end=0.5)
crossover_points.append([new_crossover[0], dynamics.get_chem_conc("E")])
new_crossover

# %% [markdown]
# ## Notice the continued - and substantial - speedup of the reaction, over the earlier runs

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %%

# %% [markdown]
# # 4. Re-start all the reactions from the same initial concentrations - except for now having a good amount of enzyme

# %%
dynamics = ReactionDynamics(chem_data=chem_data)   # A brand-new simulation  

# %%
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 5.},
                  snapshot=True)      # A good amount of enzyme `E`

# %%
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Re-take the new system to equilibrium

# %%
dynamics.single_compartment_react(reaction_duration=0.2, 
                                  initial_step=0.01, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="With a good amount of enzyme")

# %%
new_crossover = dynamics.curve_intersection("S", "P", t_start=0, t_end=0.5)
crossover_points.append([new_crossover[0], dynamics.get_chem_conc("E")])
new_crossover

# %% [markdown]
# ## Notice the continued - and substantial - speedup of the reaction, over the earlier runs

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %%

# %% [markdown]
# # 5. Re-start all the reactions from the same initial concentrations - except for now having a lot of enzyme (same as the initial [A])

# %%
dynamics = ReactionDynamics(chem_data=chem_data)   # A brand-new simulation  

# %%
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 20.},
                  snapshot=True)      # A lot of enzyme `E`

# %%
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Re-take the new system (now with a lot of enzyme) to equilibrium

# %%
dynamics.single_compartment_react(reaction_duration=0.05, 
                                  initial_step=0.005, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="With a lot of enzyme")

# %%
new_crossover = dynamics.curve_intersection("S", "P", t_start=0, t_end=0.5)
crossover_points.append([new_crossover[0], dynamics.get_chem_conc("E")])
new_crossover

# %% [markdown]
# ## Notice the continued - and substantial - speedup of the reaction, over the earlier runs

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %%

# %% [markdown]
# # 6. Re-start all the reactions from the same initial concentrations - except for now having a very large amount of enzyme (more than the initial substrate concentration [S])

# %%
dynamics = ReactionDynamics(chem_data=chem_data)   # A brand-new simulation  

# %%
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 30.},
                  snapshot=True)      # A very large amount of enzyme `E`

# %%
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Re-take the new system to equilibrium

# %%
dynamics.single_compartment_react(reaction_duration=0.02, 
                                  initial_step=0.001, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="With a very large amount of enzyme")

# %%
new_crossover = dynamics.curve_intersection("S", "P", t_start=0, t_end=0.5)
crossover_points.append([new_crossover[0], dynamics.get_chem_conc("E")])
new_crossover

# %% [markdown]
# ## Yet more speedup of the reaction, over the previous run

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %%

# %% [markdown]
# # 7. Finally, re-start all the reactions from the same initial concentrations - except for now having a huge amount of enzyme (far more than the initial [A])

# %%
dynamics = ReactionDynamics(chem_data=chem_data)   # A brand-new simulation  

# %%
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 100.},
                  snapshot=True)      # A lavish amount of enzyme `E`

# %%
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Re-take the new system to equilibrium

# %%
dynamics.single_compartment_react(reaction_duration=0.02, 
                                  initial_step=0.0005, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="With a huge amount of enzyme")

# %%
new_crossover = dynamics.curve_intersection("S", "P", t_start=0, t_end=0.5)
crossover_points.append([new_crossover[0], dynamics.get_chem_conc("E")])
new_crossover

# %% [markdown]
# ## Yet more speedup of the reaction, over the previous run

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(explain=False)

# %%

# %% [markdown]
# # 8. Finally, re-start all the reactions from the same initial concentrations - except for now having a LAVISH amount of enzyme (two orders of magnitude more than the starting substrate concentration [S])

# %%
dynamics = ReactionDynamics(chem_data=chem_data)   # A brand-new simulation  

# %%
dynamics.set_conc(conc={"S": 20., "P": 0., "E": 2000.},
                  snapshot=True)      # A lavish amount of enzyme `E`

# %%
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Re-take the new system (now with a lavish amount of enzyme) to equilibrium

# %%
dynamics.single_compartment_react(reaction_duration=0.0015, 
                                  initial_step=0.000005, variable_steps=True, explain_variable_steps=False)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_curves(chemicals=['S', 'P'],
                     colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="With a LAVISH amount of enzyme (NOT shown)")

# %% [markdown]
# _Note: NOT showing the enzyme (concentration 2,000) in the graph, to avoid squishing down the other curves!_

# %%
new_crossover = dynamics.curve_intersection("S", "P", t_start=0, t_end=0.5)
crossover_points.append([new_crossover[0], dynamics.get_chem_conc("E")])
new_crossover

# %% [markdown]
# ## Yet more speedup of the reaction, over the previous run

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%

# %% [markdown]
# ## Now, let's look at the time of the crossover points (for the [S] and [P] curves to intersect), as a function of the Enzyme concentration

# %%
crossover_points

# %% [markdown]
# #### As we previously observed, as the Enzyme concentration is increased over repeated runs, the crossover happens earlier and earlier

# %%
# Same data, as a Pandas dataframe
df = pd.DataFrame(crossover_points, columns = ['crossover time', 'E'])
df

# %%
# Let's plot just the first 6 data points, to avoid squishing the curve
px.line(data_frame=df.loc[0:5],
              x="E", y=["crossover time"],
              title="Time of crossover of [S] and [P] , as a function of Enzyme concentration (first 6 points)",
              labels={"value":"Crossover time"},
              line_shape="spline")

# %%
# Same plot, but with log scales on both axes (all data points used this time, 
# but the value E = 0 is automatically dropped by the graphic function because of the log scale)
px.line(data_frame=df,
              x="E", y=["crossover time"],
              log_x=True, log_y=True,
              title="Time of crossover of [S] and [P] , as a function of Enzyme conc (log scales on both axes)",
              labels={"value":"Crossover time"},
              line_shape="spline")

# %%