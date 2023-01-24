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
# ## A cycle of reactions A <-> B <-> C <-> A, 
# ### The "closing" of the above cycle (the "return" parth from C to A) is couple with an "energy donor" reaction:
# ### C + E_High <-> A + E_Low
# where E_High and E_Low are, respectively, the high- and low- energy molecules that drive the cycle (for example, think of ATP/ADP)
#
# All 1st-order kinetics.   
#
# LAST REVISED: Jan. 23, 2023

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from modules.reactions.reaction_data import ReactionData as chem
from modules.reactions.reaction_dynamics import ReactionDynamics
from modules.numerical.numerical import Numerical as num

import numpy as np
import plotly.express as px
from modules.visualization.graphic_log import GraphicLog

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
chem_data = chem(names=["A", "B", "C", "E_high", "E_low"])

# Reaction A <-> B, mostly in forward direction (favored energetically)
# Note: all reactions in this experiment have 1st-order kinetics for all species
chem_data.add_reaction(reactants="A", products="B",
                       forward_rate=9., reverse_rate=3.)

# Reaction B <-> C, also favored energetically
chem_data.add_reaction(reactants="B", products="C",
                       forward_rate=8., reverse_rate=4.)

# Reaction C + E_High <-> A + E_Low, also favored energetically, but kinetically slow
# Note that, thanks to the energy donation from E, we can go "upstream" from C, to the higher-energy level of "A"
chem_data.add_reaction(reactants=["C" , "E_high"], products=["A", "E_low"],
                       forward_rate=1., reverse_rate=0.2)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}  # Note the abundant energy source
initial_conc

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc=initial_conc, snapshot=True)
dynamics.describe_state()

# %% [markdown] tags=[]
# ## Start the reaction
# ### We'll split the simuation in three segments (apparent from the graph, further down):  
# Time [0-0.03] fast changes   
# Time [0.03-5.] medium changes   
# Time [5.-8.] slow changes, as we approach equilibrium

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# %% [markdown]
# # 1. For starters, we'll use FIXED time resolution, with COARSE time steps   
# (trial and error, not shown, reveals that increasing any of the time steps below, soon leads to "excessive time step" errors)

# %%
dynamics.single_compartment_react(time_step=0.0008, stop_time=0.03)
dynamics.get_history()

# %%
dynamics.single_compartment_react(time_step=0.001, stop_time=5.)
dynamics.get_history()

# %%
dynamics.single_compartment_react(time_step=0.005, stop_time=8.)

# %%
df = dynamics.get_history()
df

# %% [markdown]
# ### Notice we created 5,609 data points

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(chemicals=["E_high", "E_low"], colors=["red", "grey"])

# %%
dynamics.plot_curves(chemicals=["A", "B", "C"])

# %% [markdown]
# ### The plots has 4 distinctive intersections; locate them and save them for later comparisons across repeated runs:

# %%
run1 = []

# %%
run1.append(dynamics.curve_intersection(t_start=1., t_end=2., var1="E_high", var2="E_low"))

# %%
run1.append(dynamics.curve_intersection(t_start=2.31, t_end=2.33, var1="A", var2="B"))

# %%
run1.append(dynamics.curve_intersection(t_start=3., t_end=4., var1="A", var2="C"))

# %%
run1.append(dynamics.curve_intersection(t_start=3., t_end=4., var1="B", var2="C"))

# %%
run1

# %% [markdown]
# ### Verify the final equilibrium state:

# %%
dynamics.is_in_equilibrium()

# %%

# %%

# %% [markdown]
# # Everything below is just a repeat of the experiment, with different timesteps, for technical comparisons

# %% [markdown]
# # 2. VARIABLE time resolution, using the same primary time steps as before 
#

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)   # Note: OVER-WRITING the "dynamics" object
dynamics.set_conc(conc=initial_conc, snapshot=True) 
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# %%
dynamics.single_compartment_react(time_step=0.0008, stop_time=0.03,
                                  dynamic_steps=2, abs_fast_threshold=750.)

# %%
dynamics.single_compartment_react(time_step=0.001, stop_time=5.,
                                 dynamic_steps=2, abs_fast_threshold=250.)

# %%
dynamics.single_compartment_react(time_step=0.005, stop_time=8.,
                                  dynamic_steps=2, abs_fast_threshold=2.)

# %%
df = dynamics.get_history()
df

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(chemicals=["E_high", "E_low"], colors=["red", "grey"])

# %%
dynamics.plot_curves(chemicals=["A", "B", "C"])

# %% [markdown]
# ### The plots has 4 distinctive intersections; locate them:

# %%
run2 = []

# %%
run2.append(dynamics.curve_intersection(t_start=1., t_end=2., var1="E_high", var2="E_low"))

# %%
run2.append(dynamics.curve_intersection(t_start=2.31, t_end=2.33, var1="A", var2="B"))

# %%
run2.append(dynamics.curve_intersection(t_start=3., t_end=4., var1="A", var2="C"))

# %%
run2.append(dynamics.curve_intersection(t_start=3., t_end=4., var1="B", var2="C"))

# %%
run2

# %%

# %%

# %% [markdown]
# # 3. FIXED time resolution, using the smallest time substeps used in the variable step 2 
#

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)   # Note: OVER-WRITING the "dynamics" object
dynamics.set_conc(conc=initial_conc, snapshot=True) 
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# %% tags=[]
dynamics.single_compartment_react(time_step=0.0004, stop_time=0.03)

# %%
dynamics.single_compartment_react(time_step=0.0005, stop_time=5.)

# %%
dynamics.single_compartment_react(time_step=0.0025, stop_time=8.)

# %%
df = dynamics.get_history()
df

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(chemicals=["E_high", "E_low"], colors=["red", "grey"])

# %%
dynamics.plot_curves(chemicals=["A", "B", "C"])

# %% [markdown]
# ### The plots has 4 distinctive intersections; locate them:

# %%
run3 = []

# %%
run3.append(dynamics.curve_intersection(t_start=1., t_end=2., var1="E_high", var2="E_low"))

# %%
run3.append(dynamics.curve_intersection(t_start=2.31, t_end=2.33, var1="A", var2="B"))

# %%
run3.append(dynamics.curve_intersection(t_start=3., t_end=4., var1="A", var2="C"))

# %%
run3.append(dynamics.curve_intersection(t_start=3., t_end=4., var1="B", var2="C"))

# %%
run3

# %%

# %%

# %%

# %%
num.compare_results(run1, run3)

# %%
num.compare_results(run2, run3)
