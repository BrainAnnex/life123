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
# ## `A` up-regulates `B` , 
# ### by being *the limiting reagent* in the reaction `A + X <-> 2B` (mostly forward), where `X` is plentiful
# 1st-order kinetics.   
# If [A] is low, [B] remains low, too.  Then, if [A] goes high, then so does [B].  However, at that point, A can no longer bring B down to any substantial extent.
#
# See also the experiment "1D/reactions/up_regulation_1"
#
# LAST REVISED: Feb. 5, 2023

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.express as px
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
chem_data = chem(names=["A", "X", "B"])

# Reaction A + X <-> 2B , with 1st-order kinetics for all species
chem_data.add_reaction(reactants=[("A") , ("X")], products=[(2, "B")],
                 forward_rate=8., reverse_rate=2.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc={"A": 5., "X": 100., "B": 0.},
                  snapshot=True)      # A is scarce, X is plentiful, B is absent
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Take the initial system to equilibrium

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.0005, stop_time=0.015,
                                  dynamic_substeps=2, rel_fast_threshold=15)

df = dynamics.get_history()
df

# %% [markdown]
# A, as the scarse limiting reagent, stops the reaction.  
# When A is low, B is also low.

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=['red', 'darkorange', 'green'])

# %% [markdown] tags=[]
# # Now, let's suddenly increase [A]

# %%
dynamics.set_chem_conc(species_name="A", conc=50., snapshot=True)
dynamics.describe_state()

# %%
dynamics.get_history(tail=5)

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.0005, stop_time=0.035,
                                  dynamic_substeps=2, rel_fast_threshold=15)

df = dynamics.get_history()
df

# %%
dynamics.explain_time_advance()

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=2)

# %%
dynamics.plot_curves(colors=['red', 'darkorange', 'green'])

# %% [markdown]
# **A**, still the limiting reagent, is again stopping the reaction.  
# The (transiently) high value of [A] led to a high value of [B]

# %% [markdown] tags=[]
# # Let's again suddenly increase [A]

# %%
dynamics.set_chem_conc(species_name="A", conc=30., snapshot=True)
dynamics.describe_state()

# %%
dynamics.get_history(tail=5)

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
dynamics.single_compartment_react(time_step=0.0005, stop_time=0.070,
                                  dynamic_substeps=2, rel_fast_threshold=15)

df = dynamics.history.get()
df

# %%
dynamics.explain_time_advance()

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=2)

# %%
dynamics.plot_curves(colors=['red', 'darkorange', 'green'])

# %% [markdown]
# `A`, again the scarce limiting reagent, stops the reaction yet again.  
# And, again, the (transiently) high value of [A] up-regulated [B]  
#
# Notes:   
# `A` can up-regulate `B`, but it cannot bring it down.  
# `X` will soon need to be replenished, if `A` is to continue being the limiting reagent.

# %%
# Look up the some of the intersections of the [A] and [B] curves
dynamics.curve_intersection(t_start=0, t_end=0.01, var1="A", var2="B")

# %%
dynamics.curve_intersection(t_start=0.0151, t_end=0.02, var1="A", var2="B")

# %% [markdown]
# Note: the _curve_intersection()_ function currently cannot location the intersection at t=0.015 (the vertical rise in the red line); this issue will get addressed in future versions...

# %%
