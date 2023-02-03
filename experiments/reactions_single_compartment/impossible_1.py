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
# ## Violating the Laws of Physics for Fun and Insight!
# ###  A cascade of reactions `A <-> B <-> C` , mostly in the forward direction
# ### PART 1 : the above, together with a PHYSICALLY-IMPOSSIBLE "closing" of the above cycle with :
# #### `C <-> A`,  *ALSO* mostly in the forward direction (never mind the laws of thermodymics)!
# ### PART 2 : restoring the law of physics
#
# All 1st-order kinetics.    
#
# LAST REVISED: Feb. 2, 2023

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
chem_data = chem(names=["A", "B", "C"])

# Reaction A <-> B, mostly in forward direction (favored energetically)
# Note: all reactions in this experiment have 1st-order kinetics for all species
chem_data.add_reaction(reactants="A", products="B",
                       forward_rate=9., reverse_rate=3.)

# Reaction B <-> C, also favored energetically
chem_data.add_reaction(reactants="B", products="C",
                       forward_rate=8., reverse_rate=4.)

# %% [markdown]
# # Part 1 - "Turning off the Laws of Physics"!

# %%
# LET'S VIOLATE THE LAWS OF PHYSICS!
# Reaction C <-> A, also mostly in forward direction - MAGICALLY GOING "UPSTREAM" from C, to the higher-energy level of "A"
chem_data.add_reaction(reactants="C" , products="A",
                       forward_rate=3., reverse_rate=2.)   # PHYSICALLY IMPOSSIBLE! Future versions of Life123 may flag this!

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# # Notice the absurdity of the energy levels always going down, throughout the cycle (like in an Escher painting!)

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
initial_conc = {"A": 100., "B": 0., "C": 0.} 
initial_conc

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc(conc=initial_conc, snapshot=True)
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.001, stop_time=0.05,
                                  dynamic_steps=4, fast_threshold=80.)

dynamics.explain_time_advance()

# dynamics.get_history()

# %%
dynamics.single_compartment_react(time_step=0.005, stop_time=0.3,
                                 dynamic_steps=4, fast_threshold=100.)

dynamics.explain_time_advance()

#dynamics.get_history()

# %%
dynamics.single_compartment_react(time_step=0.01, stop_time=2.,
                                 dynamic_steps=4, fast_threshold=120.)

dynamics.explain_time_advance()

#dynamics.get_history()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(chemicals=["A", "B", "C"])

# %% [markdown]
# ### It might look like an equilibrium has been reached.  But NOT!  Verify the LACK of final equilibrium state:

# %%
dynamics.is_in_equilibrium()

# %% [markdown]
# ## Not surprisingly, none of the reactions of this physically-impossible hypothetical system are in equilibrium
# ### Even though the concentrations don't change, the cause of that is NOT equilibrium in the reactions - but rather a balancing out of consuming and replenishing across reactions
# #### Consider, for example, the concentrations of `A` at the end time, and its contributions from _individual_ reactions affecting `A`, as available from the diagnostic data:

# %%
dynamics.get_diagnostic_data(rxn_index=0, tail=1)

# %%
dynamics.get_diagnostic_data(rxn_index=2, tail=1)

# %% [markdown]
# ### Looking at the last row from each of the 2 dataframes above, one case see that, at every reaction cycle, [A] gets reduced by 0.914286 by the reaction `A <-> B`, while simultaneously getting increased by the SAME amount by the (fictional) reaction `C <-> A`.   
# ### Hence, the concentration of A remains constant - but neither reaction is in equilibrium!

# %%

# %% [markdown]
# # PART 2 - Let's restore the Laws of Physics!

# %%
chem_data.describe_reactions()

# %%
dynamics.clear_reactions()       # Let's start over with the reactions  (without affecting the data from the reactions)

# %%
# For the reactions A <-> B, and B <-> C, everything the same as before
chem_data.add_reaction(reactants="A", products="B",
                       forward_rate=9., reverse_rate=3.)

# Reaction , also favored energetically
chem_data.add_reaction(reactants="B", products="C",
                       forward_rate=8., reverse_rate=4.)

# %%
chem_data.describe_reactions()

# %%
# But for the reaction C <-> A, this time we'll "bend the knee" to the laws of thermodynamics!
# We'll use the same forward rate as before, but we'll let the reverse rate be picked by the system, 
# based of thermodynamic data consistent with the previous 2 reactions : i.e. an energy difference of -(-2,723.41 - 1,718.28) = +4,441.69
chem_data.add_reaction(reactants="C" , products="A",
                       forward_rate=3., Delta_G=4441.69)   # Notice the the positive Delta H: we're going from C, to the higher-energy level of "A"

# %%
chem_data.describe_reactions()

# %% [markdown]
# # Notice how, now that we're again following the laws of thermodynamics, the last reaction is mostly IN REVERSE (low K < 1), as it ought to be!

# %% [markdown]
# ### Now, let's continue with this "legit" set of reactions, from where we left off in our fantasy world:

# %%
dynamics.single_compartment_react(time_step=0.008, stop_time=4.,
                                 dynamic_steps=4, fast_threshold=120.)

dynamics.explain_time_advance()

#dynamics.get_history()

# %%
dynamics.plot_curves(chemicals=["A", "B", "C"])

# %% [markdown]
# ### Notice how [A] drops at time t=2, when we re-enact the Laws of Physics, because A no longer receives the extra boost from the previous mostly-forward (and thus physically-impossible given the unfavorable energy levels!) reaction `C <-> A`.   
# ### Back to the real world, that (energetically unfavored) reaction now mostly goes IN REVERSE; hence, the boost in [C] as well

# %% [markdown]
# ### Now, we have a REAL equilibrium!

# %%
dynamics.is_in_equilibrium()

# %% [markdown]
# ### The fact that individual reactions are now in actual, real equilibrium, can be easily seen from the last rows in the diagnostic data.  Notice all the delta-concentration values at the final times are virtually zero:

# %%
dynamics.get_diagnostic_data(rxn_index=0, tail=1)

# %%
dynamics.get_diagnostic_data(rxn_index=1, tail=1)

# %%
dynamics.get_diagnostic_data(rxn_index=2, tail=1)

# %%
