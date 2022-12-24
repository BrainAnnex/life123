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
# ### Adaptive variable time resolution on reaction A <-> B  between 2 species,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# Same as the experiment _"react_1"_ , but with an adaptive variable time scale
#
# LAST REVISED: Dec. 21, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from modules.reactions.reaction_data import ReactionData as chem
from modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.express as px
from modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # Initialize the System
# Specify the chemicals and the reactions

# %% tags=[]
# Specify the chemicals
chem_data = chem(names=["A", "B"])

# Reaction A <-> B , with 1st-order kinetics in both directions
chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

print("Number of reactions: ", chem_data.number_of_reactions())

# %%
chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# # Start the simulation

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)

# %%
# Initial concentrations of all the chemicals, in index order
dynamics.set_conc([10., 50.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.single_compartment_react(time_step=0.1, n_steps=11,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_step=2)      
                                  # Accepting the default:  fast_threshold=5

# %% [markdown]
# ## The argument _dynamic_step=2_ splits the time steps in 2 whenever the reaction is "fast" (as determined using fast_threshold=5)

# %%
df = dynamics.history.get()
df

# %% [markdown]
# ## Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
#
# * For example, upon completing the half step to t=0.30, i.e. **going from 0.25 to 0.30**, the last change in [A] was (21.508301 - 20.677734) = 0.830567  
# Relative to the [A] baseline of 20.677734, that change is **4.02%**.  The DEFAULT threshold for a reaction to be considered fast is 5% per full step, for any of the involved chemicals.  For a half step, that corresponds to 2.5%... and abs(4.02%) is LARGER than that.  
# The reaction is therefore marked "FAST" (as it has been so far), and the simulation then proceeds in a half step, to t=0.35

# %% [markdown]
# * (Note: at t=0, in the absence of any simulation data, ALL reactions are _assumed_ to be fast)

# %% [markdown]
# * By contrast, upon completing the half step to t=0.40, i.e. **going from 0.35 to 0.40**, the following changes occur in [A] and [B]:  

# %%
df.iloc[7]

# %%
df.iloc[8]

# %%
s_0_35 = df.iloc[7][['A', 'B']].to_numpy()
s_0_35     # Concentrations of A and B at t=0.35

# %%
s_0_40 = df.iloc[8][['A', 'B']].to_numpy()
s_0_40     # Concentrations of A and B at t=0.40

# %%
(s_0_40 - s_0_35) / s_0_35 * 100

# %% [markdown]
# The DEFAULT threshold for a reaction to be considered fast is 5% per full step, for any of the involved chemicals.  
# For a half step, that corresponds to 2.5%.   
# BOTH A's change of abs(2.11%) AND B's change of abs(-1.23%) are SMALLER than that.   
# The reaction is therefore marked "SLOW", and the simulation then proceeds in a _full step_ (i.e. a more relaxed time resolution), to t=0.50

# %% [markdown]
# ### Check the final equilibrium

# %%
dynamics.get_system_conc()

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order reactions),
#  the systems settles in the following equilibrium:
#
# [A] = 23.99316406
#  
# [B] = 36.00683594
#

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(rxn_index=0, conc=dynamics.get_conc_dict())

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="Reaction A <-> B .  Changes in concentrations with time",
              color_discrete_sequence = ['navy', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# ## Note how the left-hand side of this plot is much smoother than it was in experiment "react_1", where no adaptive time substeps were used!

# %%
