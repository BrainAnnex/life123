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
# ## Association/Dissociation reaction A + B <-> C
# #### with 1st-order kinetics for each species, taken to equilibrium
# (Adaptive variable time resolution is used)
#
# _See also the experiment "1D/reactions/reaction_4"_ 
#
# LAST REVISED: Dec. 29, 2022

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
chem_data = chem(names=["A", "B", "C"])

# Reaction A + B <-> C , with 1st-order kinetics for each species
chem_data.add_reaction(reactants=["A" , "B"], products=["C"],
                       forward_rate=5., reverse_rate=2.)

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
dynamics.set_conc([10., 50., 20.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.verbose_list = [2]
dynamics.single_compartment_react(time_step=0.004, reaction_duration=0.06,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_step=2)      
                                  # Accepting the default:  fast_threshold=5

# %% [markdown]
# ### Note: the argument _dynamic_step=2_ splits the time steps in 2 whenever the reaction is "fast" (as determined using fast_threshold=5)

# %%
df = dynamics.history.get()
df

# %% [markdown]
# ### Notice how the reaction proceeds in smaller steps in the early times, when the concentrations are changing much more rapidly

# %% [markdown]
# ## Note: "A" (now largely depleted) is largely the limiting reagent

# %% [markdown]
# ### Check the final equilibrium

# %%
dynamics.get_system_conc()

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order reactions),  
#  the systems settles in the following equilibrium:
#
# [A] = 0.29487741 , [B] = 40.29487741 , [C] = 29.70512259
#

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(rxn_index=0, conc=dynamics.get_conc_dict())

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "B", "C"], 
              title="Reaction A + B <-> C .  Changes in concentrations with time",
              color_discrete_sequence = ['red', 'violet', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
df.loc[10]

# %%
df.iloc[10]

# %%
list(df.loc[10])

# %%
list(df.loc[11])

# %%
a2 = df.loc[11].to_numpy()
a2

# %%
a1 = df.loc[10].to_numpy()
a1

# %%
a2 = df.iloc[11][['A', 'B', 'C']].to_numpy()
a2

# %%
a1 = df.iloc[10][['A', 'B', 'C']].to_numpy()
a1

# %%
a2 - a1

# %%
(a2 - a1) / a1 * 100.

# %%
abs(a2 - a1) / a1 * 100.

# %%
max(abs(a2 - a1) / a1 * 100.)

# %%
i = 0

# %%
before = df.iloc[i][['A', 'B', 'C']].to_numpy()
before

# %%
after = df.iloc[i+1][['A', 'B', 'C']].to_numpy()
after

# %%
max(abs(after - before) / before * 100.)

# %%
for i in range(22):
    print(f"---- {i} ----")
    before = df.iloc[i][['A', 'B', 'C']].to_numpy()
    print(before)
    after = df.iloc[i+1][['A', 'B', 'C']].to_numpy()
    print(after)
    print(after - before)
    print(abs(after - before) / before * 100.)
    largest_rel_change = max(abs(after - before) / before * 100.)
    print(largest_rel_change)

# %%
df

# %%
debug_df = dynamics.debug_data.get()
debug_df

# %%
i = 0

delta = debug_df.iloc[i][['A', 'B', 'C']].to_numpy()
print(delta)

baseline = df.iloc[i][['A', 'B', 'C']].to_numpy()
print(baseline)

ratio = delta / baseline * 100.
print(ratio)
print(max(abs(ratio)))

# %%
i = 1

delta = debug_df.iloc[i][['A', 'B', 'C']].to_numpy()
print(delta)

baseline = df.iloc[i][['A', 'B', 'C']].to_numpy()
print(baseline)

ratio = delta / baseline * 100.
print(ratio)
print(max(abs(ratio)))

# %%
for i in range(21):
    print(f"---- {i} ----")
    debug_time = debug_df.iloc[i]['TIME']
    print(f"debug_time: {debug_time:.5g} (Start of main t interval)")

    time_subdivision = debug_df.iloc[i]['time_subdivision']
    print(f"time_subdivision: {time_subdivision}")
    
    delta = debug_df.iloc[i][['A', 'B', 'C']].to_numpy()
    print("Delta:", delta)

    baseline = df.iloc[i][['A', 'B', 'C']].to_numpy()
    print("Baseline:", baseline)

    ratio = delta / baseline * 100.
    print("Ratio:", ratio)
    print("Max abs:", max(abs(ratio))) 
    print("Comparing the above against ", 5/time_subdivision)
    if max(abs(ratio)) > 5/time_subdivision:
        print("FAST")
    else:
        print("Slow")

# %%
