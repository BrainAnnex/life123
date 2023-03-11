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
# ### Adaptive time steps (variable time resolution) for reaction `A <-> B`,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# Same as the experiment _"react_1"_ , but with adaptive variable time steps
#
# LAST REVISED: Mar. 11, 2023

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
chem_data.add_reaction(reactants=["A"], products=["B"], 
                       forward_rate=3., reverse_rate=2.)

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
dynamics.set_conc([10., 50.])

dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(initial_step=0.1, end_time=1.2,
                                  variable_steps=True, thresholds={"low": 0.5, "high": 0.8}, 
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  )

# %% [markdown]
# ## The flag _variable_steps_ automatically adjusts up or down the time step,  whenever the changes of concentrations are, respectively, "slow" or "fast" (as determined using the specified _thresholds_ )

# %%
df = dynamics.get_history()   # The system's history, saved during the run of single_compartment_react()
df

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ## Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# ### That resulted from passing the flag _variable_steps_ to single_compartment_react()

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
(s_0_40 - s_0_35)

# %% [markdown]
# Again, the threshold we specified for a reaction to be considered fast is 100% per full step, for any of the involved chemicals.  
# For a half step, that corresponds to 50%, i.e. 0.5    
# BOTH A's change of abs(0.46) AND B's change of abs(-0.46) are SMALLER than that.   
# The reaction is therefore marked "SLOW", and the simulation then proceeds in a _full time step_ of 0.1 (i.e. a more relaxed time resolution), to t=0.50

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

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

# %% [markdown]
#

# %% [markdown]
# # Diagnostics of the run may be investigated as follows:

# %%
dynamics.get_diagnostic_conc_data()

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)      # For the 0-th reaction (the only reaction in our case)

# %%
