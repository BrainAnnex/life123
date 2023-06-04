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
# ## A simple `A <-> B` reaction between 2 species
# with 1st-order kinetics in both directions, taken to equilibrium
#
# See also the experiment _"1D/reactions/reaction_1"_ ; this is the "single-compartment" version of it.
#
# LAST REVISED: May 16, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ChemData as chem
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

# %% tags=[]
# Initialize the reaction
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

# %%
dynamics.describe_state()

# %%
dynamics.get_history()

# %% [markdown] tags=[]
# ## Start the reaction

# %%
# First step of reaction
dynamics.single_compartment_react(initial_step=0.1, n_steps=1,
                                  snapshots={"initial_caption": "first reaction step"})

# %%
dynamics.get_history()

# %%
# Numerous more fixed steps
dynamics.single_compartment_react(initial_step=0.1, n_steps=10,
                                  snapshots={"initial_caption": "2nd reaction step",
                                             "final_caption": "last reaction step"})

# %% [markdown]
# #### NOTE: for demonstration purposes, we're using FIXED time steps...  Typically, one would use the option for automated variable time steps (see experiment `react 2`)

# %%
dynamics.get_history()

# %% [markdown]
# ### Check the final equilibrium

# %%
dynamics.get_system_conc()

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1:1 stoichiometry, and the 1st order reactions), the systems settles in the following equilibrium:
#
# [A] = 23.99316406
#  
# [B] = 36.00683594
#

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown]
# ### Note that, because of the high initial concentration of B relative to A, the overall reaction has proceeded **in reverse**

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=['blue', 'orange'])

# %% [markdown]
# ### Note the raggedness of the left-side (early times) of the curves.  
# ### In experiment `react_2` this simulation gets repeated with an adaptive variable time resolution that takes smaller steps at the beginning, when the reaction is proceeding faster   
# ### By contrast, here we used FIXED steps (see below), which is generally a bad approach

# %%
dynamics.plot_curves(colors=['blue', 'orange'], show_intervals=True)

# %%
df = dynamics.get_history()

# %%
df

# %% [markdown]
# ## Now investigate A_dot, i.e. d[A]/dt

# %%
A = list(df.A)

# %%
A

# %%
len(A)

# %%
A_dot = np.gradient(A, 0.1)      # 0.1 is the constant step size

# %%
A_dot

# %%
df['A_dot'] = A_dot     # Add a column to the dataframe

# %%
df

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "A_dot"], 
              title="Changes in concentration of A with time (blue) , and changes in its rate of change (brown)",
              color_discrete_sequence = ['navy', 'brown'],
              labels={"value":"concentration (blue) /<br> concentration change per unit time (brown)", "variable":"Chemical"})
fig.show()

# %% [markdown]
# ### At t=0, [A]=10 and [A] has a high rate of change (70);  
# ### as the system approaches equilibrium, [A] approaches a value of 24, and its rate of change decays to zero.
#
# The curves are jagged because of limitations of numerically estimating derivatives, as well as the large time steps taken (especially in the early times, when there's a lot of change.)

# %%
