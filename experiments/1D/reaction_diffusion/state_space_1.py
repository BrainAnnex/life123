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
# ### One-bin A <-> 3B reaction, with 1st-order kinetics in both directions, taken to equilibrium
# ### Examine State Space trajectory, using [A] and [B] as state variables
#
# Based on reactions/reaction_2
#
# Diffusion not applicable (just 1 bin)
#
# LAST REVISED: Sep. 13, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %%
from experiments.get_notebook_info import get_notebook_basename

from modules.reactions.reaction_data import ReactionData as chem
from modules.reactions.reaction_dynamics import ReactionDynamics
from life_1D.bio_sim_1d import BioSim1D

import plotly.express as px
import plotly.graph_objects as go
from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system
chem_data = chem(names=["A", "B"])     # NOTE: Diffusion not applicable (just 1 bin)

rxn = Reactions(chem_data)

# Reaction A <-> 3B , with 1st-order kinetics in both directions
rxn.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

bio = BioSim1D(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)

bio.describe_state()

# %%
rxn.describe_reactions()

# %%
# Save the state of the concentrations of all species at bin 0
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %%
# Send the plot to the HTML log file
graph_data = rxn.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown] tags=[]
# ### To equilibrium

# %%
# Using smaller steps that in experiment reaction_2, to avoid the initial overshooting
bio.react(time_step=0.05, n_steps=10, snapshots={"frequency": 1, "sample_bin": 0})
bio.describe_state()

# %%
bio.get_history()

# %%
# Verify that the reaction has reached equilibrium
rxn.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
df = bio.get_history()

px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B"], 
              title="Reaction A <-> 3B .  Changes in [A] and [B] over time",
              color_discrete_sequence = ['navy', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"})

# %%
# Same data, but shown differently
fig0 = px.line(data_frame=bio.get_history(), x="A", y="B", 
              title="State space of reaction A <-> 3B : [A] vs. [B]",
              color_discrete_sequence = ['#C83778'],
              labels={"value":"concentration", "variable":"Chemical"})
fig0.show()

# %%
# Now show the individual data points

df['SYSTEM TIME'] = round(df['SYSTEM TIME'], 2)    # To avoid clutter from too many digits

fig1 = px.scatter(data_frame=df, x="A", y="B",
                  title="Trajectory in State space: [A] vs. [B]",
                  hover_data=['SYSTEM TIME'])

fig1.update_traces(marker={"size": 6, "color": "#2FAC74"})    # Modify the style of the dots

# Add annotations (showing the System Time value) to SOME of the points, to avoid clutter
for ind in df.index:
    label = df["SYSTEM TIME"][ind]
    if ind == 0:
        label = f"t={label}"
        
    label_x = ind*16
    label_y = 20 + ind*8    # A greater y value here means further DOWN!!
        
    if (ind <= 3) or (ind%2 == 0):
        fig1.add_annotation(x=df["A"][ind], y=df["B"][ind],
                            text=label,
                            font=dict(
                                size=10,
                                color="grey"
                            ),
                            showarrow=True, arrowhead=0, ax=label_x, ay=label_y, arrowcolor="#b0b0b0",
                            bordercolor="#c7c7c7")

fig1.show()

# %%
# Combine the two above plots, while using the layout of fig1 (which includes the title and annotations)
all_fig = go.Figure(data=fig0.data + fig1.data, layout = fig1.layout)
all_fig.show()

# %% [markdown]
# #### Note how the trajectory is progressively slowing down towards the dynamical system's "attractor" (equilibrium state of the reaction)

# %%
