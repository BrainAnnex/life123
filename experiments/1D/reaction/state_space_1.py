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
# ### One-bin A <-> 3B reaction, taken to equilibrium.  
# #### A hypothetical scenario with 1st-order kinetics in both directions.  
# ### Examine State Space trajectory, using [A] and [B] as state variables
#
# Based on experiment `1D/reaction/reaction_2`
#
# Diffusion not applicable (just 1 bin).
#
# This is the 1D version of the single-compartment reaction by the same name.

# %% [markdown]
# ### TAGS :  "reactions 1D"

# %%
LAST_REVISED = "Dec. 16, 2024"
LIFE123_VERSION = "1.0-rc.1"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename

from life123 import UniformCompartment,  BioSim1D

import plotly.express as px
import plotly.graph_objects as go
from life123 import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system.  NOTE: Diffusion not applicable (just 1 bin)
uc = UniformCompartment(names=["A", "B"])


# Reaction A <-> 3B , with 1st-order kinetics in both directions
uc.add_reaction(reactants="A", products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

uc.describe_reactions()

# %%
bio = BioSim1D(n_bins=1, reaction_handler=uc)

bio.set_uniform_concentration(species_name="A", conc=10.)
bio.set_uniform_concentration(species_name="B", conc=50.)

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %%
# Send the plot to the HTML log file
uc.plot_reaction_network("vue_cytoscape_2")

# %%

# %% [markdown] tags=[]
# ### To equilibrium

# %%
# Using smaller steps that in experiment reaction_2, to avoid the initial overshooting
bio.react(time_step=0.05, n_steps=10, snapshots={"frequency": 1, "sample_bin": 0})

# %%
bio.describe_state()

# %%
bio.get_history()

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
df = bio.get_history()

px.line(data_frame=df, x="SYSTEM TIME", y=["A", "B"], 
              title="Reaction A <-> 3B .  Changes in [A] and [B] over time",
              color_discrete_sequence = ['darkturquoise', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})

# %%

# %% [markdown]
# ## Same data, but shown differently

# %%
fig0 = px.line(data_frame=bio.get_history(), x="A", y="B", 
              title="State space of reaction A <-> 3B : [A] vs. [B]",
              color_discrete_sequence = ['#C83778'],
              labels={"value":"concentration", "variable":"Chemical"})
fig0.show()

# %%
# Now show the individual data points

df['SYSTEM TIME'] = round(df['SYSTEM TIME'], 2)    # To avoid clutter from too many digits, in the column

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
# ### Note how the trajectory is progressively slowing down towards the dynamical system's "attractor" (equilibrium state of the reaction)

# %%
