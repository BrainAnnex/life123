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
# ### `A <-> 3B` reaction, taken to equilibrium.  
# #### A hypothetical scenario with 1st-order kinetics in both directions. 
# ### Examine State Space trajectory, using [A] and [B] as state variables
#
# 1st-order kinetics in both directions
#
# See also the experiment `1D/reaction/reaction_2`
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta36)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData
from life123 import UniformCompartment

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
# Initialize the system
chem_data = ChemData(names=["A", "B"])


# Reaction A <-> 3B , with 1st-order kinetics in both directions
chem_data.add_reaction(reactants="A", products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

chem_data.describe_reactions()

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc={"A": 10., "B": 50.},
                  snapshot=True)
dynamics.describe_state()

# %%
# Send the plot to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %% [markdown] tags=[]
# ### To equilibrium

# %%
dynamics.single_compartment_react(initial_step=0.05, n_steps=10, variable_steps=False)

# %%
df = dynamics.get_history()
df

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%
dynamics.plot_history(colors=['navy', 'orange'])

# %%

# %% [markdown]
# ## Same data, but shown differently

# %%
fig0 = px.line(data_frame=dynamics.get_history(), x="A", y="B", 
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
for ind in df.index:     # for each row in the Pandas dataframe
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
# ### **Note how the trajectory is progressively slowing down towards the dynamical system's "attractor" (equilibrium state of the reaction)**

# %%
