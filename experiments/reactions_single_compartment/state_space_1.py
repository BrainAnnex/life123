# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ### `A <-> 2B` elementary reversible decomposition reaction, taken to equilibrium.  
# ### Examine State Space trajectory, using [A] and [B] as state variables
#
# See also the experiment `1D/reaction/reaction_2`

# %% [markdown]
# ### TAGS :  "uniform compartment"

# %%
LAST_REVISED = "Mar. 13, 2026"
LIFE123_VERSION = "1.0.0rc7"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import check_version, UniformCompartment, ReactionKinetics, PlotlyHelper

import plotly.express as px
import plotly.graph_objects as go

# %%
uc = UniformCompartment(names=["A", "B"])
uc.set_conc(conc={"A": 40., "B": 0.})

# %%
# Elementary reaction A <-> 2B
uc.add_reaction(reactants="A", products=[(2,"B")], kF=3., kR=2.)

uc.describe_reactions()

# %%

# %%

# %% [markdown]
# ### Simulate the reaction to equilibrium

# %%
#uc.single_compartment_react(initial_step=0.005, n_steps=12, variable_steps=False)

# %%
uc.single_compartment_react(initial_step=0.005, duration=0.07)

# %%
df = uc.get_history()
df

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium()

# %%
uc.plot_history(colors=['navy', 'orange'])

# %%

# %% [markdown]
# ## Same data, but shown differently

# %%
fig0 = PlotlyHelper.plot_pandas(df=uc.get_history(), x_var="A", fields="B", 
                                title="State space of reaction A <-> 2B : [A] vs. [B]", 
                                colors="#C83778", show=True)

# %%
fig0 = PlotlyHelper.plot_pandas(df=uc.get_history(), x_var="B", fields="A", 
                                title="State space of reaction A <-> 2B : [A] vs. [B]", 
                                colors="#C83778", show=True)

# %%
# Now show the individual data points

df['SYSTEM TIME'] = round(df['SYSTEM TIME'], 3)    # To avoid clutter from too many digits, in the column

fig1 = px.scatter(data_frame=df, x="A", y="B",
                  title="Trajectory in State space: [A] vs. [B]",
                  hover_data=['SYSTEM TIME'])

fig1.update_traces(marker={"size": 6, "color": "#2FAC74"})    # Modify the style of the dots

# Add annotations (showing the System Time value) to SOME of the points, to avoid clutter
for ind in df.index:     # for each row in the Pandas dataframe
    label = df["SYSTEM TIME"][ind]
    if ind == 0:
        label = f"t={label}"
        
    label_x = ind*(-12)
    label_y = 20 + ind*9    # A greater y value here means further DOWN!!
        
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

# %%
# Combine the two above plots, while using the layout of fig1 (which includes the title and annotations)
all_fig = go.Figure(data=fig0.data + fig1.data, layout = fig1.layout)
all_fig.show()

# %% [markdown]
# ### **Note how the trajectory is progressively slowing down towards the dynamical system's "attractor" (equilibrium state of the reaction)**

# %%
