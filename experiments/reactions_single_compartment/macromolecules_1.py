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
# ## Macromolecules : Binding Affinity and Fractional Occupancy
#
# LAST REVISED: July 14, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from src.modules.chemicals.chem_data import ChemData
from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.modules.movies.movies import MovieTabular

import numpy as np

import plotly.express as px

# %%
# Initialize the system
chem = ChemData(names=["A", "B", "C"])


# %% [markdown]
# ## Explore methods to manage the data structure for macromolecules

# %%
chem.add_macromolecules("M1")
chem.get_macromolecules()

# %%
chem.set_binding_site_affinity("M1", site_number=3,  ligand="A", Kd=1.0)
chem.set_binding_site_affinity("M1", site_number=8,  ligand="B", Kd=3.2)
chem.set_binding_site_affinity("M1", site_number=15, ligand="A", Kd=10.0)

chem.set_binding_site_affinity("M2", site_number=1, ligand="C", Kd=5.6)    # "M2" will get automatically added
chem.set_binding_site_affinity("M2", site_number=2, ligand="A", Kd=0.01)

# %%
chem.show_binding_affinities()        # Review the values we have given for the binding affinities

# %%
chem.get_binding_sites("M1")

# %%
chem.get_binding_sites_and_ligands("M1")

# %%
chem.get_binding_sites("M2")

# %%
chem.get_binding_sites_and_ligands("M2")

# %%
aff = chem.get_binding_site_affinity(macromolecule="M2", site_number=1)   # A "NamedTuple" gets returned
aff

# %%
aff.chemical

# %%
aff.Kd

# %%

# %% [markdown]
# ## Start setting up the dynamical system

# %%
dynamics = ReactionDynamics(chem_data=chem)

# %%
dynamics.set_macromolecules()      # By default, set counts to 1 for all the registered macromolecules

# %%
dynamics.describe_state()

# %% [markdown]
# ### Inspect some class attributes (not to be directly modified by the end user!)

# %%
dynamics.macro_system

# %%
dynamics.macro_system_state

# %% [markdown]
# ### Set the initial concentrations of all the ligands

# %%
dynamics.set_conc(conc={"A": 10., "B": 0., "C": 0.56})
dynamics.describe_state()

# %% [markdown]
# ### Determine and adjust the fractional occupancy of the various sites on the macromolecules, based on the current ligand concentrations

# %%
dynamics.update_occupancy()

# %%
dynamics.describe_state()

# %%
dynamics.chem_data.show_binding_affinities()        # Review the values we had given for the dissociation constants

# %% [markdown]
# #### Notes:
# **[B] = 0** => Occupancy of binding site 8 of M1 is also zero
#
# **[A] = 10.0** :   
#             * 10x the dissociation constant of A to site 3 of M1 (resulting in occupancy 0.9)  
#             * same as the dissociation constant of A to site 15 of M1 (occupancy 0.5)  
#             * 1,000x the dissociation constant of A to site 2 of M2 (occupancy almost 1, i.e. nearly saturated)
#         
#             
# **[C] = 0.56** => 1/10 of the dissociation constant of C to site 1 of M2 (occupancy 0.1)

# %%

# %% [markdown]
# ### Adjust the concentration of one ligand, [A], and update all the fractional occupancies accordingly

# %%
dynamics.set_chem_conc(conc=1000., species_name="A", snapshot=False)

# %%
dynamics.update_occupancy()

# %%
dynamics.describe_state()

# %% [markdown]
# #### Note how all the various binding sites for ligand A, across all macromolecules, now have a different value for the fractional occupancy (very close to 1 because of the large value of [A] relative to each of the dissociation constants for A.)
# The fractional occupancies for the other ligands (B and C) did not change

# %%

# %% [markdown]
# ### Sweep the values of [A] across a wide range, and compute/store how the fractional occupancies of A change

# %%
history = MovieTabular(parameter_name="[A]")  # A convenient way to store a sequence of "state snapshots" as a Pandas dataframe

# %%
print(history)

# %%
# Generate a sweep of [A] values along a log scale, from very low to very high (relative to the dissociation constants)
start = 0.001
stop = 200.
num_points = 100

log_values = np.logspace(np.log10(start), np.log10(stop), num=num_points)

print(log_values)

# %%
# Set [A] to each of the above values in turn, and determine/store the applicable fractional occupancies (for the sites where A binds)
for A_conc in log_values:
    dynamics.set_chem_conc(conc=A_conc, species_name="A", snapshot=False)
    dynamics.update_occupancy()
    history.store(A_conc, {"M1 site 3": dynamics.get_occupancy(macromolecule="M1", site_number=3), 
                           "M1 site 15": dynamics.get_occupancy(macromolecule="M1", site_number=15), 
                           "M2 site 2": dynamics.get_occupancy(macromolecule="M2", site_number=2)})

# %%
df = history.get_dataframe()
df

# %%
# Plot each of the fractional occupancies as a function of [A]

fig = px.line(data_frame=df, 
              x="[A]", y=["M2 site 2", "M1 site 3", "M1 site 15"],
              color_discrete_sequence = ["seagreen", "purple", "darkorange"],
              title="<b>Fractional Occupancy as a function of Ligand Concentration</b>",
              labels={"value":"Fractional Occupancy", "variable":"Binding site"})

fig.add_hline(y=0.5, line_width=1, line_dash="dot", line_color="gray")  # Horizontal line at 50% occupancy

fig.show()

# %%
import plotly.graph_objects as go 

# %%
# Same plot, but use a log scale for [A]

fig = px.line(data_frame=df, 
              x="[A]", y=["M2 site 2", "M1 site 3", "M1 site 15"],
              color_discrete_sequence = ["seagreen", "purple", "darkorange"],           
              log_x=True, range_x=[start,200],              
              title="<b>Fractional Occupancy as a function of Ligand Concentration</b> <br>(log plot on x-axis. Highlighting 0.1/0.5/0.9 occupancies)",
              labels={"value":"Fractional Occupancy", "variable":"Binding site"})

# Horizontal lines
fig.add_hline(y=0.1, line_width=1, line_dash="dot", line_color="gray")
fig.add_hline(y=0.5, line_width=1, line_dash="dot", line_color="gray")
fig.add_hline(y=0.9, line_width=1, line_dash="dot", line_color="gray")

# Annotations (x values adjusted for the log scale)
fig.add_annotation(x=-1.715, y=0.65, 
                   text="Dissociation constant: 0.01<br>(HIGHER Binding Affinity)", font=dict(size=12, color="seagreen"), showarrow=True, ax=-100, ay=-20)
fig.add_annotation(x=0.3, y=0.65, 
                   text="Dissociation constant: 1", font=dict(size=12, color="purple"), showarrow=True, ax=-100, ay=-20)
fig.add_annotation(x=0.6, y=0.3, 
                   text="Dissociation constant: 10<br>(LOWER Binding Affinity)", font=dict(size=12, color="darkorange"), showarrow=True, ax=100, ay=10)

fig.add_annotation(x=1.8, y=0.51, 
                   text="50% OCCUPANCY", font=dict(size=14, color="gray"), bgcolor="white", opacity=0.8, showarrow=False)

# Customize y-axis tick values
additional_y_values = [0.1, 0.5, 0.9]  # Additional values to show on y-axis
fig.update_layout(yaxis={"tickvals": list(fig.layout.yaxis.domain) + additional_y_values})
 
# Add scatter points (dots) to the plot, to highlight the 0.1/0.5/0.9 occupancies
fig.add_scatter(x=[0.01, 1, 10,     0.001, 0.1, 1.,     0.1, 10, 100 ], 
                y=[0.5, 0.5, 0.5,    0.1, 0.1, 0.1,     0.9, 0.9, 0.9], 
                mode="markers", marker={"color": "red"}, name="key points")
    
fig.show()

# %% [markdown]
# #### When the binding affinity is lower (i.e. higher Dissociation Constant, Kd, orange curve), it takes higher ligand concentrations to attain the same fractional occupancies

# %% [markdown]
# Note that fractional occupancy 0.1 occurs at ligand concentrations of 1/10 the dissociation constant (Kd);   
# occupancy 0.5 occurs at ligand concentrations equals to the dissociation constant;  
# occupancy 0.9 occurs at ligand concentrations of 10x the dissociation constant.  

# %% [markdown]
# ## The above simulation captures what's shown on Fig. 3A of 
# #### https://doi.org/10.1146/annurev-cellbio-100617-062719 
# ("Low-Affinity Binding Sites and the Transcription Factor Specificity Paradox in Eukaryotes"), a paper that guided this simulation

# %% [markdown]
# #### In upcoming versions of Life123, the fractional occupancy values will regulate the rates of reactions catalyzed by the macromolecules...

# %%
