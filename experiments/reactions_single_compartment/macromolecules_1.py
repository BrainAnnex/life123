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
# ## Macromolecules : Binding Affinity
#
# LAST REVISED: July 12, 2023

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
# ### Explore methods to manage the data structure for macromolecules

# %%
chem.add_macromolecules("M1")
chem.get_macromolecules()

# %%
chem.set_binding_site_affinity("M1", 3, "A", 1.0)
chem.set_binding_site_affinity("M1", 8, "B", 3.2)
chem.set_binding_site_affinity("M1", 15, "A", 10.0)

chem.set_binding_site_affinity("M2", 1, "C", 5.6)
chem.set_binding_site_affinity("M2", 2, "A", 0.01)

# %%
chem.show_binding_affinities()        # Review the values the had given for the binding affinities

# %%
chem.get_binding_sites("M1")

# %%
chem.get_binding_sites_and_ligands("M1")

# %%
chem.get_binding_sites("M2")

# %%
chem.get_binding_sites_and_ligands("M2")

# %%
chem.get_binding_site_affinity(macromolecule="M1", site_number=3)   # A "NamedTuple" gets returned

# %%
aff = chem.get_binding_site_affinity(macromolecule="M1", site_number=8)
aff

# %%
aff.chemical

# %%
aff.affinity

# %%

# %% [markdown]
# ### Start setting up the dynamical system

# %%
dynamics = ReactionDynamics(chem_data=chem)

# %%
dynamics.set_macromolecules()      # By default, set counts to 1 for all the registered macromolecules

# %% [markdown]
# ### Inspect some class attributes

# %%
dynamics.macro_system

# %%
dynamics.macro_system_state

# %%
dynamics.describe_state()

# %%

# %% [markdown]
# ### Set the initial concentrations of all the ligands

# %%
dynamics.set_conc(conc={"A": 10., "B": 0., "C": 0.56})
dynamics.describe_state()

# %% [markdown]
# ### Adjust the fractional occupancy of the various sites on the macromolecules, based on the current ligand concentrations

# %%
dynamics.update_occupancy()

# %%
dynamics.describe_state()

# %%
dynamics.chem_data.show_binding_affinities()        # Review the values the had given for the binding affinities

# %% [markdown]
# #### Notes:
# **[B] = 0** => Occupancy of binding site 8 on M1 is also zero
#
# **[A] = 10.0** :   
#             * 10x the binding affinity of A to site 3 on M1 (occupancy 0.9)  
#             * same as the binding affinity of A to site 15 on M1 (occupancy 0.5)  
#             * 100x the binding affinity of A to site 2 on M2 (occupancy almost 1, i.e. nearly saturated)
#         
#             
# **[C] = 0.56** => 1/10 of the binding affinity of A to site 1 on M2 (occupancy 0.1)

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
# #### Note how all the various binding sites for ligand A, across all macromolecules, now have a different value for the fractional occupancy (very close to 1 because of the large value of [A] relative to each of the binding affinities for A 

# %%

# %%
history = MovieTabular(parameter_name="[A]")

# %%
print(history)

# %%
print(history)

# %%

# %%
history.clear_dataframe()

# %%
start = 0.001
stop = 200
num_points = 100  # Number of points you want

log_values = np.logspace(np.log10(start), np.log10(stop), num=num_points)

print(log_values)

# %%
for i in log_values:
    A_conc = i
    dynamics.set_chem_conc(conc=A_conc, species_name="A", snapshot=False)
    dynamics.update_occupancy()
    history.store(A_conc, {"M1 site 3": dynamics.get_occupancy(macromolecule="M1", site_number=3), 
                           "M1 site 15": dynamics.get_occupancy(macromolecule="M1", site_number=15), 
                            "M2 site 2": dynamics.get_occupancy(macromolecule="M2", site_number=2)})

# %%
history.get_dataframe()

# %%
df = history.get_dataframe()

# %%
fig = px.line(data_frame=df, 
              x="[A]", y=["M2 site 2", "M1 site 3", "M1 site 15"],
              color_discrete_sequence = ["seagreen", "purple", "darkorange"],
              title="<b>Fractional Occupancy as a function of Ligand Concentration</b>",
              labels={"value":"Fractional Occupancy", "variable":"Binding site"})

fig.add_hline(y=0.5, line_width=1, line_dash="dot", line_color="gray")

fig.show()

# %%
fig = px.line(data_frame=df, 
                x="[A]", y=["M2 site 2", "M1 site 3", "M1 site 15"],
                log_x=True, range_x=[start,200],              
                title="Fractional Occupancy as a function of Ligand Concentration <br>(log plot on x-axis. Showing 0.1/0.5/0.9 horizontals)",
                labels={"value":"Fractional Occupancy", "variable":"Binding site"})

fig.add_hline(y=0.1, line_width=1, line_dash="dot", line_color="gray")
fig.add_hline(y=0.5, line_width=1, line_dash="dot", line_color="gray")
fig.add_hline(y=0.9, line_width=1, line_dash="dot", line_color="gray")

fig.show()

# %%
