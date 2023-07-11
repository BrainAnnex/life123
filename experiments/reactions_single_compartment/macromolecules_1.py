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
# LAST REVISED: July 10, 2023

# %% tags=[]
from src.modules.chemicals.chem_data import ChemData
from src.modules.reactions.reaction_dynamics import ReactionDynamics

# %%
# Initialize the system
chem = ChemData(names=["A", "B", "C"])


# %% [markdown]
# ### Explore methods to manage the data structure for macromolecules

# %%
chem.add_macromolecules("M1")
chem.get_macromolecules()

# %%
chem.set_binding_site_affinity("M1", 3, "A", 0.92)
chem.set_binding_site_affinity("M1", 8, "B", 1.1)
chem.set_binding_site_affinity("M1", 15, "A", 9.2)

chem.set_binding_site_affinity("M2", 1, "C", 5.6)
chem.set_binding_site_affinity("M2", 2, "A", 0.092)

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
dynamics.set_conc(conc={"A": 9.2, "B": 0., "C": 0.56})
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
# **[A] = 9.2** :   
#             * 10x the binding affinity of A to site 3 on M1 (occupancy 0.9)  
#             * same as the binding affinity of A to site 15 on M1 (occupancy 0.5)  
#             * 100x the binding affinity of A to site 2 on M2 (occupancy almost 1, i.e. nearly saturated)
#         
#             
# **[C] = 0.56** => 1/10 of the binding affinity of A to site 1 on M2 (occupancy 0.1)

# %%
