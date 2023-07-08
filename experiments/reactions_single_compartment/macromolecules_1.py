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
# LAST REVISED: July 7, 2023

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
chem.set_binding_site_affinity("M1", 3, "A", 3.7)
chem.set_binding_site_affinity("M1", 8, "B", 1.1)
chem.set_binding_site_affinity("M1", 15, "A", 9.2)

chem.set_binding_site_affinity("M2", 1, "C", 5.6)

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

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = ReactionDynamics(chem_data=chem)
dynamics.set_conc(conc={"A": 20., "B": 0., "C": 10.})
dynamics.describe_state()

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
