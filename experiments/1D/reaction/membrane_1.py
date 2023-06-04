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
# ### Exploring the data structures of MEMBRANES, and reactions in them 
# #### - with NO DIFFUSION
#
# LAST REVISED: May 28, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from src.modules.reactions.reaction_data import ChemData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

from src.life_1D.bio_sim_1d import BioSim1D

# %%
chem_data = chem(names=["A", "B", "C"])     # NOTE: Diffusion not done
bio = BioSim1D(n_bins=5, chem_data=chem_data)

bio.set_membranes(membrane_pos=[1])   # A single membrane, passing thru bin 1

bio.set_all_uniform_concentrations(conc_list=[4., 8., 12.])

bio.describe_state()

# %%
bio.set_bin_conc(bin_address=1, species_name="A", conc=10.)
bio.set_bin_conc(bin_address=1, species_name="A", conc=55., across_membrane=True)
bio.set_bin_conc(bin_address=1, species_name="B", conc=20.)
bio.set_bin_conc(bin_address=1, species_name="C", conc=30., both_sides=True)

bio.describe_state()

# %%
# Make the last bin match all the concentrations of the "post-membrane" section of bin 1
bio.set_bin_conc(bin_address=4, species_name="A", conc=55.)
bio.set_bin_conc(bin_address=4, species_name="C", conc=30.)

bio.describe_state()

# %%
# Reaction A + B <-> C , with 1st-order kinetics in both directions, mostly forward
bio.chem_data.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=8., reverse_rate=2.)

bio.chem_data.describe_reactions()

# %%
bio.react(time_step=0.002, n_steps=1)
bio.describe_state()

# %% [markdown]
# ### Note how (in the absence of diffusion, which we're neglecting) the concentrations on the "post-membrane side" of bin 1 continue to match those of bin 5

# %% [markdown]
# ## Now continue to reaction equilibrium

# %%
bio.react(time_step=0.002, n_steps=100)
bio.describe_state()

# %% [markdown]
# ### The system has now reached equilibrium
# ### in individual bins, which remain separate because we're NOT doing diffusion in this experiment

# %% [markdown]
# Verify the equilibrium in each of the active bins

# %%
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc={"A": 0.7932534523016195, "B": 4.793253452301622, "C": 15.206746547698396})
# A was largely the limiting reagent

# %%
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc={"A": 0.897094737056602, "B": 10.897094737056594, "C": 39.10290526294337})
# A was largely the limiting reagent

# %%
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc={"A": 47.20020986266437, "B": 0.20020986266437912, "C": 37.79979013733563})
# This time, with ample [A], the limiting reagent was B

# %%
