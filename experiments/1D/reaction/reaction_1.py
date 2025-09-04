# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# ### A simple `A <-> B` reaction between 2 species with initial uniform concentrations across 3 bins,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# Diffusion NOT taken into account
#
# See also the experiment `reactions_single_compartment/react_1`

# %% [markdown]
# ### TAGS :  "reactions 1D", "quick-start"

# %%
LAST_REVISED = "June 6, 2025"
LIFE123_VERSION = "1.0.0rc5"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path
  
from life123 import ChemData, BioSim1D, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# # Initialize the System

# %%
# Initialize the system
chem_data = ChemData(names=["A", "B"])              # Diffusion NOT taken into account
bio = BioSim1D(n_bins=3, chem_data=chem_data)       # We'll specify the reactions later

bio.set_uniform_concentration(chem_label="A", conc=10.)   # Same across all bins
bio.set_uniform_concentration(chem_label="B", conc=50.)   # Same across all bins

bio.describe_state()

# %%

# %% [markdown]
# ## Enable History

# %%
# Let's enable history - by default for all chemicals and all bins
bio.enable_history(take_snapshot=True)

# %%
bio.get_bin_history(bin_address=0)

# %%
bio.get_bin_history(bin_address=1)

# %%
bio.get_bin_history(bin_address=2)

# %%

# %%
# Specify the reaction
reactions = bio.get_reactions()

# Reaction A <-> B , with 1st-order kinetics in both directions
reactions.add_reaction(reactants="A", products="B", forward_rate=3., reverse_rate=2.)

reactions.describe_reactions()

# %%

# %% [markdown]
# ### First Reaction Step

# %%
# First step of reaction
bio.react(time_step=0.1, n_steps=1)
bio.describe_state()

# %% [markdown]
# NOTE: the concentration of the chemical species `A` is increasing, while that of `B` is decreasing.
# All bins have identical concentrations; so, there's no diffusion (and we're not attempting to compute it; didn't specify diffusion rates)  

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown]
# ### Several more reaction steps

# %%
# Several more steps
bio.react(time_step=0.1, n_steps=10)

bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %%
bio.get_bin_history(bin_address=2)

# %%

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order reactions),
#  the systems settles in the following equilibrium:
#

# %%
bio.bin_snapshot(bin_address=0)

# %%
# Verify that the reaction has reached equilibrium
bio.get_reaction_handler().is_in_equilibrium(conc=bio.bin_snapshot(bin_address=0))

# %%

# %% [markdown]
# ## Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title="Reaction  A <-> B . Concentrations at bin 0")

# %%
# Same plot, but with a smoothed line
bio.plot_history_single_bin(bin_address=0, 
                            title="Reaction  A <-> B . Concentrations at bin 0",
                            smoothed=True)

# %% [markdown]
# ## For more in-depth analysis of this reaction, see the experiment `reactions_single_compartment/react_1`

# %%
