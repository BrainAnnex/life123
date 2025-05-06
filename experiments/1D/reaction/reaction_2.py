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
# ### One-bin `A <-> 3B` reaction, with 1st-order kinetics in both directions,
# ### taken to equilibrium
#
# Diffusion not applicable (just 1 bin)
#
# * [First Step](#reaction_2_sec_2_first_step)
# * [Numerous more steps](#reaction_2_sec_2)
# * [Equilibrium](#reaction_2_sec_2_equilibrium)

# %% [markdown]
# ### TAGS :  "reactions 1D"

# %%
LAST_REVISED = "May 4, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData, BioSim1D, check_version

from life123 import GraphicLog

# %%
check_version(LIFE123_VERSION)

# %%

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system
chem_data = ChemData(names=["A", "B"])     # NOTE: Diffusion not applicable (just 1 bin)

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_uniform_concentration(chem_index=0, conc=10.)
bio.set_uniform_concentration(chem_index=1, conc=50.)

bio.describe_state()

# %%
# Specify the reaction
reactions = bio.get_reactions()

# Reaction A <-> 3B , with 1st-order kinetics in both directions
reactions.add_reaction(reactants="A", products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

reactions.describe_reactions()

# %%
# Send the plot of the reaction network to the HTML log file
reactions.plot_reaction_network("vue_cytoscape_2")

# %%

# %%
# Let's enable history - by default for all chemicals and all bins
bio.enable_history(take_snapshot=True, caption="Initial state")

# %%
bio.get_bin_history(bin_address=0)

# %%

# %%

# %% [markdown]
# ### <a name="reaction_2_sec_2_first_step"></a>First step

# %%
# First step
bio.react(time_step=0.1, n_steps=1)
bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown]
# ### <a name="reaction_2_sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.1, n_steps=10)

bio.describe_state()

# %% [markdown]
# ### <a name="reaction_2_sec_2_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the equilibrium:   [A] = 14.54545455 , [B] = 36.36363636

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0))

# %%
bio.get_bin_history(bin_address=0)

# %% [markdown]
# Note how the simulation initially **OVERSHOT** the equilibrium values; the first step was too large!

# %% [markdown]
# # Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title="Reaction  A <-> B . Concentrations at bin 0")

# %%
# Same plot, but with smooth line
bio.plot_history_single_bin(bin_address=0, 
                            title="Reaction  A <-> B . Concentrations at bin 0",
                            smoothed=True)

# %% [markdown]
# The early **OVERSHOOTING** of the equilibrium values shows prominently in the last plot!

# %%
