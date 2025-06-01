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
# ### One-bin `2A <-> 3B` reaction, with 1st-order kinetics in both directions, taken to equilibrium
#
# Diffusion not applicable (just 1 bin)

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

# Reaction 2A <-> 3B , with 1st-order kinetics in both directions
reactions.add_reaction(reactants=[(2,"A",1)], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)
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

# %% [markdown]
# ### <a name="sec_2_first_step"></a>First step

# %%
# First step
bio.react(time_step=0.05, n_steps=1)
bio.describe_state()

# %% [markdown]
# We're taking a smaller first step than in experiment "reaction_2", to avoid over-shooting the equilibrium value with too large a step!

# %%

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.1, n_steps=100)

bio.describe_state()

# %% [markdown]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:  [A] = 16.25 , [B] = 40.625

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
bio.get_bin_history(bin_address=0)

# %% [markdown]
# # Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title="Reaction  2A <-> 3B . Concentrations at bin 0")

# %% [markdown]
# ### Notice the *early overshoots* (the time step is too large early in the simulation!)
# Variable, adaptive time steps are explored at length in the  _"reactions_single_compartment"_ experiments

# %%
