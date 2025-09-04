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
# ## One-bin reaction `A <-> 2C + D`, with 1st-order kinetics for each species, taken to equilibrium
#
# Diffusion not applicable (just 1 bin)

# %% [markdown]
# ### TAGS :  "reactions 1D", "basic"

# %%
LAST_REVISED = "June 6, 2025"
LIFE123_VERSION = "1.0.0rc5"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename

from life123 import BioSim1D, ChemData, check_version

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
# Initialize the system; NOTE: Diffusion not applicable (just 1 bin)
chem_data = ChemData(names=["A", "C", "D"], plot_colors=['navy', 'violet', 'red'])

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_all_uniform_concentrations( [4., 7., 2.] )

bio.describe_state()

# %%
# Specify the reaction
reactions = bio.get_reactions()

# Reaction A <-> 2C + D , with 1st-order kinetics for each species
reactions.add_reaction(reactants="A", products=[(2, "C", 1) , "D"],
                       forward_rate=5., reverse_rate=2.)

reactions.describe_reactions()

# %%
# Send the plot to the HTML log file
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
# ### First step

# %%
# First step
bio.react(time_step=0.2, n_steps=1)
bio.describe_state()

# %% [markdown]
# ---   
#     Note:  the above values are quite INaccurate because of the large time step 0.2
#
#     For example, the value for the concentration of D (0.4) is a wild overshot from the initial 2.0 to the equilibrium value of 1.68941267
#        
#     A more precise calculation with bio.react(time_step=0.1, n_steps=2) gives conc_D(0.2) = 2.304
#        
#     An even more precise calculation with bio.react(time_step=0.05, n_steps=4) gives conc_D(0.2) = 1.69037202
#        
#     I.e. the system is almost at equilibrium already at t=0.2 !
# ---

# %%

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.05, n_steps=30)

bio.describe_state()

# %% [markdown]
# ## Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:  
# [A] = 4.31058733 , [C] = 6.37882534 , [D] = 1.68941267

# %%
# Verify that the reaction has reached equilibrium
bio.get_reaction_handler().is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
bio.get_bin_history(bin_address=0)

# %% [markdown]
# C and D get depleted, while A gets produced.
# A wild overshoot is present at t=0.2

# %% [markdown]
# # Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title_prefix="Reaction  `A + B <-> C`")

# %% [markdown]
# ### Notice the **wild overshoot** present at t=0.2 !  (Too large a time step, early in the reaction!)
# #### Variable, adaptive time steps are explored at length in the  _"reactions_single_compartment"_ experiments

# %%
