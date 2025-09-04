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
# ## One-bin reaction `2A <-> B`  
# ### COMPARING 1st-order and 2nd-order kinetics in *forward* direction; reverse direction always 1-st order
#
# Diffusion not applicable (just 1 bin)
#
# See also the experiment _"reactions_single_compartment/react_4"_ 

# %% [markdown]
# ### TAGS :  "reactions 1D", "under-the-hood"

# %%
LAST_REVISED = "June 5, 2025"
LIFE123_VERSION = "1.0.0rc5"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename
from life123 import BioSim1D, ChemData, check_version

from life123 import HtmlLog as log
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

# %% [markdown]
# # INITIALLY, with 1st-order kinetics in both directions

# %%
# Initialize the system; NOTE: Diffusion not applicable (just 1 bin)
chem_data = ChemData(names=["A", "B"], plot_colors=['turquoise', 'green'])

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_all_uniform_concentrations( [3., 5.] )

bio.describe_state()

# %%
# Specify the reaction
reactions = bio.get_reactions()

# Reaction  2A <-> B , FOR NOW with 1st-order kinetics in both directions
reactions.add_reaction(reactants=[(2, "A", 1)], products="B", forward_rate=5., reverse_rate=2.)

reactions.describe_reactions()

# %%
# Send a header and a plot to the HTML log file
log.write("Reaction 2A <-> B is 1st order in all species:",
          style=log.h2)

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
bio.get_reaction_handler().enable_diagnostics()   # To save diagnostic information for the simulation run, below

# %%
# First step
bio.react(time_step=0.02, n_steps=1)
bio.describe_state()

# %% [markdown]
# Small conc. changes so far:  [A] = 2.8 , [B] = 5.1

# %%
bio.get_bin_history(bin_address=0)

# %%
# Numerous more steps, to equilibrium
bio.react(time_step=0.02, n_steps=20)

bio.describe_state()

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the *1st order* reactions),
# the systems settles in the following equilibrium:  
# [A] = 2.16928427 , [B] = 5.41535786

# %%
# Verify that the reaction has reached equilibrium
bio.get_reaction_handler().is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
df = bio.get_bin_history(bin_address=0)
df

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title_prefix="Reaction `2A <-> B`")

# %% [markdown]
# A gets depleted, while B gets produced.

# %% [markdown]
# #### Let's verify that the stoichiometry is being respected

# %%
# We'll check the first two arrays of concentrations, from the run's history
arr0 = bio.get_reaction_handler().get_historical_concentrations(row=0, df=df)
arr1 = bio.get_reaction_handler().get_historical_concentrations(row=1, df=df)
arr0, arr1

# %%
bio.get_reaction_handler().get_diagnostics().stoichiometry_checker(rxn_index=0, 
                               conc_arr_before = arr0, 
                               conc_arr_after = arr1)

# %%

# %%

# %% [markdown]
# # STARTING OVER, this time with 2nd-order kinetics in the forward reaction

# %%
reactions = bio.get_reactions()

# %%
reactions.clear_reactions_data()

# %%
# Reaction  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction
reactions.add_reaction(reactants=[(2, "A", 2)], products="B", forward_rate=5., reverse_rate=2.)

# %%
reactions.describe_reactions()

# %%
# RESET the concentrations to their original values
bio.set_all_uniform_concentrations( [3., 5.] )

bio.describe_state()

# %%
# Save this manual concentration change into the ongoing history
bio.capture_snapshot(caption="RESET all concentrations to initial values")
bio.get_bin_history(bin_address=0)

# %%
# Send a header and a plot to the HTML log file
log.write("Reaction 2A <-> B is 2nd order in A, and 1st order in B:",
          style=log.h2)
reactions.plot_reaction_network("vue_cytoscape_2")

# %%

# %%
# First step
bio.react(time_step=0.02, n_steps=1)
bio.describe_state()

# %% [markdown]
# [A] = 1.6 , [B] = 5.7
# _(Contrast with the counterpart in the 1st order kinetics:  [A] = 2.8 , [B] = 5.1)_

# %%
bio.get_bin_history(bin_address=0)

# %%

# %%
# Numerous more steps
bio.react(time_step=0.02, n_steps=20)

bio.describe_state()

# %% [markdown]
# The systems settles in the following equilibrium:  [A] = 1.51554944 , [B] = 5.74222528

# %%
# Verify that the reaction has reached equilibrium
bio.get_reaction_handler().is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title_prefix="Reaction `2A <-> B` (the jump at 0.42 is the concentration reset)",
                            vertical_lines_to_add=[0.42])

# %%
df2 = bio.get_bin_history(bin_address=0)
df2

# %% [markdown]
# **Compared to first-order kinetics in A**, the (2nd order in A) reaction now takes place much more quickly, and proceeds to almost complete depletion of A

# %% [markdown]
# #### Let's verify that the stoichiometry is still being respected

# %%
# We'll check the first two arrays of concentrations, from the run's history
arr0 = bio.get_reaction_handler().get_historical_concentrations(row=22, df=df2)  # Row 22 is the conc. reset
arr1 = bio.get_reaction_handler().get_historical_concentrations(row=23, df=df2)
arr0, arr1

# %%
bio.get_reaction_handler().get_diagnostics().stoichiometry_checker(rxn_index=0, 
                                                       conc_arr_before = arr0, 
                                                       conc_arr_after = arr1)

# %%
