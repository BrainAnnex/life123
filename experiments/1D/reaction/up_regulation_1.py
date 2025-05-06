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
# ### `A` up-regulates `B` , by being *the limiting reagent* in the reaction:     
# ### `A + X <-> 2B` (mostly forward), where `X` is plentiful
# 1st-order kinetics.   
# If [A] is low, [B] remains low, too.  Then, if [A] goes high, then so does [B].  However, at that point, A can no longer bring B down to any substantial extent.
#
# **Single-bin reaction**
#
# Based on experiment `reactions_single_compartment/up_regulate_1`

# %%
LAST_REVISED = "May 5, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename

from life123 import check_version, ChemData, BioSim1D, GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system.  NOTE: Diffusion not applicable (just 1 bin)
chem_data = ChemData(names=["A", "X", "B"], plot_colors=['red', 'darkorange', 'green']) 

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_uniform_concentration(chem_label="A", conc=5.)     # Scarce
bio.set_uniform_concentration(chem_label="X", conc=100.)   # Plentiful
# Initially, no "B" is present

bio.describe_state()

# %%
reactions = bio.get_reactions()

# Reaction A + X <-> 2B , with 1st-order kinetics for all species
reactions.add_reaction(reactants=["A" , "X"], products=[(2, "B", 1)],
                       forward_rate=8., reverse_rate=2.)

reactions.describe_reactions()

# Send the plot of the reaction network to the HTML log file
reactions.plot_reaction_network("vue_cytoscape_2")

# %%

# %% [markdown]
# ### Enable History

# %%
# Let's enable history for all the chemicals
bio.enable_history(take_snapshot=True)     # Taking a snapshot to include the current initial state in the history

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown]
# ## Take the initial system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=30)
bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %% [markdown]
# A, as the scarse limiting reagent, stops the reaction.  
# When A is low, B is also low.

# %%

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 4/1 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%

# %% [markdown]
# # Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title_prefix="Reaction A + X <-> 2B")

# %%

# %% [markdown]
# # Now, let's suddenly increase [A]

# %%
bio.set_bin_conc(bin_address=0, chem_index=0, conc=50.)
bio.describe_state()

# %%
bio.capture_snapshot(caption="[A] suddenly increased externally")

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=40)
bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %% [markdown]
# A, still the limiting reagent, is again stopping the reaction.  
# The (transiently) high value of [A] led to a high value of [B]

# %%
# Verify the equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title_prefix="Reaction A + X <-> 2B")

# %% [markdown]
# `A`, still the limiting reagent, is again stopping the reaction.  
# The (transiently) high value of [A] led to a high value of [B]

# %%

# %% [markdown]
# # Let's again suddenly increase [A]

# %%
bio.set_bin_conc(bin_address=0, chem_index=0, conc=30.)
bio.describe_state()

# %%
bio.capture_snapshot(caption="[A] again suddenly increased externally")

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=70)
bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %% [markdown]
# A, again the scarse limiting reagent, stops the reaction yet again

# %%
# Verify the equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0), 
                                        tolerance=2)

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title_prefix="Reaction A + X <-> 2B")

# %% [markdown]
# `A`, again the scarse limiting reagent, stops the reaction yet again.  
# And, again, the (transiently) high value of [A] up-regulated [B]  
#
# Note: `A` can up-regulate `B`, but it cannot bring it down.  
# `X` will soon need to be replenished, if `A` is to continue being the limiting reagent.

# %% [markdown]
# # For additional exploration, see the experiment `reactions_single_compartment/up_regulate_1`

# %%
