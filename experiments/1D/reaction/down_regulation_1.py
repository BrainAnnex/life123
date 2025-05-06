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
# ### `A` down-regulates `B` , by being the *limiting reagent* in reaction `A + 2 B <-> Y` (mostly forward)
# 1st-order kinetics.   
# If [A] is low and [B] is high, then [B] remains high.  
# If [A] goes high, [B] goes low.  
# However, at that point, `A` can no longer bring `B` up to any substantial extent.
#
# *Single-bin* reaction
#
# Based on experiment `reactions_single_compartment/down_regulate_2`

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

from life123 import ChemData, UniformCompartment, BioSim1D, GraphicLog, check_version

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
# Initialize the system.  NOTE: Diffusion not applicable (just 1 bin)
chem_data = ChemData(names=["A", "B", "Y"], plot_colors=['red', 'darkorange', 'green'])

uc = UniformCompartment(chem_data=chem_data)


# Reaction A + 2 B <-> Y , with 1st-order kinetics for all chemical species
uc.add_reaction(reactants=["A", (2, "B", 1)], products="Y",
                forward_rate=8., reverse_rate=2.)

uc.describe_reactions()

# Send the plot of the reaction network to the HTML log file
uc.plot_reaction_network("vue_cytoscape_2")

# %%
bio = BioSim1D(n_bins=1, reaction_handler=uc)

bio.set_uniform_concentration(chem_label="A", conc=5.)     # Scarce
bio.set_uniform_concentration(chem_label="B", conc=100.)   # Plentiful
# Initially, no "Y" is present

bio.describe_state()

# %%

# %%
# Let's enable history - by default for all chemicals and all bins (we're only using one)
bio.enable_history(take_snapshot=True, caption="Initial setup")

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown]
# ### Take the initial system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=30) 
bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %% [markdown]
# A, as the scarse limiting reagent, stops the reaction.  
# When A is low, B is also low.

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0))

# %% [markdown]
# ## Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title="Reaction A + 2 B <-> Y . Concentrations at bin 0")

# %%

# %% [markdown]
# # Now, let's suddenly increase [A]

# %%
bio.set_bin_conc(bin_address=0, chem_index=0, conc=40.)
bio.describe_state()

# %%
bio.capture_snapshot(caption="[A] suddenly increased externally")

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=80) 
bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0), tolerance=7)

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title="Reaction A + 2 B <-> Y . Concentrations at bin 0")

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

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=70)
bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0))

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title="Reaction A + 2 B <-> Y . Concentrations at bin 0")

# %% [markdown]
# `A`, again the scarse limiting reagent, stops the reaction yet again   
#
# Note: `A` can down-regulate `B`, but it cannot bring it up.

# %% [markdown]
# # For additional exploration, see the experiment "reactions_single_compartment/down_regulate_2"

# %%
