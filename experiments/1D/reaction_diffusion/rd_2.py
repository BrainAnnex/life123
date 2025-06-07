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
# # "Crossing the border in disguise" : a chemical crosses a membrane impermeable to it,  
# #### by first undergoing a reaction converting it into a different chemical, to which the membrane is permeable,  
# #### and then re-constituting the first chemical outside of the membrane, thru the reverse reaction.
# #### Reaction-diffusion `A <-> B` in 1D, in the presence of membranes with selective permeability to passive transport.
#  
# Eventually, the reaction, the passive membrane transport and the diffusion all come to an equilibrium.

# %% [markdown]
# ### TAGS :  "reactions 1D", "diffusion 1D", "membranes 1D"

# %%
LAST_REVISED = "June 6, 2025"
LIFE123_VERSION = "1.0.0rc5"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import check_version, BioSim1D, ChemData, UniformCompartment

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Prepare the initial system, for starters with IM-PERMEABLE membranes
# with two chemicals `A` and `B`

# %%
chem_data = ChemData(diffusion_rates=[2., 2.], plot_colors=["turquoise", "green"])   # Names "A", "B" automatically assigned

bio = BioSim1D(n_bins=9, chem_data=chem_data)

# %%
reactions = bio.get_reactions()

# Reaction A <-> B , 1st-order kinetics, mostly (but not hugely) in the forward direction
reactions.add_reaction(reactants="A", products="B", forward_rate=8., reverse_rate=2.)
reactions.describe_reactions()

# %%
bio.set_bin_conc(bin_address=4, chem_label="A", conc=10.)   # Set the initial concentration of `A` in middle bin

bio.membranes().set_membranes(membranes=[ (4,5) ])    # By default impermeable

bio.describe_state()

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ## Request history-keeping for some bins

# %%
# Request to save the concentration history at the bin with the initial concentration injection, 
# and at a couple of other bins
bio.enable_history(bins=[0, 2, 4], frequency=10, take_snapshot=True)    

# %%

# %%

# %% [markdown]
# ### Advance reaction to equilibrium
# The membranes are impermeable for now; so, no transport across them occurs

# %%
delta_t = 0.002   # This will be our time "quantum" (fixed time step for reactions, passive transport and diffusion) for this experiment

# %%
bio.react_diffuse(time_step=delta_t, n_steps=350)
bio.describe_state()

# %%
bio.system_heatmaps()

# %%
bio.plot_history_single_bin(bin_address=4)

# %%
# Check the reaction equilibrium
bio.reaction_in_equilibrium(bin_address=4, rxn_index=0, explain=True) 

# %% [markdown]
#

# %%

# %% [markdown]
# # Let's siphon off the product `B` from the central bin 4, by making the membrane permeable to it   
# while still remaining impermeable to `A`.   Note the system time when this happens:

# %%
bio.get_system_time()

# %%
bio.membranes().change_permeability("B", 1.)          # Make the membrane permeable to `B` (and only to `B`!)

# %%

# %% [markdown]
# ### Advance the reaction - and now also the passive transport of `B` across the membrane, and its diffusion outside

# %%
bio.react_diffuse(time_step=delta_t, n_steps=30)
bio.describe_state()

# %%
bio.system_heatmaps()

# %% [markdown]
# Notice what's happening:  
# 1. `B` is diffusing away (because we made the membranes permeable to it)  
# 2. By Le Chatelier's principle, the reaction `A <-> B` in bin 4 is moving forward, because we're siphoning off the product `B`  
# 3. In the other bins, `A` is re-forming from `B`, from the reverse reaction

# %%

# %% [markdown]
# ### Continue advancing the reaction, passive transport of `B`, and diffusion

# %%
bio.react_diffuse(time_step=delta_t, n_steps=50)
bio.describe_state()

# %%
bio.system_heatmaps()

# %% [markdown]
# Notice how the concentration of `A` in the central bin 4 is continuing to drop

# %%
bio.plot_history_single_bin(bin_address=4, vertical_lines_to_add=[0.7], 
                            title_prefix="Dashed vertical line at time when the membranes were made permeable to `B`")

# %%

# %%
bio.react_diffuse(time_step=delta_t, n_steps=150)
bio.describe_state()

# %%
bio.system_heatmaps()

# %% [markdown]
# Notice how [A] in bin 4 is continuing to drop, as `A` keeps turns into `B`

# %%

# %%
bio.react_diffuse(time_step=delta_t, n_steps=500)
bio.describe_state()

# %%
bio.system_heatmaps()

# %%

# %%
bio.react_diffuse(time_step=delta_t, n_steps=1500)
bio.describe_state()

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ### Let's take the system to equilibrium, in reactions and diffusion

# %%
bio.react_diffuse(time_step=delta_t, n_steps=3000)
bio.describe_state()

# %%
bio.system_heatmaps()

# %%
# Check the reaction equilibrium in one of the bins (doesn't matter which one, since all bin concentrations are now the same)
bio.reaction_in_equilibrium(bin_address=0, rxn_index=0, explain=True) 

# %% [markdown]
# ### From a casual glance at the past heatmaps, it might seem that `A` has diffused to equilibrium across all bins...  
# ### But the actual mechanism was that `A` (trapped in bin 4 by a membrane impermeable to it) has transformed into `B`, which in turn has passively crossed the membrane, and then diffused across all bins, while at the same time forming `A` by the reverse reaction in all bins outside the membrane.  
# #### All said and done, `A` has "found a way, thru its disguise as `B`" to indirectly pass thru an impermeable membrane!

# %%
# Let's look at the central bin 4
bio.plot_history_single_bin(bin_address=4, vertical_lines_to_add=[0.7], 
                            title_prefix="Dashed vertical line at time when the membranes were made permeable to `B`;" \
                            "<br>notice the syphoning off of `B` and the resumed advance of the reaction that consumes `A`.")

# %%
# Let's look at a bin CLOSE to the central bin 4
bio.plot_history_single_bin(bin_address=2, vertical_lines_to_add=[0.7], 
                            title_prefix="Dashed vertical line at time when the membranes were made permeable to `B`;" \
                            "<br>notice the QUICK arrival of `B` by diffusion, and the subsequent formation of `A` by reverse reaction.")

# %%
# Let's look at a bin FAR from the central bin 4
bio.plot_history_single_bin(bin_address=0, vertical_lines_to_add=[0.7], 
                            title_prefix="Dashed vertical line at time when the membranes were made permeable to `B`;" \
                            "<br>notice the SLOW arrival of `B` by diffusion, and the subsequent formation of `A` by reverse reaction.")

# %%
