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
# # Reaction-Diffusion in 1-D: `A <-> B` in 1-D, taken to equilibrium   
# The reaction product `B` diffuses much faster than the reactant `A`.    
# Eventually, both the reaction and the diffusion come to an equilibrium.

# %% [markdown]
# ### TAGS :  "reactions 1D", "diffusion 1D"

# %%
LAST_REVISED = "May 19, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

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

# %%
# Initialize the system
chem_data = ChemData(names=["A", "B"], diffusion_rates=[1., 50.],   # `B`, the reaction produc diffuses much faster
                     plot_colors=["turquoise", "green"])                  

# %%

# %% [markdown]
# # Let's consider a 1-bin system for starter

# %%
bio = BioSim1D(n_bins=1, chem_data=chem_data)

# %%
reactions = bio.get_reactions()

# Reaction A <-> B , 1st-order kinetics, mostly (but not hugely) in the REVERSE direction
reactions.add_reaction(reactants="A", products="B", forward_rate=2., reverse_rate=8.)
reactions.describe_reactions()

# %%

# %% [markdown]
# ## Set the initial concentrations in the single bin

# %%
bio.set_bin_conc(bin_address=0, chem_label="A", conc=100.)
bio.set_bin_conc(bin_address=0, chem_label="B", conc=100.)

bio.describe_state()

# %%

# %% [markdown]
# ## Enable History

# %%
# Let's enable history for the one bin we currently have
bio.enable_history(bins=[0], take_snapshot=True)     # Taking a snapshot to include the current initial state in the history

# %%

# %%

# %% [markdown]
# ### Advance reaction to equilibrium

# %%
delta_t = 0.002   # This will be our time "quantum" for this experiment

# %%
bio.react_diffuse(time_step=delta_t, n_steps=500)
bio.describe_state()

# %%
bio.show_system_snapshot()

# %%

# %%

# %% [markdown]
# ## Equilibrium 

# %%
bio.reaction_in_equilibrium(bin_address=0, rxn_index=0, explain=True)  # Choice of bin is immaterial now, because they have all equilibrated

# %%

# %% [markdown]
# # Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0)

# %%

# %%

# %% [markdown]
# # Repeat in multi-bin system

# %%
bio = BioSim1D(n_bins=11, chem_data=chem_data)

# %%
reactions = bio.get_reactions()

# Reaction A <-> B , 1st-order kinetics, mostly (but not hugely) in the REVERSE direction
reactions.add_reaction(reactants="A", products="B", forward_rate=2., reverse_rate=8.)
reactions.describe_reactions()

# %% [markdown]
# # TIME 0 : Inject initial concentrations in the middle

# %%
bio.set_bin_conc(bin_address=5, chem_label="A", conc=100.)
bio.set_bin_conc(bin_address=5, chem_label="B", conc=100.)

bio.describe_state()

# %%
# Let's enable history for the middle bin, and one of the end ones
bio.enable_history(bins=[0, 5], frequency=2, take_snapshot=True)     # Taking a snapshot to include the current initial state in the history

# %%

# %% [markdown]
# ### Advance reaction and diffusion to equilibrium

# %%
bio.system_heatmaps()

# %%

# %%
bio.react_diffuse(time_step=0.005, n_steps=1)

# %%
bio.show_system_snapshot()

# %% [markdown]
# #### Notice how much more rapidly `B` is diffusing, relative to `A`

# %%
bio.system_heatmaps()

# %%
bio.chem_quantity(chem_label="A")

# %%
bio.chem_quantity(chem_label="B")

# %% [markdown]
# As expected from the stoichiometry, the total increase in `A` matches total decrease in `B`

# %%

# %%
bio.react_diffuse(time_step=0.005, n_steps=4)
bio.system_heatmaps()

# %%
bio.chem_quantity(chem_label="A")

# %%
bio.chem_quantity(chem_label="B")

# %%

# %%
bio.react_diffuse(time_step=0.005, n_steps=45)
bio.system_heatmaps()

# %%
bio.chem_quantity(chem_label="A")

# %%
bio.chem_quantity(chem_label="B")

# %%

# %%
bio.react_diffuse(time_step=0.005, n_steps=350)
bio.system_heatmaps()

# %%
bio.chem_quantity(chem_label="A")

# %%
bio.chem_quantity(chem_label="B")

# %%

# %%

# %%

# %%

# %%
bio.plot_history_single_bin(bin_address=5)

# %%
bio.plot_history_single_bin(bin_address=0)

# %%
