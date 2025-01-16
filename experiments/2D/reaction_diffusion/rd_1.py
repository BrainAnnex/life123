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
# ### Reaction  `A + B <-> C` in 2-D, mostly forward and with 1st-order kinetics for each species, taken to equilibrium
#
# Initial concentrations of `A` and `B` are spacially separated to the opposite ends of the system;
# as a result, no `C` is being generated.
#
# But, as soon as `A` and `B`, from their respective distant originating points at the oppposite corners, 
# diffuse into the middle - and into each other - the reaction starts,
# consuming both `A` and `B` (the forward reaction is much more substantial than the reverse one),
# until an equilibrium is reached in both diffusion and reactions.
#
# Note: This is a 2D version of the 1D experiment by the same name.

# %% [markdown]
# # TODO: 1) respect declared plot colors in the heatmaps

# %% [markdown]
# ### TAGS :  "reactions 2D", "diffusion 2D"

# %%
LAST_REVISED = "Jan. 8, 2025"
LIFE123_VERSION = "1.0.0rc2"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

from life123 import check_version, BioSim2D, ChemData, UniformCompartment

# %%
check_version(LIFE123_VERSION)

# %%

# %%
# Initialize the system
chem_data = ChemData(names=["A", "B", "C"], diffusion_rates=[50., 50., 1.],
                     plot_colors=['red', 'orange', 'green'])

uc = UniformCompartment(chem_data=chem_data)

# Reaction A + B <-> C , with 1st-order kinetics for each species; note that it's mostly in the forward direction
uc.add_reaction(reactants=["A", "B"], products="C", forward_rate=20., reverse_rate=2.)
uc.describe_reactions()

# %%
bio = BioSim2D(x_bins=7, y_bins=7, reaction_handler=uc)

# %%

# %% [markdown]
# # TIME 0 : Inject initial concentrations of `A` and `B` at opposite ends of the system

# %%
bio.set_bin_conc(bin_x = 0, bin_y = 0, chem_label="A", conc=20.)
bio.set_bin_conc(bin_x = 6, bin_y = 6, chem_label="B", conc=20.)

bio.describe_state()

# %%
bio.system_heatmaps()

# %% [markdown]
# ## Enable History

# %%
# Let's take a peek at the current concentrations of all chemicals in the bin with the initial concentration injections, as well as at the bin in the very center
bio.selected_concentrations(bins=[(0,0), (6,6), (3,3)])

# %%

# %%
# Let's enable history for those 3 bins
bio.enable_history(bins=[(0,0), (6,6), (3,3)])

# %%
# Enabling history has the effect of taking a snapshot; let's take a look at it
bio.conc_history.history.get_collection()

# %%

# %%

# %% [markdown] tags=[]
# ### Part 1 : advance to time t=0.002 (with smaller fixed steps)

# %%
bio.react_diffuse(total_duration=0.002, n_steps=10)
bio.describe_state()

# %%
bio.system_heatmaps()

# %%
# Let's take a peek at the history saved so far; we'll plot it at the end
# Each list entry is a triplet of the form (time, data, caption)     Note: this is a low-level data structure, typically not directly interacted with by the end user
bio.conc_history.history.get_collection()

# %%
# A more readable extraction of the historical data in one of the corner bins
bio.conc_history.bin_history(bin_address = (0,0))

# %%
# And in the central bin
bio.conc_history.bin_history(bin_address = (3,3))

# %%

# %%

# %% [markdown]
# ### Part 2 : continue advancing the simulation, with occasional visualization of system snapshots as heatmaps

# %%
# Continue with a few larger steps
for _ in range(5):
    bio.react_diffuse(total_duration=0.03, n_steps=50)
    fig = bio.system_heatmaps()
    fig.show()

# %%

# %%

# %% [markdown]
# ### Part 3 : advance the diffusion/reaction to equilibrium

# %%
# Continue with a more even larger steps steps
for _ in range(5):
    bio.react_diffuse(total_duration=0.3, n_steps=500)
    fig = bio.system_heatmaps()
    fig.show()

# %%

# %%
# The number of historical snapshots we have accumulated
len(bio.conc_history.history)

# %%

# %%
# Let's plot the concentration histories of all chemicals in the central bin
df = bio.conc_history.bin_history(bin_address = (3,3))
df

# %%
bio.plot_history_single_bin(bin_address = (3,3))

# %%
# And in one of the corner bins
df = bio.conc_history.bin_history(bin_address = (0,0))
df

# %%
bio.plot_history_single_bin(bin_address = (0,0))

# %%
