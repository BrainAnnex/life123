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
# # Reaction-Diffusion in 2-D :  `A + B <-> C`, taken to equilibrium
# #### Mostly forward reaction, with 1st-order kinetics for each species
#
# Initial concentrations of `A` and `B` are spatially separated to the opposite ends of the system;
# as a result, no `C` is being generated.
#
# But, as soon as `A` and `B`, from their respective distant originating points, 
# diffuse into the middle - and into each other - the reaction starts,
# consuming most of `A` and `B`,
# until an equilibrium is reached in both diffusion and reactions.
#
# Note: This is a 2D version of the 1D experiment by the same name.

# %% [markdown]
# ### TAGS :  "reactions 2D", "diffusion 2D", "quick-start"

# %%
LAST_REVISED = "Jan. 20, 2025"
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
chem_data = ChemData(names=["A", "B", "C"], diffusion_rates=[50., 50., 1.],   # `A` and `B` diffuse fast; `C` diffuses slowly
                     plot_colors=["red", "blue", "purple"])                   # Color choice is a reminder that red + blue = purple

# %%
bio = BioSim2D(x_bins=7, y_bins=7, chem_data=chem_data)

# %%
reactions = bio.get_reactions()

# Reaction A + B <-> C , with 1st-order kinetics for each species; note that it's mostly in the forward direction
reactions.add_reaction(reactants=["A", "B"], products="C", forward_rate=20., reverse_rate=2.)
reactions.describe_reactions()

# %%

# %% [markdown]
# # TIME 0 : Inject initial concentrations of `A` and `B` at opposite ends of the system

# %%
bio.set_bin_conc(bin_address = (0,0), chem_label="A", conc=20.)
bio.set_bin_conc(bin_address = (6,6), chem_label="B", conc=20.)

bio.describe_state()   # A minimalist view of all the chemical concentrations

# %%
bio.system_snapshot(chem_label="A")

# %%
bio.system_snapshot(chem_label="B")

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ## Enable History

# %%
# Let's take a peek at the current concentrations of all chemicals in the bins with the initial concentration injections, as well as at the bin in the very center
bio.selected_concentrations(bins=[(0,0), (6,6), (3,3)])

# %%
# Let's enable history for those same 3 bins
bio.enable_history(bins=[(0,0), (6,6), (3,3)], frequency=2, take_snapshot=True)     # Taking a snapshot to include the current initial state in the history

# %%

# %%

# %% [markdown] tags=[]
# ### Part 1 : advance to time t=0.002 (with smaller fixed steps)

# %%

# %%
bio.react_diffuse(total_duration=0.002, n_steps=10)
bio.describe_state()

# %%
bio.system_heatmaps()

# %%
# Let's take a peek at the history saved so far, for the bins we requested history-keeping; we'll plot it at the end
bio.conc_history.bin_history(bin_address = (0,0))

# %%
bio.conc_history.bin_history(bin_address = (6,6))  # Notice the symmetry between `B` at bin (6,6) and `A` at bin (0,0)

# %%
# And in the central bin
bio.conc_history.bin_history(bin_address = (3,3))   # Notice how `A` and `B` are just beginning to diffuse into this bin, and the reaction `A + B <-> C` is just beginning to produce `C`

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

# %% [markdown]
# ### Notice how `C` begins to get produced when `A` and `B` diffuse into each other, starting at the center bin (3,3)

# %%

# %%

# %% [markdown]
# ### Part 3 : advance the reaction/diffusion to equilibrium

# %%
# Continue with more, and larger, steps
for _ in range(4):
    bio.react_diffuse(total_duration=0.3, n_steps=300)
    fig = bio.system_heatmaps()
    fig.show()

# %%

# %% [markdown]
# ## Visualization of concentration changes with time at particular bins

# %%
bio.plot_history_single_bin(bin_address = (0,0))    # The bin with the initial injection of 'A'

# %% [markdown]
# To separate the curves for `B` and `C`, one needs to substantially magnify the y-axis, because their concentrations at that bin are very similar throughout

# %%

# %%
bio.plot_history_single_bin(bin_address = (6,6))    # The bin with the initial injection of 'B'

# %% [markdown]
# To separate the curves for `A` and `C`, one needs to substantially magnify the y-axis, because their concentrations at that bin are very similar throughout

# %%

# %%
bio.plot_history_single_bin(bin_address = (3,3))      # The midpoint bin

# %% [markdown]
# `A` and `B` are completely superposed, from the perfect symmetry

# %%
