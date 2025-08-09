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
# # Reaction-Diffusion in 1-D: `A + B <-> C` in 1-D, taken to equilibrium
# #### Mostly forward reaction, with 1st-order kinetics for each species
#
# Initial concentrations of `A` and `B` are spatially separated to the opposite ends of the system;
# as a result, no `C` is being generated.
#
# But, as soon as `A` and `B`, from their respective distant originating points at the edges, 
# diffuse into the middle - and into each other - the reaction starts,
# consuming both `A` and `B`,
# until an equilibrium is reached in both diffusion and reactions.

# %% [markdown]
# ### TAGS :  "reactions 1D", "diffusion 1D", "quick-start"

# %%
LAST_REVISED = "Aug. 7, 2025"
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

# %%
# Initialize the system
chem_data = ChemData(names=["A", "B", "C"], diffusion_rates=[50., 50., 1.],   # `A` and `B` diffuse fast; `C` diffuses slowly
                     plot_colors=["red", "blue", "purple"])                   # Color choice is a reminder that red + blue = purple

# %%
bio = BioSim1D(n_bins=7, chem_data=chem_data)

# %%
reactions = bio.get_reactions()

# Reaction A + B <-> C , with 1st-order kinetics for each species; note that it's mostly in the forward direction
reactions.add_reaction(reactants=["A", "B"], products="C", forward_rate=20., reverse_rate=2.)
reactions.describe_reactions()

# %%

# %% [markdown]
# # TIME 0 : Inject initial concentrations of `A` and `B` at opposite ends of the system

# %%
bio.set_bin_conc(bin_address=0, chem_label="A", conc=20.)
bio.set_bin_conc(bin_address=6, chem_label="B", conc=20.)

bio.describe_state()

# %%
bio.show_system_snapshot()     # A more streamlined alternate view

# %%
bio.visualize_system(title_prefix="Reaction-Diffusion A + B <-> C")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Reaction-Diffusion A + B <-> C")

# %%

# %% [markdown]
# ## Enable History

# %%
# Let's take a peek at the current concentrations of all chemicals in the bins with the initial concentration injections,
# as well as at the bin in the very center
bio.selected_concentrations(bins=[0, 6, 3])

# %%
# Let's enable history for those same 3 bins
bio.enable_history(bins=[0, 6, 3], frequency=2, take_snapshot=True)     # Taking a snapshot to include the current initial state in the history

# %%

# %%

# %% [markdown]
# ### <a name="sec_2_first_step"></a>First step : advance to time t=0.002

# %%
delta_t = 0.002   # This will be our time "quantum" for this experiment

# %%
# First step
bio.react_diffuse(time_step=delta_t, n_steps=1)
bio.describe_state()

# %%
bio.show_system_snapshot()

# %%
bio.visualize_system(title_prefix="Reaction-Diffusion A + B <-> C")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Reaction-Diffusion A + B <-> C")

# %%

# %% [markdown]
# ### <a name="sec_2"></a>Several more steps : advance to time t=0.016

# %%
# Continue with several delta_t steps
for _ in range(7):
    bio.react_diffuse(time_step=delta_t, n_steps=1)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
bio.visualize_system(title_prefix="Reaction-Diffusion A + B <-> C")   # Line curve view

# %% [markdown]
# `A` is continuing to diffuse from the left.  
# `B` is continuing to diffuse from the right.  
# They're finally beginning to overlap in the middle bins!

# %%
bio.system_heatmaps(title_prefix="Reaction-Diffusion A + B <-> C")

# %%

# %% [markdown]
# ### <a name="sec_3"></a>Several groups of longer runs : advance to time t=0.096

# %%
# Now, do several group of longer runs
for _ in range(4):
    print("\n\n+ 10 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=10)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
bio.visualize_system(title_prefix="Reaction-Diffusion A + B <-> C")   # Line curve view

# %% [markdown]
# `A` is continuing to diffuse from the left.  
# `B` is continuing to diffuse from the right.  
# The overlap of `A` and `B`, especially in the central bins, is by now extensive, and the reaction is proceeding in earnest.  
# Notice the continue symmetry about the center of the system

# %%
bio.system_heatmaps(title_prefix="Reaction-Diffusion A + B <-> C")

# %%

# %% [markdown]
# ### <a name="sec_4"></a>Advance to time t=0.336

# %%
# Continue the simulation
for _ in range(4):
    print("\n\n+++ 30 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=30)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
bio.visualize_system(title_prefix="Reaction-Diffusion A + B <-> C")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Reaction-Diffusion A + B <-> C")

# %%

# %% [markdown]
# ### <a name="sec_5"></a>Advance to time t=0.736

# %%
# Continue the simulation
for _ in range(4):
    print("\n+++++ 50 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=50)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
bio.visualize_system(title_prefix="Reaction-Diffusion A + B <-> C")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Reaction-Diffusion A + B <-> C")

# %%

# %% [markdown]
# ### <a name="sec_6"></a>Advance to time t=1.936

# %%
# Continue the simulation
for _ in range(4):
    print("\n+++++++++++++++ 150 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=150)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
bio.visualize_system(title_prefix="Reaction-Diffusion A + B <-> C")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Reaction-Diffusion A + B <-> C")

# %%

# %% [markdown]
# ### <a name="sec_7"></a>Advance to time t=5.936

# %%
# Continue the simulation
for _ in range(2):
    print("\n++++++++++ ... ++++++++++ 1,000 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=1000)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
bio.visualize_system(title_prefix="Reaction-Diffusion A + B <-> C")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Reaction-Diffusion A + B <-> C")

# %%

# %%

# %% [markdown]
# ## Equilibrium in both diffusion and reaction  
# All bins now have essentially uniform concentration, for each of the chemicals

# %% [markdown]
# In each bin at equilibrium, [A] = 0.49, [B] = 0.49, [C] = 2.37   
# Let's verify that it's consistent with our (mostly-forward) equation `A + B <-> C`:

# %%
bio.reaction_in_equilibrium(bin_address=0, rxn_index=0, explain=True)  # Choice of bin is immaterial now, because they have all equilibrated

# %%

# %% [markdown]
# **Mass conservation**: The initial total quantities of `A` and `B` were 20 each, and zero `C`.   Here's how they have changed:

# %%
20. - bio.chem_quantity(chem_label="A")  # Tot. decrease in `A`

# %%
20. - bio.chem_quantity(chem_label="B")  # Tot. decrease in `B`

# %%
bio.chem_quantity(chem_label="C")  # Tot. increase in `C`

# %% [markdown]
# All consistent with the stoichiometry of our reaction `A + B <-> C`

# %%

# %%

# %% [markdown]
# # Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0)

# %%
bio.plot_history_single_bin(bin_address=6)

# %%
bio.plot_history_single_bin(bin_address=3)

# %% [markdown]
# `A` and `B` overlap on the plot, due to the symmetry of the system.  
# Initially, in the middle bin, neither `A` nor `B` are present; over time they diffuse there... but then they react and get consumed (producing `C`), to an equilibrium value.  
# Meanwhile, `C` gradually diffuses to uniformity.

# %%
