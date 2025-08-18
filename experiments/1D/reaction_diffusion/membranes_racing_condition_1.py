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
# ## TBA
#
# **Recommended background:**  experiment `1D/reaction_diffusion/transient_getting_mopped_up`

# %% [markdown]
# ### TAGS : "reactions 1D", "diffusion 1D", "membranes 1D"

# %%
LAST_REVISED = "Aug. 15, 2025"
LIFE123_VERSION = "1.0.0rc5"       # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim1D, ChemData, Reactions, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ### Initialize the Chemical Data and the Reactions

# %%
# Initialize the chemical data
chem_data = ChemData(names=["A", "B", "C"], 
                     diffusion_rates=[100., 600., 120.],  
                     plot_colors=["red", "turquoise", "green"]) 

rxns = Reactions(chem_data=chem_data)

# Reaction A + B <-> C , with 1st-order kinetics for each species; note that it's mostly in the forward direction
# The reaction is mostly in the forward direction
rxns.add_reaction(reactants=["A", "B"], products="C", forward_rate=0.1, reverse_rate=0.02)
rxns.describe_reactions()

# %%

# %%

# %% [markdown]
# ### Initialize the 1D System, including Membranes

# %%
bio = BioSim1D(n_bins=30, chem_data=chem_data, reactions=rxns)

# %%
bio.membranes().set_membranes(membranes=[ (2, 18) ])

# %%
bio.membranes().membrane_list

# %%
# We'll use 1/2 of the diffusion rate of `A` and `B` 
# as their respective membrane permeability (by passive transport)
# `C`, by constrast, keep the default 0 permeability (i.e., can't cross membranes)
bio.membranes().change_permeability("A", 50.)
bio.membranes().change_permeability("B", 400.)

# %%

# %%

# %% [markdown]
# ### Initialize the initial concentrations

# %%
# Set up the initial bell-shape concentration of `A`, with the very narrow peak close to one end of the system,
# centered at 1/10 of the width of the system, i.e. at bin 30
bio.inject_bell_curve(chem_label="A", center=0.2, sd=0.05, max_amplitude=100., bias=0., clip=(2,17))

# %%
# Chemical `B`, by contrast, is uniformly distributed
bio.set_uniform_concentration(chem_label="B", conc=80.)

# %%

# %% [markdown]
# ### Request history-keeping for some bins

# %%
# Request to save the concentration history at the bins with the initial concentration injection, 
# and the bins at, or near, both ends of the system
bio.enable_history(bins=[0, 10, 20, 29], frequency=15, take_snapshot=True)    

# %%

# %%
df = bio.describe_state()
df

# %%
df[df.columns[2:22]]  # Zoom in where the action is

# %%
# Show as heatmap (including the membranes, shown in brown)
bio.system_heatmaps(title_prefix="Initial strong, localized transient of chemical `A` (membranes shown in brown)")

# %%
# Visualize the system state so far
bio.visualize_system(title_prefix="Initial strong, localized transient of chemical `A` (membranes shown in brown)")

# %% [markdown]
# ### The initial transient of `A` is localized within the compartment (organelle) between bins 2 and 18

# %%

# %%

# %% [markdown]
# ## Start the simulation of the reaction-diffusion

# %%
# The first round of reaction-diffusion, over a small time duration
bio.react_diffuse(total_duration=0.025, fraction_max_step=0.5)      # time_step=0.0005
bio.visualize_system(title_prefix=["The localized transient `A` starts turning into `C` by `A + B <-> C`, ",
                                   "before it can diffuse away much.  Notice the production of `C`, which can't cross the membrane"])

# %%
# SAME IN HEATMAP VIEW
bio.system_heatmaps(title_prefix=["The localized transient `A` starts turning into `C` by `A + B <-> C`, ",
                                  "before it can diffuse away much.  Notice the production of `C`, which can't cross the membrane"])

# %%

# %% [markdown]
# ### Let's continue the reaction-diffusion

# %%
bio.react_diffuse(total_duration=0.025, time_step=0.0005)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.025, time_step=0.0005)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.025, time_step=0.0005)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.05, time_step=0.0005)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.15, time_step=0.0005)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=2.7, time_step=0.0005)
bio.visualize_system()

# %%
bio.system_heatmaps()

# %%

# %%
bio.plot_history_single_bin(title_prefix="Reaction-Diffusion. Bin where the transient originates.", 
                             bin_address=10)

# %%
bio.plot_history_single_bin(title_prefix="Reaction-Diffusion. At left-most edge of system.", bin_address=0)

# %% [markdown]
# #### Only a small amount of `A` reaches bin 0 early on, and later converts to `C`

# %%
bio.plot_history_single_bin(title_prefix="Reaction-Diffusion.", bin_address=20)

# %%
bio.plot_history_single_bin(title_prefix="Reaction-Diffusion.  Faraway bin.", 
                             bin_address=29)

# %% [markdown]
# ## Virtually no `A` ever reaches the faraway bins!

# %%

# %%

# %%

# %% [markdown]
# # Now, let's repeat, starting like before but with a much higher diffusion rate of `A`

# %%
# Initial conditioned just like before

bio = BioSim1D(n_bins=30, chem_data=chem_data, reactions=rxns)

bio.membranes().set_membranes(membranes=[ (2, 18) ])

bio.inject_bell_curve(chem_label="A", center=0.2, sd=0.05, max_amplitude=100., bias=0., clip=(2,17))
bio.set_uniform_concentration(chem_label="B", conc=80.)

# %%
# 10 times more than before

chem_data.set_diffusion_rate(chem_label="A", diff_rate = 1000)     
bio.membranes().change_permeability("A", 500.)

# %%
# Request history-keeping for some bins
bio.enable_history(bins=[0, 10, 20, 29], frequency=15, take_snapshot=True)   

# %%
# Show as heatmap (including the membranes, shown in brown)
bio.system_heatmaps(title_prefix="Initial strong, localized transient of chemical `A` (membranes shown in brown)")

# %%
# Visualize the system state so far
bio.visualize_system(title_prefix="Initial strong, localized transient of chemical `A` (membranes shown in brown)")

# %% [markdown]
# ### Repeat the simulation

# %%
# The first round of REACTION-diffusion, over a small time duration
bio.react_diffuse(total_duration=0.025, time_step=0.0003)
bio.visualize_system(title_prefix=["This time, the localized transient `A` diffuses more substantially across the membrane,",
                                   "before being converted into `C`, which can't cross the membrane"])

# %%
bio.react_diffuse(total_duration=0.025, time_step=0.0002)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.025, time_step=0.0002)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.025, time_step=0.0002)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.05, time_step=0.0003)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.15, time_step=0.0003)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=2.7, time_step=0.0003)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=5, time_step=0.0003)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=10, time_step=0.0003)
bio.visualize_system()

# %%
bio.plot_history_single_bin(title_prefix="Reaction-Diffusion.  Faraway bin.", 
                             bin_address=29)

# %%

# %%
