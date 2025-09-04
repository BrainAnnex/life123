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
# ## A strong, localized transient concentration gets mopped up by a fast chemical reaction before it has a chance to diffuse far
# A sizable transient injection of a chemical `A` near one end of the system, normally diffuses to the opposite end...   
# But if it reacts on a faster timescale, virtually none of it will make it to the far end.  
# It's a *racing condition* between the rate of diffusion of `A`, and the rate at which it's consumed by the plentiful `B` in the reaction `A + B <-> C`  
# When the reaction is fast, relative to `A`'s diffusion, then only `C` (rather than `A`) will make it in substantial amounts to the far end of the system.  
#
# PART 1 explores the diffusion without the concomitant reaction.   
# PART 2 rewinds to the original initial states and carries out the reaction-diffusion.  
#
# **Recommended background:**  experiment `1D/diffusion/localized_transient`

# %% [markdown]
# ### TAGS : "reactions 1D", "diffusion 1D"

# %%
LAST_REVISED = "Aug. 28, 2025"
LIFE123_VERSION = "1.0.0rc6"       # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim1D, ChemData, ReactionRegistry, PlotlyHelper, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Initialize the Chemical Data and the Reactions; they will be used in both Part 1 and Part 2

# %%
# Initialize the chemical data
chem_data = ChemData(names=["A", "B", "C"], diffusion_rates=[100., 80., 120.],  
                     plot_colors=["red", "turquoise", "green"]) 

rxns = ReactionRegistry(chem_data=chem_data)

# Reaction A + B <-> C , with 1st-order kinetics for each species; note that it's mostly in the forward direction
rxns.add_reaction(reactants=["A", "B"], products="C", kF=0.1, kR=0.02)
rxns.describe_reactions()

# %%

# %%

# %% [markdown]
# # PART 1 - JUST DIFFUSION, no reactions

# %%
bio1 = BioSim1D(n_bins=50, chem_data=chem_data, reactions=rxns)

# %%
# Set up the initial bell-shape concentration of `A`, with the very narrow peak close to one end of the system,
# centered at 1/10 of the width of the system, i.e. at bin 30
bio1.inject_bell_curve(chem_label="A", center=0.2, sd=0.06, max_amplitude=100., bias=0.)

# %%
# Chemical `B`, by contrast, is uniformly distributed
bio1.set_uniform_concentration(chem_label="B", conc=100.)

# %%

# %%
df = bio1.describe_state()
df

# %%
df[df.columns[3:22]]  # Zoom in where the action is

# %%
tot_A = bio1.total_chem_mass(chem_label="A")    # We'll use this value for a mass-conservation check, later on
tot_A

# %%

# %% [markdown]
# ## Request history-keeping for some bins

# %%
# Request to save the concentration history at the bins with the initial concentration injection, 
# and the bins at, or near, both ends of the system
bio1.enable_history(bins=[0, 10, 20, 49], frequency=15, take_snapshot=True)    

# %%

# %%
# Visualize the system state so far
bio1.visualize_system(title_prefix="Initial strong, localized transient of chemical `A`")

# %%
# Show as heatmap
bio1.system_heatmaps(title_prefix="Initial strong, localized transient")

# %%
# The first round of diffusion, over a small time duration
bio1.diffuse(total_duration=0.05, time_step=0.002)
diff_only = bio1.visualize_system(title_prefix="Diffusion only", show=True)  # SAVE this plot

# %%
# Do several more rounds of diffusion, each over the same small time duration
for _ in range(3):
    bio1.diffuse(total_duration=0.05, time_step=0.002)
    bio1.visualize_system(title_prefix="Diffusion only.", show=True)

# %%
# Advance the diffusion
bio1.diffuse(total_duration=0.5, time_step=0.002)
bio1.visualize_system(title_prefix="Diffusion only.")

# %%
# Continue advancing the diffusion
bio1.diffuse(total_duration=0.7, time_step=0.002)
bio1.visualize_system(title_prefix="Diffusion only.")

# %%
# Final rounds of diffusion, sampled over larger time spans
for _ in range(2):
    bio1.diffuse(total_duration=1.8, time_step=0.002)
    bio1.visualize_system(title_prefix="Diffusion only.", show=True)

# %%
bio1.check_mass_conservation(expected=tot_A, chem_label="A")    # `A` diffused away, but is still present in the same amount

# %%
bio1.plot_history_single_bin(title_prefix="Diffusion only. Bin where the transient originates.", 
                             bin_address=10)

# %%
bio1.plot_history_single_bin(title_prefix="Diffusion only.", bin_address=0)

# %%
bio1.plot_history_single_bin(title_prefix="Diffusion only.", bin_address=20)

# %%
bio1.plot_history_single_bin(title_prefix="Diffusion only.  Faraway bin.", 
                             bin_address=49)

# %%

# %%

# %%

# %% [markdown]
# # PART 2 - Diffusion AND Reaction

# %%

# %%
# Set up just like before
bio2 = BioSim1D(n_bins=50, chem_data=chem_data, reactions=rxns)
bio2.inject_bell_curve(chem_label="A", center=0.2, sd=0.06, amplitude=15.25, bias=0.)
bio2.set_uniform_concentration(chem_label="B", conc=100.)

# %%
bio2.enable_history(bins=[0, 10, 20, 49], frequency=15, take_snapshot=True)

# %%
bio2.visualize_system(title_prefix="Initial strong, localized transient of chemical `A`")

# %%
# The first round of REACTION-diffusion, over a small time duration
bio2.react_diffuse(total_duration=0.05, time_step=0.002)
react_diff = bio2.visualize_system(title_prefix=["The localized transient `A` starts getting consumed by the reaction `A + B <-> C`, ",
                                                 "before it has a chance to diffuse away much.  Notice the production of `C`"])
react_diff

# %%
# SAME IN HEATMAP VIEW
bio2.system_heatmaps(title_prefix=["HEATMAP VIEW.  The initial localized concentration of `A` starts getting consumed by the reaction `A + B <-> C`, ",
                                    "before it has a chance to diffuse away much.  Notice the production of `C`"])

# %% [markdown]
#

# %% [markdown]
# ## Contrast-and-compare with the earlier diffusion-only

# %%
PlotlyHelper.combine_plots([diff_only, react_diff], 
                           title="Diffusion-only (dotted) vs. React-diffuse (solid)<br>Reaction `A + B <-> C`<br>Early system snapshot", 
                           layout_index=0, modify={0: "dot"})

# %%

# %% [markdown]
# ### Let's continue the reaction-diffusion

# %%
# Do several more rounds of reaction-diffusion, each over the same small time duration
for _ in range(3):
    bio2.react_diffuse(total_duration=0.05, time_step=0.002)
    bio2.visualize_system(title_prefix="Reaction-Diffusion.", show=True)

# %%
# Advance the reaction-diffusion
bio2.react_diffuse(total_duration=0.5, time_step=0.002)
bio2.visualize_system(title_prefix="Reaction-Diffusion.")

# %% [markdown]
# ### Notice how `A` has been largely consumed by now...  and what is diffusing away is `C` (as well as the dip in `B`) rather than `A`

# %%
# Continue advancing the reaction-diffusion
bio2.react_diffuse(total_duration=0.7, time_step=0.002)
bio2.visualize_system(title_prefix="Reaction-Diffusion.")

# %%
# Final rounds of reaction-diffusion, sampled over larger time spans
for _ in range(2):
    bio2.react_diffuse(total_duration=1.8, time_step=0.002)
    bio2.visualize_system(title_prefix="Reaction-Diffusion.", show=True)

# %%

# %%
bio2.plot_history_single_bin(title_prefix="Reaction-Diffusion. Bin where the transient originates.", 
                             bin_address=10)

# %%
bio2.plot_history_single_bin(title_prefix="Reaction-Diffusion.", bin_address=0)

# %%
bio2.plot_history_single_bin(title_prefix="Reaction-Diffusion.", bin_address=20)

# %%
bio2.plot_history_single_bin(title_prefix="Reaction-Diffusion.  Faraway bin.", 
                             bin_address=49)

# %% [markdown]
# ## Virtually no `A` ever reaches the faraway bins!

# %%
