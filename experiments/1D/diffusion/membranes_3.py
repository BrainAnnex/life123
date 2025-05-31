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
# # Membranes, with selective permeability, separate 2 previously fully-mixed chemicals
#
# #### The membranes are permeable to one of the chemicals, and impermeable to the other one 
#
# See also experiment `diffusion_3`, where a transient separation is achieved WITHOUT membranes

# %% [markdown]
# ### TAGS :  "membranes 1D", "basic", "diffusion 1D"

# %%
LAST_REVISED = "May 19, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

from life123 import BioSim1D, ChemData, PlotlyHelper, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Prepare the initial system
# with two chemicals `A` and `B`

# %%
chem_data = ChemData(diffusion_rates=[2., 2.], plot_colors=["turquoise", "green"])   # Name "A", "B" automatically assigned

bio = BioSim1D(n_bins=9, chem_data=chem_data)

# %%
bio.set_bin_conc(bin_address=4, chem_label="A", conc=10.)
bio.set_bin_conc(bin_address=4, chem_label="B", conc=10.)

bio.set_membranes(membranes=[ (4,5) ])    # By default impermeable
bio.change_permeability("B", 1.)          # Permeable to `B` (and only to `B`)

bio.describe_state()

# %%
bio.system_snapshot()

# %%
bio.system_heatmaps()

# %%



# %% [markdown]
# ## Request history-keeping for some bins

# %%
# Request to save the concentration history at the bin with the initial concentration injection, 
# and at a couple of other bins
bio.enable_history(bins=[0, 2, 4], frequency=3, take_snapshot=True)    

# %%

# %%

# %% [markdown]
# ## Now, let's start the diffusion

# %%
bio.diffuse(time_step=0.02, n_steps=1)

# %%
bio.system_snapshot()

# %%
bio.system_heatmaps()

# %%
bio.visualize_system()   # Line curve view

# %% [markdown]
# ### The concentrations of `A` and `B` start to differ
# The membranes are impermeable to `A` but permeable to `B`
#

# %%

# %% [markdown]
# ### Let's advance the diffusion

# %%
bio.diffuse(time_step=0.02, n_steps=4)

# %%
bio.system_heatmaps()

# %%
bio.visualize_system()   # Line curve view

# %%

# %%
bio.diffuse(time_step=0.02, n_steps=15)

# %%
bio.system_heatmaps()

# %%
bio.visualize_system()   # Line curve view

# %%

# %%

# %% [markdown]
# ## Finally, let's advance the diffusion to equilibrium

# %%
bio.diffuse(time_step=0.02, n_steps=600)

# %%
bio.system_heatmaps()

# %%
bio.visualize_system()   # Line curve view

# %% [markdown]
# # The membranes, with their selective permeability, have managed to separate the previously fully-mixed `A` and `B`|

# %%

# %% [markdown]
# #### Verify that the total amount of each chemicals hasn't changed; in the case of `B`, it has simply shifted across bins

# %%
# Let's verify mass conservation
bio.check_mass_conservation(chem_label="A", expected=10.)

# %%
bio.check_mass_conservation(chem_label="B", expected=10.)

# %%

# %%

# %% [markdown]
# ## Visualization of time changes at particular bins

# %% [markdown]
# #### Instead of visualizing the entire system at a moment in time, like in the previous diagrams, let's now look at the time evolution of the concentrations of `A` and `B` at selected bins (bins whose history-keeping we requested prior to running the simulation)

# %%
bio.plot_history_single_bin(bin_address=4, title_prefix="`A` and `B` get separated")

# %% [markdown]
# ## The transient separation of `A` and `B` is best seen by plotting the difference of their concentrations, as a function of time

# %%
bin4_hist = bio.conc_history.bin_history(bin_address=4)   # The bin where the initial concentration of `A` and `B` was applied
bin4_hist

# %%
PlotlyHelper.plot_curves(x=bin4_hist['SYSTEM TIME'], y = bin4_hist['A'] - bin4_hist['B'], 
                         x_label="SYSTEM TIME", y_label="[A] - [B]", 
                         colors="purple", title="The separation of `A` and `B`, as seen in central bin 4")

# %%

# %% [markdown]
# ## The separation of `A` and `B` is also seen in other bins...

# %%
bio.plot_history_single_bin(bin_address=0, title_prefix="The separation of `A` and `B` happens at other bins, too...")

# %%
bio.plot_history_single_bin(bin_address=2, title_prefix="The separation of `A` and `B` happens at other bins, too...")

# %%
