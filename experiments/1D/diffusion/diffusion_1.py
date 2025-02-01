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
# # Diffusion of 1 chemical in 1D
#
# ## An initial concentration pulse (near the left edge of the system) diffusing out towards equilibrium
#
# The system starts out with a "concentration pulse" in bin 2 (the 3rd bin from the left) - i.e. that bin is initially the only one with a non-zero concentration of the only chemical species.
# Then the system is left undisturbed, and followed to equilibrium.

# %% [markdown]
# ### TAGS :  "diffusion 1D", "quick-start"

# %%
LAST_REVISED = "Jan. 31, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

from life123 import BioSim1D, ChemData, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Prepare the initial system
# with a single non-zero bin concentration of the single chemical `A`, near the left edge of the system

# %%
chem_data = ChemData(names="A", diffusion_rates=0.1)     # If you want to assign a default color, pass arg:  plot_colors=["SOME_COLOR_NAME"]

bio = BioSim1D(n_bins=10, chem_data=chem_data)

# %%
bio.inject_conc_to_bin(bin_address=2, chem_label="A", delta_conc=10.)

bio.describe_state()

# %%
bio.system_snapshot()

# %%
bio.visualize_system(title_prefix="Diffusion")   # Line curve view

# %%
bio.system_heatmap(title_prefix="Diffusion")

# %%



# %% [markdown]
# ## Request history-keeping for some bins

# %%
# Request to save the concentration history at the bin with the initial concentration injection, 
# and the bins at the ends of the system
bio.enable_history(bins=[0, 2, 9], frequency=4, take_snapshot=True)    

# %%

# %%

# %% [markdown]
# # Initial Diffusion Step

# %%
#Advancing to time t=10, with time steps of 0.1 ...
delta_time = 10.

status = bio.diffuse(total_duration=delta_time, time_step=0.1)
print(status)

bio.describe_state(concise=True)

# %%
bio.visualize_system(title_prefix="Diffusion")   # Line curve view

# %%
bio.system_heatmap(title_prefix="Diffusion")

# %%

# %% [markdown]
# ## This is still an early stage in the diffusion process; let's advance it more...  
# (Visualization from results shown at selected times)

# %%
for i in range(50):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    if i<2 or i==6 or i>=49:
        bio.describe_state(concise=True)
    
        # Line curve view
        fig = bio.visualize_system(title_prefix="Diffusion")   # Line curve view
        fig.show()
        
        # Heatmap view
        fig = bio.system_heatmap()
        fig.show()


# %% [markdown]
# ## All bins now have essentially uniform concentration
#
# **Mass conservations**: The initial "10 units of concentration" are now uniformly spread across the 10 bins, leading to a near-constant concentration of 10/10 = **1.0**

# %%
# Mass conservation can also be verified as follows:
bio.check_mass_conservation(chem_label="A", expected=10.)

# %%

# %%

# %% [markdown]
# ## Visualization of time changes at particular bins

# %% [markdown]
# #### Instead of visualizing the entire system at a moment of time, like in the previous heatmaps, let's now look at the time evolution of the (only) chemical `A` at either of the bins whose history we requested prior to running the simulation

# %%
bio.conc_history.bin_history(bin_address=2)   # The bin where the initial concentration was applied

# %%
bio.plot_history_single_bin(bin_address=2)

# %%
bio.plot_history_single_bin(bin_address=0)   # Left "edge" of the 1D system

# %%
bio.plot_history_single_bin(bin_address=9)   # Right "edge" of the 1D system

# %%
