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
# ## A strong, localized transient becomes irrelevant with distance
#
# Diffusion of a narrow, bell-shaped, initial concentration of a single chemical.   
# With increasing distance from the location of the transient signal, hardly any change with time is detected, as the system goes to equilibrium

# %% [markdown]
# ### TAGS :  "diffusion 1D"

# %%
LAST_REVISED = "June 8, 2025"
LIFE123_VERSION = "1.0.0rc5"       # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim1D, ChemData, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %%
# Initialize the system
chem_data = ChemData(diffusion_rates=10.)    # Name "A" automatically assigned to only chemical

bio = BioSim1D(n_bins=500, chem_data=chem_data)

# %%
# Set up the initial bell-shape concentration, with the very narrow peak close to one end of the system
bio.inject_bell_curve(chem_label="A", mean=0.1, sd=0.005, amplitude=0.1, bias=10.)

# %%
bio.describe_state()

# %%
# Visualize the system state so far
bio.visualize_system(title_prefix="Initial strong, localized transient")

# %%
# Show as heatmap
bio.system_heatmaps(title_prefix="Initial strong, localized transient")

# %%

# %% [markdown]
# ## Request history-keeping for some bins

# %%
# Request to save the concentration history at the bins with the initial concentration injection, 
# and the bins at the ends of the system
bio.enable_history(bins=[0, 25, 48, 50, 100, 200, 499], frequency=25, take_snapshot=True)    

# %%

# %%

# %%
bio.diffuse(total_duration=5, time_step=0.01)

# %%
bio.visualize_system()

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ### Scrutinize the changes with time, at bins increasingly further away from the transient peak of bin 50

# %%
bio.plot_history_single_bin(bin_address=50)

# %%
bio.plot_history_single_bin(bin_address=48)

# %%
bio.plot_history_single_bin(bin_address=25)

# %%

# %% [markdown]
# ### Continue the diffusion, to equilibrium

# %%
# Do several rounds of diffusion, over relatively small times
for _ in range(5):
    bio.diffuse(total_duration=25, time_step=0.02)
    bio.visualize_system(show=True)

# %%

# %%
# Do more rounds of diffusion, over larger times
for _ in range(5):
    bio.diffuse(total_duration=400, time_step=0.025)
    bio.visualize_system(show=True)

# %%
bio.system_heatmaps()

# %% [markdown]
# #### We're now close to equilibrium

# %% [markdown]
# ### Scrutinize the changes with time, at bins increasingly further away from the transient peak of bin 50

# %%
bio.plot_history_single_bin(bin_address=50)

# %%
bio.plot_history_single_bin(bin_address=48)

# %%
bio.plot_history_single_bin(bin_address=25)

# %%
bio.plot_history_single_bin(bin_address=0)

# %%
bio.plot_history_single_bin(bin_address=100)

# %%
bio.plot_history_single_bin(bin_address=200)

# %%
bio.plot_history_single_bin(bin_address=499)  # This is at the far end of the system

# %% [markdown]
# # Notice how faraway locations barely register that the distant transient ever happened

# %%

# %%
bio.conc_history.bin_history(bin_address=0)

# %%
