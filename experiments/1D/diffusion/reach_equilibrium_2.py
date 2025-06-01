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
# ## Exploring reaching equilibrium
#
# The system starts out with a pulse in bins near the *left* and the *right* endpoints

# %% [markdown]
# ### TAGS :  "diffusion 1D", "basic"

# %%
LAST_REVISED = "May 3, 2025"
LIFE123_VERSION = "1.0.0rc3"       # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

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
chem_data = ChemData(names="A", diffusion_rates=0.1)
bio = BioSim1D(n_bins=9, chem_data=chem_data)

# Start out with a pulse in bins near the *left* and the *right* endpoints.  
# A total of 20 "units of concentration" is injected
bio.inject_conc_to_bin(chem_label="A", bin_address=2, delta_conc=10.)
bio.inject_conc_to_bin(chem_label="A", bin_address=6, delta_conc=10.)

bio.describe_state()

# %%
bio.show_system_snapshot()

# %%
# Visualize the system's initial state
bio.visualize_system(title_prefix="Diffusion")

# %%
# Show as heatmap
bio.system_heatmaps(title_prefix="Diffusion")

# %%

# %% [markdown]
# # Start the simulation steps

# %%
delta_time = 3.

# %%
for i in range(15):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    print(f"\nAfter Delta time {delta_time}.  TOTAL TIME {bio.system_time}  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)

    bio.visualize_system(title_prefix="Diffusion", show=True)
    
    bio.system_heatmaps(title_prefix="Diffusion", show=True)


# %% [markdown]
# **All cells now have essentially uniform concentration**
#
# The "20 units of concentration" are now uniformly spread across the 9 bins, leading to a near-constant concentration of 20/9 = **2.22**

# %%
bio.check_mass_conservation(expected=20., chem_label="A")

# %%
