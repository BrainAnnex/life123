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
# ## Diffusion of an initial concentration gradient of a single chemical, in the presence of a membrane in its path  
# Two scenarios: 1) membrane permeability equal to diffusion rate ; 2) much smaller

# %% [markdown]
# ### TAGS : "membranes 1D", "diffusion 1D", "basic"

# %%
LAST_REVISED = "Aug. 16, 2025"
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
chem_data = ChemData(names="A", 
                     diffusion_rates=600.)

bio = BioSim1D(n_bins=100, chem_data=chem_data)

bio.membranes().set_membranes(membranes=[ (20, 40) ])
bio.membranes().change_permeability("A", 600.)  # Unrealistically, for demonstration purposes, 
                                                # we'll use the same value as the diffusion rate in water

# %%
# Adding a gradient slanting to the left
bio.inject_gradient("A", conc_left = 10., conc_right = 60.)

# %%
bio.describe_state()

# %%
# Visualize the complete initial state
bio.visualize_system(title_prefix="Initial concentration gradient.  Membranes shown in brown")

# %%

# %%
# The first round of diffusion
bio.diffuse(total_duration=2.5, fraction_max_step=0.5, show_status=True) 
bio.visualize_system()

# %%
# Continue the diffusion
bio.diffuse(total_duration=5, fraction_max_step=0.5, show_status=True) 
bio.visualize_system()

# %%
bio.diffuse(total_duration=12.5, fraction_max_step=0.5, show_status=True) 
bio.visualize_system()

# %%

# %% [markdown]
#

# %% [markdown]
# # Part 2
# ## Identical re-run, EXCEPT that the permeability of `A` across the membrane (passive transport) is much less than its diffusion rate in water

# %%
# Initialize the system
chem_data = ChemData(names="A", 
                     diffusion_rates=600.)

bio = BioSim1D(n_bins=100, chem_data=chem_data)

bio.membranes().set_membranes(membranes=[ (20, 40) ])
bio.membranes().change_permeability("A", 30.)  # This time, MUCH smaller than before (1/20-th)

# %%
# Adding a gradient slanting to the left
bio.inject_gradient("A", conc_left = 10., conc_right = 60.)

# %%
bio.describe_state()

# %%

# %%
# Visualize the complete initial state
bio.visualize_system(title_prefix="Initial concentration gradient.  Membranes shown in brown")

# %%

# %%
# The first round of diffusion
bio.diffuse(total_duration=2.5, fraction_max_step=0.5, show_status=True) 
bio.visualize_system()

# %%
# Continue the diffusion
bio.diffuse(total_duration=5, fraction_max_step=0.5, show_status=True) 
bio.visualize_system()

# %%
bio.diffuse(total_duration=40, fraction_max_step=0.9, show_status=True) 
bio.visualize_system()

# %%
bio.diffuse(total_duration=40, fraction_max_step=0.9, show_status=True) 
bio.visualize_system()

# %%
