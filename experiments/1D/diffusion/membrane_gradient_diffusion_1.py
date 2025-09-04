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
# ## Diffusion of the initial concentration gradient of a single chemical, in the presence of membranes in its path (passive transport)   
# Three scenarios:  
# 1) no membranes  
# 2) membranes with permeability to that chemical equal its diffusion rate
# 3) membranes with much smaller permeability value

# %% [markdown]
# ### TAGS : "membranes 1D", "diffusion 1D", "basic"

# %%
LAST_REVISED = "Aug. 18, 2025"
LIFE123_VERSION = "1.0.0rc6"       # Library version this experiment is based on

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

# %% [markdown]
# # Part 1
# ## No membranes; a simple diffusion

# %%
# Initialize the system
chem_data = ChemData(names="A", 
                     diffusion_rates=600.)

bio = BioSim1D(n_bins=100, chem_data=chem_data)

# %%
# Adding a gradient slanting to the left
bio.inject_gradient("A", conc_left = 10., conc_right = 60.)
bio.describe_state()

# %%
# Visualize the complete initial state
bio.visualize_system(title_prefix="Initial concentration gradient.  NO membranes present")

# %%

# %%
# The first round of diffusion
bio.diffuse(total_duration=2.5, fraction_max_step=0.5, show_status=True) 
bio.visualize_system()

# %%
# Continue the diffusion
bio.diffuse(total_duration=5, fraction_max_step=0.5, show_status=True) 
bio.visualize_system()

# %% [markdown]
# ## Let's save the complete state of the system, to contrast it with the simulation state at the same time in Part 2

# %%
saved_state_no_membranes = bio.snapshot_system()
saved_state_no_membranes

# %%
# Further continue the diffusion
bio.diffuse(total_duration=12.5, fraction_max_step=0.9, show_status=True)     # we'll take bigger steps now
bio.visualize_system()

# %% [markdown]
# #### Almost completely equilibrated...

# %%

# %%

# %% [markdown]
# # Part 2
# ## With membranes. 
# ### The permeability of `A` across the membranes (passive transport) is EQUAL to its diffusion rate in water

# %%
# Initialize the system
chem_data = ChemData(names="A", 
                     diffusion_rates=600.)

bio = BioSim1D(n_bins=100, chem_data=chem_data)

bio.membranes().set_membranes(membranes=[ (20, 40) ])
bio.membranes().change_permeability("A", 600.)  # *********  We'll use the same value as the diffusion rate in water

# %%
# Adding a gradient slanting to the left
bio.inject_gradient("A", conc_left = 10., conc_right = 60.)
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

# %% [markdown]
# ## Does that turquoise curve (the concentrations of `A` across the length of the system) look familiar??  
# It's **absolutely identical** to what we had in part 1 (with no membranes) - as if the membranes weren't there!  
# We can can verify below that they're indeed identical, from the previously-saved system snapshot at the same time in the simulation:

# %%
bio.snapshot_system()["system_time"] - saved_state_no_membranes["system_time"]   # Indeed, the snapshots were taken at the same times

# %%
bio.snapshot_system()["concentrations"] - saved_state_no_membranes["concentrations"]

# %%

# %%
# Further continue the diffusion
bio.diffuse(total_duration=12.5, fraction_max_step=0.9, show_status=True)    # we'll take bigger steps now
bio.visualize_system()

# %% [markdown]
# #### Almost completely equilibrated...

# %%

# %%

# %% [markdown]
# # Part 3
# ## Identical re-run of Part 2, EXCEPT that the permeability of `A` across the membranes (passive transport) is now MUCH LESS than its diffusion rate in water

# %%
# Initialize the system
chem_data = ChemData(names="A", 
                     diffusion_rates=600.)

bio = BioSim1D(n_bins=100, chem_data=chem_data)

bio.membranes().set_membranes(membranes=[ (20, 40) ])
bio.membranes().change_permeability("A", 30.)  # ******** This time, MUCH smaller than before (1/20-th)

# %%
# Adding a gradient slanting to the left
bio.inject_gradient("A", conc_left = 10., conc_right = 60.)
bio.describe_state()

# %%
# Visualize the complete initial state
bio.visualize_system(title_prefix="Initial concentration gradient.  Membranes shown in brown")

# %%

# %%
# The first round of diffusion
bio.diffuse(total_duration=2.5, fraction_max_step=0.5, show_status=True) 
bio.visualize_system()

# %% [markdown]
# ## Notice the "slowdown" across the membranes, because (passive) transport across them is much slower that diffusion across the system

# %%
# Continue the diffusion
bio.diffuse(total_duration=5, fraction_max_step=0.5, show_status=True) 
bio.visualize_system()

# %%
# Further continue the diffusion
bio.diffuse(total_duration=12.5, fraction_max_step=0.9, show_status=True)   # we'll take bigger steps now
bio.visualize_system()

# %% [markdown]
# ## The system is almost completely equilibrated, but that process continues to happen at a slower rate across the membranes

# %%
