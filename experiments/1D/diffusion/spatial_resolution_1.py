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
# # Exploring the change of delta_x (spatial resolution) in diffusion accuracy.
# #### From the same initial setup, diffusion is carried out over a fixed time span,
# #### at different spatial resolutions - and then the respective results are compared

# %% [markdown]
# ### TAGS :  "diffusion 1D", "under-the-hood"

# %%
LAST_REVISED = "May 3, 2025"
LIFE123_VERSION = "1.0.0rc3"       # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim1D, ChemData, Numerical, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Prepare the initial system

# %%
chem_data = ChemData(names="A", diffusion_rates=0.1)

conc_list=[10,13,17,21,25,28,30,38,42,55,65,47,35,32,27,23,20,17,14,8,3,10,16,18,
           20,25,30,35,40,65,85,115,150,92,73,69,65,50,42,36,20,45,50,55,69,82,95,
           77,60,43,37,31,25,22,20,18,15,11,9, 8]

bio = BioSim1D(n_bins=len(conc_list), chem_data=chem_data)

bio.set_species_conc(chem_label="A", conc_list=conc_list)

bio.describe_state()

# %%
# Visualize the initial system state
bio.visualize_system(title_prefix="Diffusion")

# %%
# Show as heatmap
bio.system_heatmaps(title_prefix="Diffusion")

# %%
bio.describe_state(concise=True)

# %%

# %% [markdown]
# ### Populate the data set with more bins, using interpolated concentration values
# ### IMPORTANT: we're **NOT** changing spacial resolution here; we're just creating a less ragged dataset, as *our initial system state*

# %%
bio.smooth_spatial_resolution()
bio.describe_state()

# %%
bio.n_bins

# %%

# %% [markdown]
# # The STARTING POINT
# ### This system setup will be our starting point in exploring diffusion using different spacial resolutions

# %%
original_state = bio.save_system()    # SAVE a copy of the system state, to do multiple runs starting from it

# %%
# Line plot
bio.visualize_system(title_prefix="Diffusion")

# %%
# Show as heatmap
bio.system_heatmaps(title_prefix="Diffusion")

# %%

# %% [markdown]
# # Initial Diffusions with delta_x = 1

# %%
bio.describe_state(concise=True)   # Our initial state

# %%
bio.diffuse(total_duration=7, time_step=0.0005)
bio.describe_state(concise=True)

# %%
# SAVE the above system data (a matrix of dimension n_species x n_bins):  this is the result of diffusion with delta_x = 1
diffuse_dx_1 = bio.system

# %%
# Line plot
bio.visualize_system(title_prefix="Diffusion")

# %% [markdown]
# ### Enough time has proceeded to result in some smoothing, and non-puny changes in most values - but still nowhere near equilibrium

# %% [markdown]
# # Now restore the system to its initial (pre-diffusion) state
# ### and then perform a diffusion over the same time span, but with DOUBLE the spacial resolution
# #### delta_x will be be 1/2 instead of the original default 1

# %%
bio.restore_system(original_state)

# %%
bio.describe_state()

# %%
# Double the spacial resolution
bio.increase_spatial_resolution(2)
bio.describe_state()

# %%
# Now repeat the idential diffusion process as before, but with half the delta_x
bio.diffuse(total_duration=7, time_step=0.0005, delta_x=0.5)

# %%
# Finally, halve the resolution, to return to the original number of bins
bio.decrease_spatial_resolution(2)

# %%
bio.describe_state(concise=True)

# %%
# SAVE the above system data: this is the result of diffusion with delta_x of 1/2
diffuse_dx_1_2 = bio.system

# %% [markdown]
# ### Compare the last 2 runs (with dx=1 and dx=1/2)

# %%
Numerical .compare_states(diffuse_dx_1 , diffuse_dx_1_2, verbose=True)

# %% [markdown]
# # Again, restore the system to its initial (pre-diffusion) state
# ### and then perform a diffusion over the same time span, but with QUADRUPLE the spacial resolution
# ### delta_x will be be 1/4 instead of the original default 1

# %%
bio.restore_system(original_state)

# %%
bio.describe_state()

# %%
# Quadruple the spacial resolution
bio.increase_spatial_resolution(4)
bio.describe_state()

# %%
# Now repeat the idential diffusion process as before, but with 1/4 the delta_x
bio.diffuse(total_duration=7, time_step=0.0005, delta_x=0.25)

# %%
# Finally, reduce the resolution by a factor 4, to return to the original number of bins
bio.decrease_spatial_resolution(4)
bio.describe_state(concise=True)

# %%
# SAVE the above system data: this is the result of diffusion with delta_x of 1/4
diffuse_dx_1_4 = bio.system

# %% [markdown]
# ### Compare the latest 2 runs (with dx=1/2 and dx=1/4)

# %%
Numerical .compare_states(diffuse_dx_1_2 , diffuse_dx_1_4)

# %% [markdown]
# ### Notice how the discrepancies have gone down

# %% [markdown]
# # One last time, restore the system to its initial (pre-diffusion) state
# ### and then perform a diffusion over the same time span, but with 10x the spacial resolution
# ### delta_x will be be 1/10 instead of the original default 1

# %%
bio.restore_system(original_state)

# Increase by a factor 10 the spacial resolution
bio.increase_spatial_resolution(10)
bio.n_bins

# %%
# Now repeat the idential diffusion process as before, but with 1/10 the delta_x
bio.diffuse(total_duration=7, time_step=0.0005, delta_x=0.1)

# %%
# Finally, reduce the resolution by a factor 10, to return to the original number of bins
bio.decrease_spatial_resolution(10)
bio.describe_state(concise=True)

# %%
# SAVE the above system data: this is the result of diffusion with delta_x of 1/10
diffuse_dx_1_10 = bio.system

# %% [markdown]
# ### Again, compare the latest 2 runs (with dx=1/4 and dx=1/10)

# %%
Numerical.compare_states(diffuse_dx_1_4 , diffuse_dx_1_10)

# %% [markdown]
# ### Notice how the discrepancies have gone down even more
# ### This matches expectations that we're getting closer and closer to a "true" (very high precision) value, as we keep increasing the spacial resolution

# %%
