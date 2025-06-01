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
# ### THE DANGER OF EXCESSIVELY LARGE SINGLE TIME STEPS IN DIFFUSION
#
# ->  When the time step = (0.5 / diffusion rate),  
#     a 2-bin system equilibrates in a single step,  
#     and some 3-bin systems can over-shoot equilibrium!
#
# ->  When the time step = (0.33333 / diffusion rate),  
#     some 3-bin systems equilibrate in a single step
#
# So, (0.33333 / diffusion rate) is a - rather lax - upper bound for
# sensible single time steps!  
#
# IMPORTANT: The above is for delta_x = 1; in general, multiply by delta_x**2
#
# The **"Von Neumann stability analysis"**, which provides a slighly looser max value of time steps, is also discussed.
#
# That value of 0.33333 is saved in the Class variable "time_step_threshold"

# %% [markdown]
# ### TAGS :  "diffusion 1D", "under-the-hood"

# %%
LAST_REVISED = "May 27, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

from life123 import ChemData, BioSim1D, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Simulate a 2-bin system, with 1 chemical

# %%
chem_data = ChemData(diffusion_rates=10.)

# %%
bio = BioSim1D(n_bins=2, chem_data=chem_data)
bio.inject_conc_to_bin(bin_address=0, delta_conc=100., chem_index=0)
bio.describe_state()

# %%
bio.time_step_threshold = 0.51    # TO BYPASS THE SAFETY LIMIT 
                                  # (for the allowed maximum delta Time) that is typically in place

# When the time step is (0.5 / diffusion rate),
# a 2-bin system equilibrates in a single step!
bio.diffuse(time_step=0.05, n_steps=1)
bio.describe_state()

# %% [markdown]
# _Note the [50. 50.] values : **the two bins have equilibratedin a single simulation step!**_  
# Mighty suspicious!!

# %%

# %% [markdown]
# ## Start over with a 3-bin system

# %%
bio = BioSim1D(n_bins=3, chem_data=chem_data)
bio.inject_conc_to_bin(bin_address=1, delta_conc=100., chem_index=0)
bio.describe_state()

# %%
bio.time_step_threshold = 0.51    # TO BYPASS THE SAFETY LIMIT 
                                  # (for the allowed maximum delta Time) that is typically in place

# When the time step is (0.5 / diffusion rate),
# a 3-bin system can overshoot equilibrium!
bio.diffuse(time_step=0.05, n_steps=1)
bio.describe_state()

# %% [markdown]
# _Note the concentrations of [50.  0. 50.] : the diffusion has **over-shot equilibrium!!!**_
# #### Cleary, a very excessively large time step!

# %%

# %% [markdown]
# ## Start over again with a new 3-bin system

# %%
bio = BioSim1D(n_bins=3, chem_data=chem_data)
bio.inject_conc_to_bin(bin_address=1, delta_conc=100., chem_index=0)
bio.describe_state()

# %%
bio.time_step_threshold = 0.34    # TO slightly BYPASS THE SAFETY LIMIT 
                                  # (for the allowed maximum delta Time) that is typically in place

bio.diffuse(time_step=0.033333, n_steps=1)
bio.describe_state()

# %% [markdown]
# #### When the time step is (0.33333 / diffusion rate),
# #### a 3-bin system, configured as above, equilibrates in *a single step!*
# #### Again, an excessively large time step!  (note that we're using the default value of 1 for delta_x)

# %% [markdown]
# ## Start over again with a new 3-bin system : same as the last round, but this time use a delta_x = 10 (rather than the default 1)

# %%
bio = BioSim1D(n_bins=3, chem_data=chem_data)
bio.inject_conc_to_bin(bin_address=1, delta_conc=100., chem_index=0)
bio.describe_state()

# %%
bio.time_step_threshold = 0.34    # TO slightly BYPASS THE SAFETY LIMIT 
                                  # (for the allowed maximum delta Time) that is typically in place
    
# When the time step is (delta_x**2 * 0.33333 / diffusion rate),
# a 3-bin system, configured as above, equilibrates in a single step!


# %%
bio.diffuse(time_step=3.3333, n_steps=1, delta_x=10)
bio.describe_state()

# %% [markdown]
# ### This generalizes the previous result to any arbitrary delta_x :
# #### in this scenario, a time step of (delta_x**2 * 0.33333 / diffusion rate)
# #### lead to equilibrium in a single step!

# %%

# %% [markdown]
# # The physical intuition explored in this experiment suggests enforcing:
# #### `delta_t < delta_x**2 * 0.33333 / diffusion rate`
#
# #### Interestingly, this is only slightly stricter than the upper bound that emerges from the **"Von Neumann stability analysis"** of the diffusion equation in 1-D, which states that solutions may become unstable unless
# delta_t < delta_x**2 * 0.5 / diffusion rate
#
# ## NOTE: The method `BioSim1D.diffuse_step()` enforces the stricter upper bound

# %% [markdown]
# #### To keep in mind that the _"Von Neumann stability analysis"_ is ONLY applicable when the diffusion equation is being solved with the **"explicit Forward-Time Centered Space"** method, which is the one currently used by the `diffuse_step()` method.  
# An explanation can be found at: https://www.youtube.com/watch?v=QUiUGNwNNmo

# %%
