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
# ### Exploring effect of Time Resolution on accuracy
#
# In the examples below, the _time advance_ always remains constant,
# but the _number of steps_ used to arrive there vary

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

from life123 import BioSim1D, ChemData, check_version

# %%
check_version(LIFE123_VERSION)


# %%

# %%
def set_initial_condition():
    # Set or reset the initial concentrations
    bio.set_uniform_concentration(chem_index=0, conc=0.)
    bio.inject_conc_to_bin(bin_address=2, chem_index=0, delta_conc=10.)


# %%
chem_data = ChemData(diffusion_rates=0.1)
bio = BioSim1D(n_bins=10, chem_data=chem_data)

set_initial_condition()

bio.describe_state()

# %%
t_final = 33.3

bio.debug = True

bio.diffuse(total_duration=t_final, n_steps=10)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=20)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=30)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=50)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=100)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=1000)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=10000)

set_initial_condition()     # Reset the concentrations
bio.diffuse(total_duration=t_final, n_steps=100000)
