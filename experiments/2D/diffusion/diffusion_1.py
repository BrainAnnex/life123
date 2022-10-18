# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## An initial concentration pulse (near the left edge of the system, and halfway vertically) moving towards equilibrium
#
# The system starts out with a "concentration pulse" in just one bin - i.e. that bin is initially the only one with a non-zero concentration of the only chemical species.
# Then the system is left undisturbed, and followed to equilibrium.
#
# (Note: this is the 2D counterpart of the 1D experiment by the same name)
#
# LAST REVISED: Oct. 17, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %%
from experiments.get_notebook_info import get_notebook_basename

from life_2D.bio_sim_2d import BioSim2D

import plotly.express as px
import plotly.graph_objects as go

from modules.chemicals.chemicals import Chemicals as chem

# %%
# Prepare the initial system, with a single non-zero bin, near the left edge of the system, positioned halfway vertically
chem_data = chem(names=["A"], diffusion_rates=[0.1])
bio = BioSim2D(n_bins=(7, 10), chem_data=chem_data)

bio.inject_conc_to_bin(bin_address=(3, 2), species_index=0, delta_conc=10.)

bio.describe_state()

# %% [markdown]
# # Initial Diffusion Step

# %%
delta_time = 10.

status = bio.diffuse(total_duration=delta_time, time_step=0.1)
print("\n", status)

bio.describe_state()

# %% [markdown]
# ## This is still an early stage in the diffusion process; let's advance it more... (Visualization from results shown at selected times)

# %% tags=[]
for i in range(80):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    if i<2 or i==6 or i>=79:
        bio.describe_state()


# %% [markdown]
# ## All bins now have essentially uniform concentration
#
# The "10 units of concentration" are now uniformly spread across the 70 bins, leading to a near-constant concentration of 10/70

# %%
10/70.

# %%
