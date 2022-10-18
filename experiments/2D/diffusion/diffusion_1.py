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
chem_data = chem(names=["A"], diffusion_rates=[0.02])
bio = BioSim2D(n_bins=(5, 8), chem_data=chem_data)

bio.inject_conc_to_bin(bin_address=(2, 1), species_index=0, delta_conc=10.)

bio.describe_state()

# %%
bio.system_snapshot()

# %%
fig = px.imshow(bio.system_snapshot(), 
                title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                labels=dict(x="x (col. number)", y="y (row number)", color="Concentration"),
                text_auto=True, color_continuous_scale="gray_r")         # text_auto=’.2f’

fig.data[0].xgap=2
fig.data[0].ygap=2

fig.show()

# %% [markdown]
# # Initial Diffusion Step

# %%
delta_time = 10.

status = bio.diffuse(total_duration=delta_time, time_step=0.1)
print("\n", status)

bio.describe_state()

# %%
fig = px.imshow(bio.system_snapshot(), 
                title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                labels=dict(x="x (col. number)", y="y (row number)", color="Concentration"),
                text_auto='.2f', color_continuous_scale="gray_r")

fig.data[0].xgap=2
fig.data[0].ygap=2

fig.show()

# %% [markdown]
# ## This is still an early stage in the diffusion process; let's advance it more... (Visualization from results shown at selected times)

# %% tags=[]
for i in range(200):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    if i<2 or i==6 or i>=199:
        bio.describe_state()
        fig = px.imshow(bio.system_snapshot(), 
                title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                labels=dict(x="x (col. number)", y="y (row number)", color="Concentration"),
                text_auto='.2f', color_continuous_scale="gray_r")

        fig.data[0].xgap=2
        fig.data[0].ygap=2

        fig.show()


# %% [markdown]
# ## All bins now have essentially uniform concentration
#
# Notice the continued symmetry across the mid-row.
#
# **Mass conservations**: the "10. units of concentration" are now uniformly spread across the 40 bins, leading to a near-constant concentration of 10./40

# %%
10./40

# %%
