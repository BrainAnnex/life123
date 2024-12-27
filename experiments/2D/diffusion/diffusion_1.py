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
# ## An initial concentration pulse of a single chemical (near the left edge of the system, and halfway vertically), diffusing towards equilibrium
#
# The system starts out with a "concentration pulse" in just one bin - i.e. that bin is initially the only one with a non-zero concentration of the only chemical species.
# Then the system is left undisturbed, and followed to equilibrium.
#
# (Note: this is the 2D counterpart of the 1D experiment by the same name)

# %% [markdown]
# ### TAGS :  "diffusion 2D"

# %%
LAST_REVISED = "Dec. 26, 2024"
LIFE123_VERSION = "1.0.0rc2"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim2D, ChemData, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %%
# Prepare the initial system, with a single non-zero bin, near the left edge of the system, positioned halfway vertically
chem_data = ChemData(names="A", diffusion_rates=0.02)
bio = BioSim2D(n_bins=(5, 8), chem_data=chem_data)

bio.set_bin_conc(bin_x = 2, bin_y = 1, chem_label="A", conc=10.)

bio.describe_state()

# %%
bio.system_snapshot(chem_label="A")

# %%
bio.heatmap_single_chem(chem_label="A", title_prefix="Diffusion")

# %%

# %% [markdown]
# # Initial Diffusion Step

# %%
delta_time = 10.

status = bio.diffuse(total_duration=delta_time, time_step=0.1)
print("\n", status)

bio.describe_state()

# %%
bio.heatmap_single_chem(chem_label="A", title_prefix="Diffusion")

# %% [markdown]
# ## This is still an early stage in the diffusion process; let's advance it more... (Visualization from results shown at selected times)

# %% tags=[]
for i in range(200):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    if i<2 or i==6 or i>=199:
        bio.describe_state()
        fig = bio.heatmap_single_chem(chem_label="A", title_prefix="Diffusion", height=400)
        fig.show()


# %% [markdown]
# # All bins now have essentially uniform concentration. The diffusion has reached equilibrium
#
# Notice, throughout the simulation, the continued symmetry across the mid-row (ybin 2).
#
# **Mass conservations**: the initial "10. units of concentration" are now uniformly spread across the 40 (5x8) bins, leading to a near-constant concentration of 10./40

# %%
10./40

# %%
