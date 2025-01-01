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
# ## Initial indentical concentration pulses of 2 chemicals (starting at symmetric, almost-opposite ends of the system), diffusing towards equilibrium with identical rates.
#
# Symmetry and mass conservation is observed throughout.
#
# No reaction takes place; the system is left undisturbed, and followed to equilibrium.

# %% [markdown]
# ### TAGS :  "diffusion 2D"

# %%
LAST_REVISED = "Dec. 31, 2024"
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
# Prepare the initial system, with a single non-zero bin, near the left edge of the system, positioned halfway vertically
chem_data = ChemData(names=["A", "B"], diffusion_rates=[0.05, 0.05])
bio = BioSim2D(n_bins=(6, 10), chem_data=chem_data)

bio.set_bin_conc(bin_x = 1, bin_y = 1, chem_label="A", conc=10.)
bio.set_bin_conc(bin_x = 1, bin_y = 8, chem_label="B", conc=10.)

bio.describe_state()

# %%
bio.system_snapshot(chem_label="A")

# %%
bio.system_snapshot(chem_label="B")

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# # Initial Diffusion Step

# %%
delta_time = 2.5

# %%
status = bio.diffuse(total_duration=delta_time, time_step=0.1)
print("\n", status)

bio.system_heatmaps(title_prefix="Diffusion")

# %%
# MASS-CONSERVATION CHECK. Verify that that sum of all the entries in each of the above matrices is still the initial 10.
bio.check_mass_conservation(chem_label="A", expected=10.) \
and \
bio.check_mass_conservation(chem_label="B", expected=10.)

# %% [markdown]
# # A second step

# %%
status = bio.diffuse(total_duration=delta_time, time_step=0.1)
print("\n", status)

bio.system_heatmaps(title_prefix="Diffusion")

# %%

# %% [markdown]
# ## This is still an early stage in the diffusion process; let's advance it more... (Visualization from results shown at selected times)

# %% tags=[]
for i in range(500):
    status = bio.diffuse(total_duration=delta_time, time_step=0.2)

    if i==5 or i==50 or i==100 or i==300 or i>=499:
        fig = bio.system_heatmaps(title_prefix="Diffusion")
        fig.show()


# %% [markdown]
# # All bins now have essentially uniform concentration. The diffusion has reached equilibrium
#
# Notice, throughout the simulation, the continued symmetry of `A` and `B` across the vertical axis

# %% [markdown]
# **Mass conservations**: the initial "10. units of concentration" are now uniformly spread across the 60 (10x6) bins, leading to a near-constant concentration of 10./60

# %%
10./60

# %%
# Mass conservation for both `A` and `B` can also be verified as follows:
bio.check_mass_conservation(chem_label="A", expected=10.)

# %%
bio.check_mass_conservation(chem_label="B", expected=10.)

# %%
