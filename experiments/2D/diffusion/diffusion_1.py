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
# # Diffusion of 1 chemical in 2D
#
# ## An initial concentration pulse of a single chemical (near the left edge of the system, and halfway vertically), diffusing towards equilibrium
#
# The system starts out with a "concentration pulse" in just one bin - i.e. that bin is initially the only one with a non-zero concentration of the only chemical species.
# Then the system is left undisturbed, and followed to equilibrium.
#
# (Note: this is the 2D counterpart of the 1D experiment by the same name)

# %% [markdown]
# ### TAGS :  "diffusion 2D", "quick-start"

# %%
LAST_REVISED = "Jan. 20, 2025"
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

bio = BioSim2D(x_bins=8, y_bins=5, chem_data=chem_data)

# %%
bio.describe_state()

# %%

# %%

# %% [markdown]
# ## Request history-keeping for some bins

# %%
bio.enable_history(bins=[(1,2), (7,4)], frequency=3)    # Request to save the concentration history at those bins 
                                                        # (the one with the initial injection, and one far away in a corner)

# %%

# %% [markdown]
# ## Apply the initial concentration pulse

# %%
bio.set_bin_conc(bin_address = (1,2), chem_label="A", conc=10.)

bio.describe_state()

# %%
bio.system_snapshot(chem_label="A")

# %%
bio.heatmap_single_chem(chem_label="A", title_prefix="Diffusion")

# %%

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

# %%
# MASS-CONSERVATION CHECK. Verify that that sum of all the entries in the above matrix is still the initial 10.
bio.check_mass_conservation(chem_label="A", expected=10.)

# %%

# %% [markdown]
# ## This is still an early stage in the diffusion process; let's advance it more... (Visualization from results shown at selected times)

# %% tags=[]
for i in range(180):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    if i<2 or i==6 or i>=179:
        bio.describe_state()
        fig = bio.heatmap_single_chem(chem_label="A", title_prefix="Diffusion", height=400)
        fig.show()


# %% [markdown]
# # All bins now have essentially uniform concentration. The diffusion has reached equilibrium
#
# Notice, throughout the simulation, the continued symmetry across the mid-row (ybin 2).

# %% [markdown]
# **Mass conservations**: the initial "10. units of concentration" are now uniformly spread across the 40 (5x8) bins, leading to a near-constant concentration of 10./40

# %%
10./40

# %%
# Mass conservation can also be verified as follows:
bio.check_mass_conservation(chem_label="A", expected=10.)

# %%

# %% [markdown]
# ## Visualization of time changes at particular bins

# %% [markdown]
# #### Instead of visualizing the entire system at a moment of time, like in the previous heatmaps, let's now look at the time evolution of the (only) chemical `A` at either of the bins whose history we requested prior to running the simulation

# %%
bio.conc_history.bin_history(bin_address=(1,2))   # The bin where the initial concentration was applied

# %%
bio.plot_history_single_bin(bin_address=(1,2))

# %%
bio.plot_history_single_bin(bin_address=(7,4))   # A bin in a far-away corner from the initial concentration injection

# %%
