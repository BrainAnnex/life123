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
# # Diffusion of 2 chemicals
#
# ## Initial indentical concentration pulses of 2 chemicals (starting at symmetric, almost-opposite, ends of the system), diffusing towards equilibrium with identical rates coefficients.
#
# Symmetry and mass conservation is observed throughout.
#
# No reaction takes place; the system is left undisturbed, and followed to equilibrium.

# %% [markdown]
# ### TAGS :  "diffusion 2D", "basic"

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
# Prepare the initial system, with two concentration pulses, at almost-opposite ends of the system, with the same diffusion rate coefficient
chem_data = ChemData(names=["A", "B"], diffusion_rates=[0.05, 0.05])

bio = BioSim2D(x_bins=7, y_bins=7, chem_data=chem_data)

bio.set_bin_conc(bin_address = (1,1), chem_label="A", conc=10.)
bio.set_bin_conc(bin_address = (5,5), chem_label="B", conc=10.)

# %% [markdown]
# ### Different ways to visualize the system   
# Notice the bin with the initial concentration of `A` and the other bin with the initial concentration of `B`

# %%
bio.describe_state()

# %%
bio.system_snapshot(chem_label="A")

# %%
bio.system_snapshot(chem_label="B")

# %%
bio.system_heatmaps()

# %% [markdown]
# Notice that the bins with the initial concentration pulses lie on a diagonal from each other, and are symmetric across the center

# %%

# %%

# %% [markdown]
# ## Request history-keeping for some bins

# %%
# Request to save the concentration history of all chemicals at selected bins 
bio.enable_history(bins=[(1,1), (3,3), (1,6)], frequency=5, take_snapshot=True)    # Taking a snapshot to include the current initial state in the history

# %%

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
for i in range(300):
    status = bio.diffuse(total_duration=delta_time, time_step=0.2)

    if i==3 or i==20 or i==50 or i==100 or i==299:
        fig = bio.system_heatmaps(title_prefix="Diffusion")
        fig.show()


# %% [markdown]
# # All bins now have essentially uniform concentration. The diffusion has reached equilibrium
#
# Notice how `B`, further away from walls compared to `A`, has an opportunity to spread out more quickly, in spite of the same diffusion rate coefficient

# %% [markdown]
# **Mass conservations**: the initial "10. units of concentration" are now uniformly spread across the 49 (7x7) bins, leading to a near-constant concentration of 10./49

# %%
10./49

# %%
# Mass conservation for both `A` and `B` can also be verified as follows:
bio.check_mass_conservation(chem_label="A", expected=10.)

# %%
bio.check_mass_conservation(chem_label="B", expected=10.)

# %%

# %% [markdown]
# ## Visualization of time changes at particular bins

# %%
bio.plot_history_single_bin(bin_address=(1,1))    # The bin with the initial pulse of `A`

# %%

# %%
bio.plot_history_single_bin(bin_address=(3,3))    # The bin halfway between the initial pulses of `A` and `B`

# %% [markdown]
# Bin (3,3) is exactly halfway between the bins with the initial pulses of `A` and `B`, and because of the complete symmetry of the system, the concentrations of `A` and `B` are completely overlapping

# %%

# %%
bio.plot_history_single_bin(bin_address=((1, 6))) # Distance from (1,1), with the initial pulse of `A`, is 5;  
                                                  # distance from (5,5), with the initial pulse of `B`, is 4.12

# %% [markdown]
# `B` arrives more quickly at this bin, until `A` eventually catches up

# %%
