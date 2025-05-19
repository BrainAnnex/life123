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
# # Membranes in 1D : Data Structure
#
# ## Data Structure and Visualization of Membranes in 1D
#
# No simulations done here; for diffusion and transport across membranes, please see `membranes_2` and other experiments

# %% [markdown]
# ### TAGS :  "membranes 1D", "basic", "quick-start"

# %%
LAST_REVISED = "May 18, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

from life123 import BioSim1D, ChemData, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Prepare the initial system
# with a single non-zero bin concentration of the single chemical `A`, near the left edge of the system

# %%
chem_data = ChemData(names=["A", "B"], plot_colors=["turquoise", "green"])

bio = BioSim1D(n_bins=21, chem_data=chem_data)

# %%
bio.inject_gradient(chem_label="A", conc_left = 0., conc_right = 100.)

bio.inject_sine_conc(chem_label="B", number_cycles=2, amplitude=5., bias=10., phase=0)

bio.describe_state()

# %%
bio.system_snapshot()

# %%
bio.visualize_system()   # Line curve view

# %%
bio.system_heatmaps(text_format=".0f")

# %%



# %% [markdown]
# # Add Membranes

# %%
bio.set_membranes(membranes=[ (0, 5) ])

# %%
bio.membranes

# %%
bio.system_heatmaps(text_format=".0f")

# %%
bio.set_membranes(membranes=[ (0, 5), (10, 11), (16,21) ])    # Overwrite previous membranes

# %%
bio.system_heatmaps(text_format=".0f")

# %%
