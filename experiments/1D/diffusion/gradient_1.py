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
# # Attaining a Concentration Gradient  
# ### by continuosly injecting and draining, at opposite ends
#
# The system starts out with a uniform concentration.  
# Then identical concentrations are repeatedly *injected to the left* and *drained from the right*  
# Diffusion turn the forced concentration imbalance into a smooth gradient.

# %% [markdown]
# ### TAGS :  "diffusion 1D"

# %%
LAST_REVISED = "Apr. 29, 2025"
LIFE123_VERSION = "1.0.0rc3"       # Library version this experiment is based on

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

# %%
# Initialize the system with a uniform concentration (of the only species)
chem_data = ChemData(names="A", diffusion_rates=0.6)
bio = BioSim1D(n_bins=9, chem_data=chem_data)

bio.set_uniform_concentration(chem_label="A", conc=100.)

bio.describe_state()

# %%
bio.show_system_snapshot()

# %%
# Visualize the system's initial state
bio.visualize_system(title_prefix="Diffusion")

# %%
bio.system_heatmaps(title_prefix="Diffusion", text_format=".0f")

# %%

# %% [markdown]
# # Start the simulation steps

# %%
delta_time = 1.

# %%
for i in range(501):
    # Inject to the leftmost bin
    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=4, zero_clip = False)
    
    # Drain from the rightmost bin
    bio.inject_conc_to_bin(bin_address=8, chem_index=0, delta_conc=-4, zero_clip = False)
    
    # Note: the NET GAIN of moles of A in the system is zero!
    
    # Diffuse for the time span delta_time
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    
    if (i <= 12 and i%3 == 0) or (i%100 == 0):   # Display more frequently initially
        print()
        bio.describe_state(concise=True)
        
        # Show the system state as a line plot
        fig = bio.visualize_system(title_prefix="Diffusion")
        fig.show()
        
        # Show as heatmap
        fig = bio.system_heatmaps(title_prefix="Diffusion", text_format=".3g")
        fig.show()


# %% [markdown]
# ### By now, the gradient has stabilized with  
# ### [A] = 124.67065159 on the left and [A] = 75.32934841 on the right

# %% [markdown]
# Note: if the drain is too large, relative to the diffusion rate, the smaller concentration could "saturate" at zero

# %%
