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
# # Transient separation of 2 chemicals in 1D, due to their different diffusion rates
#
# ### Starting with identical concentrations, they initially separate 
# ### until they eventually attain identical concentrations at equilibrium

# %% [markdown]
# ### TAGS :  "diffusion 1D", "basic"

# %%
LAST_REVISED = "May 29, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim1D, ChemData, PlotlyHelper, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Prepare the initial system
# with identical bin concentration of the chemicals `A` and `B`, in the center bin.  
# `B` diffuses much faster than `A`

# %%
chem_data = ChemData(names=["A", "B"], diffusion_rates=[0.1, 2.], plot_colors=["turquoise", "green"])

bio = BioSim1D(n_bins=9, chem_data=chem_data)

# %%
bio.set_bin_conc(bin_address=4, chem_label="A", conc=10.)
bio.set_bin_conc(bin_address=4, chem_label="B", conc=10.)

bio.describe_state()

# %%
bio.system_snapshot()

# %%
bio.visualize_system(title_prefix="Diffusion.  The 2 concentrations initially perfectly overlap")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Diffusion")

# %%

# %% [markdown]
# ## Request history-keeping for some bins

# %%
# Request to save the concentration history at the bins with the initial concentration injection, 
# and the bins at the ends of the system
bio.enable_history(bins=[0, 2, 4], frequency=3, take_snapshot=True)    

# %%

# %%

# %% [markdown]
# ## Initial Diffusion Step

# %%
# Advancing to time t=1, with time steps of 0.05
bio.diffuse(total_duration=1.0, time_step=0.05)

# %%
bio.describe_state()

# %% [markdown]
# ## The 2 chemicals `A` and `B` have separated due to their different diffusion rates

# %%
bio.visualize_system(title_prefix="Diffusion")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Diffusion")

# %%

# %%

# %% [markdown]
# ## Let's advance the diffusion

# %%
# Advancing to time t=3, with time steps of 0.05
bio.diffuse(total_duration=2.0, time_step=0.05)

# %%
bio.describe_state()

# %% [markdown]
# ## `B` is closer to equilibrium, while `A` is nowhere near it!

# %%
bio.visualize_system(title_prefix="Diffusion")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Diffusion")

# %%

# %%

# %% [markdown]
# ## Let's advance the diffusion some more

# %%
bio.diffuse(total_duration=2.0, time_step=0.05)

# %%
bio.describe_state()

# %% [markdown]
# ## Notice that the concentrations of `A` and `B` in the central bin 4 are again getting closer

# %%
bio.visualize_system(title_prefix="Diffusion")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Diffusion")

# %%

# %%

# %% [markdown]
# ## Finally, let's advance the diffusion to equilibrium

# %%
bio.diffuse(total_duration=95.0, time_step=0.05)

# %%
bio.describe_state()

# %% [markdown]
# ## All bins now have essentially uniform concentration  
# ## Notice that the concentrations of `A` and `B` in the central bin 4 are again essentially identical; their separation has terminated

# %%
bio.visualize_system(title_prefix="Diffusion")   # Line curve view

# %%
bio.system_heatmaps(title_prefix="Diffusion")

# %%
# Let's verify mass conservation
bio.check_mass_conservation(chem_label="A", expected=10.)

# %%
bio.check_mass_conservation(chem_label="B", expected=10.)

# %%

# %%

# %% [markdown]
# ## Visualization of time changes at particular bins

# %% [markdown]
# #### Instead of visualizing the entire system at a moment of time, like in the previous heatmaps, let's now look at the time evolution of the  chemicals `A` and `B` at selected bins (whose history we requested prior to running the simulation)

# %%
bin4_hist = bio.conc_history.bin_history(bin_address=4)   # The bin where the initial concentration of `A` and `B` was applied
bin4_hist

# %%
bio.plot_history_single_bin(bin_address=4, title_prefix="`A` and `B` separate, until they finally rejoin")

# %% [markdown]
# ## The transient separation of `A` and `B` is best seen by plotting the difference of their concentrations, as a function of time

# %%
bin4_hist["diff"] = bin4_hist["A"] - bin4_hist["B"]
bin4_hist

# %%
PlotlyHelper.plot_pandas(df=bin4_hist, x_var="SYSTEM TIME", y_label="[A] - [B]", fields="diff",
                         colors="purple", legend_header="Chemical", 
                         title="Transient separation of `A` and `B`, as seen in central bin 4")

# %%

# %%

# %% [markdown]
# ## The transient separation of `A` and `B` is also seen in other bins...

# %%
bio.plot_history_single_bin(bin_address=0, title_prefix="The transient separation of `A` and `B` happens at other bins, too")

# %%
bio.plot_history_single_bin(bin_address=2, title_prefix="The transient separation of `A` and `B` happens at other bins, too")

# %%
