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
# ## High Frequencies in the concentration get smoothed out first by diffusion   
# ### Ultimately, the only frequency passing thru is zero! (i.e., uniform concentration at equilibrium)
# #### We explore how an initial concentration with 3 different sinusoidal frequencies fares in the course of a diffusion to equilibrium
#
# **The initial system state will consist of:**  
#     0 - A constant baseline of value 30  
#     1 - A sine wave of frequency 1 (1 cycle across the system's length), of amplitude 10  
#     2 - A sine wave of frequency 10 , of amplitude 4  
#     3 - A sine wave of frequency 40 , of amplitude 2  

# %% [markdown]
# ### TAGS :  "diffusion 1D"

# %%
LAST_REVISED = "May 2, 2025"
LIFE123_VERSION = "1.0.0rc3"       # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim1D, ChemData, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %%
# Initialize the system.  We use a RELATIVELY LARGE NUMBER OF BINS, 
# to captures the many changes in the high-frequency component
chem_data = ChemData(names="A", diffusion_rates=0.5)
bio = BioSim1D(n_bins=500, chem_data=chem_data)

# %%

# %% [markdown]
# ## PART 1 (of 3) of Initial Preparation -
# ### Start with a sinusoidal concentration, with exactly 1 cycle over the length of the system
# #### (notice the bias value of 30; it will be seen at the eventual equilibrium, at the end)

# %%
bio.inject_sine_conc(chem_label="A", amplitude=10, bias=30, number_cycles=1)

# %%
bio.show_system_snapshot()

# %%
# Visualize the system's initial state
bio.visualize_system(title_prefix="Low-frequency component of the Initial System State")

# %%
# Show as heatmap
bio.system_heatmaps(title_prefix="Low-frequency component of the Initial System State")

# %%

# %% [markdown]
# ## PART 2 (of 3) of Initial Preparation -
# ### Now add a higher-frequency component (10 cycles over the length of the system)

# %%
bio.inject_sine_conc(chem_label="A", amplitude=4, bias=0, number_cycles=10)

# %%
bio.visualize_system(title_prefix="Low- and mid-frequency components of the Initial System State", 
                     show=True)

bio.system_heatmaps(title_prefix="Low- and mid-frequency components of the Initial System State")

# %%

# %% [markdown]
# ## PART 3 (of 3) of Initial Preparation -
# ### To complete the preparation of a 3-frequency component initial state, add another, even higher, frequency component 
# #### (40 cycles over the length of the system)

# %%
bio.inject_sine_conc(chem_label="A", amplitude=2, bias=0, number_cycles=40)

# %%
bio.visualize_system(title_prefix="Initial System State with 3 superposed frequencies", 
                     show=True)

bio.system_heatmaps(title_prefix="Initial System State with 3 superposed frequencies")

# %%
# Take a look at the frequency domain of the concentration values
bio.frequency_analysis(chem_label="A")

# %%

# %% [markdown]
# # Start the diffusion steps

# %%
bio.diffuse(total_duration=10, time_step=0.1)

bio.visualize_system(title_prefix="Diffusion")   # Show as a line plot

# %% [markdown]
# ### After the initial diffusion (at t=10), the highest frequency is largely gone

# %%
# Take a look at the frequency domain of the concentration values
bio.frequency_analysis(chem_label="A")

# %%
# A lot of tiny frequency components are now present; take just the largest 4
bio.frequency_analysis(chem_label="A", n_largest=4)

# %%

# %%
# Advance the diffusion
bio.diffuse(total_duration=20, time_step=0.1)

# Show as a line plot
bio.visualize_system(title_prefix="Diffusion")

# %% [markdown]
# ### After additional diffusion (at t=30), the highest frequency (40 cycles) can no longer be visually detected.
# #### The 2 lower frequency are still recognizable
# #### The 40-cycle frequency is not even among the largest of the tiny values in the spurious frequencies of the distorted signal

# %%
bio.frequency_analysis(chem_label="A", n_largest=10)

# %%

# %%
# Advance the diffusion
bio.diffuse(total_duration=90, time_step=0.1)

# Show as a line plot
bio.visualize_system(title_prefix="Diffusion")

# %% [markdown]
# ### By now (at t=120), even the middle frequency (10 cycles) is notably attenuated

# %%
bio.frequency_analysis(chem_label="A", n_largest=3)

# %%

# %%
# Advance the diffusion
bio.diffuse(total_duration=100, time_step=0.1)

# Show as a line plot
bio.visualize_system(title_prefix="Diffusion")

# %% [markdown]
# ### With still more diffusion (at t=220), the middle frequency is only weakly visible

# %%
bio.frequency_analysis(chem_label="A", n_largest=3)

# %%

# %%
# Advance the diffusion
bio.diffuse(total_duration=180, time_step=0.1)

# Show as a line plot
bio.visualize_system(title_prefix="Diffusion")

# %% [markdown]
# ### By t=400, the middle frequency is barely noticeable.  The low frequency (1 cycle) is still prominent

# %%
bio.frequency_analysis(chem_label="A", n_largest=3)

# %%

# %%
# Advance the diffusion
bio.diffuse(total_duration=600, time_step=0.3)

# Show as a line plot
bio.visualize_system(title_prefix="Diffusion")

# %% [markdown]
# ### By t=1,000 there's no visual indication of the middle frequency (f=10)

# %%
bio.frequency_analysis(chem_label="A", n_largest=10)

# %%

# %%
# Advance the diffusion
bio.diffuse(total_duration=8000, time_step=.5)

# Show as a line plot
bio.visualize_system(title_prefix="Diffusion")

# %% [markdown]
# ### By t=9,000 even the lowest frequency has lost a major part of its sinusoidal shape

# %%
bio.frequency_analysis(chem_label="A", n_largest=2)

# %% [markdown]
# Note how the zero-frequency is now gaining over the baseline 1-cycle signal

# %%

# %%
# Advance the diffusion
bio.diffuse(total_duration=91000, time_step=.5)

# Show as a line plot
bio.visualize_system(title_prefix="Diffusion")

# %%
bio.frequency_analysis(chem_label="A", n_largest=2)

# %% [markdown]
# ### By t=100,000 the system is clearly approaching equilibrium

# %%

# %%
# Advance the diffusion
bio.diffuse(total_duration=100000, time_step=.6)

# %%
# Show as a line plot
bio.visualize_system(title_prefix="Diffusion", show=True)

bio.system_heatmaps(title_prefix="Diffusion")

# %% [markdown]
# ### By t=200,000 the system is getting close to equilibrium about the value 30, which was the original "bias" (unvarying component) of the baseline frequency; the higher-frequency signals didn't have any bias

# %%
bio.frequency_analysis(chem_label="A", n_largest=2)

# %%
