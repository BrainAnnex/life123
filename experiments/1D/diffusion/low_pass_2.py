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
# ## The attenuation over time of a single-frequency component in the initial concentration
#
# ### The initial system state is a sine wave of frequency 2 (i.e. 2 cycles across the system's length), of amplitude 10, with a baseline (bias) of 30
# ### Afterward, the process is restarted and repeated with a frequency 5 times larger  

# %% [markdown]
# ### TAGS :  "diffusion 1D"

# %%
LAST_REVISED = "June 4, 2025"
LIFE123_VERSION = "1.0.0rc6"        # Library version this experiment is based on

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

# %%
# Initialize the system.  We use a RELATIVELY LARGE NUMBER OF BINS, 
# to captures the finer changes in the frequency components of the concentration function
chem_data = ChemData(names="A", diffusion_rates=0.5)
bio = BioSim1D(n_bins=500, chem_data=chem_data)

# %%

# %% [markdown]
# ## Initial Preparation -
# ### Start with a sinusoidal concentration, with exactly 2 cycles over the length of the system
# #### (of amplitude 10 and bias value of 30)

# %%
bio.inject_sine_conc(chem_label="A", amplitude=10, bias=30, number_cycles=2)

# %%
bio.show_system_snapshot()

# %%
# Visualize the system's initial state
bio.visualize_system(title_prefix="Initial System State")

# %%
# Show as heatmap
bio.system_heatmaps(title_prefix="Initial System State (as a heatmap)")

# %%
# Take a look at the frequency domain of the concentration values
frequency_data = bio.frequency_analysis(chem_label="A")
frequency_data

# %%
# The ratio of the frequency-2 wave over the constant part (the bias)
ratio = frequency_data.loc[1, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
ratio

# %% [markdown]
# #### Start saving the value of this type of ratio

# %%
bio.save_value(data_snapshot={"ratio": ratio})
bio.get_saved_values()

# %%

# %% [markdown]
# # Start the diffusion - to time t=10

# %%
bio.diffuse(total_duration=10, n_steps=100)

# %%
# Take a look at the frequency domain of the concentration values
frequency_data = bio.frequency_analysis(chem_label="A")
frequency_data

# %%
ratio = frequency_data.loc[2, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
ratio

# %%
bio.save_value(data_snapshot={"ratio": ratio})
bio.get_saved_values()

# %%

# %% [markdown]
# ## Do 49 more rounds of diffusion - to time t = 500

# %%
for i in range(49):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(chem_label="A")
    ratio = frequency_data.loc[2, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_value(data_snapshot={"ratio": ratio})

# %%
bio.get_saved_values()

# %%
bio.visualize_system()

# %% [markdown]
# #### Note how the curve is beginning to flatten; the peaks are no longer at 20 and 40

# %%

# %% [markdown]
# ## Do 150 more rounds of diffusion - to time t = 2000

# %%
for i in range(150):
    bio.diffuse(total_duration=10, n_steps=75)    # Notice the gradual decreas of the number of intermediate steps, given the smaller gradient
    frequency_data = bio.frequency_analysis(chem_label="A")
    ratio = frequency_data.loc[2, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_value(data_snapshot={"ratio": ratio})

# %%
bio.get_saved_values()

# %%
bio.visualize_system()

# %% [markdown]
# #### Note how the curve is flatter still - and even beginning to lose shape at the boundary 

# %%

# %% [markdown]
# ## Do 800 more rounds of diffusion - to time t = 10000

# %%
for i in range(800):
    bio.diffuse(total_duration=10, n_steps=20)   # Note how we're gradually increasing the number of intermediate steps, because the gradiants are now smaller
    frequency_data = bio.frequency_analysis(chem_label="A")
    ratio = frequency_data.loc[2, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_value(data_snapshot={"ratio": ratio})

# %%
bio.get_saved_values()

# %%
bio.visualize_system()

# %% [markdown]
# #### Getting decisively flatter, as it approaches equilibrium

# %%

# %% [markdown]
# ## Now, let's look how the amplitude of the sine signal, relative to the constant bias, changed over time

# %%
PlotlyHelper.plot_pandas(df=bio.get_saved_values(), x_var="SYSTEM TIME", fields="ratio",
                         y_label="ratio",
                         title= "Component of frequency=2, relative to amplitude of constant bias",
                         colors="orange")

# %%
# Same, but wigh log scale for y-axis (and save the figure object in a variable)
fig_freq_2 = PlotlyHelper.plot_pandas(df=bio.get_saved_values(), x_var="SYSTEM TIME", fields="ratio",
                                      log_y=True, y_label="ratio",
                                      title= "Component of frequency=2, relative to amplitude of constant bias<br>(log scale in y-axis)",
                                      colors="orange")
fig_freq_2.show()

# %%

# %%

# %% [markdown]
# # Start a NEW system
# #### Everything same as before EXCEPT the frequency is now 5 times bigger (10 cycles over the length of the system)

# %%
bio = BioSim1D(n_bins=500, chem_data=chem_data)

# %%
bio.inject_sine_conc(chem_label="A", amplitude=10, bias=30, number_cycles=10)   # x5 higher frequency than before

# %%
bio.visualize_system(title_prefix="Initial System State")

# %%
# Show as heatmap
bio.system_heatmaps(title_prefix="Initial System State (as a heatmap)")

# %%
# Take a look at the frequency domain of the concentration values
frequency_data = bio.frequency_analysis(chem_label="A")
frequency_data

# %%
ratio = frequency_data.loc[1, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
ratio

# %%
bio.save_value(data_snapshot={"ratio": ratio})
bio.get_saved_values()

# %%

# %% [markdown]
# ### Start the diffusion

# %%
bio.diffuse(total_duration=10, n_steps=100)

# %%
# Take a look at the frequency domain of the concentration values
frequency_data = bio.frequency_analysis(chem_label="A")
frequency_data.loc[10]

# %%
ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
ratio

# %%
bio.save_value(data_snapshot={"ratio": ratio})
bio.get_saved_values()

# %%
bio.visualize_system()

# %% [markdown]
# ### Continue the diffusion

# %%
for i in range(4):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(chem_label="A")
    ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_value(data_snapshot={"ratio": ratio})

# %%
bio.visualize_system()

# %%
bio.get_saved_values()

# %%
for i in range(10):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(chem_label="A")
    ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_value(data_snapshot={"ratio": ratio})

# %%
bio.visualize_system()

# %%
bio.get_saved_values()

# %%
for i in range(35):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(chem_label="A")
    ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_value(data_snapshot={"ratio": ratio})

# %%
bio.visualize_system()

# %%
bio.get_saved_values()

# %%
for i in range(50):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(chem_label="A")
    ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_value(data_snapshot={"ratio": ratio})

# %%
bio.visualize_system()

# %% [markdown]
# #### **Most of the concentration graph is now flat as a board!**

# %%
bio.get_saved_values()

# %%

# %% [markdown]
# ## Like done earlier for the smaller frequency, let's look how the amplitude of the sine signal, relative to the constant bias, changed over time

# %%
PlotlyHelper.plot_pandas(df=bio.get_saved_values(), x_var="SYSTEM TIME", fields="ratio",
                         y_label="ratio",
                         title= "Component of frequency=10, relative to amplitude of constant bias",
                         colors="red")

# %%
# Same, but wigh log scale for y-axis (and save the figure object in a variable)
fig_freq_10 = PlotlyHelper.plot_pandas(df=bio.get_saved_values(), x_var="SYSTEM TIME", fields="ratio",
                                       log_y=True, y_label="ratio",
                                       title= "Component of frequency=10, relative to amplitude of constant bias<br>(log scale in y-axis)",
                                       colors="red")

fig_freq_10.show()

# %% [markdown]
# ### Resurrect the earlier plot for frequency=2 (stored in a variable)
# ### and then superpose them

# %%
fig_freq_2

# %%
# Combine the last 2 plots
PlotlyHelper.combine_plots(fig_list=[fig_freq_2, fig_freq_10],  
                           title="Superposed plots for frequency 2 and 10",
                           curve_labels=["Component of frequency=2", "Component of frequency=10"])

# %% [markdown]
# ## Note how more more substantial the attenuation of the higher frequency (red) is!

# %%
