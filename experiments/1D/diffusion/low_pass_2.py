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
# ## The attenuation over time of a single-frequency in the initial concentration
#
# ### The initial system state is a sine wave of frequency 2 (i.e. 2 cycles across the system's length), of amplitude 10, with a baseline (bias) of 30
# ### Afterward, the process is restarted and repeated with a frequency 5 times larger
#
#
# LAST REVISED: Sep. 14, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home from the folder containing this notebook

# %%
from experiments.get_notebook_info import get_notebook_basename

from life_1D.bio_sim_1d import BioSim1D

import plotly.express as px
import plotly.graph_objects as go

from modules.reactions.reaction_data import ReactionData as chem

# %%
# Initialize the system.  We use a RELATIVELY LARGE NUMBER OF BINS, 
# to captures the finer changes in the frequency components of the concentration function
chem_data = chem(names=["A"], diffusion_rates=[0.5])
bio = BioSim1D(n_bins=500, chem_data=chem_data)

# %% [markdown]
# ## Initial Preparation -
# ### Start with a sinusoidal concentration, with exactly 2 cycles over the length of the system
# #### (of amplitude 10 and bias value of 30)

# %%
bio.inject_sine_conc(species_name="A", amplitude=10, bias=30, frequency=2)

# %%
bio.show_system_snapshot()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= "Initial System State",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# Show as heatmap
fig = px.imshow(bio.system_snapshot().T, 
                title= "Initial System State (as a heatmap)", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=False, color_continuous_scale="gray_r") 

fig.data[0].xgap=0
fig.data[0].ygap=0

fig.show()

# %%
# Take a look at the frequency domain of the concentration values
frequency_data = bio.frequency_analysis(species_name="A")
frequency_data

# %%
ratio = frequency_data.loc[1, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
ratio

# %%
bio.save_snapshot(data_snapshot={"ratio": ratio})
bio.get_history()

# %% [markdown]
# # Start the diffusion - to time t=10

# %%
bio.diffuse(total_duration=10, n_steps=100)

# %%
# Take a look at the frequency domain of the concentration values
frequency_data = bio.frequency_analysis(species_name="A")
frequency_data

# %%
ratio = frequency_data.loc[2, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
ratio

# %%
bio.save_snapshot(data_snapshot={"ratio": ratio})
bio.get_history()

# %% [markdown]
# ## Do 49 more rounds of diffusion - to time t = 500

# %%
for i in range(49):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(species_name="A")
    ratio = frequency_data.loc[2, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_snapshot(data_snapshot={"ratio": ratio})

# %%
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"System State at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# #### Note how the curve is beginning to flatten; the peaks are no longer between 20 and 40

# %% [markdown]
# ## Do 150 more rounds of diffusion - to time t = 2000

# %%
for i in range(150):
    bio.diffuse(total_duration=10, n_steps=75)    # Notice the gradual decreas of the number of intermediate steps, given the smaller gradient
    frequency_data = bio.frequency_analysis(species_name="A")
    ratio = frequency_data.loc[2, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_snapshot(data_snapshot={"ratio": ratio})

# %%
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# #### Note how the curve is flatter still - and even beginning to lose shape at the boundary 

# %% [markdown]
# ## Do 800 more rounds of diffusion - to time t = 10000

# %%
for i in range(800):
    bio.diffuse(total_duration=10, n_steps=20)   # Note how we're gradually increasing the time steps, because the gradiants are now smaller
    frequency_data = bio.frequency_analysis(species_name="A")
    ratio = frequency_data.loc[2, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_snapshot(data_snapshot={"ratio": ratio})

# %%
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# #### Getting decisively flatter, as it approaches equilibrium

# %% [markdown]
# ## Now, let's look how the amplitude of the sine signal, relative to the constant bias, changed over time

# %%
fig = px.line(data_frame=bio.get_history(), y=["ratio"], 
          title= "Component of frequency=2, relative to amplitude of constant bias",
          color_discrete_sequence = ['purple'],
          labels={"value":"ratio", "variable":"Freq. 2", "index":"Time (groups of 10 units)"})
fig.show()

# %%
fig = px.line(data_frame=bio.get_history(), y=["ratio"], 
              title= "Component of frequency=2, relative to amplitude of constant bias (log scale in y-axis)",
              color_discrete_sequence = ['purple'],
              labels={"value":"ratio", "variable":"Freq. 2", "index":"Time (groups of 10 units)"},
              log_y=True)
fig.show()

# %%
saved_freq_2_df = bio.get_history().copy(deep=True)   # Save the whole dataframe, for later use

# %% [markdown]
# # Start a NEW system
# #### Everything same as before EXCEPT the frequency is now 5 times bigger (10 cycles over the length of the system)

# %%
bio = BioSim1D(n_bins=500, chem_data=chem_data)

# %%
bio.inject_sine_conc(species_name="A", amplitude=10, bias=30, frequency=10)   # x5 higher frequency than before

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= "Initial System State",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# Show as heatmap
fig = px.imshow(bio.system_snapshot().T, 
                title= "Initial System State (as a heatmap)", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=False, color_continuous_scale="gray_r") 

fig.data[0].xgap=0
fig.data[0].ygap=0

fig.show()

# %%
# Take a look at the frequency domain of the concentration values
frequency_data = bio.frequency_analysis(species_name="A")
frequency_data

# %%
ratio = frequency_data.loc[1, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
ratio

# %%
bio.save_snapshot(data_snapshot={"ratio": ratio})
bio.get_history()

# %%
bio.diffuse(total_duration=10, n_steps=100)

# %%
# Take a look at the frequency domain of the concentration values
frequency_data = bio.frequency_analysis(species_name="A")
frequency_data.loc[10]

# %%
ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
ratio

# %%
bio.save_snapshot(data_snapshot={"ratio": ratio})
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= "Initial System State",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
for i in range(4):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(species_name="A")
    ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_snapshot(data_snapshot={"ratio": ratio})

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"System State at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
bio.get_history()

# %%
for i in range(10):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(species_name="A")
    ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_snapshot(data_snapshot={"ratio": ratio})

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"System State at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
bio.get_history()

# %%
for i in range(35):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(species_name="A")
    ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_snapshot(data_snapshot={"ratio": ratio})

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"System State at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
bio.get_history()

# %%
for i in range(50):
    bio.diffuse(total_duration=10, n_steps=100)
    frequency_data = bio.frequency_analysis(species_name="A")
    ratio = frequency_data.loc[10, "Relative Amplitude"] / frequency_data.loc[0, "Relative Amplitude"]
    bio.save_snapshot(data_snapshot={"ratio": ratio})

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"System State at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# **Most of the concentration graph is now flat as a board!**

# %%
bio.get_history()

# %% [markdown]
# ## Like done earlier for the smaller frequency, let's look how the amplitude of the sine signal, relative to the constant bias, changed over time

# %%
fig = px.line(data_frame=bio.get_history(), y=["ratio"], 
          title= "Component of frequency=10, relative to amplitude of constant bias",
          color_discrete_sequence = ['green'],
          labels={"value":"ratio", "variable":"Freq. 10", "index":"Time (groups of 10 units)"})
fig.show()

# %%
fig = px.line(data_frame=bio.get_history(), y=["ratio"], 
              title= "Component of frequency=10, relative to amplitude of constant bias (log scale in y-axis)",
              color_discrete_sequence = ['green'],
              labels={"value":"ratio", "variable":"Freq. 10", "index":"Time (groups of 10 units)"},
              log_y=True)
fig.show()

# %% [markdown]
# ### Resurrect the earlier plot for frequency=2 (using a saved dataset)
# ### and then superpose them

# %%
fig_2 = px.line(data_frame=saved_freq_2_df, y=["ratio"], 
              title= "Component of frequency=2, relative to amplitude of constant bias (log scale in y-axis)",
              color_discrete_sequence = ['purple'],
              labels={"value":"ratio", "variable":"Freq. 2", "index":"Time (groups of 10 units)"},
              log_y=True)
fig_2.show()

# %%
# Combine the last 2 plots
all_fig = go.Figure(data=fig.data + fig_2.data)    # Note that the + is concatenating lists
all_fig.update_layout(title="Superposed plots for frequency 2 (purple) and 10 (green)")
all_fig.show()

# %% [markdown]
# ## Note how more more substantial the attenuation of the higher frequency (green) is!

# %%
