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
# ## High Frequencies in the concentration get smoothed out first  
# ### Ultimately, the only frequency passing thru is zero! (i.e., uniform concentration at equilibrium)
# #### We explore how an initial concentration with 3 different sinusoidal frequencies fares in the course of a diffusion to equilibrium
#
# **The initial system state will consist of:**  
#     0 - A constant baseline of value 30  
#     1 - A sine wave of frequency 1 (1 cycle across the system's length), of amplitude 10  
#     2 - A sine wave of frequency 10 , of amplitude 4  
#     3 - A sine wave of frequency 40 , of amplitude 2  
#
# LAST REVISED: Nov. 28, 2022

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
from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file
GraphicLog.config(filename=log_file,
                  components=["vue_heatmap_11", "vue_curves_3"])

# %%
# Set the heatmap parameters
heatmap_pars = {"range": [10, 50],
                "outer_width": 850, "outer_height": 150,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# Set the parameters of the line plots
lineplot_pars = {"range": [10, 50],
                "outer_width": 850, "outer_height": 250,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# %%
# Initialize the system.  We use a RELATIVELY LARGE NUMBER OF BINS, 
# to captures the many changes in the high-frequency component
chem_data = chem(names=["A"], diffusion_rates=[0.5])
bio = BioSim1D(n_bins=500, chem_data=chem_data)

# %% [markdown]
# ## PART 1 (of 3) of Initial Preparation -
# ### Start with a sinusoidal concentration, with exactly 1 cycle over the length of the system
# #### (notice the bias value of 30; it will be seen at the eventual equilibrium, at the end)

# %%
bio.inject_sine_conc(species_name="A", amplitude=10, bias=30, frequency=1)

# %%
bio.show_system_snapshot()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= "Low-frequency component of the Initial System State",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# Show as heatmap
fig = px.imshow(bio.system_snapshot().T, 
                title= "Low-frequency component of the Initial System State (as a heatmap)", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=False, color_continuous_scale="gray_r") 

fig.data[0].xgap=0
fig.data[0].ygap=0

fig.show()

# %%
# Output to the log file
log.write("Diffusion as a Low-Pass Filter", style=log.h3)

log.write(f"Low-frequency component of the Initial System State:", blanks_before=2, style=log.bold)

# Output a heatmap to the log file
bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, header=f"Time {bio.system_time} :\n", graphic_component="vue_heatmap_11")
# Output a line plot the log file
bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# %% [markdown]
# ## PART 2 (of 3) of Initial Preparation -
# ### Now add a higher-frequency component (10 cycles over the length of the system)

# %%
bio.inject_sine_conc(species_name="A", amplitude=4, bias=0, frequency=10)

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= "Low- and mid-frequency components of the Initial System State",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# Show as heatmap
fig = px.imshow(bio.system_snapshot().T, 
                title= "Low- and mid-frequency components of the Initial System State (as a heatmap)", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=False, color_continuous_scale="gray_r") 

fig.data[0].xgap=0
fig.data[0].ygap=0

fig.show()

# %% [markdown]
# ## PART 3 (of 3) of Initial Preparation -
# ### To complete the preparation of a 3-frequency component initial state, add another, even higher, frequency component 
# #### (40 cycles over the length of the system)

# %%
bio.inject_sine_conc(species_name="A", amplitude=2, bias=0, frequency=40)

fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= "Initial System State with 3 superposed frequencies",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()


# Show as heatmap
fig = px.imshow(bio.system_snapshot().T, 
                title= "Initial System State with 3 superposed frequencies (as a heatmap)", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=False, color_continuous_scale="gray_r") 

fig.data[0].xgap=0
fig.data[0].ygap=0

fig.show()

# %%
# Take a look at the frequency domain of the concentration values
bio.frequency_analysis(species_name="A")

# %% [markdown]
# # Start the diffusion steps

# %%
bio.diffuse(total_duration=10, time_step=0.1)

# Show as a line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# ### After the initial diffusion (at t=10), the highest frequency is largely gone

# %%
# Take a look at the frequency domain of the concentration values
bio.frequency_analysis(species_name="A")

# %%
# A lot of tiny frequency components are now present; take just the largest 4
bio.frequency_analysis(species_name="A", n_largest=4)

# %%
# Advance the diffusion
bio.diffuse(total_duration=20, time_step=0.1)

# Show as a line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# ### After additional diffusion (at t=30), the highest frequency (40 cycles) can no longer be visually detected.
# #### The 2 lower frequency are still recognizable
# #### The 40-cycle frequency is not even among the largest of the tiny values in the spurious frequencies of the distorted signal

# %%
bio.frequency_analysis(species_name="A", n_largest=10)

# %%
# Advance the diffusion
bio.diffuse(total_duration=90, time_step=0.1)

# Show as a line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# ### By now (at t=120), even the middle frequency (10 cycles) is notably attenuated

# %%
bio.frequency_analysis(species_name="A", n_largest=3)

# %%
# Advance the diffusion
bio.diffuse(total_duration=100, time_step=0.1)

# Show as a line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# ### With still more diffusion (at t=220), the middle frequency is only weakly visible

# %%
bio.frequency_analysis(species_name="A", n_largest=3)

# %%
# Advance the diffusion
bio.diffuse(total_duration=180, time_step=0.1)

# Show as a line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# ### By t=400, the middle frequency is barely noticeable.  The low frequency (1 cycle) is still prominent

# %%
bio.frequency_analysis(species_name="A", n_largest=3)

# %%
# Advance the diffusion
bio.diffuse(total_duration=600, time_step=0.3)

# Show as a line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# ### By t=1,000 there's no visual indication of the middle frequency (f=10)

# %%
bio.frequency_analysis(species_name="A", n_largest=10)

# %%
# Advance the diffusion
bio.diffuse(total_duration=8000, time_step=.5)

# Show as a line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# ### By t=9,000 even the lowest frequency has lost a major part of its sinusoidal shape

# %%
bio.frequency_analysis(species_name="A", n_largest=2)

# %% [markdown]
# Note how the zero-frequency is now gaining over the baseline 1-cycle signal

# %%
# Advance the diffusion
bio.diffuse(total_duration=91000, time_step=.5)

# Show as a line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
bio.frequency_analysis(species_name="A", n_largest=2)

# %% [markdown]
# ### By t=100,000 the system is clearly approaching equilibrium

# %%
# Advance the diffusion
bio.diffuse(total_duration=100000, time_step=.6)

# Show as a line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
          title= f"Diffusion. System snapshot at time t={bio.system_time}",
          color_discrete_sequence = ['red'],
          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# ### By t=200,000 the system is getting close to equilibrium about the value 30, which was the original "bias" (unvarying component) of the baseline frequency (the higher-frequency signals didn't have any bias)

# %%
bio.frequency_analysis(species_name="A", n_largest=2)
