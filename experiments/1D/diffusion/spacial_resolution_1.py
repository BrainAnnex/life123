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
# # Change of delta_x in diffusion computations
#
# LAST REVISED: Sep. 2, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %%
from experiments.get_notebook_info import get_notebook_basename

from life_1D.bio_sim_1d import BioSim1D as bio

import plotly.express as px
import plotly.graph_objects as go

from modules.chemicals.chemicals import Chemicals as chem
from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_heatmap_11", "vue_curves_3"])

# %%
# Set the heatmap parameters (for the log file)
heatmap_pars = {"range": [0, 150],
                "outer_width": 850, "outer_height": 150,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# Set the parameters of the line plots
lineplot_pars = {"range": [0, 150],
                "outer_width": 850, "outer_height": 250,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# %%
# Prepare the initial system
chem_data = chem(names=["A"], diffusion_rates=[0.1])

conc_list=[10,13,17,21,25,28,30,38,42,55,65,47,35,32,27,23,20,17,14,8,3,10,16,18,
           20,25,30,35,40,65,85,115,150,92,73,69,65,50,42,36,20,45,50,55,69,82,95,
           77,60,43,37,31,25,22,20,18,15,11,9, 8]

bio.initialize_system(n_bins=len(conc_list), chem_data=chem_data)

bio.set_species_conc(species_name="A", conc_list=conc_list)

bio.describe_state()

# %%
# Line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# Show as a heatmap
fig = px.imshow(bio.system_snapshot().T, 
                title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=False, color_continuous_scale="gray_r")
fig.data[0].xgap=2
fig.data[0].ygap=4

fig.show()

# %%
log.write("1-D diffusion of a single species, with Diffusion rate 0.1",
          style=log.h2)
log.write("Initial system state at time t=0:", blanks_before=2, style=log.bold)

# Output a heatmap to the log file
bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output a line plot the log file
bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# %%
bio.describe_state(concise=True)

# %% [markdown]
# # Populate the data set with more bins, using interpolated concentration values
# ### IMPORTANT: we're **NOT** changing spacial resolution here; we're just creating a smoother dataset, as *our initial system state*

# %%
bio.smooth_spacial_resolution()
bio.describe_state()

# %%
bio.n_bins

# %% [markdown]
# #### This system setup will be our starting point in exploring diffusion using different spacial resolutions

# %%
# Line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# Show as a heatmap
fig = px.imshow(bio.system_snapshot().T, 
                title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=False, color_continuous_scale="gray_r")
fig.data[0].xgap=2
fig.data[0].ygap=4

fig.show()

# %% [markdown]
# # Diffusions with dx = 1

# %%
bio.describe_state(concise=True)

# %%
bio.diffuse(total_duration=7, time_step=0.0005)
bio.describe_state(concise=True)

# %%
# Line plot
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# ### Enough time has proceeded some smoothing, and non-puny changes in most values - but still nowhere near equilibrium

# %%
