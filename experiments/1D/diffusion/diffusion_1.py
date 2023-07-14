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
# # An initial concentration pulse (near the left edge of the system) moving towards equilibrium
#
# The system starts out with a "concentration pulse" in bin 2 (the 3rd bin from the left) - i.e. that bin is initially the only one with a non-zero concentration of the only chemical species.
# Then the system is left undisturbed, and followed to equilibrium.
#
# LAST REVISED: July 14, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from src.life_1D.bio_sim_1d import BioSim1D

import plotly.express as px
import plotly.graph_objects as go

from src.modules.chemicals.chem_data import ChemData as chem
from src.modules.html_log.html_log import HtmlLog as log
from src.modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_heatmap_11", "vue_curves_3"])

# %%
# Set the heatmap parameters (for the log file)
heatmap_pars = {"range": [0, 2.5],
                "outer_width": 850, "outer_height": 150,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# Set the parameters of the line plots
lineplot_pars = {"range": [0, 10],
                "outer_width": 850, "outer_height": 250,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# %%
# Prepare the initial system, with a single non-zero bin, near the left edge of the system
chem_data = chem(names=["A"], diffusion_rates=[0.1])
bio = BioSim1D(n_bins=10, chem_data=chem_data)

bio.inject_conc_to_bin(bin_address=2, species_index=0, delta_conc=10.)

bio.describe_state()

# %%
bio.system_snapshot()

# %%
# Line curve view
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# ONE APPROACH TO CREATE A PLOTLY HEATMAP, using imshow() from plotly.express
fig = px.imshow(bio.system_snapshot().T, 
                title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=True, color_continuous_scale="gray_r")         # text_auto=’.2f’

fig.data[0].xgap=4
fig.data[0].ygap=4

fig.show()

# %%
# ANOTHER APPROACH TO CREATE A PLOTLY HEATMAP, using Heatmap() from plotly.graph_objects
data = go.Heatmap(z=bio.system_snapshot().T,
                    y=['A'],
                    colorscale='gray_r', colorbar={'title': 'Concentration'},
                    xgap=4, ygap=4, texttemplate = '%{z}', hovertemplate= 'Bin number: %{x}<br>Chem. species: %{y}<br>Concentration: %{z}<extra></extra>')

fig = go.Figure(data,
                layout=go.Layout(title=f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}",
                                 xaxis={'title': 'Bin number'}, yaxis={'title': 'Chem. species'}
                                )
               )
fig.show()

# %%
log.write("1-D diffusion to equilibrium of a single species, with Diffusion rate 0.1",
          style=log.h2)
log.write("Initial system state at time t=0:", blanks_before=2, style=log.bold)

# Output a heatmap to the log file
bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output a line plot the log file
bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# %% [markdown]
# # Initial Diffusion Step

# %%
log.write("Advancing to time t=10, with time steps of 0.1 ... ", blanks_before=2, newline=False)

# %%
delta_time = 10.

status = bio.diffuse(total_duration=delta_time, time_step=0.1)
print("\n", status)

log.write(f"After delta time {delta_time}.  TOTAL TIME {bio.system_time}  ({status['steps']} steps taken):")
bio.describe_state(concise=True)

# %%
# Line curve view
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# Heatmap view
fig = px.imshow(bio.system_snapshot().T, 
                title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto='.3f', color_continuous_scale="gray_r")

fig.data[0].xgap=4
fig.data[0].ygap=4

fig.show()

# %%
# Output a heatmap into the log file
bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output a line plot the log file
bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# %% [markdown]
# ## This is still an early stage in the diffusion process; let's advance it more... (Visualization from results shown at selected times)

# %% tags=[]
for i in range(50):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    print(f"\nAfter Delta time {delta_time}.  TOTAL TIME {bio.system_time}  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)

    if i<2 or i==6 or i>=49:
        # Line curve view
        fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
        fig.show()
        
        # Heatmap view
        fig = px.imshow(bio.system_snapshot().T, 
                        title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                        labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                        text_auto='.2f', color_continuous_scale="gray_r")
        fig.data[0].xgap=4
        fig.data[0].ygap=4
        fig.show()
        
        # Output a heatmap to the log file
        bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, header=f"Time {bio.system_time} :\n", graphic_component="vue_heatmap_11")
        # Output a line plot the log file
        bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")


# %% [markdown]
# ## All bins now have essentially uniform concentration
#
# **Mass conservations**: The "10 units of concentration" are now uniformly spread across the 10 bins, leading to a near-constant concentration of 10/10 = **1.0**

# %%
