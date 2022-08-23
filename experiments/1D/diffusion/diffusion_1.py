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
# **EXTRA OUTPUT (incl. graphics):** overwritten into the .htm file with the same base name.
#
# LAST REVISED: Aug. 22, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home from the folder containing this notebook

# %%
from experiments.get_notebook_info import get_notebook_basename

from life_1D.bio_sim_1d import BioSim1D as bio

import plotly.express as px

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
# Set the heatmap parameters
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
# Prepare the initial system
chem_data = chem(names=["A"], diffusion_rates=[0.1])
bio.initialize_system(n_bins=10, chem_data=chem_data)

bio.inject_conc_to_bin(bin=2, species_index=0, delta_conc=10.)

bio.describe_state()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
log.write("1-D diffusion to equilibrium of a single species, with Diffusion rate 0.1",
          style=log.h2)
log.write("Initial system state at time t=0:", blanks_before=2, style=log.bold)

# Output a heatmap to the log file
bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output a line plot the log file
bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# %%
log.write("Advancing to time t=10, with time steps of 0.1 ... ", blanks_before=2, newline=False)

# %%
delta_time = 10.

status = bio.diffuse(total_duration=delta_time, time_step=0.1)
print("\n", status)

log.write(f"After delta time {delta_time}.  TOTAL TIME {bio.system_time}  ({status['steps']} steps taken):")
bio.describe_state(concise=True)

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# Output a heatmap into the log file
bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output a line plot the log file
bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# %% [markdown]
# _Note: this is still an early stage in the diffusion process; let's advance it more... (Results at selected times plotted to log file)_
# ---

# %% tags=[]
for i in range(50):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    print(f"\nAfter Delta time {delta_time}.  TOTAL TIME {bio.system_time}  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)

    if i<2 or i==6 or i>=49:
        fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
        fig.show()
        
        # Output a heatmap to the log file
        bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, header=f"Time {bio.system_time} :\n", graphic_component="vue_heatmap_11")
        # Output a line plot the log file
        bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")


# %% [markdown]
# **All cells now have essentially uniform concentration**
#
# The "10 units of concentration" are now uniformly spread across the 10 bins, leading to a near-constant concentration of 10/10 = **1.0**

# %%
