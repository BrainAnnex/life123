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
# ## Exploring reaching equilibrium
#
# The system starts out with a pulse in bins near the *left* and the *right* endpoints
#
# LAST REVISED: May 27, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from src.life_1D.bio_sim_1d import BioSim1D

import plotly.express as px

from src.modules.reactions.reaction_data import ChemData as chem
from src.modules.html_log.html_log import HtmlLog as log
from src.modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file
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
# Initialize the system
chem_data = chem(names=["A"], diffusion_rates=[0.1])
bio = BioSim1D(n_bins=9, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=0.)

# Start out with a pulse in bins near the *left* and the *right* endpoints.  
# A total of 20 "units of concentration" is injected
bio.inject_conc_to_bin(species_index=0, bin_address=2, delta_conc=10.)
bio.inject_conc_to_bin(species_index=0, bin_address=6, delta_conc=10.)

bio.describe_state()

# %%
bio.show_system_snapshot()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# Output to the log file
log.write("1-D diffusion to equilibrium of a single species, with Diffusion rate 0.1. Time steps of 0.1",
          style=log.h2)

log.write(f"System state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output a heatmap to the log file
bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, header=f"Time {bio.system_time} :\n", graphic_component="vue_heatmap_11")
# Output a line plot the log file
bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# %% [markdown]
# # Start the simulation steps

# %%
delta_time = 3.

# %%
for i in range(15):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    print(f"\nAfter Delta time {delta_time}.  TOTAL TIME {bio.system_time}  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)

    fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= f"Diffusion. System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
    fig.show()

    #if i<2 or i==6 or i>=14:
    bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, header=f"Time {bio.system_time}\n", graphic_component="vue_heatmap_11")
    # Output a line plot the log file
    bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")


# %% [markdown] tags=[]
# **All cells now have essentially uniform concentration**
#
# The "20 units of concentration" are now uniformly spread across the 9 bins, leading to a near-constant concentration of 20/9 = **2.22**

# %%
