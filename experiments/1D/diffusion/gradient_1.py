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
# ## Attaining a Concentration Gradient  
# ### by continuosly injecting and draining, at opposite ends
#
# The system starts out with a uniform concentration.  
# Then identical concentrations are repeatedly *injected to the left* and *drained from the right*
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta34.1)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from life123 import BioSim1D

import plotly.express as px

from life123 import ChemData as chem
from life123 import HtmlLog as log
from life123 import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file
GraphicLog.config(filename=log_file,
                  components=["vue_heatmap_11", "vue_curves_3"])

# %%
# Set the heatmap parameters
heatmap_pars = {"range": [75, 125],
                "outer_width": 850, "outer_height": 150,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# Set the parameters of the line plots
lineplot_pars = {"range": [75, 125],
                "outer_width": 850, "outer_height": 250,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# %%
# Initialize the system with a uniform concentration (of the only species)
chem_data = chem(names=["A"], diffusion_rates=[0.6])
bio = BioSim1D(n_bins=9, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=100.)

bio.describe_state()

# %%
bio.show_system_snapshot()

# %%
# Visualize the system's initial state
bio.visualize_system(caption="Diffusion")

# %%
# Show as heatmap
fig = px.imshow(bio.system_snapshot().T, 
                title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=True, color_continuous_scale="gray_r")         # text_auto='.2f'

fig.data[0].xgap=4
fig.data[0].ygap=4

fig.show()

# %%
# Output to the log file
log.write("Creation of a gradient", style=log.h3)

log.write(f"System state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output a heatmap to the log file
bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, header=f"Time {bio.system_time} :\n", graphic_component="vue_heatmap_11")
# Output a line plot the log file
bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# %% [markdown]
# # Start the simulation steps

# %%
delta_time = 1.

# %% tags=[]
for i in range(501):
    # Inject to the leftmost bin
    bio.inject_conc_to_bin(bin_address=0, species_index=0, delta_conc=4, zero_clip = False)
    
    # Drain from the rightmost bin
    bio.inject_conc_to_bin(bin_address=8, species_index=0, delta_conc=-4, zero_clip = False)
    
    # Note: the NET GAIN of moles of A in the system is zero!
    
    # Diffuse for the time span delta_time
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    
    if (i <= 12 and i%3 == 0) or (i%100 == 0):   # Display more frequently initially
        print()
        bio.describe_state(concise=True)
        
        # Show the system state as a line plot
        bio.visualize_system(caption="Diffusion")
        
        # Show as heatmap
        fig = px.imshow(bio.system_snapshot().T, 
                        title= f"Diffusion. System snapshot as a heatmap at time t={bio.system_time}", 
                        labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                        text_auto='.2f', color_continuous_scale="gray_r")

        fig.data[0].xgap=4
        fig.data[0].ygap=4

        fig.show()

        # Output a heatmap the log file
        bio.single_species_heatmap(species_index=0, heatmap_pars=heatmap_pars, header=f"Time {bio.system_time}\n", graphic_component="vue_heatmap_11")
        # Output a line plot the log file
        bio.single_species_line_plot(species_index=0, plot_pars=lineplot_pars, graphic_component="vue_curves_3")


# %% [markdown]
# ### By now, the gradient has stabilized with  
# ### [A] = 124.67065159 on the left and [A] = 75.32934841 on the right

# %% [markdown]
# Note: if the drain is too large, relative to the diffusion rate, the smaller concentration could "saturate" at zero

# %%
