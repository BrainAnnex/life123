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
# ### Reaction  A + B <-> C, mostly forward and with 1st-order kinetics for each species,
# ### taken to equilibrium
#
# Initial concentrations of A and B are spacially separated to the opposite ends of the system;
# as a result, no C is being generated.
#
# But, as soon as A and B, from their respective distant originating points at the edges, 
# diffuse into the middle - and into each other - the reaction starts,
# consuming both A and B (the forward reaction is much more substantial than the reverse one),
# until an equilibrium is reached in both diffusion and reactions.
#
# A LOT of plots are sent to the log file from this experiment; the reason is to compare two
# graphic elements, "vue_curves_3" and "vue_curves_4"
#
# LAST REVISED: May 6, 2024

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from src.life_1D.bio_sim_1d import BioSim1D

import plotly.express as px
from src.modules.chemicals.chem_data import ChemData as chem
from src.modules.html_log.html_log import HtmlLog as log
from src.modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_heatmap_11", "vue_curves_3", "vue_curves_4", "vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system
chem_data = chem(names=["A", "B", "C"], diffusion_rates=[50., 50., 1.])



# Reaction A + B <-> C , with 1st-order kinetics for each species; note that it's mostly in the forward direction
chem_data.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=20., reverse_rate=2.)
bio = BioSim1D(n_bins=7, chem_data=chem_data)

# %%
bio.show_system_snapshot()

# %%
chem_data.describe_reactions()

# %%
# Send a header and a plot to the HTML log file
log.write("Reaction:  A + B <-> C",
          style=log.h2)
chem_data.plot_reaction_network("vue_cytoscape_2")

# %%
# Set the heatmap parameters
heatmap_pars = {"range": [0, 20],
                "outer_width": 850, "outer_height": 100,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# Set the parameters of the line plots (for now, same for single-curve and multiple-curves)
lineplot_pars = {"range": [0, 20],
                "outer_width": 850, "outer_height": 200,
                "margins": {"top": 30, "right": 30, "bottom": 30, "left": 55}
                }

# %% [markdown]
# # Inject initial concentrations of A and B at opposite ends of the system

# %%
bio.set_bin_conc(bin_address=0, species_index=0, conc=20.)
bio.set_bin_conc(bin_address=6, species_index=1, conc=20.)

bio.describe_state()

# %%
bio.show_system_snapshot()

# %%
# Save the state of the concentrations of all species at the middle bin
bio.add_snapshot(bio.bin_snapshot(bin_address = 3))
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A", "B", "C"], 
              title= f"A + B <-> C . System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red', 'orange', 'green'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
log.write(f"Initial system state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output to the log file a heatmap for each chemical species
for i in range(bio.n_species):
    log.write(f"{bio.chem_data.get_name(i)}:", also_print=False)
    bio.single_species_heatmap(species_index=i, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output to the log file a one-curve line plot for each chemical species
for i in range(bio.n_species):
    log.write(f"{bio.chem_data.get_name(i)}:", also_print=False)
    bio.single_species_line_plot(species_index=i, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# Output to the log file a line plot for ALL the chemicals together (same color as used for plotly elsewhere)
bio.line_plot(plot_pars=lineplot_pars, graphic_component="vue_curves_4", color_mapping={0: 'red', 1: 'orange', 2: 'green'})

# %% [markdown] tags=[]
# ### <a name="sec_2_first_step"></a>First step

# %%
delta_t = 0.002   # This will be our time "quantum" for this experiment

# %%
# First step
bio.react_diffuse(time_step=delta_t, n_steps=1)
bio.describe_state()

# %% [markdown]
# _After the first delta_t time step_:
#
#   Species 0 (A). Diff rate: 50.0. Conc:  [18.  2.  0.  0.  0.  0.  0.]
#   
#   Species 1 (B). Diff rate: 50.0. Conc:  [ 0.  0.  0.  0.  0.  2. 18.]
#   
#   Species 2 (C). Diff rate: 1.0. Conc:  [0. 0. 0. 0. 0. 0. 0.]
#

# %%
bio.show_system_snapshot()

# %%
# Save the state of the concentrations of all species at the middle bin
bio.add_snapshot(bio.bin_snapshot(bin_address = 3))
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A", "B", "C"], 
              title= f"A + B <-> C . System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red', 'orange', 'green'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
log.write(f"System state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output to the log file a heatmap for each chemical species
for i in range(bio.n_species):
    log.write(f"{bio.chem_data.get_name(i)}:", also_print=False)
    bio.single_species_heatmap(species_index=i, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output to the log file a one-curve line plot for each chemical species
for i in range(bio.n_species):
    log.write(f"{bio.chem_data.get_name(i)}:", also_print=False)
    bio.single_species_line_plot(species_index=i, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# Output to the log file a line plot for ALL the chemicals together (same color as used for plotly elsewhere)
bio.line_plot(plot_pars=lineplot_pars, graphic_component="vue_curves_4", color_mapping={0: 'red', 1: 'orange', 2: 'green'})

# %% [markdown]
# ### <a name="sec_2"></a>Several more steps

# %%
# Continue with several delta_t steps
for _ in range(7):
    bio.react_diffuse(time_step=delta_t, n_steps=1)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
# Save the state of the concentrations of all species at the middle bin
bio.add_snapshot(bio.bin_snapshot(bin_address = 3))
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A", "B", "C"], 
              title= f"A + B <-> C . System snapshot (interpolated) at time t={bio.system_time}",
              color_discrete_sequence = ['red', 'orange', 'green'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"},
              line_shape="spline")
fig.show()

# %% [markdown]
# A is continuing to diffuse from the left.  
# B is continuing to diffuse from the right.  
# They're finally beginning to overlap in the middle bin

# %%
log.write(f"System state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output to the log file a heatmap for each chemical species
for i in range(bio.n_species):
    log.write(f"{bio.chem_data.get_name(i)}:", also_print=False)
    bio.single_species_heatmap(species_index=i, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output to the log file a one-curve line plot for each chemical species
for i in range(bio.n_species):
    log.write(f"{bio.chem_data.get_name(i)}:", also_print=False)
    bio.single_species_line_plot(species_index=i, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# Output to the log file a line plot for ALL the chemicals together (same color as used for plotly elsewhere)
bio.line_plot(plot_pars=lineplot_pars, graphic_component="vue_curves_4", color_mapping={0: 'red', 1: 'orange', 2: 'green'})

# %% [markdown]
# ### <a name="sec_2"></a>Several group of longer runs

# %%
# Now, do several group of longer runs
for _ in range(4):
    print("\n\n+ 10 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=10)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
# Save the state of the concentrations of all species at the middle bin
bio.add_snapshot(bio.bin_snapshot(bin_address = 3))
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A", "B", "C"], 
              title= f"A + B <-> C . System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red', 'orange', 'green'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"},
              line_shape="spline")
fig.show()

# %% [markdown]
# A is continuing to diffuse from the left.  
# B is continuing to diffuse from the right.  
# By now, they're overlapping in the middle bin sufficiently to react and generate C

# %%
log.write(f"System state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output to the log file a heatmap for each chemical species
for i in range(bio.n_species):
    log.write(f"{bio.chem_data.get_name(i)}:", also_print=False)
    bio.single_species_heatmap(species_index=i, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output to the log file a one-curve line plot for each chemical species
for i in range(bio.n_species):
    log.write(f"{bio.chem_data.get_name(i)}:", also_print=False)
    bio.single_species_line_plot(species_index=i, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# Output to the log file a line plot for ALL the chemicals together (same color as used for plotly elsewhere)
bio.line_plot(plot_pars=lineplot_pars, graphic_component="vue_curves_4", color_mapping={0: 'red', 1: 'orange', 2: 'green'})

# %%
# Continue the simulation
for _ in range(4):
    print("\n\n+++ 30 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=30)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
# Save the state of the concentrations of all species at the middle bin
bio.add_snapshot(bio.bin_snapshot(bin_address = 3))
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A", "B", "C"], 
              title= f"A + B <-> C . System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red', 'orange', 'green'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"},
              line_shape="spline")
fig.show()

# %%
log.write(f"System state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output to the log file a heatmap for each chemical species
for i in range(3):
    bio.single_species_heatmap(species_index=i, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output to the log file a one-curve line plot for each chemical species
for i in range(3):
    bio.single_species_line_plot(species_index=i, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# Output to the log file a line plot for ALL the chemicals together (same color as used for plotly elsewhere)
bio.line_plot(plot_pars=lineplot_pars, graphic_component="vue_curves_4", color_mapping={0: 'red', 1: 'orange', 2: 'green'})

# %%
# Continue the simulation
for _ in range(4):
    print("\n+++++ 50 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=50)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
# Save the state of the concentrations of all species at the middle bin
bio.add_snapshot(bio.bin_snapshot(bin_address = 3))
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A", "B", "C"], 
              title= f"A + B <-> C . System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red', 'orange', 'green'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"},
              line_shape="spline")
fig.show()

# %%
log.write(f"System state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output to the log file a heatmap for each chemical species
for i in range(3):
    bio.single_species_heatmap(species_index=i, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output to the log file a one-curve line plot for each chemical species
for i in range(3):
    bio.single_species_line_plot(species_index=i, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# Output to the log file a line plot for ALL the chemicals together (same color as used for plotly elsewhere)
bio.line_plot(plot_pars=lineplot_pars, graphic_component="vue_curves_4", color_mapping={0: 'red', 1: 'orange', 2: 'green'})

# %%
# Continue the simulation
for _ in range(4):
    print("\n+++++++++++++++ 150 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=150)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
# Save the state of the concentrations of all species at the middle bin
bio.add_snapshot(bio.bin_snapshot(bin_address = 3))
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A", "B", "C"], 
              title= f"A + B <-> C . System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red', 'orange', 'green'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"},
              line_shape="spline")
fig.show()

# %%
log.write(f"System state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output to the log file a heatmap for each chemical species
for i in range(3):
    bio.single_species_heatmap(species_index=i, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output to the log file a one-curve line plot for each chemical species
for i in range(3):
    bio.single_species_line_plot(species_index=i, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# Output to the log file a line plot for ALL the chemicals together (same color as used for plotly elsewhere)
bio.line_plot(plot_pars=lineplot_pars, graphic_component="vue_curves_4", color_mapping={0: 'red', 1: 'orange', 2: 'green'})

# %%
# Continue the simulation
for _ in range(2):
    print("\n++++++++++ ... ++++++++++ 1,000 steps later:")
    bio.react_diffuse(time_step=delta_t, n_steps=1000)
    bio.describe_state(concise=True)

# %%
bio.show_system_snapshot()

# %%
# Save the state of the concentrations of all species at the middle bin
bio.add_snapshot(bio.bin_snapshot(bin_address = 3))
bio.get_history()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A", "B", "C"], 
              title= f"A + B <-> C . System snapshot at time t={bio.system_time}",
              color_discrete_sequence = ['red', 'orange', 'green'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"},
              line_shape="spline")
fig.show()

# %%
log.write(f"System state at time t={bio.system_time}:", blanks_before=2, style=log.bold)

# Output to the log file a heatmap for each chemical species
for i in range(3):
    bio.single_species_heatmap(species_index=i, heatmap_pars=heatmap_pars, graphic_component="vue_heatmap_11")

# Output to the log file a one-curve line plot for each chemical species
for i in range(3):
    bio.single_species_line_plot(species_index=i, plot_pars=lineplot_pars, graphic_component="vue_curves_3")

# Output to the log file a line plot for ALL the chemicals together (same color as used for plotly elsewhere)
bio.line_plot(plot_pars=lineplot_pars, graphic_component="vue_curves_4", color_mapping={0: 'red', 1: 'orange', 2: 'green'})

# %% [markdown] tags=[]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %%
# Verify equilibrium concentrations (sampled in the 1st bin; at this point, all bins have equilibrated)
A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
print(f"\nRatio of equilibrium concentrations ((C_eq) / (A_eq * B_eq)) : {(C_eq) / (A_eq * B_eq)}")
print(f"Ratio of forward/reverse rates: {chem_data.get_forward_rate(0) / chem_data.get_reverse_rate(0)}")
# Both are essentially equal, as expected

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %% tags=[]
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "C"], 
              title="Reactions:  A + B <-> C . Changes in concentrations in the MIDDLE bin",
              color_discrete_sequence = ['navy', 'cyan', 'red'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# A and B overlap on the plot, due to the symmetry of the system.  
# Initially, in the middle bin, neither A nor B are present; over time they diffuse there... but then they react and get consumed (producing C), to an equilibrium value.
# C gradually diffuses away.

# %%
