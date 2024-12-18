# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# ### A simple `A <-> B` reaction between 2 species with initial uniform concentrations across 3 bins,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# Diffusion NOT taken into account
#
# See also the experiment _"reactions_single_compartment/react_1"_ 
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta34.1)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData as chem
from life123 import BioSim1D

import plotly.express as px
from life123 import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # Initialize the System

# %% tags=[]
# Initialize the system
chem_data = chem(names=["A", "B"])              # Diffusion NOT taken into account
bio = BioSim1D(n_bins=3, chem_data=chem_data)   # We'll specify the reactions later

bio.set_uniform_concentration(species_name="A", conc=10.)
bio.set_uniform_concentration(species_name="B", conc=50.)

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0), caption="Initial state")

# %%
bio.get_history()

# %%
# Specify the reaction

# Reaction A <-> B , with 1st-order kinetics in both directions
chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

print("Number of reactions: ", chem_data.number_of_reactions())

# %%
chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %% [markdown] tags=[]
# ### <a name="sec_1"></a>First step

# %%
# First step of reaction
bio.react(time_step=0.1, n_steps=1)
bio.describe_state()

# %% [markdown]
# NOTE: the concentration of species A is increasing, while that of species B is decreasing.
# All bins have identical concentrations; so, there's no diffusion (and we're not attempting to compute it):  
# [[17. 17. 17.]  
#  [43. 43. 43.]]

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# ### <a name="sec_more_steps"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.1, n_steps=10, snapshots={"sample_bin": 0})

bio.describe_state()

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order reactions),
#  the systems settles in the following equilibrium:
#
# [A] = 23.99316406
#  
# [B] = 36.00683594
#

# %%
bio.bin_snapshot(bin_address = 0)

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0))

# %%
# Save the state of the concentrations of all species at bin 0
#bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="Changes in concentrations with time",
              color_discrete_sequence = ['navy', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
# Same plot, but with a smoothed line
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="Changes in concentrations with time (smoothed)",
              color_discrete_sequence = ['navy', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"},
              line_shape="spline")
fig.show()

# %% [markdown]
# ## For more in-depth analysis of this reaction, see the experiment _"reactions_single_compartment/react_1"_ 

# %%
