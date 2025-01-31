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
# ### One-bin `A <-> 3B` reaction, with 1st-order kinetics in both directions,
# ### taken to equilibrium
#
# Diffusion not applicable (just 1 bin)
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta34.1)
#
# * [First Step](#reaction_2_sec_2_first_step)
# * [Numerous more steps](#reaction_2_sec_2)
# * [Equilibrium](#reaction_2_sec_2_equilibrium)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData as chem
from life123 import BioSim1D

import plotly.express as px
from life123 import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system
chem_data = chem(names=["A", "B"])     # NOTE: Diffusion not applicable (just 1 bin)

# Reaction A <-> 3B , with 1st-order kinetics in both directions
chem_data.add_reaction(reactants=["A"], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)

bio.describe_state()

# %%
chem_data.describe_reactions()

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0), caption="Initial state")
bio.get_history()

# %%
# Send the plot to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %% [markdown] tags=[]
# ### <a name="reaction_2_sec_2_first_step"></a>First step

# %%
# First step
bio.react(time_step=0.1, n_steps=1, snapshots={"sample_bin": 0})
bio.describe_state()

# %% [markdown]
# _Early in the reaction :_  
# [A] = 15.   [B] = 35.

# %%
bio.get_history()

# %% [markdown]
# ### <a name="reaction_2_sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.1, n_steps=10, snapshots={"sample_bin": 0})

bio.describe_state()

# %% [markdown] tags=[]
# ### <a name="reaction_2_sec_2_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the equilibrium:   [A] = 14.54545455 , [B] = 36.36363636

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0))

# %%
bio.get_history()

# %% [markdown]
# Note how the simulation initially **OVERSHOT** the equilibrium values; the first step was too large!

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="Changes in concentrations with time",
              color_discrete_sequence = ['navy', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
# Same plot, but with smooth line
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="Changes in concentrations with time (smoothed)",
              color_discrete_sequence = ['navy', 'darkorange'],
              labels={"value":"concentration", "variable":"Chemical"},
              line_shape="spline")
fig.show()

# %% [markdown]
# The early **OVERSHOOTING** of the equilibrium values shows prominently in the last plot!

# %%
