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
# ### `A` down-regulates `B` , by being the *limiting reagent* in reaction `A + 2 B <-> Y` (mostly forward)
# 1st-order kinetics.   
# If [A] is low and [B] is high, then [B] remains high.  If [A] goes high, [B] goes low.  However, at that point, A can no longer bring B up to any substantial extent.
#
# Single-bin reaction
#
# Based on experiment "reactions_single_compartment/down_regulate_2"
#
# LAST REVISED: May 28, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ChemData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.life_1D.bio_sim_1d import BioSim1D

import plotly.express as px
from src.modules.html_log.html_log import HtmlLog as log
from src.modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system
chem_data = chem(names=["A", "B", "Y"])     # NOTE: Diffusion not applicable (just 1 bin)

# Reaction A + 2 B <-> Y , with 1st-order kinetics for all species
chem_data.add_reaction(reactants=[("A") , (2, "B")], products=[("Y")],
                       forward_rate=8., reverse_rate=2.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %%
bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_uniform_concentration(species_name="A", conc=5.)     # Scarce
bio.set_uniform_concentration(species_name="B", conc=100.)   # Plentiful
# Initially, no "Y" is present

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0 (the only bin in this system)
bio.add_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown] tags=[]
# ### Take the initial system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=30, snapshots={"frequency": 2, "sample_bin": 0})  # At every other step, take a snapshot 
                                                                                      # of all species at bin 0
bio.describe_state()
bio.get_history()

# %% [markdown]
# A, as the scarse limiting reagent, stops the reaction.  
# When A is low, B is also low.

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0))

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations (reaction A + 2 B <-> Y)",
              color_discrete_sequence = ['red', 'blue', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown] tags=[]
# # Now, let's suddenly increase [A]

# %%
bio.set_bin_conc(bin_address=0, species_index=0, conc=40.)
bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0 (the only bin in this system)
bio.add_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=80, snapshots={"frequency": 2, "sample_bin": 0})  # At every other step, take a snapshot 
                                                                                      # of all species at bin 0
bio.describe_state()
bio.get_history()

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0), tolerance=7)

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations (reaction A + 2 B <-> Y)",
              color_discrete_sequence = ['red', 'blue', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# **A**, still the limiting reagent, is again stopping the reaction.  
# The (transiently) high value of [A] led to a high value of [B]

# %% [markdown] tags=[]
# # Let's again suddenly increase [A]

# %%
bio.set_bin_conc(bin_address=0, species_index=0, conc=30.)
bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0 (the only bin in this system)
bio.add_snapshot(bio.bin_snapshot(bin_address = 0))

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=70, snapshots={"frequency": 2, "sample_bin": 0})  # At every other step, take a snapshot 
                                                                                      # of all species at bin 0
bio.describe_state()
bio.get_history()

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(conc=bio.bin_snapshot(bin_address = 0))

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "Y"], 
              title="Changes in concentrations (reaction A + 2 B <-> Y)",
              color_discrete_sequence = ['red', 'blue', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# **A**, again the scarse limiting reagent, stops the reaction yet again   
#
# Note: A can down-regulate B, but it cannot bring it up.

# %% [markdown]
# # For additional exploration, see the experiment "reactions_single_compartment/down_regulate_2"

# %%
