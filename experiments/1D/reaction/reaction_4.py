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
# ## One-bin Association/Dissociation reaction `A + B <-> C`
# ### with 1st-order kinetics for each species, taken to equilibrium
#
# Diffusion not applicable (just 1 bin)
#
# See also the experiment _"reactions_single_compartment/react_3"_ 
#
# LAST REVISED: May 28, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.life_1D.bio_sim_1d import BioSim1D

import plotly.express as px
from src.modules.visualization.graphic_log import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Specify the chemicals
chem_data = chem(names=["A", "B", "C"])     # NOTE: Diffusion not applicable (using just 1 bin)


# Reaction A + B <-> C , with 1st-order kinetics for each species
chem_data.add_reaction(reactants=["A" , "B"], products=["C"],
                       forward_rate=5., reverse_rate=2.)

chem_data.describe_reactions()

# %%
# Initialize the system
bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)
bio.set_uniform_concentration(species_index=2, conc=20.)

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0), caption="Initial state")
bio.get_history()

# %%
# Send the plot to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown] tags=[]
# ### <a name="sec_2_first_step"></a>First step

# %%
# First step
bio.react(time_step=0.002, n_steps=1)
bio.describe_state()

# %% [markdown]
# _Early in the reaction :_
# [A] = 5.08 ,  [B] = 45.08 ,  [C] = [24.92]

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.002, n_steps=29, snapshots={"sample_bin": 0})

bio.describe_state()

# %% [markdown] tags=[]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:  
# [A] = 0.29487831 , [B] = 40.29487831 , [C] = 29.70512169

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
# Save the state of the concentrations of all species at bin 0
bio.get_history()

# %% [markdown]
# ## Note: "A" (now largely depleted) is largely the limiting reagent

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "C"], 
              title="Reaction A + B <-> C .  Changes in concentrations with time",
              color_discrete_sequence = ['red', 'violet', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# ## For more in-depth analysis of this reaction, including variable time steps, see the experiment _"reactions_single_compartment/react_3"_ 

# %%
