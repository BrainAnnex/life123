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
# ### A up-regulates B , by being *the limiting reagent* in reaction A + X <-> 2B (mostly forward), where X is plentiful
# 1st-order kinetics.   
# If [A] is low, [B] remains low, too.  Then, if [A] goes high, then so does [B].  However, at that point, A can no longer bring B down to any substantial extent.
#
# Single-bin reaction
#
# LAST REVISED: Aug. 22, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

import plotly.express as px
from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system
chem_data = chem(names=["A", "X", "B"])     # NOTE: Diffusion not applicable (just 1 bin)

rxn = Reactions(chem_data)

# Reaction A + X <-> 2B , with 1st-order kinetics for all species
rxn.add_reaction(reactants=[("A") , ("X")], products=[(2, "B")],
                 forward_rate=8., reverse_rate=2.)

rxn.describe_reactions()

# Send the plot of the reaction network to the HTML log file
graph_data = rxn.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %%
bio.initialize_system(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_uniform_concentration(species_name="A", conc=5.)     # Scarce
bio.set_uniform_concentration(species_name="X", conc=100.)   # Plentiful
# Initially, no "B" is present

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0 (the only bin in this system)
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
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

# %% [markdown]
# Consistent with the 4/1 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:

# %%
# Verify that the reaction has reached equilibrium
rxn.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "X", "B"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['red', 'darkorange', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown] tags=[]
# ## Now, let's suddenly increase [A]

# %%
bio.set_bin_conc(bin=0, species_index=0, conc=50.)
bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0 (the only bin in this system)
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))
bio.get_history()

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=40, snapshots={"frequency": 2, "sample_bin": 0})  # At every other step, take a snapshot 
                                                                                      # of all species at bin 0
bio.describe_state()
bio.get_history()

# %% [markdown]
# A, still the limiting reagent, is again stopping the reaction.  
# The (transiently) high value of [A] led to a high value of [B]

# %%
# Verify the equilibrium
A_eq = bio.bin_concentration(0, 0)
X_eq = bio.bin_concentration(0, 1)
B_eq = bio.bin_concentration(0, 2)
print("Ratio of equilibrium concentrations (B_eq / (A_eq * X_eq)): ", (B_eq / (A_eq * X_eq)))
print("Ratio of forward/reverse rates: ", rxn.get_forward_rate(0) / rxn.get_back_rate(0))

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "X", "B"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['red', 'darkorange', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown] tags=[]
# ## Let's again suddenly increase [A]

# %%
bio.set_bin_conc(bin=0, species_index=0, conc=30.)
bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0 (the only bin in this system)
bio.save_snapshot(bio.bin_snapshot(bin_address = 0))

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
bio.react(time_step=0.0005, n_steps=70, snapshots={"frequency": 2, "sample_bin": 0})  # At every other step, take a snapshot 
                                                                                      # of all species at bin 0
bio.describe_state()
bio.get_history()

# %% [markdown]
# A, again the scarse limiting reagent, stops the reaction yet again

# %%
# Verify the equilibrium
A_eq = bio.bin_concentration(0, 0)
X_eq = bio.bin_concentration(0, 1)
B_eq = bio.bin_concentration(0, 2)
print("Ratio of equilibrium concentrations (B_eq / (A_eq * X_eq)): ", (B_eq / (A_eq * X_eq)))
print("Ratio of forward/reverse rates: ", rxn.get_forward_rate(0) / rxn.get_back_rate(0))

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "X", "B"], 
              title="Changes in concentrations",
              color_discrete_sequence = ['red', 'darkorange', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# Note: A can up-regulate B, but it cannot bring it down.  
# X will soon need to be replenished, if A is to continue being the limiting reagent.

# %%