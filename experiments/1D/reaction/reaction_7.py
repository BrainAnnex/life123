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
# ## One-bin  `2A <-> B` reaction,  
# ### COMPARING 1st-order and 2nd-order kinetics in *forward* direction; reverse direction 1-st order
#
# Diffusion not applicable (just 1 bin)
#
# See also the experiment _"reactions_single_compartment/react_4"_ 
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta34.1)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData as chem
from life123 import BioSim1D

import plotly.express as px
from life123 import HtmlLog as log
from life123 import GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # INITIALLY, with 1st-order kinetics in both directions

# %%
# Initialize the system
chem_data = chem(names=["A", "B"])     # NOTE: Diffusion not applicable (just 1 bin)



# Reaction  2A <-> B , FOR NOW with 1st-order kinetics in both directions
chem_data.add_reaction(reactants=[(2, "A", 1)], products=["B"], forward_rate=5., reverse_rate=2.)

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_all_uniform_concentrations( [3., 5.] )

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0), caption="Initial state")
bio.get_history()

# %%
chem_data.describe_reactions()

# %%
# Send a header and a plot to the HTML log file
log.write("Reaction 2A <-> B is 1st order in all species:",
          style=log.h2)

chem_data.plot_reaction_network("vue_cytoscape_2")

# %%
# First step
bio.react(time_step=0.02, n_steps=1, snapshots={"sample_bin": 0})
bio.describe_state()

# %% [markdown]
# Small conc. changes so far:  [A] = 2.8 , [B] = 5.1

# %%
bio.get_history()

# %%
# Numerous more steps, to equilibrium
bio.react(time_step=0.02, n_steps=20, snapshots={"sample_bin": 0})

bio.describe_state()

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the *1st order* reactions),
# the systems settles in the following equilibrium:  
# [A] = 2.16928427 , [B] = 5.41535786

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
df = bio.get_history()
df

# %% tags=[]
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="2A <-> B : changes in concentrations with time",
              color_discrete_sequence = ['navy', 'orange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# A gets depleted, while B gets produced.

# %% [markdown]
# #### Let's verify that the stoichiometry is being respected

# %%
# We'll check the first two arrays of concentrations, from the run's history
arr0 = bio.reaction_dynamics.get_historical_concentrations(row=0, df=df)
arr1 = bio.reaction_dynamics.get_historical_concentrations(row=1, df=df)
arr0, arr1

# %%
bio.reaction_dynamics.stoichiometry_checker(rxn_index=0, 
                               conc_arr_before = arr0, 
                               conc_arr_after = arr1)

# %% [markdown]
# # STARTING OVER, this time with 2nd-order kinetics in the forward reaction

# %%
bio.reaction_dynamics.clear_reactions()

# %%
# Reaction  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction
chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)

# %%
# RESET the concentrations to their original values
bio.set_all_uniform_concentrations( [3., 5.] )

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0),
                 caption = "RESET all concentrations to initial values")
bio.get_history()

# %%
chem_data.describe_reactions()

# %%
# Send a header and a plot to the HTML log file
log.write("Reaction 2A <-> B is 2nd order in A, and 1st order in B:",
          style=log.h2)
chem_data.plot_reaction_network("vue_cytoscape_2")

# %%
# First step
bio.react(time_step=0.02, n_steps=1, snapshots={"sample_bin": 0})
bio.describe_state()

# %% [markdown]
# [A] = 1.6 , [B] = 5.7
# _(Contrast with the counterpart in the 1st order kinetics:  [A] = 2.8 , [B] = 5.1)_

# %%
bio.get_history()

# %%
# Numerous more steps
bio.react(time_step=0.02, n_steps=20, snapshots={"sample_bin": 0})

bio.describe_state()

# %% [markdown]
# The systems settles in the following equilibrium:  [A] = 1.51554944 , [B] = 5.74222528

# %%
# Verify that the reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
df2 = bio.get_history()
df2

# %%
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B"], 
              title="2A <-> B : changes in concentrations (the jump at 0.42 is the concentration reset)",
              color_discrete_sequence = ['navy', 'orange'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# **Compared to first-order kinetics in A**, the (2nd order in A) reaction now takes place much more quickly, and proceeds to almost complete depletion of A

# %% [markdown]
# #### Let's verify that the stoichiometry is still being respected

# %%
# We'll check the first two arrays of concentrations, from the run's history
arr0 = bio.reaction_dynamics.get_historical_concentrations(row=22, df=df2)  # Row 22 is the conc. reset
arr1 = bio.reaction_dynamics.get_historical_concentrations(row=23, df=df2)
arr0, arr1

# %%
bio.reaction_dynamics.stoichiometry_checker(rxn_index=0, 
                               conc_arr_before = arr0, 
                               conc_arr_after = arr1)

# %%
