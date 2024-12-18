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
# ## 2 COUPLED reactions:  `A + B <-> C`  and  `C + D <-> E` , 
# ### with 1st-order kinetics for each species, taken to equilibrium
#
# Both reactions are stronger in their respective forward rates.  For the most part, "C" is produced by the 1st reaction, and consumed by the 2nd one
#
# Diffusion not applicable (just 1 bin)
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

# %%
# Initialize the system
chem_data = chem(names=["A", "B", "C", "D", "E"])     # NOTE: Diffusion not applicable (just 1 bin)

# Specify the reactions


# Reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species
chem_data.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
chem_data.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_all_uniform_concentrations( [3., 5., 1., 0.4, 0.1] )

bio.describe_state()

# %%
# Save the state of the concentrations of all species at bin 0
bio.add_snapshot(bio.bin_snapshot(bin_address = 0), caption="Initial state")
bio.get_history()

# %%
chem_data.describe_reactions()

# %%
# Send a header and a plot to the HTML log file
log.write("2 COUPLED reactions:  A + B <-> C  and  C + D <-> E",
          style=log.h2)
# Send the plot to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %% [markdown] tags=[]
# ### <a name="sec_2_first_step"></a>First step

# %%
# First step
bio.react(time_step=0.01, n_steps=1, snapshots={"sample_bin": 0})
bio.describe_state()

# %% [markdown]
# 1 bins and 5 species:
#  [[2.27 ]
#  [4.27 ]
#  [1.702]
#  [0.372]
#  [0.128]]

# %%
bio.get_history()

# %%
# Identical 2nd step
bio.react(time_step=0.01, n_steps=1, snapshots={"sample_bin": 0})
bio.describe_state()

# %% [markdown]
# 1 bins and 5 species:
#  [[1.819395  ]
#  [3.819395  ]
#  [2.10707348]
#  [0.32646848]
#  [0.17353152]]

# %%
bio.get_history()

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more identical steps, to equilibrium
bio.react(time_step=0.01, n_steps=200, snapshots={"sample_bin": 0, "frequency": 10})
bio.describe_state()

# %% [markdown]
# 1 bins and 5 species:
#  [[0.50508029]
#  [2.50508029]
#  [3.16316668]
#  [0.06824696]
#  [0.43175304]]

# %% [markdown] tags=[]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %%
# Verify that each reaction has reached equilibrium
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
bio.reaction_dynamics.is_in_equilibrium(rxn_index=1, conc=bio.bin_snapshot(bin_address = 0))

# %%
# Do a consistent check with the equilibrium concentrations:

A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
D_eq = bio.bin_concentration(0, 3)
E_eq = bio.bin_concentration(0, 4)

Rf0 = chem_data.get_forward_rate(0)
Rb0 = chem_data.get_reverse_rate(0)

Rf1 = chem_data.get_forward_rate(1)
Rb1 = chem_data.get_reverse_rate(1)

equil = -(Rf0 * A_eq * B_eq - Rf1 * C_eq * D_eq) + (Rb0 * C_eq - Rb1 * E_eq)

print("\nAt equilibrium: ", equil, " (this should be close to 0 at equilibrium)")

# %%
bio.get_history()

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %% tags=[]
fig = px.line(data_frame=bio.get_history(), x="SYSTEM TIME", y=["A", "B", "C", "D", "E"], 
              title="2 COUPLED reactions:  A + B <-> C  and  C + D <-> E . Changes in concentrations",
              color_discrete_sequence = ['navy', 'cyan', 'red', 'orange', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# A and B get consumed.  
# C gets produced by the 1st reaction more quickly than consumed by the 2nd one.
# D gets consumed, while E gets produced.

# %%
