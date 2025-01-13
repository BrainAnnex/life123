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

# %%
LAST_REVISED = "Jan. 12, 2025"
LIFE123_VERSION = "1.0.0rc2"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData, UniformCompartment, BioSim1D

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
chem_data = ChemData(names=["A", "B", "C", "D", "E"])     # NOTE: Diffusion not applicable (just 1 bin)

uc = UniformCompartment(chem_data=chem_data)

# Specify the reactions

# Reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species
uc.add_reaction(reactants=["A", "B"], products="C", forward_rate=5., reverse_rate=2.)
uc.add_reaction(reactants=["C", "D"], products="E", forward_rate=8., reverse_rate=4.)
uc.describe_reactions()

# %%

# %%
# Set up the 1-D system
bio = BioSim1D(n_bins=1, reaction_handler=uc)

# %%

# %% [markdown]
# ## Enable History

# %%
# Let's enable history - by default for all chemicals and all bins (we just have 1 bin)
bio.enable_history()

# %%

# %% [markdown]
# ### Set the initial state

# %%
bio.set_all_uniform_concentrations( [3., 5., 1., 0.4, 0.1] )

bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %%
uc.describe_reactions()

# %%
# Send a header and a plot to the HTML log file
log.write("2 COUPLED reactions:  A + B <-> C  and  C + D <-> E",
          style=log.h2)

# Send the plot to the HTML log file
uc.plot_reaction_network("vue_cytoscape_2")

# %%

# %%

# %% [markdown] tags=[]
# ### Start the reaction

# %%
# First step
bio.react(time_step=0.01, n_steps=1)

# %%
bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %%

# %%
# Identical 2nd step
bio.react(time_step=0.01, n_steps=1)

# %%
bio.describe_state()

# %%
bio.get_bin_history(bin_address=0)

# %%

# %%

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more identical steps, to equilibrium
bio.react(time_step=0.01, n_steps=200)

# %%
bio.describe_state()

# %%

# %% [markdown] tags=[]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %%
# Verify that each reaction has reached equilibrium
uc.is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
uc.is_in_equilibrium(rxn_index=1, conc=bio.bin_snapshot(bin_address = 0))

# %%
# Do a consistent check with the equilibrium concentrations:

A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
D_eq = bio.bin_concentration(0, 3)
E_eq = bio.bin_concentration(0, 4)

Rf0 = uc.get_reactions().get_forward_rate(0)
Rb0 = uc.get_reactions().get_reverse_rate(0)

Rf1 = uc.get_reactions().get_forward_rate(1)
Rb1 = uc.get_reactions().get_reverse_rate(1)

equil = -(Rf0 * A_eq * B_eq - Rf1 * C_eq * D_eq) + (Rb0 * C_eq - Rb1 * E_eq)

print("\nAt equilibrium: ", equil, " (this should be close to 0 at equilibrium)")

# %%

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown] tags=[]
# # Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title="2 COUPLED reactions:  A + B <-> C  and  C + D <-> E . Concentrations at bin 0")

# %% [markdown]
# A and B get consumed.  
# C gets produced by the 1st reaction more quickly than consumed by the 2nd one.
# D gets consumed, while E gets produced.

# %%
