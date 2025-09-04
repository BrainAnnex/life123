# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# ## One-bin reaction `2A + 5B <-> 4C + 3D`, with 1st-order kinetics for each species, taken to equilibrium
#
# Diffusion not applicable (just 1 bin)

# %% [markdown]
# ### TAGS :  "reactions 1D", "under-the-hood"

# %%
LAST_REVISED = "June 6, 2025"
LIFE123_VERSION = "1.0.0rc6"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename

from life123 import BioSim1D, ChemData, check_version

from life123 import HtmlLog as log
from life123 import GraphicLog

# %%
check_version(LIFE123_VERSION)

# %%

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system; NOTE: Diffusion not applicable (just 1 bin)
chem_data = ChemData(names=["A", "B", "C", "D"], plot_colors=['navy', 'cyan', 'red', 'orange'])    

bio = BioSim1D(n_bins=1, chem_data=chem_data)

bio.set_all_uniform_concentrations( [4., 7., 5., 2.] )

bio.describe_state()

# %%
# Specify the reaction
reactions = bio.get_reactions()

# Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
reactions.add_reaction(reactants=[(2,"A",1) , (5,"B",1)], products=[(4,"C",1) , (3,"D",1)],
                       forward_rate=5., reverse_rate=2.)

reactions.describe_reactions()

# %%
# Send a header and a plot to the HTML log file
log.write("Reaction 2 A + 5 B <-> 4 C + 3 D",
          style=log.h2)
reactions.plot_reaction_network("vue_cytoscape_2")

# %%

# %%
# Let's enable history - by default for all chemicals and all bins
bio.enable_history(take_snapshot=True, caption="Initial state")

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown]
# ### First step

# %%
bio.get_reaction_handler().enable_diagnostics()   # To save diagnostic information for the simulation run, below

# %%
# First step
bio.react(time_step=0.001, n_steps=1)
bio.describe_state()

# %% [markdown]
# _Early in the reaction :_
# [A] = 3.76 ,  [B] = 6.4 ,  [C] = 5.48 ,  [D] = 2.36

# %%
bio.get_bin_history(bin_address=0)

# %%

# %% [markdown]
# ### <a name="sec_2"></a>Numerous more steps

# %%
# Numerous more steps
bio.react(time_step=0.001, n_steps=40)

bio.describe_state()

# %% [markdown]
# ### <a name="sec_2_equilibrium"></a>Equilibrium

# %% [markdown]
# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:  
# [A] = 2.80284552 , [B] = 4.00711381 , [C] = 7.39430896 , [D] = 3.79573172

# %%
# Verify that the reaction has reached equilibrium
bio.get_reaction_handler().is_in_equilibrium(rxn_index=0, conc=bio.bin_snapshot(bin_address = 0))

# %%
df = bio.get_bin_history(bin_address=0)
df

# %% [markdown]
# A and B get depleted, while C and D get produced.
#
# **2A + 5B <-> 4C + 3D**

# %% [markdown]
# #### Let's verify that the stoichiometry is being respected

# %%
# We'll check the first two arrays of concentrations, from the run's history
arr0 = bio.get_reaction_handler().get_historical_concentrations(row=0, df=df)
arr1 = bio.get_reaction_handler().get_historical_concentrations(row=1, df=df)
arr0, arr1

# %%
# Check that the changes in the first reaction step conform to the stoichiometry
bio.get_reaction_handler().get_diagnostics().stoichiometry_checker(rxn_index=0, 
                                                                   conc_arr_before = arr0, 
                                                                   conc_arr_after = arr1)

# %%
arr1 - arr0

# %% [markdown]
# Indeed, the change in [A] is -2 x 0.12, and the change in [B] is -5 X 0.12,  
#   while the change in [C] is  4 x 0.12, and the change in [D] is  3 X 0.12

# %%
(arr1 - arr0) / 0.12    # 0.12 is the "moles of reactions" change

# %%

# %% [markdown]
# # Plots of changes of concentration with time

# %%
bio.plot_history_single_bin(bin_address=0, 
                            title_prefix="Reaction `2A + 5B <-> 4C + 3D`")

# %%
