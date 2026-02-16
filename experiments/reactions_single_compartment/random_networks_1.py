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
# ### Random Networks
#
# Early implementation

# %% [markdown]
# ### TAGS :   "uniform compartment"

# %%
LAST_REVISED = "Feb. 15, 2026"
LIFE123_VERSION = "1.0.0rc7"         # Library version this experiment is based on

# %%
#import set_path            # Using MyBinder?  Uncomment this before running the next cell!
                            # Importing this module will add the project's home directory to sys.path

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

import ipynbname

from life123 import check_version, RandomReactionNetwork, UniformCompartment

# %%
check_version(LIFE123_VERSION)    # To check compatibility

# %%
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file
                                            # IN CASE OF PROBLEMS, set manually to any desired name

# %%

# %%

# %% [markdown]
# ## Set up Random Network (no thermodynamic nor kinetic data)

# %%
net = RandomReactionNetwork(n_chems=4, n_rxns=6)

# %%
net.chem_data.get_all_labels()

# %% [markdown]
# ### The above chemical labels were auto-generated

# %%
net.registry.number_of_reactions()

# %%
for i in range(net.registry.number_of_reactions()):
    rxn = net.registry.get_reaction(i)
    print(f"({i}) {rxn.describe(concise=False)}")

# %%

# %%
# Instantiate the simulator with our chemicals and reactions
uc = UniformCompartment(reactions=net.registry)

# %%

# %%
# Send a plot of the network of reactions to the HTML log file
uc.plot_reaction_network(log_file=log_file)
