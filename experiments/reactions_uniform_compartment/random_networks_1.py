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
# ## Random Reaction Networks
#
# #### Early implementation of random networks of reactions that are **thermodynamically possible**,
# #### and **biologically plausible**

# %% [markdown]
# ### TAGS :   "uniform compartment"

# %%
LAST_REVISED = "Oct. 19, 2026"
LIFE123_VERSION = "1.0.0rc9"         # Library version this experiment is based on

# %%
#import set_path            # Using MyBinder?  Uncomment this before running the next cell!
                            # Importing this module will add the project's home directory to sys.path

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

import ipynbname
from IPython.display import IFrame

from life123 import check_version, RandomReactionNetwork, ReactionSimulator

# %%
check_version(LIFE123_VERSION)    # To check compatibility

# %%
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file
                                            # IN CASE OF PROBLEMS, set manually to any desired name
log_file

# %%

# %%

# %% [markdown]
# ## Set up Random Network

# %%
net = RandomReactionNetwork(n_species=4, n_rxns=7, seed=47242315, thermodynamic_ruggedness=0.2)

# %%
net.species_data.get_all_species_ids()

# %% [markdown]
# #### The above species id's were auto-generated

# %%

# %%
rxns = net.get_reaction_data()          # Object of type "ReactionRegistry"

# %%
rxns.number_of_reactions()

# %%
rxns.describe_reactions(detail="minimal")

# %%
rxns.describe_reactions(detail="summary")

# %% [markdown]
# ## NOTE ON UNITS: _delta_H_ and _delta_G_ are in **kJ/mol**  ;  _delta_S_ are in **J/(mol·K)**

# %%

# %% [markdown]
# #### Investigate the random values that were assigned to the reactions (in a manner to be consistent with thermodynamics, for example across groups of reactions):

# %%
H_list = [rxn.source_object.thermodynamics.delta_H
          for rxn in net.get_all_random_reactions()]    # rxn is an object of type "SimulationReaction"
H_list   # kJ/mol

# %%
S_list = [rxn.source_object.thermodynamics.delta_S
          for rxn in net.get_all_random_reactions()]    # rxn is an object of type "SimulationReaction"
S_list   # J/(mol·K)

# %%
G_list = [rxn.source_object.thermodynamics.delta_G
          for rxn in net.get_all_random_reactions()]    # rxn is an object of type "SimulationReaction"
G_list   # kJ/mol

# %%
# Verify G_list
[H_list[i] - (25+273.15) * S_list[i] / 1000 for i in range(len(H_list))]

# %%

# %%

# %%
# Instantiate the simulator with our species and reactions
sim = ReactionSimulator(reaction_registry=rxns)

# %%

# %%
# Send a plot of the network of reactions to the HTML log file
sim.plot_reaction_network(log_file=log_file)

# %%
IFrame(log_file, width=1200, height=700)         # You may also open the above log file in a browser

# %%
