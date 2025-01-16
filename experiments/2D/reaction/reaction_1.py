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
# ## `A <-> B` reaction, with 1st-order kinetics in both directions,
# ### taken to equilibrium
#
# Diffusion NOT done
# %% [markdown]
# ### TAGS :  "reactions 2D"

# %%
LAST_REVISED = "Dec. 16, 2024"
LIFE123_VERSION = "1.0-rc.1"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename

from life123 import UniformCompartment, BioSim2D, GraphicLog

# %%
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file
# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system.  NOTE: Diffusion not done
uc = UniformCompartment(names=["A", "B"])

# Reaction A <-> B , with 1st-order kinetics in both directions
uc.add_reaction(reactants="A", products="B", forward_rate=3., reverse_rate=2.)

uc.describe_reactions()

# %%
# Send the plot to the HTML log file
uc.plot_reaction_network("vue_cytoscape_2")

# %%
bio = BioSim2D(x_bins=3, y_bins=4, reaction_handler=uc)

bio.set_bin_conc_all_species(bin_x=0, bin_y=0, conc_list=[10.,50.])
bio.set_bin_conc_all_species(bin_x=0, bin_y=1, conc_list=[20.,35.])
bio.set_bin_conc_all_species(bin_x=2, bin_y=3, conc_list=[5.,100.])

bio.describe_state()

# %%

# %% [markdown] tags=[]
# ## First step

# %%
# First step (NOTE: here we're using a lower-level function that doesn't update the system state;
#                   it only computes the delta_reactions array)
bio.reaction_step(delta_time=0.1)
print("bio.delta_reactions:\n", bio.delta_reactions)

# %%
bio.system += bio.delta_reactions       # Matrix operation to update all the concentrations
bio.system_time += 0.1

bio.describe_state()

# %% [markdown] tags=[]
# ## Second step

# %%
# NOTE: now, we're using a highel-level function that also updates the system state
bio.react(time_step=0.1, n_steps=1)
bio.describe_state()

# %% [markdown]
# ## Many more steps, to equilibrium

# %%
bio.react(time_step=0.1, n_steps=8)
bio.describe_state()

# %%
bio.react(time_step=0.1, n_steps=10)
bio.describe_state()

# %% [markdown]
# ### The system has now reached equilibrium
# ### in individual bins, which remain separate because we're NOT doing diffusion in this experiment

# %% [markdown]
# Verify the equilibrium in each of the active bins

# %%
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc={"A": 23.99998665, "B": 36.00001335})

# %%
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc={"A": 21.99999809, "B": 33.00000191})

# %%
bio.reaction_dynamics.is_in_equilibrium(rxn_index=0, conc={"A": 41.99996471, "B": 63.00003529}, explain=False)

# %%
