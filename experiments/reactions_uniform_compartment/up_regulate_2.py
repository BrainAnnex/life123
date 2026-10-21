# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# ## `U` ("Up-regulator") up-regulates `X` , by sharing a reaction product `D` ("Drain") across 2 separate reactions:   
# ### `U <-> 2 D` and `X <-> D` (both mostly forward)
#
# 1st-order kinetics throughout; i.e. `U <-> 2 D` is a HYPOTHETICAL reaction with first-order kinetic in all species.  
# This experiment also highlights how to use **custom** kinetic reaction models, since one of our reactions doesn't follow mass-action kinetics.
#
# Invoking [Le Chatelier's principle](https://www.chemguide.co.uk/physical/equilibria/lechatelier.html), it can be seen that, starting from equilibrium, when [U] goes up, so does [D]; and when [D] goes up, so does [X].   
# Conversely, when [U] goes down, so does [D]; and when [D] goes down, so does [X].     

# %%
LAST_REVISED = "Oct. 20, 2026"
LIFE123_VERSION = "1.0.0rc9"     # Library version this experiment is based on

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

from life123 import check_version, ReactionSimulator
from life123.kinetic_models import Custom_Model

# %%
check_version(LIFE123_VERSION)    # To check compatibility

# %%
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file
                                            # IN CASE OF PROBLEMS, set manually to any desired name
log_file

# %%

# %% [markdown]
# ### Initialize the system

# %%
# Instantiate the simulator and specify the reaction
sim = ReactionSimulator(preset="fast", 
                        enable_diagnostics=True) 

# (HYPOTHETICAL!) Reaction U <-> 2D , with custom kinetics: 1st-order for all species
sim.add_reaction(reactants="U", products=[(2, "D")],
                 reaction_model="custom",
                 kinetic_parameters={"kF": 8., "kR": 2.,
                                     "rate_function": Custom_Model.kinetic_rate_first_order})

# Elementary Reaction X <-> D
sim.add_reaction(reactants="X", products="D",
                 reaction_model="mass action",
                 kinetic_parameters={"kF": 6., "kR": 3.})

# %%
sim.describe_reactions(detail="summary")

# %%
# Send a plot of the network of reactions to the HTML log file
sim.plot_reaction_network(log_file=log_file)

# %%
IFrame(log_file, width=1200, height=800)         # You may also open the above log file in a browser

# %%

# %%

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
sim.set_conc(conc={"U": 50., "X": 100.})
sim.describe_state()

# %%

# %% [markdown]
# # 1. Take the initial system to equilibrium

# %%
sim.single_compartment_react(initial_step=0.03, target_end_time=0.5)

# %%
sim.get_history()

# %%
sim.diagnostics.conc_data.explain_time_advance()

# %% [markdown]
# ## Plots of changes of concentration with time

# %%
# Assign specific colors to each specie
sim.species_registry.update(species_id="U", plot_color="red")
sim.species_registry.update(species_id="D", plot_color="gray")
sim.species_registry.update(species_id="X", plot_color="green")

# %%
sim.plot_history()

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %%
# Verify that the reaction has reached equilibrium
sim.is_in_equilibrium(rxn_index=1)   # Only 'mass action' reaction models are currently supported

# %%

# %%

# %% [markdown]
# # 2. Now, let's suddenly increase [U]

# %%
sim.set_single_conc(species_id="U", conc=70.)
sim.describe_state()

# %%
sim.get_history(tail=3)

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
sim.single_compartment_react(initial_step=0.03, target_end_time=1)

# %%
#sim.get_history()
#sim.explain_time_advance()

# %%
sim.plot_history()  # The colors set earlier will be used

# %% [markdown]
# ### The (transiently) high value of [U] led to an increase in [X]

# %%
# Verify that the reaction has reached equilibrium
sim.is_in_equilibrium(rxn_index=1, tolerance=2)

# %%

# %%

# %% [markdown]
# # 3. Let's again suddenly increase [U]

# %%
sim.set_single_conc(species_id="U", conc=100.)
sim.describe_state()

# %%
sim.get_history(tail=3)

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
sim.single_compartment_react(initial_step=0.03, target_end_time=1.6)

# %%
#sim.get_history()
#sim.explain_time_advance()

# %%
sim.plot_history()

# %% [markdown]
# ### The (transiently) high value of [U] again led to an increase in [X]

# %%
# Verify that the reaction has reached equilibrium
sim.is_in_equilibrium(rxn_index=1)       # Only 'mass action' reaction models are currently supported

# %%

# %%

# %% [markdown]
# # 4. Now, instead, let's DECREASE [U]

# %%
sim.set_single_conc(species_id="U", conc=5.)
sim.describe_state()

# %%
sim.get_history(tail=3)

# %% [markdown]
# ### Take the system to equilibrium

# %%
sim.single_compartment_react(initial_step=0.03, target_end_time=2.3)

# %%
#sim.get_history()
#sim.explain_time_advance()

# %%
sim.plot_history()

# %% [markdown]
# ### The (transiently) LOW value of [U] led to an DECREASE in [X]

# %%
# Verify that the reaction has reached equilibrium
sim.is_in_equilibrium(rxn_index=1, explain=False)   # Only 'mass action' reaction models are currently supported

# %%
