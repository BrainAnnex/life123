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
# ## `A` down-regulates `B` , 
# ### by being the *limiting reagent* in reaction `A + 2 B <-> Y` (mostly forward)
# 1st-order kinetics.   
# If [A] is low and [B] is high, then [B] remains high.  If [A] goes high, [B] goes low.  However, at that point, A can no longer bring B up to any substantial extent.
#
# See also 1D/reactions/down_regulation_1
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta36)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData as chem
from life123 import UniformCompartment

from life123 import GraphicLog

# %% tags=[]
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# ### Initialize the system

# %%
# Initialize the system
chem_data = chem(names=["A", "B", "Y"])

# Reaction A + 2 B <-> Y , with 1st-order kinetics for all species
chem_data.add_reaction(reactants=["A" , (2, "B", 1)], products="Y",
                       forward_rate=8., reverse_rate=2.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = UniformCompartment(chem_data=chem_data, preset="fast")
dynamics.set_conc(conc={"A": 5., "B": 100.},
                  snapshot=True)      # A is scarce, B is plentiful, Y is absent
dynamics.describe_state()

# %%

# %%

# %% [markdown] tags=[]
# # 1. Take the initial system to equilibrium

# %%
# All of these settings are currently close to the default values... but subject to change; set for repeatability

dynamics.single_compartment_react(initial_step=0.0005, duration=0.015,
                                  variable_steps=True)

# %% [markdown]
# A, as the scarse limiting reagent, stops the reaction.  
# When A is low, B is also low.

# %%
dynamics.plot_history(colors=['red', 'darkorange', 'green'],
                      title="Changes in concentrations (reaction A + 2 B <-> Y)")

# %%
dynamics.get_history()

# %%
dynamics.explain_time_advance(sys_history=dynamics.get_history())

# %% [markdown]
# #### Equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%

# %%

# %% [markdown] tags=[]
# # 2. Now, let's suddenly increase [A]

# %%
dynamics.set_single_conc(species_name="A", conc=40., snapshot=True)
dynamics.describe_state()

# %%
dynamics.history.get_dataframe(tail=5)

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
dynamics.single_compartment_react(initial_step=0.0005, target_end_time=0.055,
                                  variable_steps=True)

# %%
dynamics.plot_history(colors=['red', 'darkorange', 'green'],
                      title="Changes in concentrations (reaction A + 2 B <-> Y)")

# %% [markdown]
# **A**, still the limiting reagent, is again stopping the reaction.  
# The (transiently) high value of [A] led to a high value of [B]

# %%
dynamics.get_history()

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%

# %%

# %% [markdown] tags=[]
# # 3. Let's again suddenly increase [A]

# %%
dynamics.set_single_conc(species_name="A", conc=30., snapshot=True)
dynamics.describe_state()

# %%
dynamics.history.get_dataframe(tail=5)

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
dynamics.single_compartment_react(initial_step=0.001, target_end_time=0.09,
                                  variable_steps=True)

# %%
dynamics.plot_history(colors=['red', 'darkorange', 'green'],
                      title="Changes in concentrations (reaction A + 2 B <-> Y)")

# %%
dynamics.get_history()

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium(tolerance=2)

# %% [markdown]
# **A**, again the scarse limiting reagent, stops the reaction yet again

# %%

# %%

# %% [markdown]
# # 4. A can down-regulate B, but it cannot bring it up to any significant amount
# #### Even if A is completely taken out (i.e., [A] set to 0), [B] can only slightly increase, from the reverse reaction ("Le Chatelier's principle".)  
# Let's try it:

# %%
dynamics.set_single_conc(species_name="A", conc=0., snapshot=True)   # Completely eliminate A
dynamics.describe_state()

# %%
dynamics.single_compartment_react(initial_step=0.001, target_end_time=0.16,
                                  variable_steps=True)


# %%
dynamics.plot_history(colors=['red', 'darkorange', 'green'],
                      title="Changes in concentrations (reaction A + 2 B <-> Y)")

# %%
dynamics.get_history()

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown]
# #### As expected, even the complete withdrawal of A (red), brings about only a modest increase of B's concentration, from the reverse reaction (i.e. [B] slightly increases at the expense of [Y].)  
# #### The change is modest because our  reaction A + 2 B <-> Y is mostly in the forward direction (K = 4)
# *Le Chatelier's principle* in action: "A change in one of the variables that describe a system at equilibrium produces a shift in the position of the equilibrium that counteracts the effect of this change."

# %%
