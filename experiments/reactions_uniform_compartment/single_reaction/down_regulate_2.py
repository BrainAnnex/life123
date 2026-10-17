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
# ## `A` down-regulates `B` , 
# ### by being the *limiting reagent* in a HYPOTHETICAL reaction `A + 2 B <-> Y` with first-order kinetic in all species (mostly forward)   
# If [A] is low and [B] is high, then [B] remains high.  If [A] goes high, [B] goes low.  However, at that point, A can no longer bring B up to any substantial extent.
#
# This experiment also highlights how to use **custom** kinetic reaction models, since our hypothetical reaction doesn't follow mass-action kinetics.

# %% [markdown]
# ### TAGS : "uniform compartment", "chemistry"

# %%
LAST_REVISED = "Oct. 16, 2026"
LIFE123_VERSION = "1.0.0rc9"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import check_version, ReactionSimulator, ReactionKinetics, PlotlyHelper
from life123.kinetic_models import Custom_Model

# %%
#from experiments.get_notebook_info import get_notebook_basename

#from life123 import SpeciesRegistry as chem
#from life123 import UniformCompartment

#from life123 import GraphicLog

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ### Initialize the system

# %%
# Instantiate the simulator and specify the reaction
uc = ReactionSimulator(preset="fast", 
                       enable_diagnostics=True)     # We call this "uc" because it's a "uniform compartment"

# Reaction A + 2 B <-> Y , with custom kinetics: 1st-order for all species
uc.add_reaction(reactants=["A" , (2, "B")], products="Y",
                reaction_model="custom",
                kinetic_parameters={"kF": 8., "kR": 2.,
                                    "rate_function": Custom_Model.kinetic_rate_first_order})

# %%
uc.describe_reactions()

# %%
# Set initial concentrations of all the chemicals
uc.set_conc({"A": 5., "B": 100.})   # A is scarce, B is plentiful, Y is absent
uc.describe_state()

# %%

# %%

# %% [markdown]
# # 1. Take the initial system to equilibrium

# %%
uc.single_compartment_react(initial_step=0.0005, duration=0.015)

# %% [markdown]
# A, as the scarse limiting reagent, stops the reaction.  
# When A is low, B is also low.

# %%

# %% [markdown]
# ### Visualize

# %%
# Assign specific colors to each specie
uc.species_registry.update(species_id="A", plot_color="red")
uc.species_registry.update(species_id="B", plot_color="darkorange")
uc.species_registry.update(species_id="Y", plot_color="green")

# %%
uc.species_registry.as_dataframe()

# %%
uc.plot_history(title="Changes in concentrations (reaction A + 2 B <-> Y)")

# %%
uc.get_history()

# %%
uc.diagnostics.conc_data.explain_time_advance()

# %% [markdown]
# #### Equilibrium

# %%

# %%

# %% [markdown]
# # 2. Now, let's suddenly increase [A]

# %%
uc.set_single_conc(species_id="A", conc=40.)
uc.describe_state()

# %%
uc.get_history(tail=5)

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
uc.single_compartment_react(initial_step=0.0005, target_end_time=0.055)

# %%
uc.plot_history(title="Changes in concentrations (reaction A + 2 B <-> Y)")  # The colors set earlier will be used

# %% [markdown]
# **A**, still the limiting reagent, is again stopping the reaction.  
# The (transiently) high value of [A] led to a high value of [B]

# %%
uc.get_history()

# %%

# %%

# %% [markdown]
# # 3. Let's again suddenly increase [A]

# %%
uc.set_single_conc(species_id="A", conc=30.)
uc.describe_state()

# %%
uc.get_history(tail=5)

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
uc.single_compartment_react(initial_step=0.001, target_end_time=0.09)

# %%
uc.plot_history(title="Changes in concentrations (reaction A + 2 B <-> Y)")

# %%
uc.get_history()

# %% [markdown]
# **A**, again the scarse limiting reagent, stops the reaction yet again

# %%

# %%

# %% [markdown]
# # 4. A can down-regulate B, but it cannot bring it up to any significant amount
# #### Even if A is completely taken out (i.e., [A] set to 0), [B] can only slightly increase, from the reverse reaction ("Le Chatelier's principle".)  
# Let's try it:

# %%
uc.set_single_conc(species_id="A", conc=0., snapshot=True)   # Completely eliminate A
uc.describe_state()

# %%
uc.single_compartment_react(initial_step=0.001, target_end_time=0.16,
                                  variable_steps=True)


# %%
uc.plot_history(title="Changes in concentrations (reaction A + 2 B <-> Y)")

# %%
uc.get_history()

# %% [markdown]
# #### As expected, even the complete withdrawal of A (red), brings about only a modest increase of B's concentration, from the reverse reaction (i.e. [B] slightly increases at the expense of [Y].)  
# #### The change is modest because our  reaction A + 2 B <-> Y is mostly in the forward direction (K = 4)
# *Le Chatelier's principle* in action: "A change in one of the variables that describe a system at equilibrium produces a shift in the position of the equilibrium that counteracts the effect of this change."

# %%
