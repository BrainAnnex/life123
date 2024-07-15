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
# ## Coupled pair of reactions: `A <-> B` , and  `A + E <-> B + E`
# A direct reaction and the same reaction, catalyzed
# ### Enzyme `E` initially zero, and then suddenly added mid-reaction
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta36)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from life123 import ChemData
from life123 import UniformCompartment

# %%
# Initialize the system
chem_data = ChemData(names=["A", "B", "E"])

# Reaction A <-> B , with 1st-order kinetics, favorable thermodynamics in the forward direction, 
# and a forward rate that is much slower than it would be with the enzyme - as seen in the next reaction, below
chem_data.add_reaction(reactants="A", products="B",
                       forward_rate=1., delta_G=-3989.73)

# Reaction A + E <-> B + E , with 1st-order kinetics, and a forward rate that is faster than it was without the enzyme
# Thermodynamically, there's no change from the reaction without the enzyme
chem_data.add_reaction(reactants=["A", "E"], products=["B", "E"],
                       forward_rate=10., delta_G=-3989.73)

chem_data.describe_reactions()     # Notice how the enzyme `E` is noted in the printout below

# %% [markdown]
# ### Set the initial concentrations of all the chemicals - starting with no enzyme

# %%
dynamics = UniformCompartment(chem_data=chem_data, preset="mid")
dynamics.set_conc(conc={"A": 20.},
                  snapshot=True)      # Initially, no enzyme `E`
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Advance the reactions (for now without enzyme), but stopping well before equilibrium

# %%
dynamics.enable_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# Perform the reactions
dynamics.single_compartment_react(duration=0.25,
                                  initial_step=0.05, variable_steps=True)

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="WITH zero enzyme")

# %% [markdown]
# ### The reactions, lacking enzyme, are proceeding slowly towards equilibrium, just like the reaction that was discussed in part 1 of the experiment "enzyme_1"

# %% [markdown]
# # Now suddently add a lot of enzyme

# %%
dynamics.set_single_conc(30., species_name="E", snapshot=True)    # Plenty of enzyme `E`

# %%
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Take the system to equilibrium

# %%
dynamics.single_compartment_react(duration=0.04,
                                  initial_step=0.005, variable_steps=True)

# %% [markdown]
# #### Note how the (proposed) initial step - in spite of having been reduced from the earlier run - is now appearing _large_, given the much-faster reaction dynamics.  However, the variable-step engine intercepts and automatically corrects the problem!

# %%
#dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=['darkorange', 'green', 'violet'], show_intervals=True, title_prefix="WITH enzyme added mid-reaction")

# %% [markdown]
# ## Notice the dramatic acceleration of the reaction as soon as the enzyme `E` is added at t = 0.275!
# The reactions simulator automatically switches to small time steps is in order to handle the rapid amount of change

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%
dynamics.get_history()

# %%
