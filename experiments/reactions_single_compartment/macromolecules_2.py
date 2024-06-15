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
# ## Macromolecules : Binding Affinity and Fractional Occupancy. Regulation of the rates of reactions catalyzed by the macromolecule
#
# ### Reaction `A <-> B` catalyzed when ligand `L` binds to `site 1` of macromolecule `M1`
#
# # IN-PROGRESS
#
# In Part 1 we consider the un-catalyzed reaction `A <-> B`  
# In Part 2 we'll see what happens when catalysis is added
#
# LAST REVISED: June 14, 2024 (using v. 1.0 beta33)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from src.modules.chemicals.chem_data import ChemData
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import plotly.express as px

# %%

# %% [markdown]
# # PART 1 - Consider just the un-catalyzed reaction `A <-> B` by itself

# %%
# Initialize the system
chem1 = ChemData(names=["A", "B"])

# Reaction A <-> B , without catalysis (slow forward rate, relative to what we'll see in Part 2)
chem1.add_reaction(reactants="A", products="B",
                  forward_rate=1. , delta_G= -5000)

chem1.describe_reactions()


# %%
# Set the initial concentrations of all the chemicals
dynamics1 = ReactionDynamics(chem_data=chem1)

dynamics1.set_conc(conc={"A": 100., "B": 20.},
                  snapshot=True)

dynamics1.describe_state()

# %%
# Take the system to equilibrium
dynamics1.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics1.set_thresholds(norm="norm_A", low=0.6, high=1.0, abort=1.44)
dynamics1.set_thresholds(norm="norm_B", low=0.04, high=0.6, abort=1.5)
dynamics1.set_step_factors(upshift=1.2, downshift=0.7, abort=0.4)
dynamics1.set_error_step_factor(0.3)

dynamics1.single_compartment_react(initial_step=0.02, reaction_duration=4.0,
                                  variable_steps=True, explain_variable_steps=False)

# %%
#dynamics1.explain_time_advance()

# %%
dynamics1.plot_history(colors=['darkorange', 'green'], show_intervals=True, title_prefix="WITHOUT catalysis")

# %%
# Verify that the reaction has reached equilibrium
dynamics1.is_in_equilibrium(tolerance=5)

# %%

# %%

# %% [markdown]
# # PART 2

# %%
# Initialize the system
chem2 = ChemData(names=["A", "B", "L"])

chem2.add_macromolecules("M1")

chem2.set_binding_site_affinity("M1", site_number=1, ligand="L", Kd=5.0)


# %%
chem2.show_binding_affinities()        # Review the values we have given for the dissociation constants

# %%

# %% [markdown]
# ## Define 2 reactions, one with and one without catalysis

# %%
# Reaction A <-> B , without catalysis (same as he had in part 1)
chem2.add_reaction(reactants="A", products="B",
                  forward_rate=1. , delta_G= -5000)

# Reaction A <-> B , WITH catalysis (fast forward rate)
#rxn = chem2.add_reaction(reactants="A", products="B",
#                        forward_rate=10. , delta_G= -5000)

#rxn.set_macro_enzyme(macromolecule="M1", site_number=1)

chem2.describe_reactions()

# %% [markdown]
# #### Notice the thermodynamics is (as it should be!) the same - and hence the same equilibrium constants - but the kinetics are very different

# %%

# %% [markdown]
# ### Start with no macromolecule ligand `L` present; hence, only the slow reaction is in effect

# %%
# Set the initial concentrations of all the chemicals, including the macromolecule

dynamics2 = ReactionDynamics(chem_data=chem2)
dynamics2.set_conc(conc={"A": 100., "B": 20., "L": 0.},
                  snapshot=True)      # The macromolecule ligand L is absent

dynamics2.set_macromolecules()      # By default, set counts to 1 for all the registered macromolecules

dynamics2.describe_state()

# %%

# %% [markdown] tags=[]
# ### Take the initial system to equilibrium

# %%
dynamics2.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics2.set_thresholds(norm="norm_A", low=0.6, high=1.0, abort=1.44)
dynamics2.set_thresholds(norm="norm_B", low=0.04, high=0.6, abort=1.5)
dynamics2.set_step_factors(upshift=1.2, downshift=0.7, abort=0.4)
dynamics2.set_error_step_factor(0.3)

dynamics2.single_compartment_react(initial_step=0.02, reaction_duration=4.0,
                                  variable_steps=True, explain_variable_steps=False)

# %%
dynamics2.explain_time_advance()

# %%
dynamics2.plot_history(colors=['darkorange', 'green', 'darkblue'], show_intervals=True, title_prefix="WITHOUT catalysis")

# %%
# Verify that the reaction has reached equilibrium
dynamics2.is_in_equilibrium(tolerance=5)

# %%
