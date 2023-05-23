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
# ## Association/Dissociation reaction `2A <-> C`
# #### with 2nd-order kinetics for `A`,  
# #### and 1-st order kinetics for `C`
#
# Taken to equilibrium.  (Adaptive variable time teps are used)
#
# _See also the experiment "1D/reactions/reaction_7"_ 
#
# LAST REVISED: May 22, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.express as px
from src.modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # Initialize the System
# Specify the chemicals and the reactions

# %% tags=[]
# Specify the chemicals
chem_data = chem(names=["A", "C"])

# Reaction 2A <-> C , with 2nd-order kinetics for A, and 1st-order kinetics for C
chem_data.add_reaction(reactants=[(2, "A", 2)], products=["C"],
                       forward_rate=3., reverse_rate=2.)   
# Note: the first 2 in (2, "A", 2) is the stoichiometry coefficient, while the other one is the order

chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# # Start the simulation

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)

# %%
# Initial concentrations of all the chemicals, in index order
dynamics.set_conc([200., 40.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.get_history()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# All of these settings are slight variations of the currently the default values... but subject to change
dynamics.set_thresholds(norm="norm_A", low=1.0, high=2.0, abort=3.5)
dynamics.set_thresholds(norm="norm_B", low=0.1, high=0.5, abort=3.0)
dynamics.set_step_factors(abort=0.5, downshift=0.5, upshift=1.5)

dynamics.single_compartment_react(time_step=0.002, reaction_duration=0.04,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True, explain_variable_steps=False)

# %% [markdown]
# ### Note how the (tentative) original time step that we provide, 0.002, turned out to be so large that the simulation backtracks several times, because of "hard" aborts (negative concentrations) or "soft" aborts (concentration changes surpassing the thresholds we provided)

# %%
dynamics.get_history()

# %%
dynamics.plot_curves(colors=['red', 'green'],
                     title="Reaction 2A <-> C  (2nd order in A).  Changes in concentrations with time")

# %% [markdown]
# ## Note: "A" (now largely depleted) is the limiting reagent

# %% [markdown]
# #### Let's take a look at time t=0.002, which in our simulation run had proposed as the first step:

# %%
# Locate the value closest to the original time step we had requested
dynamics.get_history(t=0.002)

# %% [markdown]
# ### Because of the very large changes happening between t=0 and 0.002, the simulation automatically slowed down and opted to actually take 74 steps in lieu of the 1 step we had proposed.
# The number of variable steps actually taken can be modulated by changing the "norm thresholds" and "step factors" that we can optionally specify

# %% [markdown]
# ### Notice how, later in the simulation, the step sizes get BIGGER than the 0.002 we had originally proposed:

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Notice how the reaction proceeds in far-smaller steps in the early times, when the concentrations are changing much more rapidly

# %%
# Let's look at the first two arrays of concentrations, from the run's history
arr0 = dynamics.get_historical_concentrations(0)   # The initial concentrations
arr1 = dynamics.get_historical_concentrations(1)   # After the first actual step
arr0, arr1

# %%
# Let's verify that the reaction's stoichiometry is being respected
dynamics.stoichiometry_checker(rxn_index=0, 
                               conc_arr_before = arr0, 
                               conc_arr_after = arr1)

# %% [markdown]
# #### Indeed, it can be easy checked that the drop in [A] is twice the increase in [C], as dictated by the stoichiometry

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0, head=15)    # Easily seen in the diagnostic data

# %% [markdown]
# ### From the diagnostic data, it can be seen that the first step had several false starts - and was automatically repeatedly shrunk - but finally happened.  `Delta A` indeed equals - 2 * `Delta C`, satisfying the stoichiometry

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%
dynamics.stoichiometry_checker_entire_run()

# %% [markdown] tags=[]
# ## Display the variable time steps

# %%
dynamics.plot_curves(colors=['red', 'green'], show_intervals=True,
                     title="Reaction 2A <-> C  (2nd order in A).  Changes in concentrations with time")

# %% [markdown]
# ### The intersection of the two lines may be found as follows:

# %%
dynamics.curve_intersection('A', 'C', t_start=0, t_end=0.01)

# %%

# %% [markdown]
# #### For diagnostic insight, uncomment the following lines:

# %%
#dynamics.get_diagnostic_decisions_data()

#dynamics.get_diagnostic_rxn_data(rxn_index=0)

#dynamics.get_diagnostic_conc_data()

# %%
