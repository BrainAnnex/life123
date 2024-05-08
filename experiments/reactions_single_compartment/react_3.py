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
# ## Association/Dissociation reaction `A + B <-> C`
# #### with 1st-order kinetics for each species, taken to equilibrium.
# #### Exploration of debugging and diagnostics options
# (Adaptive variable time steps are used)
#
# _See also the experiment "1D/reactions/reaction_4"_  
#
# LAST REVISED: May 5, 2024

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%

# %% [markdown]
# # Initialize the System
# Specify the chemicals, the reactions, and the initial concentrations

# %%
# Instantiate the simulator and specify the chemicals
dynamics = ReactionDynamics(names=["A", "B", "C"])

# %%
# Reaction A + B <-> C , with 1st-order kinetics for each species
dynamics.add_reaction(reactants=["A" , "B"], products="C",
                      forward_rate=5., reverse_rate=2.)

dynamics.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
dynamics.plot_reaction_network("vue_cytoscape_2")

# %%
# Set the initial concentrations of all the chemicals, in their index order
dynamics.set_conc([10., 50., 20.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.get_history()

# %%

# %% [markdown]
# ### Sneak preview of eventual equilibrum:
# we can preview the final equilibrium concentrations without actually running the simulation

# %%
dynamics.find_equilibrium_conc(rxn_index=0)    # This is an EXACT solution

# %% [markdown]
# The reaction will proceed forward, with `A` and `B` being consumed, and `C` being produced

# %%

# %% [markdown] tags=[]
# # Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# For repeatibility, we avoid the defaults, and instead specify a particular group of preset parameters 
# applicable to the adaptive time steps.
# Here we use a "fast" heuristic: advance quickly thru time
dynamics.use_adaptive_preset(preset="fast")

dynamics.single_compartment_react(initial_step=0.004, reaction_duration=0.06,
                                  variable_steps=True, explain_variable_steps=False,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"})

# %%
dynamics.get_history()

# %% [markdown]
# ## Note: "A" (now largely depleted) is the limiting reagent

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Notice how the reaction proceeds in smaller steps in the early times, when the concentrations are changing much more rapidly.
# #### The argument argument _variable_steps=True_ dynamically adjusts the initial_step (which is initially found to be too large, leading to some backtracking)

# %%
dynamics.plot_step_sizes(show_intervals=True)

# %% [markdown] tags=[]
# ## Plots changes of concentration with time

# %%
dynamics.plot_history(colors=['red', 'violet', 'green'], show_intervals=True)

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown]
# Compare with the values we saw earlier for the exact solution of the equilibrium values:   
# {'A': 0.2948774087575341, 'B': 40.294877408757536, 'C': 29.705122591242464}  
#
# It's instructive to compare the exact values with the last few points from the simulation:  

# %%
dynamics.get_history(tail=3)

# %% [markdown]
# The 2nd-to-last simulation point, rather than the last one, is actually closer to the exact equilibrium values.  
# That's because by that time the variable steps are getting so large that they introduce some error.  
# If we were to run the simulation longer (not shown), we'd see the variable steps continuing to grow, and then suddenly being reduced;  
# then continued cycles of growth and reduction ("hitting the brakes whenever getting too fast")

# %%

# %% [markdown]
# # Everthing below is just for diagnostic insight 
# ## into the adaptive variable time steps  
# This information is available because we made a call to `dynamics.set_diagnostics()` prior to running the simulation

# %%
dynamics.get_diagnostic_decisions_data()   # diagnostic data about concentration changes at every step - EVEN aborted ones

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)   # diagnostic run data of the requested SINGLE reaction

# %%
dynamics.get_diagnostic_conc_data()   # diagnostic concentration data saved during the run, regardless of how much history we requested to save

# %%
