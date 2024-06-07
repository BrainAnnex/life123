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
# ### Adaptive time steps (variable time resolution) for reaction `A <-> B`,
# with 1st-order kinetics in both directions, taken to equilibrium
#
# This is a repeat of the experiment _"react_2_a"_ , but with adaptive variable time steps
#
# LAST REVISED: May 6, 2024

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_dynamics import UniformCompartment
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
# # PART 1 - RUN THE SIMULATION

# %% [markdown]
# ### Initialize the System

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = UniformCompartment(names=["A", "B"])

# Reaction A <-> B , with 1st-order kinetics in both directions
dynamics.add_reaction(reactants=["A"], products=["B"], 
                       forward_rate=3., reverse_rate=2.)

dynamics.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
dynamics.plot_reaction_network("vue_cytoscape_2")

# %%
# Set the initial concentrations of all the chemicals, in their index order
dynamics.set_conc([10., 50.])

dynamics.describe_state()

# %%
dynamics.get_history()

# %%
dynamics.set_diagnostics()      # To save diagnostic information about the call to single_compartment_react()
                                # Useful for insight into the inner workings of the simulation

# %%
# For experiment repeatability, we specify a particular group of preset parameters applicable to the adaptive time steps
dynamics.use_adaptive_preset(preset="mid")       # A "middle-of-the road" heuristic: somewhat "conservative" but not overly so

dynamics.show_adaptive_parameters()              # Details may vary across different versions of Life123

# %% [markdown] tags=[]
# ### Run the reaction

# %%
dynamics.single_compartment_react(initial_step=0.1, target_end_time=1.2,
                                  variable_steps=True, explain_variable_steps=False,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"}
                                  )

# %% [markdown]
# ## The flag _variable_steps_ automatically adjusts up or down the time step

# %%
history = dynamics.get_history()   # The system's history, saved during the run of single_compartment_react()
history

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# ### That resulted from passing the flag _variable_steps=True_ to single_compartment_react()

# %%

# %% [markdown] tags=[]
# # PART 2 - Scrutinizing some instances of step-size changes

# %% [markdown]
# The Delta-concentration values for all the individual reaction time steps, as contribued by a single reaction, may all be inspected at once from the diagnostic data:

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)    # For the 0-th reaction (the only reaction in our case)

# %% [markdown]
# ### Note that diagnostic data with the DELTA Concentrations - above and below - also record the values that were considered (but not actually used) during ABORTED steps

# %% [markdown]
# In the examples below, we'll re-compute Delta values for individual steps, directly from the system history.

# %% [markdown] tags=[]
# ### Example 1: **very early in the run**    

# %%
history[1:4]

# %%
delta_concentrations = dynamics.extract_delta_concentrations(history, 1, 2, ['A', 'B'])
delta_concentrations

# %% [markdown]
# As expected by the 1:1 stoichiometry, delta_A = - delta_B

# %%
# Get all the concentrations at the start of the above steps
baseline_conc = dynamics.get_historical_concentrations(row=1)
#dynamics.get_historical_concentrations(t=0.016)   # Alternate way
baseline_conc

# %%
# Computes some measures of how large delta_concentrations is, and propose a course of action
dynamics.adjust_speed(delta_conc=delta_concentrations, baseline_conc=baseline_conc)  

# %% [markdown]
# #### The above analysis indicates that the time step is just about right, and the simulations should STAY on that course : that's based on the shown computed norms (indicating the extent of the change taking place.)  
# Indeed, the simulator maintains the same time step :

# %%
original_step = history["SYSTEM TIME"][2] - history["SYSTEM TIME"][1]
original_step

# %%
next_step = history["SYSTEM TIME"][3] - history["SYSTEM TIME"][2]
next_step

# %%
next_step / original_step

# %% [markdown] tags=[]
# ### Example 2: **very late in the run**    

# %%
history[17:20]

# %%
delta_concentrations = dynamics.extract_delta_concentrations(history, 17, 18, ['A', 'B'])
delta_concentrations

# %% [markdown]
# #### Notice the far less change now that the system is approaching equilibrium

# %%
# Get all the concentrations at the start of the above steps
baseline_conc = dynamics.get_historical_concentrations(row=17)
#dynamics.get_historical_concentrations(t=0.850186)   # Alternate way
baseline_conc

# %%
# Computes a measure of how large delta_concentrations is, and propose a course of action
dynamics.adjust_speed(delta_conc=delta_concentrations, baseline_conc=baseline_conc)  

# %% [markdown]
# #### The above analysis indicates that the time step is on the "LOW" side, and the simulations should increase it by a factor 1.2 : again, that's based on the shown computed norms (indicating the extent of the change taking place.)  
# Indeed, the simulator increases the time step x1.2:

# %%
original_step = history["SYSTEM TIME"][18] - history["SYSTEM TIME"][17]
original_step

# %%
next_step = history["SYSTEM TIME"][19] - history["SYSTEM TIME"][18]
next_step

# %%
next_step / original_step

# %% [markdown]
# Where does that x1.2 factor come from?  It's one of the parameters that we passed to the simulator; they can be seen as follows:

# %%
dynamics.show_adaptive_parameters()

# %% [markdown]
# **1.2** is stored as the "step factor" (for the time steps to take) in case an 'upshift' (in step size) is the decided course of action

# %%

# %% [markdown]
# ## Diagnostics of the run may be investigated as follows:  
# _(note - this is possible because we make a call to set_diagnostics() prior to running the simulation)_

# %%
dynamics.get_diagnostic_conc_data()   # This will be complete, even if we only saved part of the history during the run

# %%
dynamics.get_diagnostic_decisions_data()

# %%

# %% [markdown]
# # PART 3 - Analyze the reaction Dynamics

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# ### Plots of changes of concentration with time

# %%
dynamics.plot_history(colors=['blue', 'orange'], show_intervals=True)

# %% [markdown]
# ## Note how the left-hand side of this plot is much smoother than it was in experiment `react_2_a`, where no adaptive time steps were used!

# %% [markdown]
# #### Compare the above with the fixed step sizes of experiment `react_2_a`    
# To see the sizes of the steps taken:

# %%
dynamics.plot_step_sizes(show_intervals=True)

# %%
