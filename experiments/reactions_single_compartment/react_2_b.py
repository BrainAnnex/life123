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
# with 1st-order kinetics in both directions, taken to equilibrium.
#
# This is a repeat of the experiment _"react_2_a"_ , but with **adaptive variable time steps** 
# and the use of **diagnostic tools** for insight into the details of the simulation.  

# %%
LAST_REVISED = "July 22, 2024"
LIFE123_VERSION = "1.0.0.beta.37"  # Version this experiment is based on

# %% tags=[]
# If can't find module, do: 1) import sys  2) sys.path.append("full path of folder containing modules")
    
from experiments.get_notebook_info import get_notebook_basename

from life123 import check_version, UniformCompartment, GraphicLog

# %%
check_version(LIFE123_VERSION)

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%

# %%

# %% [markdown]
# # PART 1 - RUN THE SIMULATION

# %% [markdown]
# ### Initialize the System

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = UniformCompartment(names=["A", "B"], preset="mid")

# Reaction A <-> B , with 1st-order kinetics in both directions
dynamics.add_reaction(reactants="A", products="B", 
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
dynamics.enable_diagnostics()   # To save diagnostic information about the call to single_compartment_react()
                                # Useful for insight into the inner workings of the simulation

# %%
# For experiment repeatability, we specified, when instantiating the "UniformCompartment" class, 
# a particular preset applicable to the adaptive time steps; 
# that preset assigned the following values
dynamics.adaptive_steps.show_adaptive_parameters() 

# %% [markdown] tags=[]
# ## Run the reaction   
# #### Passing True to _variable_steps_ automatically adjusts up or down the time steps

# %%
dynamics.single_compartment_react(initial_step=0.1, target_end_time=1.2,
                                  variable_steps=True,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"}
                                  )

# %%
history = dynamics.get_history()   # The system's history, saved during the run of single_compartment_react()
history

# %% [markdown] tags=[]
# ## Notice how the reaction proceeds in smaller steps in the early times, when [A] and [B] are changing much more rapidly
# #### That resulted from passing the flag _variable_steps=True_ to single_compartment_react()

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %%

# %%

# %% [markdown]
# # PART 2 - Visualize the Results

# %% [markdown] tags=[]
# ### Plots of changes of concentration with time

# %%
dynamics.plot_history(colors=['darkturquoise', 'orange'], show_intervals=True)

# %% [markdown]
# ## Note how the left-hand side of this plot is much smoother than it was in experiment `react_2_a`, where no adaptive time steps were used!

# %% [markdown]
# #### Compare the above with the fixed step sizes of experiment `react_2_a`    
# To see the sizes of the steps taken:

# %%
dynamics.plot_step_sizes(show_intervals=True)

# %%

# %%

# %% [markdown] tags=[]
# # PART 2 - Scrutinizing the inner workings of the step-size changes

# %% [markdown]
# NOTE: this part is NOT meant for typically end users.  It's for debugging, and for anyone interested in taking an "under the hood" look

# %%
diagnostics = dynamics.diagnostics      # Available because we turned on diagnostics

# %% [markdown]
# The "Diagnostics" object contains a treasure trove of data and methods yo get insights into the inner workings of the simulation.   
# Diagnostic data was saved because of the call to `enable_diagnostics()` prior to running the reaction

# %%
type(diagnostics)

# %%
# Let's revisit the variables steps taken
dynamics.diagnostics.explain_time_advance()

# %%

# %% [markdown]
# #### The Delta-concentration values for all the individual reaction time steps, as contribued by a single reaction, may be inspected from the diagnostic data:

# %%
diagnostics.get_diagnostic_rxn_data(rxn_index=0)    # For the 0-th reaction (the only reaction in our case)

# %% [markdown]
# ### Note that diagnostic data with the DELTA Concentrations - in the above listing - also records the values that were considered (but not actually used) during *ABORTED* steps.  For example, in steps 0-2, above, the START_TIME remains the same, as the time_step gets progressively reduced until the resulting changes are deemed acceptable.

# %% [markdown]
# Line 2, above, shows that the concentration of the product [B], which was set by us to an initial value of 50, gets affected by a **rate of change of -70**, sustained over a delta_time of 0.016 .  
# The reaction rate is negative because _the product is decreasing._  
# The new value for [B] is:  

# %%
50. - 70. * 0.016

# %% [markdown]
# That's indeed the value we saw in the history of the product [B] at time t=0.016 :

# %%
dynamics.get_history(t=0.016)

# %%

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
dynamics.adaptive_steps.adjust_timestep(delta_conc=delta_concentrations, baseline_conc=baseline_conc,
                                        n_chems=2, indexes_of_active_chemicals=dynamics.chem_data.indexes_of_active_chemicals())  

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

# %%

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
dynamics.adaptive_steps.adjust_timestep(delta_conc=delta_concentrations, baseline_conc=baseline_conc,
                                        n_chems=2, indexes_of_active_chemicals=dynamics.chem_data.indexes_of_active_chemicals())  

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
dynamics.adaptive_steps.show_adaptive_parameters()

# %% [markdown]
# **1.2** is stored as the "step factor" (for the time steps to take) in case an _'upshift'_ (in step size) is the decided course of action

# %%

# %% [markdown]
# ## Diagnostics of the run may be investigated as follows:  
# _(note - this is possible because we make a call to set_diagnostics() prior to running the simulation)_

# %%
diagnostics.get_diagnostic_conc_data()   # This will be complete, even if we only saved part of the history during the run

# %%
diagnostics.get_diagnostic_decisions_data()

# %%
