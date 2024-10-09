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
# ## The reversible Synthesis/Decomposition reaction `A + B <-> C`
# #### with 1st-order kinetics for each species, taken to equilibrium.
# #### Exploration of debugging and diagnostics options
# (Adaptive variable time steps are used)
#
# _See also the experiment "1D/reactions/reaction_4"_  

# %% [markdown]
# ### TAGS :  "uniform compartment", "under-the-hood"

# %%
LAST_REVISED = "Sep. 8, 2024"
LIFE123_VERSION = "1.0.0.beta.38"      # Version this experiment is based on

# %%
#import set_path            # Using MyBinder?  Uncomment this before running the next cell!
                            # Importing this local file will add the project's home directory to sys.path

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

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

# %% [markdown]
# # Initialize the System
# Specify the chemicals, the reactions, and the initial concentrations

# %%
# Instantiate the simulator and specify the chemicals
uc = UniformCompartment(preset="fast")

# %%
# Reaction A + B <-> C , with 1st-order kinetics for each species
uc.add_reaction(reactants=["A" , "B"], products="C",
                forward_rate=5., reverse_rate=2.)

uc.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
uc.plot_reaction_network("vue_cytoscape_2")

# %%
# Set the initial concentrations of all the chemicals
uc.set_conc({"A": 10., "B": 50., "C": 20.}, snapshot=True)

# %%
uc.describe_state()

# %%
uc.get_history()

# %%

# %% [markdown]
# ### Sneak preview of eventual equilibrum:
# we can preview the final equilibrium concentrations without actually running the simulation

# %%
uc.find_equilibrium_conc(rxn_index=0)    # This is an EXACT solution

# %% [markdown]
# The reaction will proceed forward, with `A` and `B` being consumed, and `C` being produced

# %%

# %% [markdown] tags=[]
# # Run the reaction

# %%
uc.enable_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

uc.single_compartment_react(initial_step=0.004, duration=0.06,
                                  variable_steps=True,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"})

# %%
uc.get_history()

# %% [markdown]
# ## Note: "A" (now largely depleted) is the limiting reagent

# %% tags=[]
uc.diagnostics.explain_time_advance()

# %% [markdown]
# ### Notice how the reaction proceeds in smaller steps in the early times, when the concentrations are changing much more rapidly.
# #### The argument argument _variable_steps=True_ dynamically adjusts the initial_step (which is initially found to be too large, leading to some backtracking)

# %%
uc.plot_step_sizes(show_intervals=True)

# %% [markdown] tags=[]
# ## Plots changes of concentration with time

# %%
uc.plot_history(colors=['red', 'darkorange', 'green'], show_intervals=True)

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium()

# %% [markdown]
# Compare with the values we saw earlier for the exact solution of the equilibrium values:   
# {'A': 0.2948774087575341, 'B': 40.294877408757536, 'C': 29.705122591242464}  
#
# It's instructive to compare the exact values with the last few points from the simulation:  

# %%
uc.get_history(tail=3)

# %% [markdown]
# The 2nd-to-last simulation point, rather than the last one, is actually closer to the exact equilibrium values.  
# That's because by that time the variable steps are getting so large that they introduce some error.  
# If we were to run the simulation longer (not shown), we'd see the variable steps continuing to grow, and then suddenly being reduced;  
# then continued cycles of growth and reduction ("hitting the brakes whenever getting too fast")

# %%

# %%

# %% [markdown]
# # _Everthing below is just for diagnostic insight_ 
# ## _into the adaptive variable time steps_   
# This information is available because we made a call to `dynamics.set_diagnostics()` prior to running the simulation

# %%
uc.diagnostics.get_diagnostic_decisions_data()   # diagnostic data about concentration changes at every step - EVEN aborted ones

# %%
uc.diagnostics.get_rxn_data(rxn_index=0)   # diagnostic run data of the requested SINGLE reaction

# %%
uc.diagnostics.get_diagnostic_conc_data()   # diagnostic concentration data saved during the run, regardless of how much history we requested to save

# %%
