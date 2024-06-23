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
# LAST REVISED: June 23, 2024 (using v. 1.0 beta36)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from life123 import UniformCompartment
from life123 import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%

# %% [markdown]
# # PART 1 - The Simulation

# %% [markdown]
# # Initialize the System
# Specify the chemicals, the reactions, and the initial state

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = UniformCompartment(names=["A", "C"], preset="fast")

# Reaction 2A <-> C , with 2nd-order kinetics for A, and 1st-order kinetics for C
dynamics.add_reaction(reactants=[(2, "A")], products="C",
                      forward_rate=3., reverse_rate=2.)   
# Note: the reaction order for a chemical defaults to its stoichiometry coefficient; 
#       to specify it explicitly, pass it as 3rd term in tuple:  (2, "A", 2)

dynamics.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
dynamics.plot_reaction_network("vue_cytoscape_2")

# %%
# Initial concentrations of all the chemicals, in their index order
dynamics.set_conc([200., 40.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.get_history()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(initial_step=0.002, duration=0.03,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True)

# %% [markdown]
# ### Note how the (tentative) original time step that we provide, 0.002, turned out to be so large that the simulation backtracks several times, because of "hard" aborts (negative concentrations) or "soft" aborts (concentration changes surpassing the thresholds we provided)
#
# #### For example, the first step was automatically reduced from 0.002 to 0.000008

# %%
dynamics.get_history()

# %%
dynamics.plot_history(colors=['darkturquoise', 'green'],
                      title="Reaction 2A <-> C  (2nd order in A).  Changes in concentrations with time")

# %% [markdown]
# ### Note: "A" (now largely depleted) is the limiting reagent

# %%

# %%

# %% [markdown]
# # PART 2 - Analysis and Validation

# %% [markdown]
# #### Let's take a look at time t=0.002, which in our simulation run had proposed as the first step:

# %%
# Locate the value closest to the original time step we had requested
dynamics.get_history(t=0.002)

# %% [markdown]
# ### Because of the very large changes happening between t=0 and 0.002, the simulation automatically slowed down and opted to actually take 80 steps in lieu of the 1 step we had (rather optimistically!) proposed  
# The number of variable steps actually taken can be modulated by changing the _preset_ passed when the `UniformCompartment` class is first instantiated - or, alternatively, using calls to `use_adaptive_preset()`. For finer control, advanced users may tweak internal parameters such as "norm thresholds" and "step factors"

# %% [markdown]
# ### Notice how, late in the simulation, the step sizes get BIGGER than the 0.002 we had originally proposed:

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### One can see how the reaction proceeds in far-smaller steps in the early times, when the concentrations are changing much more rapidly

# %%
# Let's look at the first two arrays of concentrations, from the run's history
arr0 = dynamics.get_historical_concentrations(0)   # The initial concentrations
arr1 = dynamics.get_historical_concentrations(1)   # After the first actual simulation step
arr0, arr1

# %%
# Let's verify that the reaction's stoichiometry is being respected
dynamics.stoichiometry_checker(rxn_index=0, 
                               conc_arr_before = arr0, 
                               conc_arr_after = arr1)

# %% [markdown]
# #### Indeed, it can be easy checked that the drop in [A] is twice the increase in [C], as dictated by the stoichiometry.
# The diagnostic data, enabled by our earlier call to `set_diagnostics()`, makes it convenient to check

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0, head=15)

# %% [markdown]
# ### From the diagnostic data, it can be seen that the first step had several false starts - and the time was automatically repeatedly shrunk - but finally happened.  `Delta A` indeed equals - 2 * `Delta C`, satisfying the stoichiometry

# %%
dynamics.stoichiometry_checker_entire_run()

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown] tags=[]
# ## Display the variable time steps

# %%
dynamics.plot_history(colors=['darkturquoise', 'green'], show_intervals=True,
                      title="Reaction 2A <-> C  (2nd order in A).  Concentrations changes")

# %% [markdown]
# ### The intersection of the two lines may be found as follows:

# %%
dynamics.curve_intersect('A', 'C', t_start=0, t_end=0.01)

# %%

# %% [markdown]
# #### For additional diagnostic insight:
# `norm_A` and `norm_B` are computed quantities that are used to guide the adaptive time steps

# %%
dynamics.get_diagnostic_decisions_data()

# %%
