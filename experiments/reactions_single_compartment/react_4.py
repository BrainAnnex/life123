# ---
# jupyter:
#   jupytext:
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
# ## Elementary Association/Dissociation reaction `2A <-> C`
# #### Standard mass-action kinetics with 2nd-order for `A`,  
# #### and 1-st order for `C`
#
# Taken to equilibrium.  (Adaptive variable time teps are used)
#
# _See also the experiment "1D/reactions/reaction_7"_ 

# %% [markdown]
# ### TAGS :  "basic", "under-the-hood", "uniform compartment"

# %%
LAST_REVISED = "Mar. 12, 2026"
LIFE123_VERSION = "1.0.0rc7"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import numpy as np
import ipynbname

from life123 import check_version, UniformCompartment

# %%
check_version(LIFE123_VERSION)

# %%
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file
                                            # IN CASE OF PROBLEMS, set manually to any desired name

# %%

# %%

# %% [markdown]
# # PART 1 - The Simulation

# %% [markdown]
# ### Initialize the System
# Specify the chemicals, the reactions, and the initial state

# %%
# Instantiate the simulator and specify the chemicals
# The diagnostics will be give insight into the inner workings of the simulation
uc = UniformCompartment(names=["A", "C"], preset="fast", enable_diagnostics=True)

# Reaction 2A <-> C , with 2nd-order kinetics for A, and 1st-order kinetics for C
uc.add_reaction(reactants=[(2, "A")], products="C", kF=3., kR=2.)

uc.describe_reactions()

# %%

# %%
# Send a plot of the network of reactions to the HTML log file
uc.plot_reaction_network(log_file=log_file)

# %%

# %%
# Initial concentrations of all the chemicals
uc.set_conc({"A": 200., "C": 40.})
uc.describe_state()

# %%
uc.get_history()

# %%

# %% [markdown]
# ### Let's look ahead at the final equilibrium point

# %%
uc.find_equilibrium_conc(rxn_index=0)    # This is an EXACT equilibrium solution, 
                                         # for 1 reaction (the only reaction)

# %%

# %%

# %% [markdown]
# ## Run the reaction

# %%
uc.single_compartment_react(initial_step=0.002, duration=0.03)    # Variable steps is the default

# %% [markdown]
# ### Note how the (tentative) original time step that we provide, 0.002, turned out to be so large that the simulation backtracks several times, because of "hard" aborts (negative concentrations) or "soft" aborts (concentration changes surpassing the thresholds we provided)
#
# #### For example, if we check the history, we can see that the first step was automatically reduced from 0.002 to 0.000008

# %%
uc.get_history()

# %%
uc.plot_history(colors=['darkturquoise', 'green'],
                title="Reaction 2A <-> C  (2nd order in A).  Changes in concentrations with time")

# %% [markdown]
# ### Note: "A" (now largely depleted) is the limiting reagent

# %% [markdown]
# The values indeed appear to converge to the exact analytic solution for the equilibrium, from earlier:   {'A': 9.49568869375716, 'C': 135.2521556531214}

# %%

# %%

# %% [markdown]
# # PART 2 - Analysis and Validation

# %% [markdown]
# #### Let's take a look at time t=0.002, which in our simulation run had proposed as the first step:

# %%
# Locate the value closest to the original time step we had requested
uc.get_history(t=0.002)

# %% [markdown]
# ### Because of the very large changes happening between t=0 and 0.002, the simulation automatically slowed down and opted to actually take 80 steps in lieu of the 1 step we had (rather optimistically!) proposed  
# The number of variable steps actually taken can be modulated by changing the _preset_ passed when the `UniformCompartment` class is first instantiated - or, alternatively, using calls to `use_adaptive_preset()`. For finer control, advanced users may tweak internal parameters such as "norm thresholds" and "step factors"

# %% [markdown]
# ### Step sizes can work both ways!  Notice how, late in the simulation, the step sizes get BIGGER than the 0.002 we had originally proposed:

# %%
uc.get_diagnostics().explain_time_advance()

# %% [markdown]
# ### One can see how the reaction proceeds in far-smaller steps in the early times, when the concentrations are changing much more rapidly

# %%
# Let's look at the first two arrays of concentrations, from the run's history
arr0 = uc.get_historical_concentrations(0)   # The initial concentrations
arr1 = uc.get_historical_concentrations(1)   # After the first actual simulation step
arr0, arr1

# %%
# Let's verify that the reaction's stoichiometry is being respected
uc.get_diagnostics().stoichiometry_checker(rxn_index=0, 
                                           conc_arr_before = arr0, 
                                           conc_arr_after = arr1)

# %% [markdown]
# #### Indeed, it can be easy checked that the drop in [A] is twice the increase in [C], as dictated by the stoichiometry.
# The diagnostic data, enabled by our earlier call to `set_diagnostics()`, makes it convenient to check

# %%
uc.get_diagnostics().get_rxn_data(rxn_index=0, head=15)

# %% [markdown]
# ### From the diagnostic data, it can be seen that the first step had several false starts - and the time was automatically repeatedly shrunk - but finally happened.  `Delta A` indeed equals - 2 * `Delta C`, satisfying the stoichiometry

# %%
uc.get_diagnostics().stoichiometry_checker_entire_run()

# %%

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium()

# %% [markdown]
# ## Display the variable time steps

# %%
uc.plot_history(colors=['darkturquoise', 'green'], show_intervals=True,
                title="Reaction 2A <-> C  (2nd order in A).  Concentrations changes")

# %% [markdown]
# ### The intersection of the two lines may be found as follows:

# %%
uc.curve_intersect('A', 'C', t_start=0, t_end=0.01)

# %%

# %% [markdown]
# #### For additional diagnostic insight:
# `norm_A` and `norm_B` are computed quantities that are used to guide the adaptive time steps

# %%
uc.get_diagnostics().get_decisions_data()

# %%
