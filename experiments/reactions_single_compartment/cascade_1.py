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
# ## 2 COUPLED reactions of different speeds, forming a "cascade":  
# ### `A <-> B` (fast) and `B <-> C` (slow)
# Taken to equilibrium. Both reactions are mostly forward. All 1st order.  
# **The concentration of the intermediate product B manifests 1 oscillation (transient "overshoot")**
#
# Adaptive variable time resolution is used, with extensive diagnostics, 
# and finally compared to a new run using fixed time intervals, with the same initial data.
#
# In part2, some diagnostic insight is explored.  
#
# In part3, two identical runs ("adaptive variable steps" vs. "fixed small steps") are compared. 

# %% [markdown]
# ## Bathtub analogy:
# A is initially full, while B and C are empty.  
# Tubs are progressively lower (reactions are mostly forward.)  
# A BIG pipe connects A and B: fast kinetics.  A small pipe connects B and C: slow kinetics. 
#
# INTUITION: B, unable to quickly drain into C while at the same time being blasted by a hefty inflow from A,  
# will experience a transient surge, in excess of its final equilibrium level.
#
# * **[Compare with the final reaction plot (the orange line is B)](#cascade_1_plot)**

# %% [markdown]
# ![2 Coupled Reactions](../../docs/2_coupled_reactions.png)

# %%
LAST_REVISED = "Oct. 11, 2024"
LIFE123_VERSION = "1.0.0.beta.39"   # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

import ipynbname

from life123 import check_version, ChemData, UniformCompartment, PlotlyHelper, GraphicLog

# %%
check_version(LIFE123_VERSION)

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # Initialize the System
# Specify the chemicals and the reactions

# %% tags=[]
# Specify the chemicals and the reactions; this data structure will get re-used in 
# the various simulations below
chem = ChemData()

# Reaction A <-> B (fast)
chem.add_reaction(reactants="A", products="B",
                       forward_rate=64., reverse_rate=8.) 

# Reaction B <-> C (slow)
chem.add_reaction(reactants="B", products="C",
                       forward_rate=12., reverse_rate=2.) 

print("Number of reactions: ", chem.number_of_reactions())

# %%
chem.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
chem.plot_reaction_network("vue_cytoscape_2")

# %%

# %% [markdown]
# ## Run the simulation

# %%
dynamics = UniformCompartment(chem_data=chem, preset="fast")

dynamics.set_conc({"A": 50.}, snapshot=True) # Set the initial concentrations of all the chemicals, in their index order
dynamics.describe_state()

# %%
dynamics.get_history()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.enable_diagnostics()    # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(initial_step=0.02, duration=0.4,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True)

# %% [markdown]
# ### <a name="cascade_1_plot"> Plots of changes of concentration with time</a>
# Notice the variable time steps (vertical dashed lines)

# %%
dynamics.plot_history(title="Coupled reactions A <-> B and B <-> C",
                      colors=['darkturquoise', 'orange', 'green'], show_intervals=True)

# %%
dynamics.curve_intersect("A", "B", t_start=0, t_end=0.05)

# %%
dynamics.curve_intersect("A", "C", t_start=0, t_end=0.05)

# %%
dynamics.curve_intersect("B", "C", t_start=0.05, t_end=0.1)

# %%
dynamics.get_history()

# %%
dynamics.diagnostics.explain_time_advance()

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown]
# ### Let's look at the final concentrations of `A` and `C` (i.e., the reactant and product of the composite reaction)

# %%
A_final = dynamics.get_chem_conc("A")
A_final

# %%
C_final = dynamics.get_chem_conc("C")
C_final

# %% [markdown]
# #### Their ratio:

# %%
C_final / A_final

# %% [markdown]
# ### As expected the equilibrium constant for the overall reaction `A <-> C` (approx. 48) is indeed the product of the equilibrium constants of the two elementary reactions (K = 8 and K = 6, respectively) that we saw earlier.

# %%

# %%

# %% [markdown]
# # PART 2 - DIAGNOSTIC INSIGHT
#
# Perform some verification

# %% [markdown]
# ### Take a peek at the diagnostic data saved during the earlier reaction simulation

# %%
# Concentration increments due to reaction 0 (A <-> B)
# Note that [C] is not affected
dynamics.diagnostics.get_rxn_data(rxn_index=0)

# %%
# Concentration increments due to reaction 1 (B <-> C) 
# Also notice that the 0-th row from the A <-> B reaction isn't seen here (start time 0 and step 0.02), 
# because that step was aborted early on, BEFORE even getting to THIS reaction
dynamics.diagnostics.get_rxn_data(rxn_index=1)

# %%

# %%

# %% [markdown] tags=[]
# # PART 3 : Re-run with very small constant steps, and compare with original run

# %% [markdown]
# We'll use **constant steps of size 0.0005** , which is 1/4 of the smallest steps (the "substep" size) previously used in the variable-step run

# %%
dynamics2 = UniformCompartment(chem_data=chem)  # Re-use the same chemicals and reactions of the previous simulation

# %% tags=[]
dynamics2.set_conc({"A": 50.}, snapshot=True)

# %%
# Notice that we're using FIXED steps this time
dynamics2.single_compartment_react(initial_step=0.0005, duration=0.4,
                                   variable_steps=False,
                                   snapshots={"initial_caption": "1st reaction step",
                                              "final_caption": "last reaction step"},
                                   )

# %%
dynamics2.plot_history(title="Coupled reactions A <-> B and B <-> C , re-run with CONSTANT STEPS",
                       colors=['darkturquoise', 'orange', 'green'], show_intervals=True)

# %% [markdown]
# _(Notice that the vertical steps are now equally spaced - and that there are so many of them that we're only showing 1 every 6)_

# %%
dynamics2.curve_intersect(t_start=0, t_end=0.05, chem1="A", chem2="B")

# %%
dynamics2.curve_intersect(t_start=0, t_end=0.05, chem1="A", chem2="C")

# %%
dynamics2.curve_intersect(t_start=0.05, t_end=0.1, chem1="B", chem2="C")

# %%
df2 = dynamics2.get_history()
df2

# %% [markdown]
# ## Notice that we now did 800 steps - vs. the 48 of the earlier variable-resolution run!

# %% [markdown]
# ## Let's compare some entries with the coarser previous variable-time run

# %% [markdown]
# #### Let's compare the plots of [B] from the earlier (variable-step) run, and the latest (high-precision, fixed-step) one:

# %%
# Earlier run (using variable time steps)
fig1 = dynamics.plot_history(chemicals='B', colors='orange', title="Adaptive variable-step run", 
                             show=True)

# %%
# Latest run (high-precision result from fine fixed-resolution run)
fig2 = dynamics2.plot_history(chemicals='B', colors=['violet'], title="Fine fixed-step run", 
                              show=True)

# %%
PlotlyHelper.combine_plots(fig_list=[fig1, fig2], title="The 2 runs, contrasted together", 
                           curve_labels=["B (adaptive variable steps)", "B (fixed small steps)"])

# %% [markdown]
# #### They overlap fairly well!  The 800 fixed-timestep points vs. the 48 adaptable variable-timestep ones.
# The adaptive algorithms avoided 752 extra steps of limited benefit...

# %%
