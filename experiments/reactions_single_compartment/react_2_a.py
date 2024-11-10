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
# ## An `A <-> B` reaction  
# with 1st-order kinetics in both directions, taken to equilibrium,
# using a simple, **coarse fixed-timestep simulation.**  
#
# In Part 2, below, perform some analysis of the results: in particular, examine the **reaction rates**. 
#
# (See also the experiment _"1D/reactions/reaction_1"_ for a multi-compartment version)  
#
# #### This experiment gets repeated in _"react_2_b"_ , with a more sophisticated approach, 
# #### involving adaptive variable time steps

# %%
LAST_REVISED = "Nov. 9, 2024"
LIFE123_VERSION = "1.0.0.rc.0"      # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import numpy as np
import ipynbname

from life123 import check_version, UniformCompartment, PlotlyHelper, GraphicLog

# %%
check_version(LIFE123_VERSION)

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file
                                            # IN CASE OF PROBLEMS, set manually to any desired name

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%

# %%

# %% [markdown]
# # PART 1 - simulation of the reaction

# %% tags=[]
# Instantiate the simulator and specify the chemicals
uc = UniformCompartment()

# Reaction A <-> B , with 1st-order kinetics in both directions
uc.add_reaction(reactants="A", products="B", 
                forward_rate=3., reverse_rate=2.)

print("Number of reactions: ", uc.number_of_reactions())

# %%
uc.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
uc.plot_reaction_network("vue_cytoscape_2")


# %%

# %%
# Initial concentrations of all the chemicals
uc.set_conc({"A": 10., "B": 50.})

# %%
uc.describe_state()

# %%
uc.get_history()

# %% [markdown]
# ### Test your intuition: 
# #### given that this reaction operates mostly in the forward direction (kF = 3 , kR = 2 , K = 1.5), 
# #### do you think that A will be consumed and B will be produced??
# We can take a sneak preview at the final equilibrium concentrations without actually running the simulation:

# %%
uc.find_equilibrium_conc(rxn_index=0)    # This is an EXACT equilibrium solution, 
                                         # for 1 reaction (the only reaction)

# %% [markdown]
# #### The reaction will actually proceed IN REVERSE, because of the large initial concentration of B (which we had set to 50), relative to the small initial concentration of A (10)
# Now, let's see the reaction in action!

# %%

# %% [markdown] tags=[]
# ### Run the reaction

# %%
# First step of reaction
uc.single_compartment_react(initial_step=0.1, n_steps=1, variable_steps=False) # NOT using variable steps!

# %%
uc.get_history()

# %% [markdown]
# We can already see the reaction proceeding in reverse...

# %%
# Numerous more fixed steps
uc.single_compartment_react(initial_step=0.1, n_steps=10, variable_steps=False, 
                            snapshots={"initial_caption": "2nd round of simulation"})

# %%
uc.get_history()

# %% [markdown]
# ## NOTE: for demonstration purposes, we're using FIXED time steps...  
# ## Typically, one would use the option for adaptive variable time steps (see experiment `react_2_b`)

# %% [markdown]
# ### Check the final equilibrium

# %%
uc.get_system_conc()   # The current concentrations, in the order the chemicals were added 

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order of the reactions), the systems settles in the following equilibrium:
#
# [A] = 23.99316406
#  
# [B] = 36.00683594
#

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium()

# %% [markdown]
# ### As noted earlier, because of the high initial concentration of B relative to A, the overall reaction has proceeded IN REVERSE

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
uc.plot_history(colors=['darkturquoise', 'green'], show_intervals=True)

# %% [markdown]
# ### Note the raggedness of the left-side (early times) of the curves.  
# ### In experiment `react_2_b` this simulation gets repeated with an _adaptive variable time resolution_ that takes smaller steps at the beginning, when the reaction is proceeding faster   
# #### By contrast, here we used _FIXED_ time steps (shown in dashed lines), which generally gives poor results, unless taking a very large number of very small steps!  In particular, the early steps we took were too large, and the later steps were unnecessarily small

# %%

# %%

# %% [markdown]
# # PART 2 - Now investigate A_dot, i.e. d[A]/dt

# %% [markdown]
# There's actually no need to compute this; it is automatically saved during the reaction simulations

# %%
rates_df = uc.get_rate_history()
rates_df

# %% [markdown]
# Reaction rates refer to products; since A is a reactant (in reaction 0, our only reaction), we need to flip its sign

# %%
rates_df['A_dot'] = - rates_df['rxn0_rate']    # Add a column to the Pandas dataframe
rates_df

# %%
uc.get_history()         # Revisited from earlier

# %% [markdown]
# Notice that we lack a rate for the last time value, in the above table, because no reaction simulation starting at that time has been performed

# %%
p1 = uc.plot_history(chemicals="A", colors="darkturquoise")   # The plot of [A] from the system history

# %%
p2 = PlotlyHelper.plot_pandas(df=rates_df, x_var="SYSTEM TIME", fields="A_dot", colors="brown")   # The plot of A_dot, from rates_df

# %%
PlotlyHelper.combine_plots([p1, p2], 
                           title="Concentration of A with time, and its rate of change (A_dot)",
                           y_label="[A] (turquoise) /<br> A_dot (brown)",
                           legend_title="Plot",
                           curve_labels=["A", "A_dot"])

# %% [markdown]
# ### At t=0 :  
# [A]=10 and [A] has a high rate of change (70)
# ### As the system approaches equilibrium :  
# [A] approaches a value of 24, and its rate of change decays to zero.

# %% [markdown]
# #### **NOTE:** The above curves are jagged because of _the large time steps taken_ (especially in the early times, when there's a lot of change.)  
# ## In experiment `react_2_b`, we revisit the same reaction using a better approach that employs **_adaptive variable time steps_**.

# %%
