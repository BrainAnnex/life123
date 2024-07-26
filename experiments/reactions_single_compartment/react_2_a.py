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
# Afterwards, perform some analysis of the results
#
# See also the experiment _"1D/reactions/reaction_1"_ for a multi-compartment version.  
#
# #### This experiment gets repeated in _"react_2_b"_ , with a more sophisticated approach, 
# #### involving adaptive variable time steps.

# %%
LAST_REVISED = "July 24, 2024"
LIFE123_VERSION = "1.0.0.beta.37"    # Version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

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
# # Initialize the System

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = UniformCompartment()

# Reaction A <-> B , with 1st-order kinetics in both directions
dynamics.add_reaction(reactants="A", products="B", 
                      forward_rate=3., reverse_rate=2.)

print("Number of reactions: ", dynamics.number_of_reactions())

# %%
dynamics.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
dynamics.plot_reaction_network("vue_cytoscape_2")


# %%
# Initial concentrations of all the chemicals
dynamics.set_conc({"A": 10., "B": 50.})

# %%
dynamics.describe_state()

# %%
dynamics.get_history()

# %% [markdown]
# ### Test your intuition: 
# #### given that this reaction operates mostly in the forward direction (kF = 3 , kR = 2 , K = 1.5), 
# #### do you think that A will be consumed and B will be produced??
# We can take a sneak preview at the final equilibrium concentrations without actually running the simulation:

# %%
dynamics.find_equilibrium_conc(rxn_index=0)    # This is an EXACT solution

# %% [markdown]
# #### The reaction will actually proceed IN REVERSE, because of the large initial concentration of B (which we had set to 50), relative to the small initial concentration of A (10)
# Now, let's see the reaction in action!

# %%

# %% [markdown] tags=[]
# ### Run the reaction

# %%
# First step of reaction
dynamics.single_compartment_react(initial_step=0.1, n_steps=1, variable_steps=False, 
                                  snapshots={"initial_caption": "first reaction step"})

# %%
dynamics.get_history()

# %% [markdown]
# We can already see the reaction proceeding in reverse...

# %%
# Numerous more fixed steps
dynamics.single_compartment_react(initial_step=0.1, n_steps=10, variable_steps=False, 
                                  snapshots={"initial_caption": "2nd reaction step",
                                             "final_caption": "last reaction step"})

# %%
dynamics.get_history()

# %% [markdown]
# ## NOTE: for demonstration purposes, we're using FIXED time steps...  
# ## Typically, one would use the option for adaptive variable time steps (see experiment `react_2_b`)

# %% [markdown]
# ### Check the final equilibrium

# %%
dynamics.get_system_conc()   # The current concentrations, in the order the chemicals were added 

# %% [markdown]
# NOTE: Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order of the reactions), the systems settles in the following equilibrium:
#
# [A] = 23.99316406
#  
# [B] = 36.00683594
#

# %%
# Verify that the reaction has reached equilibrium
dynamics.is_in_equilibrium()

# %% [markdown]
# ### As noted earlier, because of the high initial concentration of B relative to A, the overall reaction has proceeded IN REVERSE

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_history(colors=['darkturquoise', 'orange'])

# %% [markdown]
# ### Note the raggedness of the left-side (early times) of the curves.  
# ### In experiment `react_2_b` this simulation gets repeated with an _adaptive variable time resolution_ that takes smaller steps at the beginning, when the reaction is proceeding faster   
# ### By contrast, here we used _FIXED_ time steps (shown below), which generally gives poor results, unless taking a very large number of very small steps!

# %%
dynamics.plot_history(colors=['darkturquoise', 'orange'], show_intervals=True)

# %%
df = dynamics.get_history()         # Revisited from earlier
df

# %%

# %% [markdown]
# ## PART 2 - Now investigate A_dot, i.e. d[A]/dt

# %% [markdown]
# NOTE: there's actually no need to compute this; it can be automatically saved during the reaction, as demonstrated in experiment `react_2_b`

# %%
A = list(df.A)
A

# %%
len(A)

# %%
A_dot = np.gradient(A, 0.1)      # 0.1 is the constant step size

# %%
A_dot

# %%
df['A_dot'] = A_dot     # Add a column to the Pandas dataframe

# %%
df

# %%
dynamics.plot_history(chemicals=["A", "A_dot"], colors=['darkturquoise', 'brown'], 
                      ylabel="concentration (darkturquoise) /<br> concentration change per unit time (brown)",
                      title="Concentration of A with time (darkturquoise), and its rate of change (brown)")

# %% [markdown]
# ### At t=0 :  
# [A]=10 and [A] has a high rate of change (70)
# ### As the system approaches equilibrium :  
# [A] approaches a value of 24, and its rate of change decays to zero.

# %% [markdown]
# #### **NOTE:** The curves are jagged because of limitations of numerically estimating derivatives, as well as _the large time steps taken_ (especially in the early times, when there's a lot of change.)  
# ## In experiment "react_2_b", we revisit the same reaction using a better approach that employs **_adaptive variable time steps_** , and also automatically saves the reaction rates.

# %%
