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
# ## A simple `A <-> B` reaction whose rate constants are to estimated from the time evolution of [A] and [B] (values given on a variable time grid.)
#
# Assume the reaction is known to be 1st order (won't verify that.)  
#
# LAST REVISED: Nov. 19, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.chemicals.chem_data import ChemData
from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.modules.visualization.plotly_helper import PlotlyHelper

import numpy as np

# %% [markdown]
# # PART 1 - We'll generate the time evolution of [A] and [B] by simulation a reaction of known rate constants...
# ## but pretend you don't see this section!

# %% tags=[]
# Specify the chemicals
chem_data = ChemData(names=["A", "B"])

# Reaction A <-> B
chem_data.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=12., reverse_rate=2.) 
 
chem_data.describe_reactions()

# %% [markdown]
# ### Run the simulation

# %%
dynamics = ReactionDynamics(chem_data=chem_data)

# %%
dynamics.set_conc([40., 10.], snapshot=True)  # Set the initial concentrations of all the chemicals, in their index order
dynamics.describe_state()

# %%
dynamics.set_diagnostics()         # To save diagnostic information about the call to single_compartment_react()

# These settings can be tweaked to make the time resolution finer or coarser.  
# The values are the current defaults, but setting them here for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
dynamics.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4)
dynamics.set_error_step_factor(0.25)

dynamics.single_compartment_react(initial_step=0.01, reaction_duration=0.5,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True, explain_variable_steps=False)

# %% [markdown]
# ### <a name="cascade_1_plot"> Plots of changes of concentration with time</a>
# Notice the variable time steps (vertical dashed lines)

# %%
dynamics.plot_history(title="Reaction A <-> B",
                      colors=['blue', 'green'], show_intervals=True)

# %%

# %% [markdown]
# # PART 2 - This is the starting point of fitting the data from part 1.  
# ## We're given the data of the above curves, and we want to estimate the rate constants (forward and reverse) of the reaction `A <-> B`

# %% [markdown]
# Let's start by taking stock of the actual data (saved during the simulation of part 1):

# %%
df = dynamics.get_history()
df

# %% [markdown]
# #### Let's extract some columns, as Numpy arrays:

# %%
t_arr = df["SYSTEM TIME"].to_numpy()   # The independent variable : Time

# %%
A_conc = df["A"].to_numpy()

# %%
B_conc = df["B"].to_numpy()

# %% [markdown]
# #### Let's verify that the stoichiometry is satified.  From `A <-> B` we can infer that any drop in [A] corresponds to an equal increase in [B].   Their sum will remain constants:

# %%
A_conc + B_conc

# %% [markdown]
# ### Just as expected!  We'll call that constant value, "TOT_conc"

# %%
TOT_conc = 50.

# %%
# Incidentally, there's a function to verify that the stoichiometry of a single reaction holds true across the entire simulation run (overkill in this case!)
dynamics.stoichiometry_checker_entire_run() 

# %%

# %% [markdown]
# ### Now, let's investigate the rates of change of [A] and [B]

# %%
# The rate of change of [A] with time
Deriv_A = np.gradient(A_conc, t_arr, edge_order=2)

# The rate of change of [B] with time
Deriv_B = np.gradient(B_conc, t_arr, edge_order=2)

# %%
# As expected from the stoichiometry, the two derivatives are opposites: when [A] increases by a certain amount, [B] decreases by that same amount
Deriv_A + Deriv_B

# %%

# %%
PlotlyHelper.plot_curves(x=t_arr, y=[Deriv_A , Deriv_B], title="d/dt A(t) and d/dt B(t) as a function of time",
                         xlabel="t", ylabel="Time derivatives", curve_labels=["A'(t)", "B'(t)"], legend_title="Derivative", colors=['aqua', 'greenyellow'])

# %% [markdown]
# ### Now, let's determine what kF and kR rate constants for `A <-> B` will yield the above data
# (Of course, we know the answer because we ran the simulation.  We'll double-check, at the very end, that we indeed get it back)

# %% [markdown]
# The rate of change of [B] is the difference of the forward rate (producing `B`) and reverse rate (consuming it):  
# B'(t) = kF * A(t) - kR * B(t)   &nbsp; &nbsp; &nbsp;  **(Eqn. 1)**  
#   
# We also know that A(t) + B(t) = TOT_conc ,  
# i.e.  A(t) = TOT_conc - B(t)  
#
# Replacing A(t) in Eqn. 1:  
# B'(t) = kF * {TOT_conc - B(t)} - kR * B(t) = kF * TOT_conc - kF * B(t) - kR * B(t)  
#
# in conclusion:   B'(t) = {kF * TOT_conc} + {-kF -kR} * B(t)     &nbsp; &nbsp; &nbsp;  **(Eqn. 2)**  

# %% [markdown]
# We already have, from our data, B'(t) as the Numpy array `Deriv_B` , and we also have B(t) as the Numpy array `B_conc`  
# Let's verify that they indeed have a linear relationship:

# %%
PlotlyHelper.plot_curves(x=B_conc, y=Deriv_B, title="d/dt B(t) as a function of B(t)",
                         xlabel="B(t)", ylabel="B'(t)", colors="green")

# %% [markdown]
# As expected, it appears to be a straight line, and the rate of change is highest when [B] is small, which happens early in the reaction.  
#
# If we fit a linear model (least-square fit straight line), we can estimate  B'(t) = a + b * B(t) , for some numbers a and b.  
# I.e. **we want to fit: Y = a + b * X , for some numbers a and b**  
# where Y is `Deriv_B` and X is `B_conc`, the Numpy arrays we computed earlier:

# %%
Y = Deriv_B    # The dependent variable
Y

# %%
X = B_conc    # The independent variable
X

# %% [markdown]
# Let's do the least-square fit:

# %%
M = np.vstack([np.ones(len(Y)), X]).T
# N is an nx2 matrix , where n is the number of data points.  The 2nd column are the values of X
M[:10, :]   # Show the first 10 rows

# %%
a, b = np.linalg.lstsq(M, Y, rcond=None)[0]  # Carry out the least-square fit
a, b

# %% [markdown]
# #### Visually verify the least-square fit

# %%
PlotlyHelper.plot_curves(x=B_conc, y=[Deriv_B , a + b*B_conc], title="d/dt B(t) as a function of B(t), alongside its least-square fit",
                         xlabel="B(t)", ylabel="B'(t)", curve_labels=["B'(t)", "Linear Fit"], legend_title="Curve vs Fit:", colors=['green', 'red'])

# %% [markdown]
# _Virtually indistinguishable lines!_

# %% [markdown]
# Replacing the estimates for a and b into Eqn. 2, above, and using the known value TOT_conc = 50.:  
#
# ```
# kF * TOT_conc = a  
# -kF - kR = b  
# ```
# i.e.
#
# ```
# kF = a / TOT_conc
# kR = -kF - b = -(a / TOT_conc) - b
# ```
# using the values: 

# %%
kF = a / TOT_conc
kF

# %%
kR = -(a / TOT_conc) - b
kR

# %% [markdown]
# ### The values we estimated for kF and kR from the data, are in good agreement with the values we used in the simulation to get that data, respectively 12 and 2 (see PART 1, above)  
# Note that our data set is quite skimpy in the number of points:

# %%
len(B_conc)

# %% [markdown]
# and that it uses a `variable` grid, with more points where there's more change, such as in the early times:

# %%
t_arr  # Time points in our data set

# %%
np.gradient(t_arr)

# %% [markdown]
# Best seen in the plot repeated from PART 1:

# %%
dynamics.plot_history(title="Reaction A <-> B",
                      colors=['blue', 'green'], show_intervals=True)

# %%
