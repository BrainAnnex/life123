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
# ## A simple `A <-> B` reaction whose rate constants are to be estimated 
# ## from a given time evolution of [A] and [B]  
# ### (values given on a *variable-time* grid.)
#
# Assume the reaction is known to be 1st order (won't verify that.)  
#
# In PART 1, a time evolution of [A] and [B] is obtained by simulation  
# In PART 2, the time functions generated in Part 1 are taken as a _starting point,_ to estimate the rate constants of `A <-> B`  
# In PART 3, we'll repeat what we did in Part 2, but this time showing the full details of how the answer is arrived at
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta36)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from life123 import UniformCompartment
from life123.visualization.plotly_helper import PlotlyHelper

import numpy as np

# %%

# %% [markdown]
# # PART 1 - We'll generate the time evolution of [A] and [B] by simulating a reaction of KNOWN rate constants...
# ## but pretend you don't see this section! (because we later want to estimate those rate constants)

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = UniformCompartment(names=["A", "B"], preset="mid")

# Reaction A <-> B
dynamics.add_reaction(reactants="A", products="B",
                      forward_rate=12., reverse_rate=2.) 
 
dynamics.describe_reactions()

# %% [markdown]
# ### Run the simulation

# %%
dynamics.set_conc([40., 10.], snapshot=True)  # Set the initial concentrations of all the chemicals, in their index order
dynamics.describe_state()

# %%
dynamics.enable_diagnostics()         # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(initial_step=0.01, duration=0.5,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True)

# %% [markdown]
# ### <a name="cascade_1_plot"> Plots of changes of concentration with time</a>
# Notice the variable time steps (vertical dashed lines)

# %%
dynamics.plot_history(title="Reaction A <-> B",
                      colors=['darkturquoise', 'green'], show_intervals=True)

# %%

# %%

# %% [markdown]
# # PART 2 - This is the starting point of fitting the data from part 1.  
# ### We're given the data of the above curves - i.e. the system history, and we want to estimate the rate constants (forward and reverse) of the reaction `A <-> B`

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
# ### Here, we take the easy way out, using a specialized Life123 function!
# (in Part 3, we'll do a step-by-step derivation, to see how it works)

# %%
dynamics.estimate_rate_constants(t=t_arr, reactant_conc=A_conc, product_conc=B_conc, reactant_name="A", product_name="B")

# %% [markdown]
# ### The least-square fit is good...  and the values estimated from the data for kF and kR are in good agreement with the values we used in the simulation to get that data, respectively 12 and 2 (see PART 1, above)  

# %% [markdown]
# Note that our data set is quite skimpy in the number of points:

# %%
len(B_conc)

# %% [markdown]
# and that it uses a _variable_ grid, with more points where there's more change, such as in the early times:

# %%
t_arr  # Time points in our data set

# %%
np.gradient(t_arr)

# %% [markdown]
# #### The variable time grid, and the skimpy number of data points, are best seen in the plot that was shown at the end of PART 1

# %%

# %%

# %% [markdown]
# # PART 3 - investigate how the `estimate_rate_constants()` function used in part 2 works  
# #### Again, the starting point are the time evolutions of [A] and [B] , that is the system history that was given to us

# %% [markdown]
# Let's revisit the Numpy arrays that we had set up at the beginning of Part 2

# %%
t_arr    # The independent variable : Time

# %%
A_conc

# %%
B_conc

# %% [markdown]
# #### Let's verify that the stoichiometry is satified.  From the reaction `A <-> B` we can infer that any drop in [A] corresponds to an equal increase in [B].   Their sum will remain constants:

# %%
A_conc + B_conc

# %% [markdown]
# ### Just as expected!  We'll call that constant value, "TOT_conc"

# %%
TOT_conc = 50.

# %%
# Incidentally, there's a function to verify that the stoichiometry of a single reaction holds true across the entire simulation run 
# (overkill in this case!)
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

# %% tags=[]
PlotlyHelper.plot_curves(x=t_arr, y=[Deriv_A , Deriv_B], title="d/dt A(t) and d/dt B(t) as a function of time",
                         xlabel="t", ylabel="Time derivatives", curve_labels=["A'(t)", "B'(t)"], 
                         legend_title="Derivative", colors=['aqua', 'greenyellow'])

# %% [markdown]
# ### Now, let's determine what kF and kR rate constants for `A <-> B` will yield the above data

# %% [markdown]
# Assuming that A <-> B is an elementary chemical reaction (i.e. occuring in a single step)  
# OR THAT IT CAN BE APPROXIMATED AS ONE, 
# then the rate of change of the reaction product [B] is the difference of the forward rate (producing `B`) and the reverse rate (consuming it):  
#
# `B'(t) = kF * A(t) - kR * B(t)`   &nbsp; &nbsp; &nbsp;  **(Eqn. 1)**  
#   
# We also know that A(t) + B(t) = TOT_conc (a CONSTANT),  i.e.  
# `B(t) = TOT_conc - A(t)`    &nbsp; &nbsp; &nbsp;  **(Eqn. 2)**  
#
# Replacing B(t) from Eqn. 2 into Eqn. 1:
#
# `B'(t) = kF * A(t) -  kR * [TOT_conc - A(t)]`
#
# Simplifying and rearranging: 
#
# `B'(t) = kF * A(t) - kR * TOT_conc + kR * A(t)`
#
# `B'(t) = - kR * TOT_conc + kF * A(t)  + kR * A(t)`
#
# `B'(t) = [- kR * TOT_conc] + [kF + kR] * A(t)`     &nbsp; &nbsp; &nbsp;  **(Eqn. 3)** 
#
# `TOT_conc` is a known constant; `kF` and `kR` are the rate constants that we are trying to estimate.  
#
# **If we can do a satisfactory Least Square Fit to express `B'(t)` as a linear function of `A(t)`**, as:
#
# `B'(t) = a + b * A(t)` , for some constants a, b
#
# then, comparing with Eqn. 3, we get the following system of equations:
#
# * `- kR * TOT_conc = a`  
#
# * `kF + kR = b`
#
# which can be immediately solved as:
#
# * `kR = - a / TOT_conc`    &nbsp; &nbsp; &nbsp;  **(Eqn. 4)**
#
# * `kF = b - kR` 

# %% [markdown]
# Let's carry it out!  First, let's verify that `B'(t)` is indeed a linear function of `A(t)`.  
# We already have, from our data, B'(t) as the Numpy array `Deriv_B` , and we also have A(t) as the Numpy array `A_conc`  

# %%
PlotlyHelper.plot_curves(x=A_conc, y=Deriv_B, title="d/dt B(t) as a function of A(t)",
                         xlabel="A(t)", ylabel="B'(t)", colors="green")

# %% [markdown]
# As expected, it appears to be a straight line, and the rate of change in the product B is higher when the concentration of the reactant A is larger.  
#
# If we fit a linear model (least-square fit straight line), we can estimate  B'(t) = a + b * A(t) , for some numbers a and b.  
# I.e. **we want to fit: Y = a + b * X , for some numbers a and b**  
# where Y is `Deriv_B` and X is `A_conc`, the Numpy arrays we computed earlier:

# %%
Y = Deriv_B    # The dependent variable (the time-gradient of the PRODUCT concentrations)
Y

# %%
X = A_conc    # The independent variable (the REACTANT concentrations)
X

# %% [markdown]
# #### Let's do the least-square fit:

# %%
M = np.vstack([np.ones(len(Y)), X]).T
# M is an nx2 matrix , where n is the number of data points.  
# The 2nd column contains the values of X

M[:10, :]   # Show the first 10 rows

# %%
a, b = np.linalg.lstsq(M, Y, rcond=None)[0]  # Carry out the least-square fit  as: Y = a + b X
a, b

# %% [markdown]
# #### Visually verify the least-square fit

# %%
# Plot both Y and its least-square fit, as functions of X

PlotlyHelper.plot_curves(x=A_conc, y=[Deriv_B , a + b*A_conc], 
                         title="d/dt B(t) as a function of A(t), alongside its least-square fit",
                         xlabel="A(t)", ylabel="B'(t)", 
                         curve_labels=["B'(t)", "Linear Fit"], legend_title="Curve vs Fit:", colors=['green', 'red'])

# %% [markdown]
# _Virtually indistinguishable lines!  And the same plot we saw in Part 2!_

# %% [markdown]
# Finally, from equations 4, repeated here:
#
# * `kR = - a / TOT_conc`
#
# * `kF = b - kR` 

# %%
kR = - a / TOT_conc
kR

# %%
kF = b - kR
kF

# %% [markdown]
# #### We just obtained the same values of the estimated kF and kR as were computed by a call to `estimate_rate_constants()` in Part 2

# %%
