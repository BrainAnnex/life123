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
# ## A simple `A <-> B` reaction whose rate constants are to be estimated from a given time evolution of [A] and [B]  
# ### (values given on a *variable-time* grid.)
#
# Assume the reaction is known to be 1st order (won't verify that.)  
#
# In PART 1, a time evolution of [A] and [B], with known rate constants, is obtained by simulation  
#
# In PART 2, the time evolutions generated in Part 1 are taken as a _starting point,_ to estimate the rate constants of `A <-> B`  
#
# In PART 3, we'll repeat what we did in Part 2, but this time showing the full details of how the answer is arrived at

# %% [markdown]
# ### TAGS :  "numerical", "uniform compartment", "under-the-hood"

# %%
LAST_REVISED = "Nov. 12, 2024"
LIFE123_VERSION = "1.0.0.rc.0"      # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import ipynbname
import numpy as np

from life123 import check_version, UniformCompartment, ReactionKinetics, PlotlyHelper, Numerical


# %%
check_version(LIFE123_VERSION)

# %%

# %%

# %% [markdown]
# # PART 1 - We'll generate the time evolution of [A] and [B] by simulating a reaction of KNOWN rate constants...
# ## but pretend you don't see this section! (because we later want to estimate those rate constants)

# %% tags=[]
# Instantiate the simulator and specify the accuracy preset
uc = UniformCompartment(preset="mid", enable_diagnostics=True)

# Reaction A <-> B (mostly in the forward direction)
uc.add_reaction(reactants="A", products="B",
                forward_rate=12., reverse_rate=2.) 
 
uc.describe_reactions()

# %% [markdown]
# ### Run the simulation

# %%
uc.set_conc({"A": 40., "B": 10.})  # Set the initial concentrations
uc.describe_state()

# %%
uc.single_compartment_react(initial_step=0.01, duration=0.5,
                            variable_steps=True)

# %% [markdown]
# ### <a name="cascade_1_plot"> Plots of changes of concentration with time</a>

# %%
uc.plot_history(colors=['darkturquoise', 'green'], show_intervals=True)

# %% [markdown]
# Notice the variable time steps (vertical dashed lines), more frequent when there's more change

# %%

# %%

# %% [markdown]
# # PART 2 - This is the starting point of fitting the data from part 1  
# ### We're given the data of the above curves - i.e. the system history, and we want to estimate the rate constants (forward and reverse) of the reaction `A <-> B`

# %% [markdown]
# Let's start by taking stock of the actual data (saved during the simulation of part 1):

# %%
df = uc.get_history() 
df

# %% [markdown]
# The reaction is mostly forward; the reactant `A` gets consumed, while the product `B` gets produced

# %% [markdown]
# #### Let's extract some columns, as Numpy arrays:

# %%
t_arr = df["SYSTEM TIME"].to_numpy()   # The independent variable : Time

# %%
A_conc = df["A"].to_numpy()

# %%
B_conc = df["B"].to_numpy()

# %% [markdown]
# ### **Here, we take the easy way out,** using a specialized Life123 function!
# (in Part 3, we'll do a step-by-step derivation, to see how it works)

# %%
ReactionKinetics.estimate_rate_constants_simple(t=t_arr,
                                                A_conc=A_conc, B_conc=B_conc,
                                                reactant_name="A", product_name="B")

# %% [markdown]
# ### The least-square fit is good...  and the values estimated from the data for kF and kR are in good agreement with the values we used in the simulation to get that data, respectively 12 and 2 (see PART 1, above)  

# %% [markdown]
# Note that our data set is _quite skimpy_ in the number of points:

# %%
len(B_conc)

# %% [markdown]
# and that it uses a _variable_ grid, with more points where there's more change, such as in the early times:

# %%
t_arr  # Time points in our data set

# %%
np.gradient(t_arr)   # Notice how it gets larger in later times, as bigger steps get taken

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
# #### Just as expected!

# %%
# Incidentally, there's a function to verify that the stoichiometry 
# of a single reaction holds true across the entire simulation run 
# (overkill in this case!)
uc.get_diagnostics().stoichiometry_checker_entire_run() 

# %%

# %% [markdown]
# ### Now, let's investigate the rates of change of [A] and [B]

# %%
# The rate of change of [A] with time
Deriv_A = np.gradient(A_conc, t_arr, edge_order=2)

# The rate of change of [B] with time
Deriv_B = np.gradient(B_conc, t_arr, edge_order=2)

# %%
# As expected from the stoichiometry, the two derivatives are opposites: 
# when [A] increases by a certain amount, [B] decreases by that same amount
Deriv_A + Deriv_B   # Will be very close to zero throughout

# %% tags=[]
PlotlyHelper.plot_curves(x=t_arr, y=[Deriv_A , Deriv_B], title="d/dt A(t) and d/dt B(t) as a function of time",
                         x_label="t", y_label="Time derivatives", curve_labels=["A'(t)", "B'(t)"],
                         legend_title="Derivative", colors=['aqua', 'greenyellow'])

# %% [markdown]
# The rate of changes of both [A] and [B] get smaller as the reaction marches towards equilibrium

# %%

# %% [markdown]
# ### Now, let's determine what kF and kR rate constants for `A <-> B` will yield the above data

# %% [markdown]
# Assuming that `A <-> B` is an elementary chemical reaction (i.e. occuring in a single step),  
# OR THAT IT CAN BE APPROXIMATED AS ONE, 
# then the rate of change of the reaction product concentration `B(t)` is the difference of the forward rate (producing `B`) and the reverse rate (consuming it):  
#
# `B'(t) = kF * A(t) - kR * B(t)`
#   
# We can re-write it as:   
# `B'(t) = kF * {A(t)} + kR * {- B(t)}`    &nbsp; &nbsp; &nbsp;  **(Eqn. 1)**  
#
# `A(t)`, `B(t)` are given to us; `B'(t)` is a gradient we already computed numerically; `kF` and `kR` are the rate constants that we are trying to estimate.  
#
# **If we can do a satisfactory Least Square Fit to express `B'(t)` as a linear function of `{A(t)}` and `{- B(t)}`**, as:
#
# `B'(t) = a * {A(t)} + b * {- B(t)}` , for some constants a and b,  
#
# then, comparing with Eqn. 1, it immediately follows that  
# * `kF = a`  
#
# * `kR = b`

# %% [markdown]
# Let's carry it out!  First, let's verify that `B'(t)` is indeed a linear function of `A(t)`.  
# We already have, from our data, B'(t) as the Numpy array `Deriv_B` , and we also have A(t) as the Numpy arrays `A_conc`

# %%
fig_side = PlotlyHelper.plot_curves(x=A_conc, y=Deriv_B, title="B'(t) as a function of A(t)",
                                    x_label="A(t)", y_label="B'(t)", colors="purple")
fig_side

# %% [markdown]
# As expected, it appears to be a straight line (green), and the rate of change in the product B is higher when the concentration of the reactant A is larger.  

# %% [markdown]
# Since `A(t) + B(t)` is a constant, then `B'(t)`, just shown to be a linear function of `A(t)`, will also be a linear function of `B(t)`:

# %%
PlotlyHelper.plot_curves(x=B_conc, y=Deriv_B, title="B'(t) as a function of B(t)",
                         x_label="B(t)", y_label="B'(t)", colors="gray")

# %% [markdown]
# #### Let's do the least-square fit we had set out to do: `B'(t) = a * {A(t)} + b * {- B(t)}` , for some a, b

# %%
a, b = Numerical.two_vector_least_square(V = A_conc, W = -B_conc, Y = Deriv_B)
a, b

# %% [markdown]
# #### **Voila', those are, respectively, our estimated kF and kR!**

# %% [markdown]
# #### We just obtained the same values of the estimated `kF` and `kR` as were computed by a call to `estimate_rate_constants()` in Part 2

# %% [markdown]
# #### Visually verify the least-square fit:

# %%
# Plot B'(t) and its least-square approx
fig_main = \
PlotlyHelper.plot_curves(x=t_arr, y= [Deriv_B, a * A_conc - b * B_conc],
                         title="d/dt B(t) and its least-square fit",
                         x_label="t", y_label="d/dt  B(t)",
                         colors=['green', 'red'],
                         legend_title="Curves",
                         curve_labels=["d/dt B(t) : exact", "d/dt B(t) : least-square fit"])
fig_main

# %% [markdown]
# _Virtually indistinguishable lines!  And the same plot as we saw in Part 2!_

# %%
