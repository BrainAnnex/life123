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
# ## A simple `A + B <-> C` reaction whose rate constants are to be estimated from a given time evolution of [A] and [B]  
# ### (values given on a *variable-time* grid.)
#
# Assume the reaction is known to be 1st order (won't verify that.)  
#
# This is the counterpart of experiment `mystery_reaction_1` for a more complex reaction; we'll also proceed faster.  
#
# In PART 1, a time evolution of [A] and [B] is obtained by simulation  
#
# In PART 2, the time evolutions generated in Part 1 are taken as a _starting point,_ to estimate the rate constants of `A <-> B`  
#
# In PART 3, we'll repeat what we did in Part 2, but this time showing the full details of how the answer is arrived at

# %% [markdown]
# ### TAGS :  "numerical", "uniform compartment", "under-the-hood"

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

import ipynbname
import numpy as np

from life123 import check_version, UniformCompartment, PlotlyHelper, Numerical


# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# # PART 1 - We'll generate the time evolution of [A] and [B] by simulating a reaction of KNOWN rate constants...
# ## but pretend you don't see this section! (because we later want to estimate those rate constants)

# %% tags=[]
# Instantiate the simulator and specify the accuracy preset
dynamics = UniformCompartment(preset="mid")

# Reaction A <-> B (mostly in the forward direction)
dynamics.add_reaction(reactants=["A", "B"], products="C",
                      forward_rate=5., reverse_rate=1.) 
 
dynamics.describe_reactions()

# %% [markdown]
# ### Run the simulation

# %%
dynamics.set_conc({"A": 40., "B": 20.,  "C": 5.}, snapshot=True)  # Set the initial concentrations
dynamics.describe_state()

# %%
dynamics.enable_diagnostics()         # To save diagnostic information for the simulation run, below

dynamics.single_compartment_react(initial_step=0.001, duration=0.05,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True)

# %% [markdown]
# ### <a name="cascade_1_plot"> Plots of changes of concentration with time</a>

# %%
dynamics.plot_history(colors=['darkturquoise', 'violet', 'green'], show_intervals=True)

# %% [markdown]
# Notice the variable time steps (vertical dashed lines), more frequent when there's more change

# %%

# %%

# %% [markdown]
# # PART 2 - This is the starting point of fitting the data from part 1  
# ### We're given the data of the above curves - i.e. the system history, and we want to estimate the rate constants (forward and reverse) of the reaction `A + B <-> C`

# %% [markdown]
# Let's start by taking stock of the actual data (saved during the simulation of part 1):

# %%
df = dynamics.get_history() 
df

# %% [markdown]
# The reaction is mostly forward; the reactants `A` and `B` gets consumed, while the product `C` gets produced

# %% [markdown]
# #### Let's extract some columns, as Numpy arrays:

# %%
t_arr = df["SYSTEM TIME"].to_numpy()   # The independent variable : Time

# %%
A_conc = df["A"].to_numpy()

# %%
B_conc = df["B"].to_numpy()

# %%
C_conc = df["C"].to_numpy()

# %%

# %% [markdown]
# ### **Here, we take the easy way out,** using a specialized Life123 function!
# (in Part 3, we'll do a step-by-step derivation, to see how it works)

# %%
dynamics.estimate_rate_constants_synthesis(t=t_arr,
                                           A_conc=A_conc, B_conc=B_conc, C_conc=C_conc,
                                           reactants=["A", "B"], product="C")

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
# #### Again, the starting point are the time evolutions of [A], [B] and [C] , that is the system history that was given to us

# %% [markdown]
# Let's revisit the Numpy arrays that we had set up at the beginning of Part 2

# %%
t_arr    # The independent variable : Time

# %%
A_conc

# %%
B_conc

# %%
C_conc

# %% [markdown]
# #### Let's verify that the stoichiometry is satified.  From the reaction `A + B <-> C` we can infer that and any drop in [A] corresponds to an equal drop in [B] , and an equal increase in [C]  

# %%
A_conc - B_conc    # Constant differences will show that they change in lockstep

# %%
A_conc + C_conc    # Constant sums will show that any drop in [A] corresponds to an equal increase in [C]

# %%

# %%
# Incidentally, there's a function to verify that the stoichiometry 
# of a single reaction holds true across the entire simulation run 
# (overkill in this case!)
dynamics.diagnostics.stoichiometry_checker_entire_run() 

# %%

# %% [markdown]
# ### Now, let's investigate the rates of change of [A], [B] and [C]

# %%
# The rate of change of [A] with time
Deriv_A = np.gradient(A_conc, t_arr, edge_order=2)

# The rate of change of [B] with time
Deriv_B = np.gradient(B_conc, t_arr, edge_order=2)

# The rate of change of [C] with time
Deriv_C = np.gradient(C_conc, t_arr, edge_order=2)

# %%
# As expected from the stoichiometry, the two derivatives are opposites: 
# when [A] increases by a certain amount, [C] decreases by that same amount
Deriv_A + Deriv_C   # Will be very close to zero throughout

# %%
# Similarly:
Deriv_B + Deriv_C   # Will be very close to zero throughout

# %%
# As expected from the stoichiometry, [A] and [B] vary in unison in time
Deriv_A - Deriv_B   # Will be very close to zero throughout

# %% tags=[]
PlotlyHelper.plot_curves(x=t_arr, y=[Deriv_A , Deriv_C], title="d/dt A(t) and d/dt C(t) as a function of time",
                         xlabel="t", ylabel="Time derivatives", curve_labels=["A'(t)", "C'(t)"], 
                         legend_title="Derivative", colors=['aqua', 'greenyellow'])

# %% [markdown]
# The rate of changes of both [A] and [C] get smaller as the reaction marches towards equilibrium.  
# B'(t) not shown, because essentially identical to A'(t)

# %%

# %% [markdown]
# ### Now, let's determine what kF and kR rate constants for `A + B <-> C` will yield the above data

# %% [markdown]
# Assuming that A + B <-> C is an elementary chemical reaction (i.e. occuring in a single step)  
# OR THAT IT CAN BE APPROXIMATED AS ONE, 
# then the rate of change of the reaction product [C] is the difference of the forward rate (producing `C`) and the reverse rate (consuming it):  
#
# `C'(t) = kF * A(t) * B(t) - kR * C(t)`
#   
# We can re-write it as:.  
# `C'(t) = kF * {A(t) * B(t)} + kR * {- C(t)}`    &nbsp; &nbsp; &nbsp;  **(Eqn. 1)**  
#
# `A(t)` and `B(t)` are given to us; `B'(t)` is a gradient we already computed numerically; `kF` and `kR` are the rate constants that we are trying to estimate.  
#
# **If we can do a satisfactory Least Square Fit to express `C'(t)` as a linear function of `{A(t) * B(t)}` and `{- C(t)}`**, as:
#
# `C'(t) = a * {A(t) * B(t)} + b * {- C(t)}` , for some constants a and b,  
#
# then, comparing with Eqn. 1, it immediately follows that  
# * `kF = a`  
#
# * `kR = b`  
#

# %% [markdown]
# #### Let's carry out the least-square fit:

# %%
a, b = Numerical.two_vector_least_square(V = A_conc * B_conc, W = - C_conc, Y = Deriv_C)
a, b

# %% [markdown]
# #### **Voila', those are, respectively, our estimated kF and kR!**

# %% [markdown]
# #### We just obtained the same values of the estimated kF and kR as were computed by a call to `estimate_rate_constants()` in Part 2

# %%
