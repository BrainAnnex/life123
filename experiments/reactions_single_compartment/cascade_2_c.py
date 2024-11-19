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
# ## A complex (multistep) reaction `A <-> C` derived from 2 coupled elementary reactions:  
# ## `A <-> B` (fast) and `B <-> C` (slow)  
# A repeat of experiment `cascade_2_b`, but with the fast/slow elementary reactions in reverse order.
#
# In PART 1, a time evolution of [A], [B] and [C] is obtained by simulation  
# In PART 2, the time functions generated in Part 1 are taken as a _starting point,_ to explore how to model the composite reaction `A <-> C`   
#
# **Background**: please see experiment `cascade_2_b`  

# %%
LAST_REVISED = "Nov. 12, 2024"
LIFE123_VERSION = "1.0.0.rc.0"      # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename

from life123 import UniformCompartment, ReactionKinetics, PlotlyHelper

# %%

# %% [markdown]
# # PART 1 - We'll generate the time evolution of [A] and [C] by simulating coupled elementary reactions of KNOWN rate constants...
# ## but pretend you don't see this section! (because we later assume that those time evolutions are just GIVEN to us)

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = UniformCompartment(names=["A", "B", "C"], preset="mid")

# Reaction A <-> B (much faster, and with a much larger K)
dynamics.add_reaction(reactants="A", products="B",
                       forward_rate=80., reverse_rate=0.1)   # <===== SPEEDS of 2 reactions ARE REVERSED, relative to experiment `cascade_2_b`

# Reaction B <-> C (much slower, and with a much smaller K)
dynamics.add_reaction(reactants="B", products="C",                     
                      forward_rate=8., reverse_rate=2.)      # <===== SPEEDS of 2 reactions ARE REVERSED, relative to experiment `cascade_2_b`
                                   
dynamics.describe_reactions()

# %% [markdown]
# ### Run the simulation

# %%
dynamics.set_conc({"A": 50.}, snapshot=True)
dynamics.describe_state()

# %%
dynamics.single_compartment_react(initial_step=0.01, duration=0.8,
                                  variable_steps=True)

# %%
dynamics.plot_history(colors=['darkturquoise', 'orange', 'green'], show_intervals=True)

# %%
dynamics.is_in_equilibrium(tolerance=3)

# %% [markdown]
# #### A profoundly different plot than we got in experiment `cascade_2_b`

# %%

# %%

# %% [markdown]
# # PART 2 - This is the starting point of fitting the data from part 1.  
# ### We're given the data of the time evolution of `A` and `C`, and we want to try to model the complex reaction `A <-> C`

# %% [markdown]
# Let's start by taking stock of the actual data (saved during the simulation of part 1):

# %%
df = dynamics.get_history(columns=["SYSTEM TIME", "A", "C", "caption"])    # We're NOT given the intermediary B
df

# %% [markdown]
# ## Column B is NOT given to us.  For example, `B` might be an intermediary we can't measure.  
# #### Only [A] and [C] are given to us, on some variably-spaced time grid

# %% [markdown]
# #### Let's extract some columns, as Numpy arrays:

# %%
t_arr = df["SYSTEM TIME"].to_numpy()   # The independent variable : Time

# %%
A_conc = df["A"].to_numpy()

# %%
C_conc = df["C"].to_numpy()

# %% [markdown]
# #### If the composite reaction `A <-> C` could be modeled as an elementary reaction, we'd expect the rate of change of [C] to be proportional to [A]  
# Let's see what happens if we try to do such a linear fit!

# %%
ReactionKinetics.estimate_rate_constants_simple(t=t_arr, A_conc=A_conc, B_conc=C_conc, reactant_name="A", product_name="C")

# %% [markdown]
# ### A terrible fit!  But it looks like we'll do much better if split it into 3 portions, one where A(t) ranges from 0 to about 0.04, one from about 0.04 to 5, and one from about 5 to 50   
#
# This is a larger numer of splits than we did in experiments `cascade_2_a` and `cascade_2_b`
#
# Let's visually locate at what times those [A] values occur:

# %%
dynamics.plot_history(chemicals='A', colors='darkturquoise', range_x=[0, 0.15],
                      vertical_lines_to_add=[0.028, 0.1], title="[A] as a function of time")

# %%
dynamics.get_history(columns=["SYSTEM TIME", "A"], t_start=0.02, t_end=0.03) 

# %% [markdown]
# [A] assumes the value 5 around **t=0.028**

# %%
dynamics.get_history(columns=["SYSTEM TIME", "A"], t_start=0.09, t_end=0.11) 

# %% [markdown]
# [A] assumes the value 0.05 around **t=0.1**

# %% [markdown]
# ### Let's split the `A_conc` and `C_conc` arrays we extracted earlier (with the entire time evolution of, respectively, [A] and [C]) into 3 parts:  
# 1) points numbered 0-78  
# 2) points 79-128   
# 2) points 128-end

# %%
A_conc_early = A_conc[:79]
A_conc_mid = A_conc[79:129]
A_conc_late = A_conc[129:]

C_conc_early = C_conc[:79]
C_conc_mid = C_conc[79:129]
C_conc_late = C_conc[129:]

t_arr_early = t_arr[:79]
t_arr_mid = t_arr[79:129]
t_arr_late = t_arr[129:]

# %%

# %% [markdown]
# ## I. Let's start with the EARLY region, when t < 0.028

# %%
ReactionKinetics.estimate_rate_constants_simple(t=t_arr_early, A_conc=A_conc_early, B_conc=C_conc_early,
                                                reactant_name="A", product_name="C")

# %% [markdown]
# Just as we saw in experiment `cascade_2_b`, trying to fit an elementary reaction to that region leads to a **negative** reverse rate constant!  
# We won't discuss this part any further.

# %%

# %% [markdown]
# ## II. Let's now consider the MID region (not seen in earlier experiments in this series), when 0.028 < t < 0.1

# %%
ReactionKinetics.estimate_rate_constants_simple(t=t_arr_mid, A_conc=A_conc_mid, B_conc=C_conc_mid,
                                                reactant_name="A", product_name="C")

# %% [markdown]
# For this region, too, trying to fit an elementary reaction to that region leads to a **negative** reverse rate!  
# We won't discuss this part any further.

# %%

# %% [markdown]
# ## III. And now let's consider the LATE region, when t > 0.1

# %%
ReactionKinetics.estimate_rate_constants_simple(t=t_arr_late, A_conc=A_conc_late, B_conc=C_conc_late,
                                                reactant_name="A", product_name="C")

# %% [markdown]
# This time we have an adequate linear fit AND positive rate constants, but a **huge value of kF** relative to that of any of the elementary reactions!  
#
# Let's see the time evolution again, but just for `A` and `C`:

# %%
dynamics.plot_history(chemicals=['A', 'C'], colors=['darkturquoise',  'green'], range_x=[0, 0.25],
                      vertical_lines_to_add=[0.028, 0.1], title='A <-> C compound reaction')

# %% [markdown]
# For t > 0.1, the product [C] appears to have a lot of response to nearly non-existing changes in [A] !  
#
# That doesn't seem like a good fit to an elementary reaction...

# %% [markdown]
# #### This time, modeling the compound reaction `A <-> C` as an elementary reaction is not a good fit at any time

# %% [markdown]
# ### EXPANDED conclusion :  

# %% [markdown]
# ### While it's a well-known Chemistry notion that the slower reaction is the rate-determining step in a chain, we saw in this experiment that the complex reaction CANNOT always be modeled, not even roughly, with the rate constants of the slower reaction.  

# %% [markdown]
# We observed this in a scenario where the slower reaction is the later one.

# %%
