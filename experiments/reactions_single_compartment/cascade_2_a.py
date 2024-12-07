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
# ## A complex (composite/multistep) reaction `A <-> C` derived from 2 coupled elementary reactions: 
# ## `A <-> B` and `B <-> C`  
# We are given the time evolution of the complex reaction,  
# and want to determine whether it can be modeled as an elementary reaction.  
#
# In PART 1, a time evolution of [A], [B] and [C] is obtained by simulation  
# In PART 2, the time functions generated in Part 1 are taken as a _starting point,_ to explore how to model the composite reaction `A <-> C`  
#
# **Background**: please see experiments `cascade_1` and `mystery_reaction_1`

# %%
LAST_REVISED = "Dec. 6, 2024"
LIFE123_VERSION = "1.0-rc.1"     # Library version this experiment is based on

# %%
#import set_path                 # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import check_version, ChemData, UniformCompartment, ReactionKinetics, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %%

# %%

# %% [markdown]
# # PART 1 - We'll generate the time evolution of [A] and [C] by simulating coupled elementary reactions of KNOWN rate constants...
# ## but pretend you don't see this section! (because we later assume that those time evolutions are just GIVEN to us)

# %% tags=[]
# Instantiate the simulator and specify the chemicals
chem_data = ChemData()
#names=["A", "B", "C"], plot_colors=["darkturquoise", "orange", "green"]

chem_data.add_chemical(name="A", plot_color="darkturquoise")
chem_data.add_chemical(name="B", plot_color="orange")
chem_data.add_chemical(name="C", plot_color="green");

# %% tags=[]
# Here we use the "mid" preset for the variable steps, a compromise between speed and accuracy
uc = UniformCompartment(chem_data=chem_data, preset="mid")

# Reaction A <-> B (slower, and with a smaller K)
uc.add_reaction(reactants="A", products="B",
                forward_rate=8., reverse_rate=2.) 

# Reaction B <-> C (faster, and with a larger K)
uc.add_reaction(reactants="B", products="C",
                forward_rate=12., reverse_rate=1.) 
                                   
uc.describe_reactions()

# %% [markdown]
# ### Run the simulation

# %%
uc.set_conc({"A": 50.})  # Set the initial concentrations
uc.describe_state()

# %%
uc.single_compartment_react(initial_step=0.01, duration=0.8)

# %%
uc.plot_history(show_intervals=True)

# %%
uc.is_in_equilibrium(tolerance=12)

# %%

# %%

# %% [markdown]
# # PART 2 - This is the starting point of fitting the data from part 1.  
# ### We're given the data of the time evolution of `A` and `C`, and we want to try to model the complex reaction `A <-> C`

# %% [markdown]
# Let's start by taking stock of the actual data (saved during the simulation of part 1):

# %%
# For this analysis, we're NOT given the intermediary B
df = uc.get_history(columns=["SYSTEM TIME", "A", "C", "caption"])
df

# %% [markdown]
# ## Column B is NOT given to us.  For example, `B` might be an intermediary we can't measure.  
# #### Only [A] and [C] are given to us, on a variably-spaced time grid

# %% [markdown]
# #### Let's extract some columns, as Numpy arrays:

# %%
t_arr = df["SYSTEM TIME"].to_numpy()   # The independent variable : Time

# %%
A_conc = df["A"].to_numpy()

# %%
C_conc = df["C"].to_numpy()

# %% [markdown]
# #### If the composite reaction `A <-> C` could be modeled as an elementary reaction, we'd expect the rate of change of [C] to be `C'(t) = kF * A(t) - kR * C(t)` , for some values `kF` and `kR`  
# Let's see what happens if we try a fit to determine values for `kF` and `kR`!

# %%
ReactionKinetics.estimate_rate_constants_simple(t=t_arr, A_conc=A_conc, B_conc=C_conc, reactant_name="A", product_name="C")

# %% [markdown]
# ### The least-square fit is awful : the complex reaction `A <-> C` doesn't seem to be amenable to being modeled as an elementary reaction with some suitable rate constants
# Probably not too surprising given our "secret" knowledge from Part 1 that the complex reaction originates from 2 elementary reactions where the intermediate product builds up at one point

# %% [markdown]
# ### A glance at the above diagram reveals much-better linear fits, if split into 2 portions, one where A(t) ranges from 0 to about 24, and one from about 24 to 50   
# Indeed, revisiting the early portion of the time plot from Part 1, one can see an inflection in the [C] green curve roughly around time t=0.1, which is when [A] is around 24 (turquoise).  That's shortly after the peak of the mystery intermediate B (orange).    
#
# We'll pick time **t=0.1** as the divider between the 2 domains of the `A <-> C` time evolution that we want to model. 

# %%
uc.plot_history(range_x=[0, 0.4],
                vertical_lines_to_add=[0.1])

# %% [markdown]
# #### Let's locate where the t = 0.1 point occurs in the data

# %%
uc.get_history(t=0.1)

# %% [markdown]
# ### Let's split the `A_conc` and `C_conc` arrays we extracted earlier (with the entire time evolution of, respectively, [A] and [C]) into two parts:  
# 1) points numbered 0 thru 47   
# 2) points 48 - end

# %%
A_conc_early = A_conc[:48]
A_conc_early

# %%
A_conc_late = A_conc[48:]
A_conc_late

# %%
len(A_conc_early) + len(A_conc_late) - len(A_conc)    # Double-check we got all points

# %%
# Same for [C] and for t_arr
C_conc_early = C_conc[:48]
C_conc_late = C_conc[48:]

t_arr_early = t_arr[:48]
t_arr_late = t_arr[48:]

# %%

# %% [markdown]
# ## I. Let's start with the EARLY region, when t < 0.1

# %%
ReactionKinetics.estimate_rate_constants_simple(t=t_arr_early, A_conc=A_conc_early, B_conc=C_conc_early,
                                        reactant_name="A", product_name="C")

# %% [markdown]
# While the linear fits is a lot better than before - nonetheless trying to fit an elementary reaction to that region leads to a **negative** reverse rate constant!  
# It's no surprise that an elementary reaction is NOT a good fit, if one observes what happens to the time evolution of the concentrations.  Repeating the earlier plot, but only showing `A` and `C` (i.e. hiding the intermediary `B`):

# %%
uc.plot_history(range_x=[0, 0.4], vertical_lines_to_add=[0.1],
                chemicals=['A', 'C'], title="Changes in concentration for `A <-> C`")

# %% [markdown]
# In the zone to the left of the vertical dashed line:  
# when the reactant `A` is plentiful, the rate of change (slope of the greeen line) of the product `C` is low - and vice versa.  
# Does that look like an elementary reaction in its kinetics?  Nope!

# %%

# %% [markdown]
# ## II. And now let's consider the LATE region, when t > 0.1

# %%
ReactionKinetics.estimate_rate_constants_simple(t=t_arr_late, A_conc=A_conc_late, B_conc=C_conc_late,
                                                reactant_name="A", product_name="C")

# %% [markdown]
# This time we have an adequate linear fit AND meaningful rate constants : kF of about 8 and kR of about 0.  Do those numbers sound familiar?  A definite resemblance to the kF=8, kR=2 of the SLOWER elementary reaction `A <-> B`!  
#
# #### The slower `A <-> B` reaction dominates the kinetics, from about t=0.1 on  
#
# Let's see the graph again:

# %%
uc.plot_history(range_x=[0, 0.4],
                vertical_lines_to_add=[0.1])

# %% [markdown]
# A possible conclusion to draw is that, in this case:  
# 1) the earlier part of the complex (compound) reaction `A <-> C` cannot be modeled by an elementary reaction
# 2) while the later part can indeed be modeled by a 1st order elementary reaction, with kinetics similar to the slower `A <-> B` reaction

# %% [markdown]
# **The intuition:** imagine `A <-> B <-> C` as a supply line.  
# `A <-> B` is slow, but the moment something arrives in B, it's very quickly moved to C.  
# The slow link (`A <-> B`) largely determines the kinetics of the supply line.

# %% [markdown]
# ### While it's a well-known Chemistry notion that the slower reaction is the rate-determining step in a chain, we saw in this experiment  that **the complex reaction could be roughly modeled with the rate constants of the slower reaction ONLY AFTER SOME TIME**.  

# %% [markdown]
# If we were interested in early transients (for example, if diffusion quickly intervened), we couldn't use that model.

# %% [markdown]
# #### Is that surprising?  At early times, compare the inflection of the final product, `C`, of the composite reaction vs. the inflection of the product of a simple reaction (such as B in experiment `react1`, both appearing in green.)

# %%

# %% [markdown]
# #### NEXT: in the continuation experiment, `cascade_2_b`, we explore the scenario where the 2 elementary reactions are much more different from each other

# %%
