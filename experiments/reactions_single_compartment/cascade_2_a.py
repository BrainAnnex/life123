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
# ## `A <-> B` and `B <-> C`  
# We are given the time evolution of the complex reaction,  
# and want to determine whether it can be modeled as an elementary reaction.  
#
# In PART 1, a time evolution of [A], [B] and [C] is obtained by simulation  
# In PART 2, the time functions generated in Part 1 are taken as a _starting point,_ to explore how to model the composite reaction `A <-> C`  
#
# **Background**: please see experiments `cascade_1` and `mystery_reaction_1`
#
# LAST REVISED: June 14, 2024 (using v. 1.0 beta33)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.uniform_compartment import UniformCompartment
from src.modules.visualization.plotly_helper import PlotlyHelper

# %%

# %%

# %% [markdown]
# # PART 1 - We'll generate the time evolution of [A] and [C] by simulating coupled elementary reactions of KNOWN rate constants...
# ## but pretend you don't see this section! (because we later assume that those time evolutions are just GIVEN to us)

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = UniformCompartment(names=["A", "B", "C"], preset="mid")

# Reaction A <-> B (slower, and with a smaller K)
dynamics.add_reaction(reactants="A", products="B",
                      forward_rate=8., reverse_rate=2.) 

# Reaction B <-> C (faster, and with a larger K)
dynamics.add_reaction(reactants="B", products="C",
                      forward_rate=12., reverse_rate=1.) 
                                   
dynamics.describe_reactions()

# %% [markdown]
# ### Run the simulation

# %%
dynamics.set_conc({"A": 50.}, snapshot=True)  # Set the initial concentrations of all the chemicals, in their index order
dynamics.describe_state()

# %%
dynamics.single_compartment_react(initial_step=0.01, duration=0.8,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True)

# %%
dynamics.plot_history(colors=['darkturquoise', 'orange', 'green'], show_intervals=True)

# %%
dynamics.is_in_equilibrium(tolerance=15)

# %%

# %%

# %% [markdown]
# # PART 2 - This is the starting point of fitting the data from part 1.  
# ### We're given the data of the time evolution of `A` and `C`, and we want to try to model the complex reaction `A <-> C`

# %% [markdown]
# Let's start by taking stock of the actual data (saved during the simulation of part 1):

# %%
df = dynamics.get_history(columns=["SYSTEM TIME", "A", "C", "caption"])   # We're NOT given the intermediary B
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
dynamics.estimate_rate_constants(t=t_arr, reactant_conc=A_conc, product_conc=C_conc, reactant_name="A", product_name="C")

# %% [markdown]
# ### The least-square fit is awful : the complex reaction `A <-> C` doesn't seem to be amenable to being modeled as a simple reaction with some suitable rate constants
# Probably not too surprising given our "secret" knowledge from Part 1 that the complex reaction originates from 2 elementary reactions where one doesn't dominate the other one in terms of reaction kinetics

# %% [markdown]
# ### A glance at the above diagram reveals much-better linear fits, if split into 2 portions, one where A(t) ranges from 0 to about 24, and one from about 24 to 50   
# Indeed, revisiting the early portion of the time plot from Part 1, one can see an inflection in the [C] green curve roughly around time t=0.1, which is when [A] is around 24 (blue).  That's shortly after the peak of the mystery intermediate B (orange).    
#
# We'll pick time **t=0.1** as the divider between the 2 domains of the `A <-> C` time evolution that we want to model. 

# %%
dynamics.plot_history(colors=['darkturquoise', 'orange', 'green'], xrange=[0, 0.4], 
                      vertical_lines=[0.1])

# %% [markdown]
# #### Let's locate where the t = 0.1 point occurs in the data

# %%
dynamics.get_history(t=0.1)

# %% [markdown]
# ### Let's split the `A_conc` and `C_conc` arrays we extracted earlier (with the entire time evolution of, respectively, [A] and [C]) into two parts:  
# 1) points numbered 0-47   
# 2) points 48-end

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

# %% [markdown]
# ### I. Let's start with the EARLY region, when t < 0.1

# %%
dynamics.estimate_rate_constants(t=t_arr_early, reactant_conc=A_conc_early, product_conc=C_conc_early, 
                                 reactant_name="A", product_name="C")

# %% [markdown]
# Trying to fit an elementary reaction to that region leads to a **negative** reverse rate constant!  
# It's not surprise that an elementary reaction is a good fit, if one observes what happens to the time evolution of the concentrations.  Repeating the earlier plot, but only showing `A` and `C` (i.e. hiding the intermediary `B`):

# %%
dynamics.plot_history(colors=['darkturquoise', 'green'], xrange=[0, 0.4], vertical_lines=[0.1], 
                     chemicals=['A', 'C'], title="Changes in concentration for `A <-> C`")

# %% [markdown]
# In the zone to the left of the vertical dashed line:  
# when the reactant `A` is plentiful, the rate of change (gradient) of the product `C` is low - and vice versa.  
# Does that look like an elementary reaction in its kinetics?  Nope!

# %% [markdown]
# ### II. And now let's consider the LATE region, when t > 0.1

# %%
dynamics.estimate_rate_constants(t=t_arr_late, reactant_conc=A_conc_late, product_conc=C_conc_late, 
                                 reactant_name="A", product_name="C")

# %% [markdown]
# This time we have an adequate linear fit AND meaningful rate constants : kF of about 8 and kR of about 0.  Do those numbers sound familiar?  A definite resemblance to the kF=8, kR=2 of the SLOWER elementary reaction `A <-> B`!  
#
# #### The slower `A <-> B` reaction dominates the kinetics from about t=0.1  
#
# Let's see the graph again:

# %%
dynamics.plot_history(colors=['darkturquoise', 'orange', 'green'], xrange=[0, 0.4], 
                      vertical_lines=[0.1])

# %% [markdown]
# A possible conclusion to draw is that, in this case, the earlier part of the complex (compound) reaction `A <-> C` cannot be modeled by an elementary reaction, while the later part can indeed be modeled by a 1st order elementary reaction, with kinetics similar to the slower `A <-> B` reaction

# %% [markdown]
# **The intuition:** imagine `A <-> B <-> C` as a supply line.  
# `A <-> B` is slow, but the moment something arrives in B, it's very quickly moved to C.  
# The slow link (`A <-> B`) largely determines the kinetics of the supply line.

# %% [markdown]
# ### While it's a well-known Chemistry notion that the slower reaction is the rate-determining step in a chain, we saw in this experiment  that the complex reaction could be roughly modeled with the rate constants of the slower reaction only after some time.  

# %% [markdown]
# If we were interested in early transients (for example, if diffusion quickly intervened), we couldn't use that model.

# %% [markdown]
# #### In the continuation experiment, `cascade_2_b`, we explore the scenario where the 2 elementary reactions are much more different from each other

# %%
