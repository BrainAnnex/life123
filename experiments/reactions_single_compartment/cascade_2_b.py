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
# ## `A <-> B` (slow) and `B <-> C` (fast)  
# A repeat of experiment `cascade_2_a`, with more DISSIMILAR elementary reactions.
#
# In PART 1, a time evolution of [A], [B] and [C] is obtained by simulation   
# In PART 2, the time functions generated in Part 1 are taken as a _starting point,_ to explore how to model the composite reaction `A <-> C`   
#
# **Background**: please see experiment `cascade_2_a`  
#
# LAST REVISED: Dec. 6, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_dynamics import ReactionDynamics
from src.modules.visualization.plotly_helper import PlotlyHelper

# %%

# %% [markdown]
# # PART 1 - We'll generate the time evolution of [A] and [C] by simulating coupled elementary reactions of KNOWN rate constants...
# ## but pretend you don't see this section! (because we later assume that those time evolutions are just GIVEN to us)

# %% tags=[]
# Instantiate the simulator and specify the chemicals
dynamics = ReactionDynamics(names=["A", "B", "C"])

# Reaction A <-> B (much slower, and with a much smaller K)
dynamics.add_reaction(reactants="A", products="B",
                      forward_rate=8., reverse_rate=2.) 

# Reaction B <-> C (much faster, and with a much larger K)
dynamics.add_reaction(reactants="B", products="C",
                      forward_rate=80., reverse_rate=0.1)   # <===== THIS IS THE KEY DIFFERENCE FROM THE EARLIER EXPERIMENT `cascade_2_a`
                                   
dynamics.describe_reactions()

# %% [markdown]
# ### Run the simulation

# %%
dynamics.set_conc([50., 0., 0.], snapshot=True)  # Set the initial concentrations of all the chemicals, in their index order
dynamics.describe_state()

# %%
# These settings can be tweaked to make the time resolution finer or coarser.  
# Here we use a "mid" heuristic: neither too fast nor too prudent
dynamics.use_adaptive_preset(preset="mid")

dynamics.single_compartment_react(initial_step=0.01, reaction_duration=0.8,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True, explain_variable_steps=False)

# %%
dynamics.plot_history(colors=['blue', 'orange', 'green'], show_intervals=True)

# %%
dynamics.is_in_equilibrium(tolerance=15)

# %% [markdown]
# We didn't quite advance to equilibrium this time.  A separate run (not shown) displayed far-better values if we had progressed the simulation just a little more in time. 

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
dynamics.estimate_rate_constants(t=t_arr, reactant_conc=A_conc, product_conc=C_conc, reactant_name="A", product_name="C")

# %% [markdown]
# ### The least-square fit is awful : the complex reaction `A <-> C` doesn't seem to be amenable to being modeled as a simple reaction with some suitable rate constants

# %% [markdown]
# ### But it looks like we'll do much better if splitting into 2 portions, one where A(t) ranges from 0 to about 40, and one from about 40 to 50   
# Indeed, revisiting the early portion of the time plot from Part 1, one can see an inflection in the [C] green curve roughly around time t=0.028, which is when [A] is around 40 (blue).  That's around the peak of the mystery intermediate B (orange).
#
# We'll pick time **t=0.028** as the divider between the 2 domains of the `A <-> C` time evolution that we want to model. 
#
# Note that this is a much smaller time than we saw in experiment `cascade_2_a`

# %%
dynamics.plot_history(colors=['blue', 'orange', 'green'], xrange=[0, 0.4], 
                      vertical_lines=[0.028])

# %% [markdown]
# #### Let's locate where the t = 0.028 point occurs in the data

# %%
dynamics.get_history(t=0.028)

# %% [markdown]
# ### Let's split the `A_conc` and `C_conc` arrays we extracted earlier (with the entire time evolution of, respectively, [A] and [C]) into two parts:  
# 1) points numbered 0-19  
# 2) points 19-end

# %%
A_conc_early = A_conc[:20]
A_conc_late = A_conc[20:]

C_conc_early = C_conc[:20]
C_conc_late = C_conc[20:]

t_arr_early = t_arr[:20]
t_arr_late = t_arr[20:]

# %% [markdown]
# ### I. Let's start with the EARLY region, when t < 0.028

# %%
dynamics.estimate_rate_constants(t=t_arr_early, reactant_conc=A_conc_early, product_conc=C_conc_early, 
                                 reactant_name="A", product_name="C")

# %% [markdown]
# Just as we saw in experiment `cascade_2_a`, trying to fit an elementary reaction to that region leads to a **negative** reverse rate constant!   
# This time, we won't discuss this part any further.

# %% [markdown]
# ### II. And now let's consider the LATE region, when t > 0.028

# %%
dynamics.estimate_rate_constants(t=t_arr_late, reactant_conc=A_conc_late, product_conc=C_conc_late, 
                                 reactant_name="A", product_name="C")

# %% [markdown]
# This time we have an adequate linear fit AND meaningful rate constants : kF of about 8 and kR of about 0.  Do those numbers sound familiar?  A definite resemblance to the kF=8, kR=2 of the SLOWER elementary reaction `A <-> B`!  
#
# #### The slower `A <-> B` reaction dominates the kinetics from about t=0.028 on    
#
# Let's see the graph again:

# %%
dynamics.plot_history(colors=['blue', 'orange', 'green'], xrange=[0, 0.4], 
                      vertical_lines=[0.028])

# %% [markdown]
# Just as we concluded in experiment `cascade_2_a`, the earlier part of the complex (compound) reaction `A <-> C` cannot be modeled by an elementary reaction, while the later part can indeed be modeled by a 1st order elementary reaction, with kinetics similar to the slower `A <-> B` reaction.  This time, with a greater disparity between the two elementary reaction, the transition happens much sooner.

# %% [markdown]
# ### While it's a well-known Chemistry notion that the slower reaction is the rate-determining step in a chain, we saw in this experiment, and in the previous one, that the complex reaction could be roughly modeled with the rate constants of the slower reaction only after some time - especially if the 2 elementary reactions are relatively similar (as in the previous experiment).  

# %% [markdown]
# If we were interested in early transients (for example, if diffusion quickly intervened), we couldn't use that model.

# %% [markdown]
# #### In the continuation experiment, `cascade_2_c`, we explore the scenario where the 2 elementary reactions are reversed in their relative speeds

# %%