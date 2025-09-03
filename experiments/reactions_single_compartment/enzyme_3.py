# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # **Enzyme Kinetics** : 
# #### An _accurate numerical solution_ of an enzymatic reaction `E + S <-> ES -> E + P` - compared to the **Michaelis-Menten** model approximation and to the alternative **Morrison** model, using real-life kinetic parameters. 
#
# #### Scenario with _small amount of Enzyme_, relative to the initial Substrate concentration.
#
# #### Unlike in experiment `enzyme_1_a`, we'll use data from a reaction that **VIOLATES the customary Michaelis-Menten assumption** that the rate constants satisfy `k1_reverse >> k2_forward`

# %% [markdown]
# ### THE REACTION:  
# the enzyme `Aminopeptidase` with the substrate `Leu-Ala-DED`,  
# and the initial concentration values choosen below, all satisfy the customary Michaelis-Menten assumptions that  `[E] << [S]`    BUT the reaction rate constants DON'T satisfy `k1_reverse >> k2_forward`
#
# For this reaction: k1_forward = 160 , k1_reverse = 0.089 , k2_forward = 0.58 
#
# Source of kinetic parameters:  *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% [markdown]
# ### TAGS :  "uniform compartment", "chemistry", "numerical", "enzymes"

# %%
LAST_REVISED = "Sep. 2, 2025"
LIFE123_VERSION = "1.0.0rc5"         # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import pandas as pd

from life123 import check_version, ChemData, UniformCompartment, ReactionEnzyme, PlotlyHelper

# %%
check_version(LIFE123_VERSION)    # To check compatibility

# %%

# %%

# %% [markdown]
# # PART 1. Accurate numerical solution

# %%
chem_data = ChemData(names=["P", "ES"], plot_colors=["green", "red"])

# %%
# Our Enzyme
chem_data.add_chemical(name="Aminopeptidase", label="E", plot_color="violet")

# Our Substrate
chem_data.add_chemical(name="Leu-Ala-DED", label="S", plot_color="darkturquoise")

chem_data.all_chemicals()

# %%

# %% [markdown]
# ### Specify the Kinetic Parameters

# %% [markdown]
# Source: *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% [markdown]
# ### Notice that the reaction rate constants below DON'T satisfy the customary `k1_reverse >> k2_forward`

# %%
# Here we use the "slow" preset for the variable steps, a conservative option prioritizing accuracy over speed
uc = UniformCompartment(chem_data=chem_data, preset="slow")

# %%
# Enzymatic reaction `E + S <-> ES -> E + P` 
uc.add_reaction(reactants="S", products="P", enzyme="E",
                k1_F=160., k1_R=0.089, k2_F=0.58)

uc.describe_reactions()

# %%

# %%

# %% [markdown]
# ## Set the initial concentrations of all the chemicals
# ### in the scenario we're exploring, there's LITTLE ENZYME relative to initial substrate concentration,  
# just like we did in experiment `enzyme_1_a`

# %%
S0 = 20.
E0 = 1.

# %%
uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%

# %% [markdown]
# #### Simulate the very early part of the reaction

# %%
# Perform the reactions  
# (Note: by default, concentration and rate history is kept for each step; we'll later reduce this frequency)
uc.single_compartment_react(duration=0.0015, initial_step=0.00001)

# %%
uc.plot_history(show_intervals=True, 
                title_prefix="Small amout of E relative to S(0)")

# %%
# Highlight a detail about the initial buildup of ES
uc.plot_history(title_prefix="DETAIL at early times",
                range_y=[0, 1])

# %% [markdown]
# ### Notice the EXTREMELY fast kinetics involving the production of the intermediate `ES`

# %% [markdown]
# If one wishes to exam the history (of concentration or reaction-rate data) in tabular form, it can be done as follows:

# %%
uc.get_history()

# %%
uc.get_rate_history()

# %% [markdown]
# Note that rates are associated with the START times of the intervals.  So, at the end time of the simulation, there's no rate.

# %%

# %%

# %% [markdown]
# #### Advance the reactions to equilibrium

# %%
# From now on, the data from 1 of every 10 computation steps will be saved in the history (for later use in plots, etc). 
# Until this point, we've used the default of 1
uc.enable_history(frequency=10)

# Continue the reactions
uc.single_compartment_react(duration=40., initial_step=0.00001)

# %%
uc.plot_history(title_prefix="Small amout of E relative to S(0)")

# %%
# Highlight a detail about the initial buildup of ES
uc.plot_history(title_prefix="DETAIL of ES and E",
                range_y=[0, 1.5])

# %% [markdown]
# #### Notice how the bound enzyme `ES` (red) quickly builds up at the very beginning, from time 0 to roughly 0.01 ... and that, in the longer term, the enzyme returns to its unbound state `E`

# %%

# %% [markdown]
# ### What is the initial rate of production of the final reaction product `P`?   
# One could take the numerical derivative (gradient) of the time values of [P] - but no need to!  **Reaction rates are computed in the course of the simulation, and stored in a rate-history dataframe** (as we also saw earlier)

# %%
rates = uc.get_rate_history()
rates

# %% [markdown]
# Note that initially we were saving historical data at every simulation step, but later switched to every 10 steps

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time",
                         x_var="SYSTEM TIME", fields="rxn0_rate_2", 
                         x_label="time", y_label="dP/dt")

# %% [markdown]
# ### Notice how very different is it from the Reaction rate over time seen in experiment `enzyme_1_a`

# %%
# A closer peek at its very quick buildup
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time (DETAIL at early times)",
                         x_var="SYSTEM TIME", fields="rxn0_rate_2", 
                         range_x=[0, 0.0015])

# %% [markdown]
# As we saw earlier, the time it took for `ES` to build up was about 0.001  
# **Ignoring the very brief initial transient phase, the reaction rate starts at about 0.57, and stays at that value until the substrate gets depleted, arount t=33**

# %%

# %%

# %% [markdown]
# # PART 2. Comparing the results to the Michaelis-Menten model approximation

# %% [markdown]
# ### Let's compute some parameters used by the Michaelis-Menten model   
# for background reference, see:  https://vallance.chem.ox.ac.uk/pdfs/KineticsLectureNotes.pdf (p. 20)

# %%
rxn = uc.get_single_reaction(0)

# %%
rxn.kM          #  For the data in this experiment, it comes out to (0.089 + 0.58) / 160.

# %%
rxn.kcat

# %%
vmax = rxn.compute_vmax(E_tot=E0)          # kcat * E0
vmax

# %%
initial_rxn_rate = rxn.compute_rate(S_conc=S0)    #  (vmax * S0) / (kM + S0)
initial_rxn_rate

# %% [markdown]
# Pretty close the value of about 0.57 we saw earlier, after the initial transient...   

# %%

# %% [markdown]
# ### Now, let's look at the reaction rate that produces `P` as a function of [S]; we'll compare what we computed earlier vs. what is given by the approximation of the **Michaelis-Menten model**

# %% [markdown]
# First, we'll merge the concentration history and and the rate history into a single dataframe `df`

# %%
df = uc.add_rate_to_conc_history(rate_name="rxn0_rate_2", new_rate_name="P_rate")
df

# %%
# Let's add a column with the rate estimated by the Michaelis-Menten model
df["Michaelis_rate"] = rxn.compute_rate(S_conc=df["S"])
df

# %%
# Let's see how our computed rate compares with the approximations from the Michaelis-Menten model
PlotlyHelper.plot_pandas(df=df, x_var="S", fields=["P_rate", "Michaelis_rate"],
                         title="Reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt", legend_header="Rates",
                         vertical_lines_to_add=19, colors=["blue", "yellow"])

# %% [markdown]
# Let's recall that our reactions started with [S]=20  
# Virtually overlapped plots, except at the very right! (very early times, when [S] is greater than about 19)

# %% [markdown]
# The overlap at the very left (small [S], late in the reaction) has a few issues, but not too bad; let's magnify the previous plot:

# %%
PlotlyHelper.plot_pandas(df=df, x_var="S", fields=["P_rate", "Michaelis_rate"],
                         title="EARLY DETAIL of Reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt", legend_header="Rates",
                         vertical_lines_to_add=19, colors=["blue", "yellow"],
                         range_x=[0, .2])

# %% [markdown]
# ## In this scenario, the Michaelis-Menten model is remarkable accurate EXCEPT at the very early times of the reaction (the large values of `S`), which is the brief initial transient phase when `ES` builds up from zero.  
# ### Some discrepancy also at the end of the reaction (when `S` drops to zero)

# %% [markdown]
# #### This is as expected, because the kinetic parameters of this reaction satisfy the key Michaelis-Menten assumptions that `[E] << [S]` , but NOT the other assumption that the rates satisfy `k1_reverse >> k2_forward`

# %%

# %%

# %% [markdown]
# # PART 3. Comparing the results to the Morrison model

# %% [markdown]
# #### Following section 7.1 of _"Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023_, we'll test out an the alternative **Morrison** approach, which is expected to perform better than the **Michaelis-Menten** model when the Enzyme concentration isn't so small.  But, in this scenario, we only have small amounts of enzyme - and we'll see below that no significant improvement is gained by switching model.

# %%
df["Morrison_rate"] = rxn.compute_rate_morrison(E_tot=E0,
                                                S_tot=df["S"] + df["ES"])
df

# %%
PlotlyHelper.plot_pandas(df=df, x_var="S", fields=["P_rate", "Michaelis_rate", "Morrison_rate"],
                         title="Reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt", legend_header="Rates",
                         vertical_lines_to_add=19,
                         colors=["blue", "yellow", "orange"])

# %%
# Detail at the very left
PlotlyHelper.plot_pandas(df=df, x_var="S", fields=["P_rate", "Michaelis_rate", "Morrison_rate"],
                         title="EARLY DETAIL of reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt", legend_header="Rates",
                         vertical_lines_to_add=19,
                         colors=["blue", "yellow", "orange"],
                         range_x=[0, 0.02])

# %% [markdown]
# # In this scenario, the Morrison model appears better than the Michaelis-Menten one at late reaction times (small `S`, close to zero) - and virtually identical elsewhere 

# %%
