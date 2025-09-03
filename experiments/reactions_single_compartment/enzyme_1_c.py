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
# #### An _accurate numerical solution_ of an enzymatic reaction `E + S <-> ES -> E + P` - compared to the **Michaelis-Menten** model approximation and to the alternative **Morrison** model.
#
# #### Same as experiments `enzyme_1_a` and `enzyme_1_b` but with _a lavish amount of Enzyme_ relative to the initial Substrate concentration.  
# Experiment `enzyme_1_b` had a ten-fold increase in enzyme, relative to `enzyme_1_a`.  Here, we do _another x10 increase!_
#
# #### Real-life kinetic parameters.  
#
# Please refer to `enzyme_1_a` and `enzyme_1_b` for more details.

# %% [markdown]
# #### THE REACTION:  
# the enzyme `Adenosine deaminase` with the substrate `2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine` (same as in experiment `enzyme_1`)   
# satisfies the customary Michaelis-Menten assumptions that `k1_reverse >> k2_forward`
#
# However, the initial concentration values choosen below, do NOT satisfy the expected Michaelis-Menten assumptions that `[E] << [S]` 
#
# Source of kinetic parameters:  *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% [markdown]
# ### TAGS :  "uniform compartment", "chemistry", "numerical", "enzymes"

# %%
LAST_REVISED = "Sep. 2, 2025"
LIFE123_VERSION = "1.0.0rc5"         # Library version this experiment is based on

# %%
#import set_path                     # Using MyBinder?  Uncomment this before running the next cell!

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
chem_data.add_chemical(name="Adenosine deaminase", label="E", plot_color="violet") 

# Our Substrate
chem_data.add_chemical(name="2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine", label="S", plot_color="darkturquoise")

chem_data.all_chemicals()

# %% [markdown]
# ### Specify the Kinetic Parameters   
# This part is identical to experiments `enzyme_1` and `enzyme_2`

# %% [markdown]
# Source: *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %%
# Here we use the "slower" preset for the variable steps, a conservative option prioritizing accuracy over speed
uc = UniformCompartment(chem_data=chem_data, preset="slower")

# %%
# Enzymatic reaction `E + S <-> ES -> E + P` 
uc.add_reaction(reactants="S", products="P", enzyme="E",
                k1_F=18., k1_R=100., k2_F=49.)

uc.describe_reactions()

# %%

# %%

# %% [markdown]
# ## Set the initial concentrations of all the chemicals
# ### Unlike in experiments `enzyme_1_a` and `enzyme_1_b`, the reaction is "awash" in enzyme, relative to initial substrate concentration

# %%
S0 = 20.
E0 = 100.        # A LAVISH AMOUNT!  10 times more than in experiment `enzyme_2` - and x100 as much as experiment`enzyme_1`

# %%
uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %% [markdown]
# #### Advance the reactions to equilibrium

# %%
# Perform the reactions
uc.single_compartment_react(duration=0.05, initial_step=0.01)  # A much-shorter duration than in previous experiments in the series

# %%
uc.plot_history(show_intervals=True, 
                title_prefix="Almost as much E as S(0)")

# %% [markdown]
# The general shapes of the various curves is reminiscent of their counterparts in experiment `enzyme_1`, but notice:   
#
# 1. The product `P` is being produced on a much-faster timescale than before
# 2. The maximum amount of the bound enzyme `ES` is no longer puny; in fact, it's a good fraction of the initial [S]
# 3. The buildup of `ES`, at the very beginning, no longer looks instantaneous relatively to the main reaction: by the time `ES` reaches its maximum value, a non-trivial amount of product `P` has been produced
#
# In short, the "transient" is not so "transient"!

# %%

# %% [markdown]
# ### What is the rate of production of the final reaction product `P`?   

# %%
# We specify the reaction that involves `P`
rates = uc.add_rate_to_conc_history(rate_name="rxn0_rate_2")
rates

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time",
                         x_var="SYSTEM TIME", fields="rxn0_rate_2", 
                         x_label="time", y_label="dP/dt")

# %% [markdown]
# **The initial transient phase is no longer miniscule in duration**

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
rxn.kM          #  For the data in this experiment, it comes out to (49. + 100.) / 18.

# %%
rxn.kcat

# %%
vmax = rxn.compute_vmax(E_tot=E0)          # kcat * E0
vmax

# %%
initial_rxn_rate = rxn.compute_rate(S_conc=S0)    #  (vmax * S0) / (kM + S0)
initial_rxn_rate

# %% [markdown]
# #### Both `vmax` and the initial reaction rate are 10x what they were in experiment `enzyme_1_a`, because we started with an `E0` 10 times as large

# %%

# %% [markdown]
# ### Now, let's look at rate as a function of [S]; we'll compare what we computed earlier vs. as given by the approximation of the Michaelis-Menten model

# %%
# Let's add a column with the rate estimated by the Michaelis-Menten model
rates["Michaelis_rate"] = rxn.compute_rate(S_conc=rates["S"])
rates

# %%
# Let's see how our computed rate compares with the approximations from the Michaelis-Menten model
PlotlyHelper.plot_pandas(df=rates, x_var="S", fields=["rxn0_rate_2", "Michaelis_rate"],
                         title="Reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt", legend_header="Rates",
                         vertical_lines_to_add=1.65,
                         colors=["blue", "yellow"])

# %% [markdown]
# Let's recall that our reactions started out with [S]=20  
# The curve overlap is still passable at later times (left part of graph, when [S] drops below about 12), but is rather bad at earlier times (when [S] > 12)  
# ### We no longer have the striking overlap we had in experiment `enzyme_1_a`

# %% [markdown]
# ## In this scenario, with a lot more initial enzyme `E` than in experiment `enzyme_1_a`, the Michaelis-Menten model suffers from poor accuracy for an non-trivial early time, because the transient early phase (when `ES` builds up from zero) is no longer very brief

# %% [markdown]
# #### This is as expected, because we substantially deviated from the Michaelis-Menten assumptions that `[E] << [S]`

# %%

# %%

# %% [markdown]
# # PART 3. Comparing the results to the Morrison model

# %% [markdown]
# #### Following section 7.1 of _"Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023_ , we'll test out an the alternative **Morrison** approach, which is expected to perform better than the **Michaelis-Menten** model when the Enzyme concentration isn't so small

# %%
rates["Morrison_rate"] = rxn.compute_rate_morrison(E_tot=E0,
                                                   S_tot=rates["S"] + rates["ES"])
rates

# %%
PlotlyHelper.plot_pandas(df=rates, x_var="S", fields=["rxn0_rate_2", "Michaelis_rate", "Morrison_rate"],
                         title="Reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt", legend_header="Rates",
                         vertical_lines_to_add=1.65,
                         colors=["blue", "yellow", "orange"])

# %% [markdown]
# ## The **Morrison** model is **VASTLY better** than the Michaelis-Menten one in this *high-enzyme scenario*!

# %%
