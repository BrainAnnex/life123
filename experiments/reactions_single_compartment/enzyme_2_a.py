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
# # IN-PROGRESS
#
# ## **Enzyme Kinetics** : 
#
# #### Our model: `E + S <-> ES` (with kinetic parameters _k1_forward_ and _k1_reverse_), and  `ES -> E + P`  (_k2_forward_)  
#
# #### In experiment `enzyme_1_a`, we were given `k1_forward`, `k1_reverse` and `k2_forward`...  But what to do if we're **just given `kM` and `kcat`** ?  
#
# In Part 1, we'll "cheat" and use the actual value of _k1_forward_  
# In Part 2, we'll explore what happens if, lacking an actual value, we **under-estimate** _k1_forward_  
# In Part 3, we'll explore what happens if we **over-estimate** _k1_forward_
#
# Background: please see experiment `enzyme_1_a`

# %% [markdown]
# #### THE REACTION:  
# the enzyme `Adenosinedeaminase` with the substrate `2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine`,  
# and the initial concentration values choosen below, all satisfy the customary Michaelis-Menten assumptions that  
# `[E] << [S]` and that the reaction rate constants satisfy `k1_reverse >> k2_forward`
#
# Source of kinetic parameters:  *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% [markdown]
# ### TAGS :  "uniform compartment", "chemistry", "numerical", "enzymes"

# %%
LAST_REVISED = "Dec. 9, 2024"
LIFE123_VERSION = "1.0.0.rc.1"      # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import ipynbname
import pandas as pd

from life123 import check_version, ChemData, UniformCompartment, ReactionEnz, GraphicLog, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %%

# %%

# %% [markdown]
# ## Assume we're only given values for kM and kcat
# What values of k1_forward, k1_reverse and k2_forward are compatible with them?

# %%
# The following values are taken from experiment `enzyme_1`
kM = 8.27777777777777
kcat = 49

# %%
# k2_forward equals kcat ; not much to say here
k2_forward = kcat

# %% [markdown]
# By contrast kM = (k2_forward + k1_reverse) / k1_forward  
#
# We are given kM and k2_forward, which is the same value as the given kcat.   
#
# But how to solve for k1_forward and k1_reverse??  We have just 1 equation and 2 variables!  The system of equations is "underdetermined" : what can we do?   
#
# We'll venture some guesses for k1_forward, and compute the corresponding k1_reverse that will satisfy the above equation.   
#
# **k1_reverse = kM * k1_forward - kcat**  (kcat has same value as k2_forward)

# %%
# Example, using the k1_forward=18. from experiment `enzyme_1`
k1_forward = 18.
k1_reverse = kM * k1_forward - kcat 
k1_reverse

# %% [markdown]
# Naturally, we're getting the same value we had for k1_reverse in that experiment.  But what if k1_forward is quite different?

# %%
import numpy as np


# %%
def compute_k1_reverse(kM, kcat, k1_forward, verbose=True):
    k1_reverse = kM * k1_forward - kcat 
    if verbose:
        if np.allclose(k1_reverse, 0):
            print (f"k1_reverse: {k1_reverse} , K = INFINITE")
        else:
            K = k1_forward / k1_reverse
            print(f"k1_reverse: {k1_reverse} , K = {K}")
        
    return k1_reverse


# %%
compute_k1_reverse(kM, kcat, k1_forward=18.)

# %%
compute_k1_reverse(kM, kcat, k1_forward=0.1)

# %% [markdown]
# #### The negative values aren't good!  We must add a constraint that:   
# kM * k1_forward - kcat > 0  
# i.e. k1_forward > kcat/kM
#
# In other words, `k1_forward_min = kcat/kM`

# %%
k1_forward_min = kcat / kM
k1_forward_min

# %%
compute_k1_reverse(kM, kcat, k1_forward=k1_forward_min)

# %%
compute_k1_reverse(kM, kcat, k1_forward=6.)

# %%
compute_k1_reverse(kM, kcat, k1_forward=10.)

# %%
compute_k1_reverse(kM, kcat, k1_forward=40.)

# %% [markdown]
# The fixed point k1_forward = k1_reverse occurs when  k1_forward = kM * k1_forward - kcat , i.e.   
# kM * k1_forward - k1_forward = kcat  
# k1_forward * (kM - 1) = kcat  
# k1_forward = kcat / (kM - 1)
# Notice that the above makes sense only if kM > 1

# %%
# Our kM is indeed > 1, so all good:

k1_forward_fixed_pnt = kcat / (kM - 1)
k1_forward_fixed_pnt

# %%
compute_k1_reverse(kM, kcat, k1_forward=6.7328244274809235)

# %% [markdown]
# Indeed, a fixed point

# %%
k1_forward_choices = np.linspace(5.92, 40., 2000)    # Grid creation

# %%
k1_reverse_choices = compute_k1_reverse(kM, kcat, k1_forward=k1_forward_choices, verbose=False)

# %%
import plotly.express as px

# %%
px.line(x=k1_forward_choices, y=k1_reverse_choices, title="k1_forward vs. k1_reverse")

# %%
K_choices = k1_forward_choices / k1_reverse_choices

# %%
px.line(x=k1_forward_choices, y=np.log(K_choices), title="log(K) vs. k1_forward")

# %%

# %% [markdown]
# ## In the scenario we're exploring, there's LITTLE ENZYME relative to initial substrate concentration

# %%
S0 = 20.
E0 = 1.

# %%

# %%

# %% [markdown]
# # 1. Using the known, exact values for `k1_forward` and `k1_reverse`  
#
# Source: *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %%
k1_forward = 18.
k1_reverse = 100.

# %%
chem_data = ChemData(names=["P", "ES"])

# %%
# Our Enzyme
chem_data.add_chemical(name="Adenosinedeaminase", label="E") 

# Our Substrate
chem_data.add_chemical(name="2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine", label="S");

# %%
chem_data.all_chemicals()

# %% [markdown]
# ### Specify the Kinetic Parameters

# %% tags=[]
# Elementary Reaction E + S <-> ES
chem_data.add_reaction(reactants=["E", "S"], products=["ES"],
                       forward_rate=k1_forward, reverse_rate=k1_reverse) 

# Elementary Reaction ES <-> E + P 
chem_data.add_reaction(reactants=["ES"], products=["E", "P"],
                       forward_rate=kcat, reverse_rate=0)

chem_data.describe_reactions()

# %%

# %%
# Here we use the "slower" preset for the variable steps, a conservative option prioritizing accuracy over speed
uc = UniformCompartment(chem_data=chem_data, preset="slower")
uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%
uc.enable_diagnostics()   # To save diagnostic information about the simulation - in particular, the REACTION RATES

# %%
# Perform the reactions
uc.single_compartment_react(duration=1.2, initial_step=0.05)

# %%
uc.plot_history(colors=['cyan', 'green', 'violet', 'red'], show_intervals=True, 
                title_prefix="Using exact values for k1_forward and k1_reverse")

# %% [markdown]
# ### What is the initial rate of production of the final reaction product `P`?   
# One could take the numerical derivative (gradient) of the time values of [P] - but no need to!  **Reaction rates are computed in the course of the simulation, and stored alongside other diagnostic data**, provided that diagnostics were enabled (as we did indeed enable)

# %%
history_with_rates_exact = uc.get_diagnostics().get_system_history_with_rxn_rates(rxn_index=1)   # We specify the reaction that involves `P`
history_with_rates_exact

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=history_with_rates_exact, 
                         title="Reaction rate, dP/dt, over time",
                         x_var="TIME", fields="rate", 
                         x_label="time", y_label="dP/dt")

# %%

# %%

# %% [markdown]
# # 2. Guessing a SMALLER value of `k1_forward` (UNDER-ESTIMATED)

# %%
chem_data = ChemData(names=["P", "ES"])

# %%
# Our Enzyme
chem_data.add_chemical(name="Adenosinedeaminase", label="E") 

# Our Substrate
chem_data.add_chemical(name="2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine", label="S");

# %%
k1_forward = 6.5    # Close to the smallest possible value

# %%
k1_reverse = compute_k1_reverse(kM, kcat, k1_forward=k1_forward)
k1_reverse

# %% tags=[]
# Reaction E + S <-> ES , with 1st-order kinetics, 
# and a forward rate that is much faster than its revers one
chem_data.add_reaction(reactants=["E", "S"], products=["ES"],
                       forward_rate=k1_forward, reverse_rate=k1_reverse) 

# Reaction ES <-> E + P , with 1st-order kinetics
chem_data.add_reaction(reactants=["ES"], products=["E", "P"],
                       forward_rate=49.,reverse_rate=0)

chem_data.describe_reactions()

# %%
# Here we use the "slower" preset for the variable steps, a conservative option prioritizing accuracy over speed
uc = UniformCompartment(chem_data=chem_data, preset="slower")
uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%
uc.enable_diagnostics()   # To save diagnostic information about the simulation - in particular, the REACTION RATES

# %%
# Perform the reactions
uc.single_compartment_react(duration=1.2, initial_step=0.05)

# %%
uc.plot_history(colors=['cyan', 'green', 'violet', 'red'], show_intervals=True, 
                title_prefix="Using smaller value for k1_forward")

# %% [markdown]
# At first glance, not too different from before, overall, but let's take a closer look

# %%
history_with_rates_under = uc.get_diagnostics().get_system_history_with_rxn_rates(rxn_index=1)   # We specify the reaction that involves `P`
history_with_rates_under

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=history_with_rates_under, 
                         title="Reaction rate, dP/dt, over time",
                         x_var="TIME", fields="rate", 
                         x_label="time", y_label="dP/dt")

# %%

# %%

# %% [markdown]
# # 3. Guessing a LARGER value of `k1_forward` (OVER-ESTIMATED)

# %%
chem_data = ChemData(names=["P", "ES"])

# %%
# Our Enzyme
chem_data.add_chemical(name="Adenosinedeaminase", label="E") 

# Our Substrate
chem_data.add_chemical(name="2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine", label="S");

# %%
k1_forward = 40.    # Well above the known value of 18.

# %%
k1_reverse = compute_k1_reverse(kM, kcat, k1_forward=k1_forward)
k1_reverse

# %% tags=[]
# Reaction E + S <-> ES , with 1st-order kinetics, 
# and a forward rate that is much faster than its revers one
chem_data.add_reaction(reactants=["E", "S"], products=["ES"],
                       forward_rate=k1_forward, reverse_rate=k1_reverse) 

# Reaction ES <-> E + P , with 1st-order kinetics
chem_data.add_reaction(reactants=["ES"], products=["E", "P"],
                       forward_rate=49.,reverse_rate=0)

chem_data.describe_reactions()

# %%
# Here we use the "slower" preset for the variable steps, a conservative option prioritizing accuracy over speed
uc = UniformCompartment(chem_data=chem_data, preset="slower")
uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%
uc.enable_diagnostics()   # To save diagnostic information about the simulation - in particular, the REACTION RATES

# %%
# Perform the reactions
uc.single_compartment_react(duration=1.2, initial_step=0.05)

# %%
uc.plot_history(colors=['cyan', 'green', 'violet', 'red'], show_intervals=True, 
                title_prefix="Using smaller value for k1_forward")

# %% [markdown]
# At first glance, not too different from before, overall, but let's take a closer look

# %%
history_with_rates_over = uc.get_diagnostics().get_system_history_with_rxn_rates(rxn_index=1)   # We specify the reaction that involves `P`
history_with_rates_over

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=history_with_rates_over, 
                         title="Reaction rate, dP/dt, over time",
                         x_var="TIME", fields="rate", 
                         x_label="time", y_label="dP/dt")

# %%

# %%

# %% [markdown]
# # Comparing the 3 scenarios

# %% [markdown]
# #### The concentrations of the intermediate `ES` over time

# %%
es_under = PlotlyHelper.plot_pandas(df=history_with_rates_under, 
                         title="underestimated k1_forward",
                         x_var="TIME", fields="ES", 
                         x_label="time", y_label="ES", colors = "yellow")

# %%
es_exact = PlotlyHelper.plot_pandas(df=history_with_rates_exact, 
                         title="exact k1_forward",
                         x_var="TIME", fields="ES", 
                         x_label="time", y_label="ES", colors = "green")

# %%
es_over = PlotlyHelper.plot_pandas(df=history_with_rates_over, 
                         title="overestimated k1_forward",
                         x_var="TIME", fields="ES", 
                         x_label="time", y_label="ES", colors = "purple")

 # %%
 PlotlyHelper.combine_plots([es_under, es_exact, es_over], 
                            title="ES")

# %% [markdown]
# #### The concentrations of the product `P` over time

# %%
p_under = PlotlyHelper.plot_pandas(df=history_with_rates_under, 
                         title="underestimated k1_forward",
                         x_var="TIME", fields="P", 
                         x_label="time", y_label="P", colors = "yellow")

# %%
p_exact = PlotlyHelper.plot_pandas(df=history_with_rates_exact, 
                         title="exact k1_forward",
                         x_var="TIME", fields="P", 
                         x_label="time", y_label="P", colors = "cyan")

# %%
p_over = PlotlyHelper.plot_pandas(df=history_with_rates_over, 
                         title="overestimated k1_forward",
                         x_var="TIME", fields="P", 
                         x_label="time", y_label="P", colors = "purple")

 # %%
 PlotlyHelper.combine_plots([p_under, p_exact, p_over], 
                            title="P")

# %%

# %% [markdown]
# #### The Reaction Rate over time

# %%
r1_under = PlotlyHelper.plot_pandas(df=history_with_rates_under, 
                         title="underestimated k1_forward",
                         x_var="TIME", fields="rate", 
                         x_label="time", y_label="dP/dt", colors = "yellow")

# %%
r1_exact = PlotlyHelper.plot_pandas(df=history_with_rates_exact, 
                         title="exact k1_forward",
                         x_var="TIME", fields="rate", 
                         x_label="time", y_label="dP/dt", colors = "blue")

# %%
r1_over = PlotlyHelper.plot_pandas(df=history_with_rates_over, 
                         title="overestimated k1_forward",
                         x_var="TIME", fields="rate", 
                         x_label="time", y_label="dP/dt", colors = "purple")

 # %%
 PlotlyHelper.combine_plots([r1_under, r1_exact, r1_over], 
                            title="Reaction Rate")

# %%
