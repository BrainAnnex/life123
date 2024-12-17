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
LAST_REVISED = "Dec. 15, 2024"
LIFE123_VERSION = "1.0-rc.1"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import check_version, ChemData, Reactions, UniformCompartment, ReactionEnz, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %%

# %%

# %% [markdown]
# ## Assume we're only given values for `kM` and `kcat`

# %%
# The following values are taken from experiment `enzyme_1`
kM = 8.27777777777777
kcat = 49

# %% [markdown]
# What values of k1_forward, k1_reverse and k2_forward are compatible with them?  
# That question was explored in experiment `enzyme2_a`.  Here we'll explore the dynamic consequences of selecting various combinations of values.

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
chem_data = ChemData(names=["P", "ES"], plot_colors=["green", "red"])

# %%
# Our Enzyme
chem_data.add_chemical(name="Adenosinedeaminase", label="E", plot_color="violet") 

# Our Substrate
chem_data.add_chemical(name="2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine", label="S", plot_color="darkturquoise")

chem_data.all_chemicals()

# %% [markdown]
# ### Specify the Kinetic Parameters

# %%
rxns = Reactions(chem_data=chem_data)

# %% tags=[]
# Elementary Reaction E + S <-> ES
rxns.add_reaction(reactants=["E", "S"], products=["ES"],
                  forward_rate=k1_forward, reverse_rate=k1_reverse) 

# Elementary Reaction ES <-> E + P 
rxns.add_reaction(reactants=["ES"], products=["E", "P"],
                  forward_rate=kcat, reverse_rate=0)

rxns.describe_reactions()

# %%

# %%
# Here we use the "slower" preset for the variable steps, a conservative option prioritizing accuracy over speed
uc = UniformCompartment(reactions=rxns, preset="slower")

uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%
# Perform the reactions
uc.single_compartment_react(duration=1.2, initial_step=0.05)

# %%
history_exact = uc.get_history()    # Let's save the historical data for a comparison at the very end

# %%
uc.plot_history(show_intervals=True, 
                title_prefix="Using exact values for k1_forward and k1_reverse")

# %% [markdown]
# ### What is the initial rate of production of the final reaction product `P`?   
# One could take the numerical derivative (gradient) of the time values of [P] - but no need to!  **Reaction rates are computed in the course of the simulation, and stored in a rate-history dataframe**

# %%
rates_exact = uc.get_rate_history()
rates_exact

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=rates_exact, colors="blue",
                         title="Reaction rate, dP/dt, over time",
                         x_var="SYSTEM TIME", fields="rxn1_rate", 
                         x_label="time", y_label="dP/dt")

# %%

# %%

# %% [markdown]
# # 2. Guessing a SMALLER value of `k1_forward` (UNDER-ESTIMATED)

# %%
enz = ReactionEnz()

# %%
enz.min_k1_forward(kM=kM, kcat=kcat)   # The smallest possible value

# %%
k1_forward = 6.5    # Close to the smallest possible value

# %%
k1_reverse = enz.compute_k1_reverse(kM, kcat, k1_forward=k1_forward)
k1_reverse

# %%

# %%
chem_data.all_chemicals()   # Nothing changed here; same chemicals as before

# %%

# %%
rxns_underest = Reactions(chem_data=chem_data)   # New set of reactions

# %% tags=[]
# Reaction E + S <-> ES , with 1st-order kinetics, 
# and a forward rate that is much faster than its revers one
rxns_underest.add_reaction(reactants=["E", "S"], products=["ES"],
                           forward_rate=k1_forward, reverse_rate=k1_reverse) 

# Reaction ES <-> E + P , with 1st-order kinetics
rxns_underest.add_reaction(reactants=["ES"], products=["E", "P"],
                           forward_rate=kcat, reverse_rate=0)

rxns_underest.describe_reactions()

# %%
# Here we use the "slower" preset for the variable steps, a conservative option prioritizing accuracy over speed
uc = UniformCompartment(reactions=rxns_underest, preset="slower")

uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%
# Perform the reactions
uc.single_compartment_react(duration=1.2, initial_step=0.05)

# %%
history_underest = uc.get_history()    # Let's save the historical data for a comparison at the very end

# %%
uc.plot_history(show_intervals=True, 
                title_prefix="Using SMALLER (under-estimated) value for k1_forward")

# %% [markdown]
# At first glance, not too different from before, overall, but let's take a closer look

# %%
rates_underest = uc.get_rate_history()
rates_underest

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=rates_underest, colors="blue",
                         title="Reaction rate, dP/dt, over time",
                         x_var="SYSTEM TIME", fields="rxn1_rate", 
                         x_label="time", y_label="dP/dt")

# %%

# %%

# %% [markdown]
# # 3. Guessing a LARGER value of `k1_forward` (OVER-ESTIMATED)

# %%
k1_forward = 40.    # Well above the known value of 18.

# %%
k1_reverse = enz.compute_k1_reverse(kM, kcat, k1_forward=k1_forward)
k1_reverse

# %%

# %%
chem_data.all_chemicals()   # Nothing changed here; same chemicals as before

# %%

# %%
rxns_overest = Reactions(chem_data=chem_data)   # New set of reactions

# %% tags=[]
# Reaction E + S <-> ES , with 1st-order kinetics, 
# and a forward rate that is much faster than its revers one
rxns_overest.add_reaction(reactants=["E", "S"], products=["ES"],
                          forward_rate=k1_forward, reverse_rate=k1_reverse) 

# Reaction ES <-> E + P , with 1st-order kinetics
rxns_overest.add_reaction(reactants=["ES"], products=["E", "P"],
                          forward_rate=kcat, reverse_rate=0)

rxns_overest.describe_reactions()

# %%
# Here we use the "slower" preset for the variable steps, a conservative option prioritizing accuracy over speed
uc = UniformCompartment(reactions=rxns_overest, preset="slower")

uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%
# Perform the reactions
uc.single_compartment_react(duration=1.2, initial_step=0.05)

# %%
history_overest = uc.get_history()    # Let's save the historical data for a comparison at the very end

# %%
uc.plot_history(show_intervals=True, 
                title_prefix="Using LARGER (over-estimated) value for k1_forward")

# %% [markdown]
# At first glance, not too different from before, overall, but let's take a closer look

# %%
rates_overest = uc.get_rate_history()
rates_overest

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=rates_overest, colors="blue",
                         title="Reaction rate, dP/dt, over time",
                         x_var="SYSTEM TIME", fields="rxn1_rate", 
                         x_label="time", y_label="dP/dt")

# %%

# %%
history_overest = uc.get_history()
history_overest

# %% [markdown]
# # Comparing the 3 scenarios

# %% [markdown]
# #### The concentrations of the intermediate `ES` over time

# %%
es_under = PlotlyHelper.plot_pandas(df=history_underest, 
                         title="underestimated k1_forward",
                         x_var="SYSTEM TIME", fields="ES", 
                         x_label="time", y_label="ES", colors = "yellow")

# %%
es_exact = PlotlyHelper.plot_pandas(df=history_exact, 
                         title="exact k1_forward",
                         x_var="SYSTEM TIME", fields="ES", 
                         x_label="time", y_label="ES", colors = "green")

# %%
es_over = PlotlyHelper.plot_pandas(df=history_overest, 
                         title="overestimated k1_forward",
                         x_var="SYSTEM TIME", fields="ES", 
                         x_label="time", y_label="ES", colors = "purple")

 # %%
 PlotlyHelper.combine_plots([es_under, es_exact, es_over], 
                            title="ES")

# %% [markdown]
# #### The concentrations of the product `P` over time

# %%
p_under = PlotlyHelper.plot_pandas(df=history_underest, 
                         title="underestimated k1_forward",
                         x_var="SYSTEM TIME", fields="P", 
                         x_label="time", y_label="P", colors = "yellow")

# %%
p_exact = PlotlyHelper.plot_pandas(df=history_exact, 
                         title="exact k1_forward",
                         x_var="SYSTEM TIME", fields="P", 
                         x_label="time", y_label="P", colors = "cyan")

# %%
p_over = PlotlyHelper.plot_pandas(df=history_overest, 
                         title="overestimated k1_forward",
                         x_var="SYSTEM TIME", fields="P", 
                         x_label="time", y_label="P", colors = "purple")

 # %%
 PlotlyHelper.combine_plots([p_under, p_exact, p_over], 
                            title="P")

# %%

# %% [markdown]
# #### The Reaction Rate over time

# %%
r1_under = PlotlyHelper.plot_pandas(df=rates_underest, 
                         title="underestimated k1_forward",
                         x_var="SYSTEM TIME", fields="rxn1_rate", 
                         x_label="time", y_label="dP/dt", colors = "yellow")

# %%
r1_exact = PlotlyHelper.plot_pandas(df=rates_exact, 
                         title="exact k1_forward",
                         x_var="SYSTEM TIME", fields="rxn1_rate", 
                         x_label="time", y_label="dP/dt", colors = "blue")

# %%
r1_over = PlotlyHelper.plot_pandas(df=rates_overest, 
                         title="overestimated k1_forward",
                         x_var="SYSTEM TIME", fields="rxn1_rate", 
                         x_label="time", y_label="dP/dt", colors = "purple")

 # %%
 PlotlyHelper.combine_plots([r1_under, r1_exact, r1_over], 
                            title="Reaction Rate")

# %%
