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
# ## IN-PROGRESS
#
# ## **Enzyme Kinetics** : 
# #### _Accurate numerical solution_ - compared to the **Michaelis-Menten** model approximation and to the alternative **Morrison** model.
#
# #### Two Coupled Reactions: `E + S <-> ES`, and  `ES -> E + P` , using real-life kinetic parameters.  
# #### Scenario with _small amount of Enzyme_, relative to the initial Substrate concentration.
#
# #### Unlike in experiment `enzyme_1_a`, we'll use data from a reaction that VIOLATES the customary Michaelis-Menten assumption that the rate constants satisfy `k1_reverse >> k2_forward`

# %% [markdown]
# #### THE REACTION:  
# the enzyme `Aminopeptidase` with the substrate `Leu-Ala-DED`,  
# and the initial concentration values choosen below, all satisfy the customary Michaelis-Menten assumptions that  
# `[E] << [S]` BUT the reaction rate constants DON'T satisfy `k1_reverse >> k2_forward`
#
# For this reaction: k1_forward = 160 , k1_reverse = 0.089 , k2_forward = 0.58 
#
# Source of kinetic parameters:  *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% [markdown]
# ### TAGS :  "uniform compartment", "chemistry", "numerical", "enzymes"

# %%
LAST_REVISED = "Nov. 7, 2024"
LIFE123_VERSION = "1.0.0.rc.0"      # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import ipynbname
import pandas as pd

from life123 import check_version, ChemData, UniformCompartment, ReactionEnz, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %% tags=[]
# Initialize logging (for the system_state)
csv_log_file = ipynbname.name() + "_system_log.csv"   # Use the notebook base filename 
                                                      # IN CASE OF PROBLEMS, set manually to any desired name

# %%

# %%

# %% [markdown]
# # 1. Accurate numerical solution

# %%
chem_data = ChemData(names=["P", "ES"])

# %%
# Our Enzyme
chem_data.add_chemical(name="Aminopeptidase", label="E") 

# Our Substrate
chem_data.add_chemical(name="Leu-Ala-DED", label="S");

# %%
chem_data.all_chemicals()

# %% [markdown]
# ### Specify the Kinetic Parameters

# %% [markdown]
# Source: *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% tags=[]
# Reaction E + S <-> ES , with 1st-order kinetics, 
chem_data.add_reaction(reactants=["E", "S"], products=["ES"],
                       forward_rate=160., reverse_rate=0.089) 

# Reaction ES <-> E + P , with 1st-order kinetics, ignoring the reverse reaction
chem_data.add_reaction(reactants=["ES"], products=["E", "P"],
                       forward_rate=0.58, reverse_rate=0)  

chem_data.describe_reactions()

# %%

# %%

# %% [markdown]
# ## Set the initial concentrations of all the chemicals
# ### in the scenario we're exploring, there's LITTLE ENZYME relative to initial substrate concentration

# %%
S0 = 20.
E0 = 1.

# %%

# %%
uc = UniformCompartment(chem_data=chem_data, preset="slow")

# %%
uc.start_csv_log(csv_log_file)

# %%
uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%
uc.enable_diagnostics()   # To save diagnostic information about the simulation - in particular, the REACTION RATES

# %%

# %% [markdown] tags=[]
# #### Simulate the very early part of the reaction

# %%
# Perform the reactions
uc.single_compartment_react(duration=0.0015, initial_step=0.00001)

# %%
uc.plot_history(colors=['green', 'red', 'violet', 'darkturquoise'], show_intervals=True, 
                title_prefix="Small amout of E relative to S(0)")

# %%
# Highlight a detail about the initial buildup of ES
uc.plot_history(colors=['green', 'red', 'violet', 'darkturquoise'], title_prefix="DETAIL at early times",
                range_y=[0, 1])

# %%
uc.get_history()

# %%

# %%

# %% [markdown] tags=[]
# #### Advance the reactions to equilibrium

# %%
# Continue the reactions

try:
    uc.single_compartment_react(duration=1., initial_step=0.00001, 
                                snapshots={"frequency": 2})

except KeyboardInterrupt:
    print("\n*** KeyboardInterrupt exception caught")

# %%
print(uc.system_time)

# %%
uc.get_history()

# %%
uc.plot_history(colors=['green', 'red', 'violet', 'darkturquoise'], 
                title_prefix="Small amout of E relative to S(0)")

# %%
# Highlight a detail about the initial buildup of ES
uc.plot_history(colors=['green', 'red', 'violet', 'darkturquoise'], title_prefix="DETAIL at early times",
                range_y=[0, 1.5])

# %%
# Highlight a detail about the initial buildup of ES
uc.plot_history(colors=['green', 'red', 'violet', 'darkturquoise'], title_prefix="DETAIL at early times",
                range_x=[0, 0.001], range_y=[0, 1.5])

# %% [markdown]
# #### Notice how the bound enzyme `ES` (red) quickly builds up at the very beginning, from time 0 to roughly 0.01 ... and that, in the longer term, the enzyme returns to its unbound state `E`

# %%
# Perform the reactions
uc.single_compartment_react(duration=5., initial_step=0.001)

# %%
uc.adaptive_steps.show_adaptive_parameters()

# %%
uc.adaptive_steps.use_adaptive_preset("slow")

# %%
uc.adaptive_steps.show_adaptive_parameters()

# %%
# Perform the reactions
uc.single_compartment_react(duration=5., initial_step=0.001)

# %%

# %% [markdown]
# ### What is the initial rate of production of the final reaction product `P`?   
# One could take the numerical derivative (gradient) of the time values of [P] - but no need to!  **Reaction rates are computed in the course of the simulation, and stored alongside other diagnostic data**, provided that diagnostics were enabled (as we did indeed enable)

# %%
rates = uc.get_diagnostics().get_system_history_with_rxn_rates(rxn_index=1)   # We specify the reaction that involves `P`
rates

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time",
                         x_var="TIME", fields="rate", 
                         x_label="time", y_label="dP/dt")

# %%
# A closer peek at it maximum value
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time (DETAIL at early times)",
                         x_var="TIME", fields="rate", 
                         x_label="time", y_label="dP/dt",
                         range_x=[0,0.05], range_y=[33., 34.5])

# %% [markdown]
# As we saw earlier, the time it took for `ES` to build up was about 0.01  
# **Ignoring the brief initial transient phase, the reaction rate starts at about 34.1**

# %%

# %%

# %% [markdown]
# # 2. Comparing the results to the Michaelis-Menten model approximation

# %% [markdown]
# ### Let's compute some parameters used by the Michaelis-Menten model   
# for background reference, see:  https://vallance.chem.ox.ac.uk/pdfs/KineticsLectureNotes.pdf (p. 20)

# %%
rxn = ReactionEnz(enzyme="E", substrate="S", product="P",
                  k1_F=chem_data.get_forward_rate(0), k1_R=chem_data.get_reverse_rate(0), 
                  k2_F=chem_data.get_forward_rate(1))

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
# Not too far from the value of about 34.1 we saw earlier...   
# To keep in mind that the initial part of the reaction is affected by transients

# %%

# %% [markdown]
# ### Now, let's look at rate as a function of [S]; we'll compare what we computed earlier vs. what is given by the approximation of the Michaelis-Menten model

# %%
# Let's add a column with the rate estimated by the Michaelis-Menten model
rates["Michaelis_rate"] = rxn.compute_rate(S_conc=rates["S"])
rates

# %%
# Let's see how our computed rate compares with the approximations from the Michaelis-Menten model
PlotlyHelper.plot_pandas(df=rates, x_var="S", fields=["rate", "Michaelis_rate"],
                         title="Reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt", legend_header="Rates",
                         vertical_lines_to_add=18.9)

# %% [markdown]
# Let's recall that our reactions started with [S]=20  
# Virtually overlapped plots, except at the very right! (very early times, when [S] is greater than about 19)

# %% [markdown]
# ## In this scenario, the Michaelis-Menten model is remarkable accurate EXCEPT at the very early times of the reaction (the large values of `S`), which is the brief initial transient phase when `ES` builds up from zero

# %% [markdown]
# #### This is as expected, because the kinetic parameters of this reaction satisfy the Michaelis-Menten assumptions that `[E] << [S]` and that the rates satisfy `k1_reverse >> k2_forward`

# %%

# %%

# %% [markdown] tags=[]
# # 3. Comparing the results to the Morrison model

# %% [markdown]
# #### Following section 7.1 of _"Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023_ , we'll test out an the alternative **Morrison** approach, which is expected to perform better than the **Michaelis-Menten** model when the Enzyme concentration isn't so small.  But, in this scenario, we only have small amounts of enzyme - and we'll see below that no significant improvement is gained by switching model.

# %%
rates["Morrison_rate"] = rxn.compute_rate_morrison(E_tot=E0,
                                                   S_tot=rates["S"] + rates["ES"])
rates

# %%
PlotlyHelper.plot_pandas(df=rates, x_var="S", fields=["rate", "Michaelis_rate", "Morrison_rate"],
                         title="Reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt", legend_header="Rates",
                         vertical_lines_to_add=18.9,
                         colors=["green", "blue", "purple"])

# %% [markdown]
# ## The Morrison model appears to be only very mildly better than the Michaelis-Menten one at the early reaction times - and virtually identical later on.  
# Then why even bother with it?  In the next two experiments, `enzyme_2` and `enzyme_3`, we'll see how the Morrison model becomes progressively better than Michaelis-Menten at increasingly higher amounts of enzyme.

# %%
