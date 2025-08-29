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
# #### _Accurate numerical solution_ -  which gets compared to the **Michaelis-Menten** model approximation and to the alternative **Morrison** model.
#
# #### Two Coupled Reactions: `E + S <-> ES`, and  `ES -> E + P` , using real-life kinetic parameters.  
# #### Scenario with **_small amount of Enzyme_**, relative to the initial Substrate concentration.
# In the follow-up experiments `enzyme_1_b` and `enzyme_1_c`, the enzyme concentration gets progressively increased.
#
# Two other reactions being assumed negligible, and not included, are :  
# 1) the reverse direction of the 2nd reaction  
# 2) the non-catalyzed E + S <-> P reaction

# %% [markdown]
# ### THE REACTION:  
# the enzyme `Adenosine deaminase` with the substrate `2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine`,  
# and the initial concentration values choosen below, all satisfy the customary Michaelis-Menten assumptions that  `[E] << [S]` and that the reaction rate constants satisfy `k1_reverse >> k2_forward`  
#
# For this reaction: k1_forward = 18, k1_reverse = 100, k2_forward = 49  
#
# Source of kinetic parameters:  *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% [markdown]
# ### TAGS :  "uniform compartment", "chemistry", "numerical", "enzymes"

# %%
LAST_REVISED = "Jan. 12, 2025"
LIFE123_VERSION = "1.0.0rc2"         # Library version this experiment is based on

# %%
#import set_path                     # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import ipynbname
import pandas as pd

from life123 import check_version, ChemData, UniformCompartment, ReactionEnzyme, GraphicLog, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %%
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file
                                            # IN CASE OF PROBLEMS, set manually to any desired name

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%

# %%

# %% [markdown]
# # 1. Accurate numerical solution

# %%
chem_data = ChemData(names=["P", "ES"], plot_colors=["green", "red"])

# %%
# Our Enzyme
chem_data.add_chemical(name="Adenosine deaminase", label="E", plot_color="violet") 

# Our Substrate
chem_data.add_chemical(name="2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine", label="S", plot_color="darkturquoise")

chem_data.all_chemicals()

# %%

# %% [markdown]
# ### Specify the Kinetic Parameters

# %% [markdown]
# Source: *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %%
# Here we use the "slower" preset for the variable steps, a very conservative option prioritizing accuracy over speed
uc = UniformCompartment(chem_data=chem_data, preset="slower")

# %%
# Reaction E + S <-> ES , with 1st-order kinetics, 
uc.add_reaction(reactants=["E", "S"], products="ES",
                forward_rate=18., reverse_rate=100.) 

# Reaction ES <-> E + P , with 1st-order kinetics, ignoring the reverse reaction
uc.add_reaction(reactants="ES", products=["E", "P"],
                forward_rate=49., reverse_rate=0) 

uc.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
uc.plot_reaction_network("vue_cytoscape_2")

# %%

# %%

# %% [markdown]
# ## Set the initial concentrations of all the chemicals
# ### in the scenario we're exploring, there's LITTLE ENZYME relative to initial substrate concentration  
# In the follow-up experiments `enzyme_1_b` and `enzyme_1_c`, the enzyme concentration gets progressively increased.

# %%
S0 = 20.
E0 = 1.

# %%
uc.set_conc(conc={"S": S0, "E": E0})      # SMALL ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %% [markdown]
# #### Advance the reactions to equilibrium

# %%
# Perform the reactions
uc.single_compartment_react(duration=1.5, initial_step=0.05)

# %%
uc.plot_history(show_intervals=True, 
                title_prefix="Small amout of E relative to S(0)")

# %%
# Highlight a detail about the initial buildup of ES
uc.plot_history(title_prefix="DETAIL at early times",
                range_x=[0, 0.5], range_y=[0, 2.])

# %% [markdown]
# #### Notice how the bound enzyme `ES` (red) quickly builds up at the very beginning, from time 0 to roughly 0.01 ... and that, in the longer term, the enzyme returns to its unbound state `E`

# %%

# %% [markdown]
# ### What is the initial rate of production of the final reaction product `P`?   
# One could take the numerical derivative (gradient) of the time values of [P] - but no need to!  **Reaction rates are computed in the course of the simulation, and stored in a rate-history dataframe**

# %%
rates = uc.get_rate_history()   # We'll be interested in rxn1_rate (the reaction that leads to `P`)
rates

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time",
                         x_var="SYSTEM TIME", fields="rxn1_rate", 
                         x_label="time", y_label="dP/dt")

# %%
# A closer peek at its maximum value
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time (DETAIL at early times)",
                         x_var="SYSTEM TIME", fields="rxn1_rate", 
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
rxn = ReactionEnzyme(enzyme="E", substrate="S", product="P",
                     k1_F=uc.get_reactions().get_forward_rate(0),
                     k1_R=uc.get_reactions().get_reverse_rate(0),
                     k2_F=uc.get_reactions().get_forward_rate(1))

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
# ### Now, let's look at the reaction rate that produces `P` as a function of [S]; we'll compare what we computed earlier vs. what is given by the approximation of the **Michaelis-Menten model**

# %% [markdown]
# First, we'll merge the concentration history and and the rate history into a single dataframe `df`

# %%
df = uc.add_rate_to_conc_history(rate_name="rxn1_rate", new_rate_name="P_rate")
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
                         vertical_lines_to_add=18.9, colors=["blue", "yellow"])

# %% [markdown]
# Let's recall that our reactions started with [S]=20  
# Virtually overlapped plots, except at the very right! (very early times, when [S] is greater than about 19)

# %% [markdown]
# ## In this scenario, the Michaelis-Menten model is remarkable accurate EXCEPT at the very early times of the reaction (the large values of `S`), which is the brief initial transient phase when `ES` builds up from zero

# %% [markdown]
# #### This is as expected, because the kinetic parameters of this reaction satisfy the Michaelis-Menten assumptions that `[E] << [S]` and that the rates satisfy `k1_reverse >> k2_forward`

# %%

# %%

# %% [markdown]
# # 3. Comparing the results to the Morrison model

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
                         vertical_lines_to_add=18.9,
                         colors=["blue", "yellow", "orange"])

# %%
# Let's take a closer look at the rightmost portion
PlotlyHelper.plot_pandas(df=df, x_var="S", fields=["P_rate", "Michaelis_rate", "Morrison_rate"],
                         title="Reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt", legend_header="Rates",
                         vertical_lines_to_add=18.9,
                         colors=["blue", "yellow", "orange"],
                         range_x=[17, 20], range_y=[32,35])

# %% [markdown]
# # In this scenario, the Morrison model appears to be only very mildly better than the Michaelis-Menten one at the early reaction times (large [S], close to the initial value of 20.) - and virtually identical everywhere else.  
# Then why even bother with it?  In the next two experiments, `enzyme_2` and `enzyme_3`, we'll see how the Morrison model becomes progressively better than Michaelis-Menten at increasingly higher amounts of enzyme.

# %%
