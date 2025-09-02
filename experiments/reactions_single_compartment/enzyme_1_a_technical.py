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
# # **Enzyme Kinetics**  
#
# ### This is a re-visit of the first part of experiment `enzyme_1_a` 
# #### An under-the-hood look at how the enzymatic reaction is internally modeled with 2 coupled elementary reactions:  
# `E + S <-> ES` (reversible), and  `ES -> E + P` (irreversible).
#
# The results are identical, because the algorithm (including the choice of adaptive variable time steps) is the same. The only difference here is that we aren't making use of the _native support_ for enzymatic reactions. 
#
# Two other reactions being assumed negligible, and not included, are :  
# 1) the reverse direction of the 2nd reaction  
# 2) the non-catalyzed S <-> P reaction

# %% [markdown]
# ### THE REACTION:  
# the enzyme `Adenosine deaminase` with the substrate `2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine`,  
# and the initial concentration values choosen below, all satisfy the customary Michaelis-Menten assumptions that  `[E] << [S]` and that the reaction rate constants satisfy `k1_reverse >> k2_forward`  
#
# For this reaction: k1_forward = 18, k1_reverse = 100, k2_forward = 49  
#
# Source of kinetic parameters:  *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% [markdown]
# ### TAGS :  "uniform compartment", "chemistry", "numerical", "enzymes", "technical"

# %%
LAST_REVISED = "Sep. 2, 2025"
LIFE123_VERSION = "1.0.0rc5"         # Library version this experiment is based on

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
check_version(LIFE123_VERSION)    # To check compatibility

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
# # Revisit Part 1 of experiment `enzyme_1_a`, but WITHOUT making use of native support for enzymatic reactions

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

# %%
# Here we use the "slower" preset for the variable steps, a very conservative option prioritizing accuracy over speed
uc = UniformCompartment(chem_data=chem_data, preset="slower")

# %% [markdown]
# ### Here's where we start to diverge from experiment `enzyme_1_a` ...

# %%
# Reaction E + S <-> ES , with 1st-order kinetics, 
uc.add_reaction(reactants=["E", "S"], products="ES",
                kF=18., kR=100.) 

# Reaction ES <-> E + P , with 1st-order kinetics, ignoring the reverse reaction
uc.add_reaction(reactants="ES", products=["E", "P"],
                kF=49., kR=0) 

uc.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
# NOTE: 2 reactions, instead of the 1 reaction of experiment `enzyme_1_a` will now show up
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

# %% [markdown]
# ALL THE RESULTS ARE IDENTICAL TO WHAT WE HAD IN EXPERIMENT `enzyme_1_a`
#
# _(For additional analysis, see that experiment)_

# %%
