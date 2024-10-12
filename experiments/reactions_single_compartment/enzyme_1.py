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
# ## Enzyme Kinetics : accurate numerical solution compared to the Michaelis-Menten model approximation
#
# #### Two Coupled Reactions: `E + S <-> ES*`, and  `ES* -> E + P` , using real-life kinetic parameters.  
# #### Scenario with small amount of Enzyme, relative to the initial Substrate concentration.
#
# Two other reactions being assumed negligible, and not included, are :  
# 1) the reverse direction of the 2nd reaction  
# 2) the non-catalyzed E + S <-> P reaction

# %% [markdown]
# The reaction of the enzyme `Adenosinedeaminase` with the substrate `2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine`,  
# and the initial concentration values choosen below, all satisfy the customary Michaelis-Menten assumptions that  
# `[E] << [S]` and that the reaction rate constants satisfy `k1_reverse >> k2_forward`
#
# Source of kinetic parameters:  *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% [markdown]
# ### TAGS :  "uniform compartment", "chemistry", "numerical", "enzymes"

# %%
LAST_REVISED = "Oct. 11, 2024"
LIFE123_VERSION = "1.0.0.beta.39"   # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import ipynbname
import pandas as pd

from life123 import check_version, ChemData, UniformCompartment, GraphicLog, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %% tags=[]
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
chem_data = ChemData(names=["P", "ES*"])

# %%
# Our Enzyme
chem_data.add_chemical(name="Adenosinedeaminase", label="E") 

# Our Substrate
chem_data.add_chemical(name="2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine", label="S");

# %%
chem_data.all_chemicals()

# %%

# %% [markdown]
# ### Specify the Kinetic Parameters

# %% [markdown]
# Source: *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% tags=[]
# Reaction E + S <-> ES* , with 1st-order kinetics, 
# and a forward rate that is much faster than its revers one
chem_data.add_reaction(reactants=["E", "S"], products=["ES*"],
                       forward_rate=18., reverse_rate=100.) 

# Reaction ES* <-> E + P , with 1st-order kinetics
chem_data.add_reaction(reactants=["ES*"], products=["E", "P"],
                forward_rate=49.,reverse_rate=0)

chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %%

# %%

# %% [markdown]
# ## Set the initial concentrations of all the chemicals
# ### in the scenario we're exploring, there's little enzyme relative to initial substrate concentration

# %%
S0 = 20.
E0 = 1.

# %%
# Here we use the "slower" preset for the variable steps, a conservative option prioritizing accuracy over speed
uc = UniformCompartment(chem_data=chem_data, preset="slower")
uc.set_conc(conc={"S": S0, "E": E0})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%
uc.enable_diagnostics()   # To save diagnostic information about the simulation - in particular, the REACTION RATES

# %% [markdown] tags=[]
# #### Advance the reactions to equilibrium

# %%
# Perform the reactions
uc.single_compartment_react(duration=1.5, initial_step=0.05)

# %%
uc.plot_history(colors=['cyan', 'green', 'violet', 'red'], show_intervals=True, title_prefix="Small amout of E relative to S(0)")

# %%
uc.plot_history(colors=['cyan', 'green', 'violet', 'red'], title_prefix="DETAIL at early times",
                range_x=[0, 0.5], range_y=[0, 2.])

# %% [markdown]
# Notice how the bound enzyme `ES*` quickly builds up at the very beginning, from time 0 to roughly 0.01 ... and that, in the longer term, the enzyme returns to its unbound state `E`

# %%

# %% [markdown]
# ### What is the initial rate of production of the final reaction product `P`?   
# One could take the numerical derivative (gradient) of the time values of [P] - but no need to!  **Reaction rates are computed in the course of the simulation, and stored alongside other diagnostic data**, if diagnostics were enabled (as we did indeed enable)

# %%
rates = uc.get_diagnostics().get_rxn_rates(rxn_index=1)   # We specify the reaction that involves `P`
rates

# %%
# Let's take a look at how the reaction rate varies with time
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time",
                         x_var="START_TIME", fields="rate", 
                         x_label="time", y_label="dP/dt")

# %%
# A closer peek at it maximum value
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time (DETAIL at early times)",
                         x_var="START_TIME", fields="rate", 
                         x_label="time", y_label="dP/dt",
                         range_x=[0,0.05], range_y=[33., 34.5])

# %% [markdown]
# As we saw earlier, the time it took for `ES*` to build up was about 0.01  
# **Ignoring the brief initial transient phase, the reaction rate starts at about 34.1**

# %%

# %% [markdown]
# # 2. Comparing the results to the Michaelis-Menten model approximation

# %% [markdown]
# ### Let's compute some parameters used by the Michaelis-Menten model   
# for background reference, see:  https://vallance.chem.ox.ac.uk/pdfs/KineticsLectureNotes.pdf (p. 20)

# %%
kM = (chem_data.get_forward_rate(1) + chem_data.get_reverse_rate(0)) / chem_data.get_forward_rate(0)
kM      # (49. + 100.) / 18. in this experiment

# %%
kcat = chem_data.get_forward_rate(1)
kcat

# %%
vmax = kcat * E0
vmax

# %%
initial_rxn_rate = (vmax * S0) / (kM + S0)
initial_rxn_rate

# %% [markdown]
# Not too far from the value of about 34.1 we saw earlier...   
# To keep in mind that the initial part of the reaction is affected by transients

# %% [markdown]
# #### Now, let's look at rate as a function of [S]; we'll compare what we computed earlier vs. as given by the approximation of the Michaelis-Menten model

# %%

# %%
# Let's retrieve all the values of `S` at the various simulation time
hist = uc.get_history()[:-1]    # Dropping the last row, because no rate information is known about the next simulation step not taken!
hist

# %%
# Now, let's put together in the same table (Pandas dataframe) the `S` values from above, and the `rates` we extracted earlier

assert len(hist) == len(rates)   # They'd better match up!

rate_and_substrate = pd.DataFrame({
    "SYSTEM TIME": hist["SYSTEM TIME"],   # We don't actually need the time, but just for clarity
    "S": hist["S"],
    "computed_rate": rates["rate"]
})

rate_and_substrate

# %%
# Let's a column with the rate estimated by the Michaelis-Menten model
rate_and_substrate["Michaelis_rate"] = vmax * rate_and_substrate["S"] / (kM + rate_and_substrate["S"])
rate_and_substrate

# %%
# Let's see how our computed rate compares with the approximations from the Michaelis-Menten model
PlotlyHelper.plot_pandas(df=rate_and_substrate, x_var="S", fields=["computed_rate", "Michaelis_rate"],
                         title="Reaction rate, dP/dt, as a function of Substrate concentration",
                         y_label="dP/dt")

# %% [markdown]
# Virtually overlapped plots, except at the very right!

# %% [markdown]
# ## In this scenario, the Michaelis-Menten model is remarkable accurate EXCEPT at the very early times of the reaction (the large values of `S`), which is the brief initial transient phase when `ES*` builds up from zero

# %% [markdown]
# #### This is as expected, because the kinetic parameters of this reaction satisfy the customary Michaelis-Menten assumptions that `[E] << [S]` and that the rates satisfy `k1_reverse >> k2_forward`

# %%
