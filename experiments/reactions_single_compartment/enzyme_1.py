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
# ## Enzyme Kinetics not employing the Michaelis-Menten simplified model
#
# #### 2 Coupled Reactions: `E + S <-> ES*`, and  `ES* -> E + P` , with real-life kinetic parameters  
#
# Two reactions being assumed negligible, and not included, are :  
# 1) the reverse direction of the 2nd reaction  
# 2) the non-catalyzed E + S <-> P reaction

# %% [markdown]
# #### Source of kinetic parameters:  page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023

# %%
LAST_REVISED = "Oct. 9, 2024"
LIFE123_VERSION = "1.0.0.beta.39"    # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import ipynbname
import numpy as np

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
chem_data = ChemData(names=["P", "ES*"])

# %%
chem_data.add_chemical(name="Adenosinedeaminase", label="E")     # Our Enzyme

# %%
chem_data.add_chemical(name="2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine", label="S")     # Our Substrate

# %%
chem_data.all_chemicals()

# %%

# %% [markdown]
# ### Specify the Kinetic Parameters

# %% tags=[]
# Reaction E + S <-> ES* , with 1st-order kinetics, and a forward rate that is much faster than it was without the enzyme
chem_data.add_reaction(reactants=["E", "S"], products=["ES*"],
                forward_rate=18., reverse_rate=100.) 

# Reaction ES* <-> E + P , with 1st-order kinetics, and a forward rate that is much faster than it was without the enzyme
# Thermodynamically, the total energy change of this reaction and the previous one adds up to the same value as the reaction without the enzyme (-3989.73)
chem_data.add_reaction(reactants=["ES*"], products=["E", "P"],
                forward_rate=49.,reverse_rate=0)

chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %%

# %%

# %% [markdown]
# # 1. Set the initial concentrations of all the chemicals - little enzyme

# %%
# Here we use the "TBA" preset for the variable steps, a compromise between speed and accuracy
uc = UniformCompartment(chem_data=chem_data, preset="slower")
uc.set_conc(conc={"S": 20., "E": 1.})      # Small ampount of enzyme `E`, relative to substrate `S`
uc.describe_state()

# %%
uc.enable_diagnostics()   # To save diagnostic information about the call to single_compartment_react()
                          # Useful for insight into the inner workings of the simulation

# %% [markdown] tags=[]
# ### Advance the reactions to equilibrium

# %%
# Perform the reactions
uc.single_compartment_react(duration=1.5, initial_step=0.05)

# %%
uc.plot_history(colors=['cyan', 'green', 'violet', 'red'], show_intervals=True, title_prefix="Small amout of E relative to S")

# %%
uc.plot_history(colors=['cyan', 'green', 'violet', 'red'], title_prefix="DETAIL at early times",
                range_x=[0, 0.8], range_y=[0, 2.])

# %%

# %% [markdown]
# ### What is the initial rate of production of the final reaction product P?   
# One could take the numerical derivative (gradient) of the time values of [P] - but no need to!  Reaction rates are computed in the course of the simulation, and stored alongside other diagnostic data, if diagnostics were enabled (as we did indeed)

# %%
rates = uc.get_diagnostics().get_rxn_rates(rxn_index=1)
rates

# %%
rates[:310]

# %%
rates[311:330]

# %%
PlotlyHelper.plot_pandas(df=rates, 
                         title="Reaction rate, dP/dt, over time",
                         x_var="START_TIME", fields="rate", 
                         x_label="time", y_label="dP/dt")

# %%
PlotlyHelper.plot_pandas(df=diagnostic_df, 
                         title="Reaction rate, dP/dt, over time (DETAIL at early times)",
                         x_var="START_TIME", fields="rate", 
                         x_label="time", y_label="dP/dt",
                         range_x=[0,0.05], range_y=[33., 34.5])

# %%

# %% [markdown]
# ## Let's compute the Michaelis constant kM:

# %%
kM = (49. + 100.) / 18.
kM

# %%
kcat = 49.

# %%
vmax = kcat * 1.
vmax

# %%
rate = (vmax * 20.) / (kM + 20.)
rate

# %%

# %%
hist = uc.get_history()[:-1]    # Dropping the last row, because no rate information is know about the next simulation step not taken!
hist

# %%
assert len(hist) == len(rates)   # They'd better match up!

# %%
