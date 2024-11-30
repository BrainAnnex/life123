# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# ## Enzyme Kinetics in a NON Michaelis-Menten modality
# #### 3 Coupled Reactions: `S <-> P` , `E + S <-> ES`, and  `ES <-> E + P`     
# A direct reaction and the same reaction, catalyzed by an enzyme `E` and showing the intermediate state.  
# Re-run from same initial concentrations of S ("Substrate") and P ("Product"), for various concentations of the enzyme `E`: from zero to hugely abundant 
# ### We'll REJECT the customary Michaelis-Menten assumptions that `[E] << [S]` and that the rate constants satisfy `k1_reverse >> k2_forward` !   
# #### We'll **repeat runs with increasingly higher initial concentration of Enzyme**, and explore exotic scenarios with lavish amount of enzyme, leading to diminishing (though fast-produced!) products,  and a buildup of the (not-so-transient!) `ES` intermediate

# %% [markdown]
# ### TAGS :  "uniform compartment", "chemistry", "enzymes"

# %%
LAST_REVISED = "Nov. 29, 2024"
LIFE123_VERSION = "1.0-rc.1"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from experiments.get_notebook_info import get_notebook_basename

from life123 import check_version, ChemData, UniformCompartment, CollectionTabular, GraphicLog

import pandas as pd

# %%
check_version(LIFE123_VERSION)

# %%

# %% tags=[]
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%
# Initialize the system
chem_data = ChemData()

chem_data.add_chemical(name="S",  plot_color="cyan")
chem_data.add_chemical(name="E",  plot_color="violet")
chem_data.add_chemical(name="ES", plot_color="red")
chem_data.add_chemical(name="P",  plot_color="green")

                                                     
# Reaction S <-> P , with 1st-order kinetics, favorable thermodynamics in the forward direction, 
# and a forward rate that is much slower than it would be with the enzyme - as seen in the next reaction, below
chem_data.add_reaction(reactants="S", products="P",
                       forward_rate=1., delta_G=-3989.73)

     
# Reaction E + S <-> ES , with 1st-order kinetics, and a forward rate that is much faster than it was without the enzyme
# Thermodynamically, the forward direction is at a disadvantage (higher energy state) because of the activation barrier in forming the transient state ES
chem_data.add_reaction(reactants=["E", "S"], products="ES",
                       forward_rate=100., delta_G=2000)                           
                                                      
# Reaction ES <-> E + P , with 1st-order kinetics, and a forward rate that is much faster than it was without the enzyme
# Thermodynamically, the total energy change of this reaction and the previous one adds up to the same value as the reaction without the enzyme (-3989.73)
chem_data.add_reaction(reactants="ES", products=["E", "P"],
                       forward_rate=200., delta_G=-5989.73)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %% [markdown]
# Note that `E` is not labeled as an "enzyme" because it doesn't appear as a catalyst in any of the registered reactions; it only becomes an enzyme in the context of the _compound_ reaction from (2) and (3)

# %%

# %%

# %% [markdown]
# # 1. Set the initial concentrations of all the chemicals - starting with no enzyme

# %%
uc = UniformCompartment(chem_data=chem_data, preset="mid")
uc.set_conc(conc={"S": 20.})      # Initially, no enzyme `E`
uc.describe_state()

# %% [markdown] tags=[]
# ### Advance the reactions (for now without enzyme) to equilibrium

# %%
#uc.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# Perform the reactions
uc.single_compartment_react(duration=4.0,
                            initial_step=0.1)

# %%
uc.plot_history(show_intervals=True, title_prefix="With ZERO enzyme")

# %% [markdown]
# ### The reactions, lacking enzyme, are proceeding slowly towards equilibrium, just like the reaction that was discussed in part 1 of the experiment "enzyme_1"

# %%
# Verify that the reactions have reached equilibrium
uc.is_in_equilibrium(tolerance=2)

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=1.0)

# %%
uc.get_history(columns=['SYSTEM TIME', 'E', 'ES', 'P'], tail=1)  # Last point in the simulation

# %%
P_equil = uc.get_chem_conc("P")
P_equil

# %%
P_70_threshold = P_equil * 0.70
P_70_threshold

# %%
uc.reach_threshold(chem="P", threshold=P_70_threshold)

# %%

# %%

# %% [markdown]
# # 2. Re-start all the reactions from the same initial concentrations - except for now having a tiny amount of enzyme (two orders of magnitude less than the starting [S])

# %%
E_init = 0.2     # A tiny bit of enzyme `E`: 1/100 of the initial [S]

# %%
uc = UniformCompartment(chem_data=chem_data, preset="slower")      # A brand-new simulation, with the same chemicals and reactions as before

uc.set_conc(conc={"S": 20., "E": E_init})
uc.describe_state()

# %%
# Perform the reactions (The duration of the run was manually adjusted for optimal visibility)
uc.single_compartment_react(duration=1.3,
                            initial_step=0.00005)

# %%
# Verify that the reactions have reached equilibrium
uc.is_in_equilibrium(verbose=False)

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=1.0)

# %%
# Early part of the reaction
uc.plot_history(title_prefix=f"Early times when E0 = {E_init}", range_x=[0, 0.4])

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=0.4)

# %% [markdown]
# ### Now, let's investigate E and ES in the very early times

# %% [markdown]
# ## Notice how even a tiny amount of enzyme (1/100 of the initial [S])  makes a very pronounced difference!

# %%
# The very early part of the reaction
uc.plot_history(chemicals=['E', 'ES', 'P'],
                title_prefix=f"Detail when E0 = {E_init}", range_x=[0, 0.002], range_y=[0, 0.2])

# %% [markdown]
# ### Notice how, with this small initial concentration of [E], the timescale of [E] and [ES] is vastly faster than that of [P] and [S]

# %%
# The full reaction of E and ES
uc.plot_history(chemicals=['E', 'ES'], show_intervals=True, 
                title_prefix=f"E and ES when E0 = {E_init}")

# %% [markdown]
# Notice how at every onset of instability in [E] or [ES], the adaptive time steps shrink down

# %%
uc.get_history(columns=['SYSTEM TIME', 'E', 'ES', 'P'], tail=1)  # Last point in the simulation

# %% [markdown]
# Interestingly, most of the inital [E] of 0.2 is now, at equilibrium, stored as [ES]=0.119; the energy of the "activation barrier" from E + S to ES might be unrealistically low (2000 Joules).  Zooming in on the very earl part of the plot:

# %%
P_equil = uc.get_chem_conc("P")
P_equil

# %%
P_70_threshold = P_equil * 0.70
P_70_threshold

# %%
uc.reach_threshold(chem="P", threshold=P_70_threshold)

# %%

# %%

# %% [markdown] tags=[]
# # 3. Keep increasing the initial [E]

# %%
E_init = 1.0

# %%
uc = UniformCompartment(chem_data=chem_data, preset="slower")      # A brand-new simulation, with the same chemicals and reactions as before

uc.set_conc(conc={"S": 20., "E": E_init})     
uc.describe_state()

# %%
# Perform the reactions (The duration of the run was manually adjusted for optimal visibility)
uc.single_compartment_react(duration=0.4, 
                            initial_step=0.00005)

# %%
# Verify that the reactions have reached equilibrium
uc.is_in_equilibrium(verbose=False)

# %%
uc.plot_history(show_intervals=True, 
                title_prefix=f"E0 = {E_init}")

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=0.4)

# %% [markdown]
# The timescale of [S] and [P] continues to become faster with a higher initial [E]

# %%
# The very early part of the reactions
uc.plot_history(chemicals=['E', 'ES', 'P'], 
                title_prefix=f"Detail when E0 = {E_init}", range_x=[0, 0.002], range_y=[0, 1.])

# %%
# The full reaction of E and ES
uc.plot_history(chemicals=['E', 'ES'], show_intervals=True, 
                title_prefix=f"E and ES when E0 = {E_init}")

# %%
uc.get_history(columns=['SYSTEM TIME', 'E', 'ES', 'P'], tail=1)  # Last point in the simulation

# %%
P_equil = uc.get_chem_conc("P")
P_equil

# %% [markdown]
# #### P_equil continues to decrease with increasing initial [E]

# %%
P_70_threshold = P_equil * 0.70
P_70_threshold

# %%
uc.reach_threshold(chem="P", threshold=P_70_threshold)

# %% [markdown]
# #### The time at which we reach the 70% threshold of the equilibrium value of P continues to decrease with increasing initial [E]

# %%

# %%

# %% [markdown] tags=[]
# # 4. Keep increasing the initial [E]

# %%
E_init = 2.0      # 1/10 of the initial [S]

# %%
uc = UniformCompartment(chem_data=chem_data, preset="slower")      # A brand-new simulation, with the same chemicals and reactions as before

uc.set_conc(conc={"S": 20., "E": E_init})     
uc.describe_state()

# %%
# Perform the reactions (The duration of the run was manually adjusted for optimal visibility)
uc.single_compartment_react(duration=0.2,
                            initial_step=0.00005)

# %%
# Verify that the reactions have reached equilibrium
uc.is_in_equilibrium(verbose=False)

# %%
uc.plot_history(show_intervals=True, 
                title_prefix=f"E0 = {E_init}")

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=0.4)

# %% [markdown]
# The timescale of [S] and [P] continues to become faster with a higher initial [E]

# %%
#The very early part of the reaction
uc.plot_history(chemicals=['E', 'ES', 'P'],
                title_prefix=f"Detail when E0 = {E_init}", range_x=[0, 0.002], range_y=[0, 2.])

# %%
# Show the full reaction of E and ES
uc.plot_history(chemicals=['E', 'ES'], show_intervals=True, 
                title_prefix=f"E and ES when E0 = {E_init}")

# %%
uc.get_history(columns=['SYSTEM TIME', 'E', 'ES', 'P'], tail=1)  # Last point in the simulation

# %%
P_equil = uc.get_chem_conc("P")
P_equil

# %% [markdown]
# #### P_equil continues to decrease with increasing initial [E]

# %%
P_70_threshold = P_equil * 0.70
P_70_threshold

# %%
uc.reach_threshold(chem="P", threshold=P_70_threshold)

# %% [markdown]
# #### The time at which we reach the 70% threshold of the equilibrium value of P continues to decrease with increasing initial [E]

# %%

# %%

# %% [markdown] tags=[]
# # 5. Keep increasing the initial [E]

# %%
E_init = 10.0      # 1/2 of the initial [S]

# %%
uc = UniformCompartment(chem_data=chem_data, preset="slower")      # A brand-new simulation, with the same chemicals and reactions as before

uc.set_conc(conc={"S": 20., "E": E_init})     
uc.describe_state()

# %%
# Perform the reactions (The duration of the run was manually adjusted for optimal visibility)
uc.single_compartment_react(duration=0.05,
                            initial_step=0.00005)

# %%
# Verify that the reactions have reached equilibrium
uc.is_in_equilibrium(verbose=False)

# %%
uc.plot_history(show_intervals=True, 
                title_prefix=f"E0 = {E_init}")

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=0.05)

# %% [markdown]
# #### The timescale of [S] and [P] continues to become faster with a higher initial [E] -- **AND THEY'RE NOW APPROACHING THE TIMESCALES OF E AND ES**

# %%
# The very early part of the reaction
uc.plot_history(chemicals=['E', 'ES', 'P'], 
                title_prefix=f"Detail when E0 = {E_init}", range_x=[0, 0.002], range_y=[0, 10.])

# %%
# The full reaction of E and ES
uc.plot_history(chemicals=['E', 'ES'], show_intervals=True, 
                title_prefix=f"E and ES when E0 = {E_init}")

# %% [markdown]
# #### Notice that at these higher initial concentrations of [E], we're now beginning to see overshoots in [E] and [ES]

# %%
uc.get_history(columns=['SYSTEM TIME', 'E', 'ES', 'P'], tail=1)  # Last point in the simulation

# %%
P_equil = uc.get_chem_conc("P")
P_equil

# %% [markdown]
# #### P_equil continues to decrease with increasing initial [E]

# %%
P_70_threshold = P_equil * 0.70
P_70_threshold

# %%
uc.reach_threshold(chem="P", threshold=P_70_threshold)

# %% [markdown]
# #### The time at which we reach the 70% threshold of the equilibrium value of P continues to decrease with increasing initial [E]

# %%

# %%

# %% [markdown] tags=[]
# # 6. Keep increasing the initial [E]

# %%
E_init = 20.0      # Same as the initial [S]

# %%
uc = UniformCompartment(chem_data=chem_data, preset="slower")      # A brand-new simulation, with the same chemicals and reactions as before

uc.set_conc(conc={"S": 20., "E": E_init})     
uc.describe_state()

# %%
# Perform the reactions (The duration of the run was manually adjusted for optimal visibility)
uc.single_compartment_react(duration=0.02, 
                            initial_step=0.00005)

# %%
# Verify that the reactions have reached equilibrium
uc.is_in_equilibrium(verbose=False)

# %%
uc.plot_history(show_intervals=True, 
                title_prefix=f"E0 = {E_init}")

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=0.02)

# %% [markdown]
# #### The timescale of [S] and [P] continues to become faster with a higher initial [E] -- **AND THEY'RE NOW GETTING CLOSE THE TIMESCALES OF E AND ES**

# %%
# The very early part of the reaction
uc.plot_history(title_prefix=f"Detail when E0 = {E_init}", range_x=[0, 0.002], range_y=[0, 20.])

# %%
# The full reaction of E and ES
uc.plot_history(chemicals=['E', 'ES'], show_intervals=True, 
                title_prefix=f"E and ES when E0 = {E_init}")

# %% [markdown]
# #### At these higher initial concentrations of [E], we're now beginning to see overshoots in [E] and [ES] that overlap less

# %%
uc.get_history(columns=['SYSTEM TIME', 'E', 'ES', 'P'], tail=1)  # Last point in the simulation

# %%
P_equil = uc.get_chem_conc("P")
P_equil

# %% [markdown]
# #### P_equil continues to decrease with increasing initial [E]

# %%
P_70_threshold = P_equil * 0.70
P_70_threshold

# %%
uc.reach_threshold(chem="P", threshold=P_70_threshold)

# %% [markdown] tags=[]
# #### The time at which we reach the 70% threshold of the equilibrium value of P continues to decrease with increasing initial [E]

# %%

# %%

# %% [markdown] tags=[]
# # 7. Keep increasing the initial [E]

# %%
E_init = 30.0      # 50% higher than the initial [S]

# %%
uc = UniformCompartment(chem_data=chem_data, preset="slower")      # A brand-new simulation, with the same chemicals and reactions as before

uc.set_conc(conc={"S": 20., "E": E_init})     
uc.describe_state()

# %%
# Perform the reactions (The duration of the run was manually adjusted for optimal visibility)
uc.single_compartment_react(duration=0.01, 
                            initial_step=0.00005)

# %%
# Verify that the reactions have reached equilibrium
uc.is_in_equilibrium(verbose=False)

# %%
uc.plot_history(show_intervals=True, 
                title_prefix=f"E0 = {E_init}")

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=0.01)

# %% [markdown]
# #### The timescale of [S] and [P] continues to become faster with a higher initial [E] -- **AND THEY'RE NOW COMPARABLE THE TIMESCALES OF E AND ES**

# %%
#The very early part of the reaction
uc.plot_history(title_prefix=f"Detail when E0 = {E_init}", range_x=[0, 0.005], range_y=[0, 30.])

# %%
# The full reaction of E and ES
uc.plot_history(chemicals=['E', 'ES'], show_intervals=True, 
                title_prefix=f"E and ES when E0 = {E_init}")

# %% [markdown]
# #### At these high initial concentrations of [E], we're now beginning to see [E] and [ES] no longer overlapping; the overshoot is still present in both

# %%
uc.get_history(columns=['SYSTEM TIME', 'E', 'ES', 'P'], tail=1)  # Last point in the simulation

# %%
P_equil = uc.get_chem_conc("P")
P_equil

# %% [markdown]
# #### P_equil continues to decrease with increasing initial [E]

# %%
P_70_threshold = P_equil * 0.70
P_70_threshold

# %%
uc.reach_threshold(chem="P", threshold=P_70_threshold)

# %% [markdown] tags=[]
# #### The time at which we reach the 70% threshold of the equilibrium value of P continues to decrease with increasing initial [E]

# %%

# %%

# %% [markdown] tags=[]
# # 8. Keep increasing the initial [E]

# %%
E_init = 60.0      # Triple the initial [S]

# %%
uc = UniformCompartment(chem_data=chem_data, preset="slower")      # A brand-new simulation, with the same chemicals and reactions as before

uc.set_conc(conc={"S": 20., "E": E_init})     
uc.describe_state()

# %%
# Perform the reactions (The duration of the run was manually adjusted for optimal visibility)
uc.single_compartment_react(duration=0.005,
                            initial_step=0.00005)

# %%
# Verify that the reactions have reached equilibrium
uc.is_in_equilibrium(verbose=False)

# %%
uc.plot_history(show_intervals=True, 
                title_prefix=f"E0 = {E_init}")

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=0.005)

# %% [markdown]
# #### The timescales of [S] and [P] continues to become faster with a higher initial [E] -- **AND NOW THEY REMAIN COMPARABLE THE TIMESCALES OF E AND ES**

# %% [markdown]
# #### At these high initial concentrations of [E],  [E] and [ES] no longer overlapping, and are now getting further apart

# %%
uc.get_history(columns=['SYSTEM TIME', 'E', 'ES', 'P'], tail=1)  # Last point in the simulation

# %%
P_equil = uc.get_chem_conc("P")
P_equil

# %% [markdown]
# #### P_equil continues to decrease with increasing initial [E]

# %%
P_70_threshold = P_equil * 0.70
P_70_threshold

# %%
uc.reach_threshold(chem="P", threshold=P_70_threshold)

# %% [markdown] tags=[]
# #### The time at which we reach the 70% threshold of the equilibrium value of P continues to decrease with increasing initial [E]

# %%

# %%

# %% [markdown] tags=[]
# # 9. Keep increasing the initial [E]

# %%
E_init = 100.0      # Quintuple the initial [S]

# %%
uc = UniformCompartment(chem_data=chem_data, preset="slower")      # A brand-new simulation, with the same chemicals and reactions as before

uc.set_conc(conc={"S": 20., "E": E_init})     
uc.describe_state()

# %%
# Perform the reactions (The duration of the run was manually adjusted for optimal visibility)
uc.single_compartment_react(duration=0.003,
                            initial_step=0.00005)

# %%
# Verify that the reactions have reached equilibrium
uc.is_in_equilibrium(verbose=False)

# %%
uc.plot_history(show_intervals=True, 
                title_prefix=f"E0 = {E_init}")

# %%
# Locate the intersection of the curves for [S] and [P]:
uc.curve_intersect("S", "P", t_start=0, t_end=0.005)

# %% [markdown]
# #### The timescales of [S] and [P] continues to become faster with a higher initial [E] -- **AND REMAIN COMPARABLE TO THE TIMESCALES OF E AND ES**

# %% [markdown]
# #### At these high initial concentrations of [E], the time curves of [E] and [ES] no longer overlapping, and are now getting further apart

# %%
uc.get_history(columns=['SYSTEM TIME', 'E', 'ES', 'P'], tail=1)  # Last point in the simulation

# %%
P_equil = uc.get_chem_conc("P")
P_equil

# %% [markdown]
# #### P_equil continues to decrease with increasing initial [E]

# %%
P_70_threshold = P_equil * 0.70
P_70_threshold

# %%
uc.reach_threshold(chem="P", threshold=P_70_threshold)

# %%

# %% [markdown]
# # CONCLUSION:  
# ### as the initial concentration of the Enzyme E increases from zero to lavish values much larger that the concentration of the Substrate S, the reaction keeps reaching an equilibrium faster and faster...  
# #### BUT the equilibrium concentration of the Product P steadily _reduces_ - and ES, far from being transient, actually builds up.   
# #### The **strange world of departing from the customary Michaelis-Menten assumptions** that `[E] << [S]` and that the rates satisfy `k1_reverse >> k2_forward` !

# %%
uc.get_history(head=1)  # First point in the simulation

# %%
uc.get_history(tail=1)  # Last point in the simulation

# %% [markdown]
# ### When the initial [E] is quite high, relatively little of the reactant S goes into making the product P ; most of S binds with the abundant E, to produce a lasting ES 
# The following manual stoichiometry check illustrates it.

# %% [markdown]
# #### Manual overall STOICHIOMETRY CHECK:

# %%
delta_S = 0.46627 - 20.0
delta_S

# %%
delta_E = 82.779099 - 100.0
delta_E

# %%
delta_E = 82.779099 - 100.0
delta_E

# %%
delta_ES = 17.220901 - 0
delta_ES

# %%
chem_data.describe_reactions()

# %%

# %% [markdown]
# 19.53373 units of `S` are consumed : a meager 2.312829 of that goes into the production of `P` (rxn 0)  ; the remainder of Delta S is:

# %%
19.53373 - 2.312829

# %% [markdown]
# That extra Delta S of 17.2209 combines with an equal amount of `E` to produce the same amount, 17.2209, of ES (rxn 1)  
#
# `E`, now depleted by that amount of 17.2209, reaches the value:

# %%
100 - 17.2209

# %%
# Also review that the chemical equilibrium holds, with the final simulation values
uc.is_in_equilibrium()

# %%
