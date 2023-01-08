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
# ## 2 COUPLED reactions of different speeds:  
# ### A <-> B (fast) and B <-> C (slow)
# All 1st order. Taken to equilibrium.   
# The concentration of the intermediate product B manifests 1 oscillation (transient "overshoot")
#
# (Adaptive variable time resolution is used)
#
# LAST REVISED: Jan. 8, 2023

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from modules.reactions.reaction_data import ReactionData as chem
from modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.express as px
from modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # Initialize the System
# Specify the chemicals and the reactions

# %% tags=[]
# Specify the chemicals
chem_data = chem(names=["A", "B", "C"])

# Reaction A <-> B (fast)
chem_data.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=64., reverse_rate=8.) 

# Reaction B <-> C (slow)
chem_data.add_reaction(reactants=["B"], products=["C"],
                       forward_rate=12., reverse_rate=2.) 

print("Number of reactions: ", chem_data.number_of_reactions())

# %%
chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# # Start the simulation

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics.set_conc([50., 0, 0.], snapshot=True)

# %%
dynamics.describe_state()

# %%
dynamics.history.get()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()
#dynamics.verbose_list = [1, 2, 3]      # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

# The changes of concentrations vary very rapidly early on; so, we'll be using dynamic_step=4 , i.e. increase time resolution
# by x4 initially, as long as the reaction remains "fast" (based on a threshold of 5% change)
dynamics.single_compartment_react(time_step=0.02, reaction_duration=0.4,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_step=10, fast_threshold=15)      

# %% [markdown]
# ### Note: the argument _dynamic_step=10_ splits the time steps in 10 for any reactions that are "fast-changing" (as determined using _fast_threshold=15_ )

# %%
df = dynamics.history.get()
df

# %%
# Let's expand the last part
df.loc[60:]

# %% [markdown]
# ### Notice:
# * the reaction proceeds in smaller steps in the earlier times (until t=0.160, in line 80), when the concentrations are changing much more rapidly 
#
# * between lines 70 and 80, only rection #1 is regarded as fast-changing (based on the fast_threshold we specified in the _simulation run_); previously, both reactions were regarded as fast-changing
#
# * "fast-changing" and "slow-changing" is NOT the same thing as "fast" and "slow" reaction kinetics.  For example, reaction #1, though it has much slower kinetics than reaction #0, involves large relative concentration changes because [C] is small
#
# * after step 80, both reactions are regarded as slow-changing, and no more intermediate steps are used

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium(tolerance=0.2)

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
fig = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["A", "B", "C"], 
              title="Coupled reactions A <-> B and B <-> C",
              color_discrete_sequence = ['blue', 'red', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %% [markdown]
# #### For diagnostic insight, uncomment the following lines:

# %%
dynamics.diagnostic_data_baselines.get()

# %%
#dynamics.diagnose_variable_time_steps()

# %%
dynamics.get_diagnostic_data(rxn_index=0)

# %%
dynamics.get_diagnostic_data(rxn_index=0).loc[60:]

# %%
dynamics.get_diagnostic_data(rxn_index=1)

# %%
dynamics.get_diagnostic_data(rxn_index=1).loc[60:]

# %% [markdown]
# ## Perform some verification

# %%
# Verify that the stoichiometry is respected at every reaction step/substep (NOTE: it requires earlier activation of saving diagnostic data)
dynamics.stoichiometry_checker_entire_run()

# %%

# %%

# %%
chemical_list = dynamics.reaction_data.get_all_names()
chemical_list

# %%
chemical_delta_list = dynamics.delta_names()
chemical_delta_list

# %%
row_baseline = 1
row_0 = 1
row_1 = 1

# %%
df_row = dynamics.diagnostic_data_baselines.get().loc[row_baseline]
df_row

# %%
conc_arr_before = df_row[chemical_list].to_numpy()
conc_arr_before

# %%
# For reaction 0
delta_0 = dynamics.get_diagnostic_data(rxn_index=0).loc[row_0][chemical_delta_list].to_numpy()
delta_0

# %%
# For reaction 1
delta_1 = dynamics.get_diagnostic_data(rxn_index=1).loc[row_1][chemical_delta_list].to_numpy()
delta_1

# %%
conc_after = conc_arr_before + delta_0 + delta_1
conc_after

# %%
dynamics.diagnostic_data_baselines.get().loc[row_baseline+1]

# %%
import numpy as np

# %%
dynamics.get_diagnostic_data(rxn_index=0).loc[59:72]

# %%
dynamics.get_diagnostic_data(rxn_index=1).loc[59:72]


# %%
def foo(row_baseline):
    chemical_list = dynamics.reaction_data.get_all_names()  # EXAMPLE: ["A", "B", "C"]
    chemical_delta_list = dynamics.delta_names()  # EXAMPLE: ["Delta A", "Delta B", "Delta C"]
    
    row_0 = row_baseline
    row_1 = row_baseline
    conc_arr_before = dynamics.diagnostic_data_baselines.get().loc[row_baseline][chemical_list].to_numpy()
    
    # For reaction 0
    delta_0 = dynamics.get_diagnostic_data(rxn_index=0).loc[row_0][chemical_delta_list].to_numpy()
    delta_0
    
    # For reaction 1
    delta_1 = dynamics.get_diagnostic_data(rxn_index=1).loc[row_1][chemical_delta_list].to_numpy()
    delta_1
    
    conc_after = conc_arr_before + delta_0 + delta_1
    print(conc_after)

    next_system_state = dynamics.diagnostic_data_baselines.get().loc[row_baseline+1][chemical_list].to_numpy()
    print(next_system_state)
    
    #return conc_after, next_system_state
    
    return np.allclose(conc_after.astype(float), next_system_state.astype(float))
    


# %%

# %%
foo(69)

# %%
chemical_list = dynamics.reaction_data.get_all_names()  # EXAMPLE: ["A", "B", "C"]
chemical_delta_list = dynamics.delta_names()  # EXAMPLE: ["Delta A", "Delta B", "Delta C"]

row_baseline = 70
row_0 = row_baseline
row_1 = row_baseline
conc_arr_before = dynamics.diagnostic_data_baselines.get().loc[row_baseline][chemical_list].to_numpy()

# For reaction 0
delta_0 = 0
#delta_0 = dynamics.get_diagnostic_data(rxn_index=0).loc[row_0][chemical_delta_list].to_numpy()
delta_0

# For reaction 1
delta_1 = dynamics.get_diagnostic_data(rxn_index=1).loc[row_1][chemical_delta_list].to_numpy()
delta_1

conc_after = conc_arr_before + delta_0 + delta_1
print(conc_after)

next_system_state = dynamics.diagnostic_data_baselines.get().loc[row_baseline+1][chemical_list].to_numpy()
print(next_system_state)


# %%
row_baseline += 1
row_1 += 1

# %%
conc_arr_before = dynamics.diagnostic_data_baselines.get().loc[row_baseline][chemical_list].to_numpy()
print(conc_arr_before)

# For reaction 0
delta_0 = 0
#delta_0 = dynamics.get_diagnostic_data(rxn_index=0).loc[row_0][chemical_delta_list].to_numpy()
delta_0

# For reaction 1
delta_1 = dynamics.get_diagnostic_data(rxn_index=1).loc[row_1][chemical_delta_list].to_numpy()
delta_1

conc_after = conc_arr_before + delta_0 + delta_1
print(conc_after)

next_system_state = dynamics.diagnostic_data_baselines.get().loc[row_baseline+1][chemical_list].to_numpy()
print(next_system_state)

print(np.allclose(conc_after.astype(float), next_system_state.astype(float)))


# %%
row_baseline += 1
row_1 += 1

# %%
conc_arr_before = dynamics.diagnostic_data_baselines.get().loc[row_baseline][chemical_list].to_numpy()
print(conc_arr_before)

# For reaction 0
delta_0 = 0
#delta_0 = dynamics.get_diagnostic_data(rxn_index=0).loc[row_0][chemical_delta_list].to_numpy()
delta_0

# For reaction 1
delta_1 = dynamics.get_diagnostic_data(rxn_index=1).loc[row_1][chemical_delta_list].to_numpy()
delta_1

conc_after = conc_arr_before + delta_0 + delta_1
print(conc_after)

next_system_state = dynamics.diagnostic_data_baselines.get().loc[row_baseline+1][chemical_list].to_numpy()
print(next_system_state)

print(np.allclose(conc_after.astype(float), next_system_state.astype(float)))


# %%
def f(fast_list):
    conc_arr_before = dynamics.diagnostic_data_baselines.get().loc[row_baseline][chemical_list].to_numpy()
    print(conc_arr_before)

    # For reaction 0
    delta_0 = 0
    #delta_0 = dynamics.get_diagnostic_data(rxn_index=0).loc[row_0][chemical_delta_list].to_numpy()
    delta_0

    # For reaction 1
    delta_1 = dynamics.get_diagnostic_data(rxn_index=1).loc[row_1][chemical_delta_list].to_numpy()
    delta_1

    conc_after = conc_arr_before + delta_0 + delta_1
    print(conc_after)

    next_system_state = dynamics.diagnostic_data_baselines.get().loc[row_baseline+1][chemical_list].to_numpy()
    print(next_system_state)

    print(np.allclose(conc_after.astype(float), next_system_state.astype(float))) 


# %%
row_baseline += 1
row_1 += 1

# %%
f([1])


# %%
def ff(fast_list):
    conc_arr_before = dynamics.diagnostic_data_baselines.get().loc[row_baseline][chemical_list].to_numpy()
    print(conc_arr_before)

    # For reaction 0
    if (0 in fast_list):
        delta_0 = dynamics.get_diagnostic_data(rxn_index=0).loc[row_0][chemical_delta_list].to_numpy()
    else:
        delta_0 = 0   # TODO: this will change in a future version!
        
    # For reaction 1    
    if (1 in fast_list):
        delta_1 = dynamics.get_diagnostic_data(rxn_index=1).loc[row_1][chemical_delta_list].to_numpy()
    else:
        delta_1 = 0 
        
    conc_after = conc_arr_before + delta_0 + delta_1
    print(conc_after)

    next_system_state = dynamics.diagnostic_data_baselines.get().loc[row_baseline+1][chemical_list].to_numpy()
    print(next_system_state)

    print(np.allclose(conc_after.astype(float), next_system_state.astype(float))) 


# %%
row_baseline += 1
row_1 += 1

# %%
ff([1])


# %%
def g(fast_list, row_baseline, row_list):
    print("ROW of baseline data: ", row_baseline)
    print("row_list: ", row_list)
    
    conc_arr_before = dynamics.diagnostic_data_baselines.get().loc[row_baseline][chemical_list].to_numpy().astype(float)
    print("baseline concentration: ", conc_arr_before)

    delta_cumulative = np.zeros(dynamics.reaction_data.number_of_chemicals(), 
                                dtype=float)  # One element per chemical species
    
    # For each reaction
    for rxn_index in range(dynamics.reaction_data.number_of_reactions()):
        if (rxn_index in fast_list):
            row = row_list[rxn_index]
            delta_rxn = dynamics.get_diagnostic_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy().astype(float)
        else:
            delta_rxn = np.zeros(dynamics.reaction_data.number_of_chemicals(), 
                                dtype=float)   # TODO: this will change in a future version!
        
        print(f"For rxn {rxn_index}: delta_rxn = {delta_rxn}")
        delta_cumulative += delta_rxn
              
    print("delta_cumulative: ", delta_cumulative)
    
    conc_after = conc_arr_before + delta_cumulative
    print("updated concentration: ", conc_after)

    next_system_state = dynamics.diagnostic_data_baselines.get().loc[row_baseline+1][chemical_list].to_numpy()
    print("concentration from the system state: ", next_system_state)

    print(np.allclose(conc_after.astype(float), next_system_state.astype(float))) 


# %%
row_baseline += 1
row_1 += 1

# %%
g([1])

# %%
delta_cumulative = np.zeros(dynamics.reaction_data.number_of_chemicals(), 
                                dtype=float)
delta_cumulative

# %%
ff([1])

# %%
row_list = [70,76]

# %%
g(fast_list=[1], row_baseline=row_baseline, row_list=row_list)

# %%
row_baseline += 1
row_list[1] += 1
row_list

# %%
g(fast_list=[1], row_baseline=row_baseline, row_list=row_list)


# %%

# %%

# %%

# %%

# %%

# %%

# %%

# %%
def good(active_list, row_baseline, row_list):
    print("ROW of baseline data: ", row_baseline)
    print("row_list: ", row_list)
    
    conc_arr_before = dynamics.diagnostic_data_baselines.get().loc[row_baseline][chemical_list].to_numpy().astype(float)
    print("baseline concentration: ", conc_arr_before)

    delta_cumulative = np.zeros(dynamics.reaction_data.number_of_chemicals(), 
                                dtype=float)  # One element per chemical species
    
    # For each reaction
    for rxn_index in range(dynamics.reaction_data.number_of_reactions()):
        if (rxn_index in active_list):
            row = row_list[rxn_index]
            delta_rxn = dynamics.get_diagnostic_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy().astype(float)
        else:
            delta_rxn = np.zeros(dynamics.reaction_data.number_of_chemicals(), 
                                dtype=float)   # TODO: this will change in a future version!
        
        print(f"For rxn {rxn_index}: delta_rxn = {delta_rxn}")
        delta_cumulative += delta_rxn
              
    print("delta_cumulative: ", delta_cumulative)
    
    conc_after = conc_arr_before + delta_cumulative
    print("updated concentration: ", conc_after)

    next_system_state = dynamics.diagnostic_data_baselines.get().loc[row_baseline+1][chemical_list].to_numpy()
    print("concentration from the system state: ", next_system_state)

    print(np.allclose(conc_after.astype(float), next_system_state.astype(float))) 


# %%
row_baseline = 68
row_list = [68, 68]
active_list = [0, 1]

# %%
good(active_list=active_list, row_baseline=row_baseline, row_list=row_list)

# %%
dynamics.get_diagnostic_data(rxn_index=0).loc[row_list[0]]['time_subdivision']

# %%
dynamics.get_diagnostic_data(rxn_index=0).loc[row_list[0]]['substep']

# %%
dynamics.get_diagnostic_data(rxn_index=1).loc[row_list[1]]['time_subdivision']

# %%
dynamics.get_diagnostic_data(rxn_index=1).loc[row_list[1]]['substep']

# %%
row_baseline += 1

# %%
if dynamics.get_diagnostic_data(rxn_index=0).loc[row_list[0]]['time_subdivision'] == dynamics.get_diagnostic_data(rxn_index=1).loc[row_list[1]]['time_subdivision']:
    row_list[0]  += 1
    row_list[1]  += 1
else:
    print("TBA")

# %%
good(active_list=active_list, row_baseline=row_baseline, row_list=row_list)

# %%
row_baseline += 1

if dynamics.get_diagnostic_data(rxn_index=0).loc[row_list[0]]['time_subdivision'] == dynamics.get_diagnostic_data(rxn_index=1).loc[row_list[1]]['time_subdivision']:
    row_list[0]  += 1
    row_list[1]  += 1
else:
    print("TBA")

# %%
row_list

# %%
time_sub0 = dynamics.get_diagnostic_data(rxn_index=0).loc[row_list[0]]['time_subdivision']
time_sub0

# %%
time_sub1 = dynamics.get_diagnostic_data(rxn_index=1).loc[row_list[1]]['time_subdivision']
time_sub1

# %%
if time_sub0 > time_sub1:
    active_list = [0]
elif time_sub1 > time_sub0:
    active_list = [1]

# %%
active_list

# %%
#active_list = [1]   # ************** HOW TO DETECT THIS???

# %%
good(active_list=active_list, row_baseline=row_baseline, row_list=row_list)

# %%
row_baseline += 1
for i in active_list:
    row_list[i] += 1

# %%
good(active_list=active_list, row_baseline=row_baseline, row_list=row_list)

# %%
row_baseline += 1
for i in active_list:
    row_list[i] += 1

# %%
good(active_list=active_list, row_baseline=row_baseline, row_list=row_list)

# %%
row_baseline += 1
for i in active_list:
    row_list[i] += 1

good(active_list=active_list, row_baseline=row_baseline, row_list=row_list)    

# %%
row_baseline += 1
for i in active_list:
    row_list[i] += 1

# %%
row_list

# %%
row_baseline

# %%
substep1 = dynamics.get_diagnostic_data(rxn_index=1).loc[row_list[1]]['substep']
substep1

# %%
time_sub1

# %%
if substep1 == time_sub1 - 1:
    active_list = [0,1]

# %%
active_list

# %%
#active_list = [0,1]      # ************** HOW TO DETECT THIS???

# %%
good(active_list=active_list, row_baseline=row_baseline, row_list=row_list) 

# %%
row_baseline += 1
for i in active_list:
    row_list[i] += 1

# %%
good(active_list=active_list, row_baseline=row_baseline, row_list=row_list) 

# %%
row_baseline += 1
for i in active_list:
    row_list[i] += 1

# %%
good(active_list=active_list, row_baseline=row_baseline, row_list=row_list) 

# %%
row_baseline += 1
for i in active_list:
    row_list[i] += 1

# %%
good(active_list=active_list, row_baseline=row_baseline, row_list=row_list) 

# %%
