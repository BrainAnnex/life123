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
# ## `A <-> B` , downregulated by the "shunt" (coupled reaction) `A <-> S`
# ### _Kinetic_ advantage (downregulation in early phase) vs. _Thermodynamic_ advantage (long-term downregulation) 
#
# **[Scenario 1](#down_regulate_1_scenario_1)** : No downregulation on `A <-> B `
#
# **[Scenario 2](#down_regulate_1_scenario_2)** : The shunt (`A <-> S`) has a *kinetic* advantage but *thermodynamic* DIS-advantage compared to `A <-> B `  
# (i.e. `A <-> S` is fast, but energetically unfavored) 
#
# **[Scenario 3](#down_regulate_1_scenario_3)** : The shunt (`A <-> S`) is has a *kinetic* DIS-advantage but a *thermodynamic* advantage compared to `A <-> B`     
# (i.e. `A <-> S` is slow, but energetically favored)  
#
# All reactions 1st order, mostly forward.  Taken to equilibrium.
#
# LAST REVISED: Feb. 4, 2023

# %% [markdown]
# ## Bathtub analogy:
# A is initially full, while B and S are empty.  
# If the "shunt" S is present, scenario 2 corresponds to a large pipe and a small elevation change...  
# while scenario 3 corresponds to a narrow pipe and a large elevation change.

# %% [markdown]
# ![Downregulated by shunt](../../docs/down_regulate_1.png)

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(2)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import numpy as np
import plotly.express as px
from src.modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# # <a name="down_regulate_1_scenario_1"></a>  Scenario 1: A <-> B in the absence of the 2nd reaction

# %% [markdown]
# ### Initialize the System
# Specify the chemicals and the reaction

# %% tags=[]
# Specify the chemicals
chem_data = chem(names=["A", "B"])

# Reaction A <-> B
chem_data.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=30., reverse_rate=5.)

chem_data.describe_reactions()

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)
dynamics.set_conc([50., 0.], snapshot=True)
dynamics.describe_state()

# %% [markdown]
# ### Run the reaction

# %%
dynamics.set_diagnostics()          # To save diagnostic information about the call to single_compartment_react()
#dynamics.verbose_list = [1, 2, 3]   # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

# The changes of concentrations vary very rapidly early on; 
# so, we'll be using the dynamic_substeps option to increase time resolution,
# as long as the reaction remains "fast" (based on a threshold of % change, as specified by fast_threshold)
dynamics.single_compartment_react(time_step=0.001, reaction_duration=0.05,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_substeps=4, rel_fast_threshold=25)

# %%
df_iterm = dynamics.get_history()
df_iterm

# %%
# Continue running the reaction at lover resolution
dynamics.single_compartment_react(time_step=0.002, reaction_duration=0.25,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_substeps=4, fast_threshold=10)

# %%
df = dynamics.get_history()
df

# %%
dynamics.explain_time_advance()

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium()

# %%
dynamics.plot_curves(colors=["blue", "green"], title="Single reaction A <-> B (no downregulation)")

# %%

# %% [markdown]
# # <a name="down_regulate_1_scenario_2"></a>  Scenario 2: 
# ### downregulated by shunt: 
# ### kinetically fast,   
# ### but with thermodynamical dis-advantage (i.e. energetically un-favored)

# %% tags=[]
# Register the new chemical ("S")
chem_data.add_chemical("S")

# Add the reaction A <-> S (fast shunt, poor thermodynical energetic advantage)
chem_data.add_reaction(reactants=["A"], products=["S"],
                       forward_rate=150., reverse_rate=100.) 

chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics = ReactionDynamics(reaction_data=chem_data)   # Notice we're over-writing the earlier "dynamics" object
dynamics.set_conc([50., 0, 0.], snapshot=True)
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Run the reaction

# %%
dynamics.set_diagnostics()         # To save diagnostic information about the call to single_compartment_react()
#dynamics.verbose_list = [1, 2, 3]  # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

# The changes of concentrations vary very rapidly early on; 
# so, we'll be using the dynamic_substeps option to increase time resolution,
# as long as the reaction remains "fast" (based on a threshold of % change, as specified by fast_threshold)
dynamics.single_compartment_react(time_step=0.001, reaction_duration=0.05,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_substeps=4, rel_fast_threshold=10)

# %%
# Continue running the reaction at lover resolution
dynamics.single_compartment_react(time_step=0.002, reaction_duration=0.25,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  dynamic_substeps=4, rel_fast_threshold=10)

# %%
df = dynamics.get_history()
df

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium(tolerance=3)

# %% [markdown] tags=[]
# ## <a name="react_5_plot"></a>Plots of changes of concentration with time

# %%
dynamics.plot_curves(colors=["blue", "green", "red"], 
                     title="Coupled reactions A <-> B and A <-> S (fast but disadvantaged energetically)")

# %% [markdown]
# ### Notice how the "alternate (shunt) path" of the reaction, i.e. S (red)   
# ### has a FAST START (fast kinetics),
# ### but EVENTUALLY PETERS OUT (energy dis-advantage)

# %%
# Verify that the stoichiometry is respected at every reaction step/substep (NOTE: it requires earlier activation of saving diagnostic data)
dynamics.stoichiometry_checker_entire_run()

# %%

# %% [markdown]
# # <a name="down_regulate_1_scenario_3"></a> Scenario 3: 
# ### downregulated by shunt:   
# ### kinetically slow,   
# ### but with thermodynamical advantage (i.e. energetically favored)

# %%
# Specify the chemicals  (notice that we're starting with new objects)
chem_data3 = chem(names=["A", "B", "S"])

# Reaction A <-> B (as before)
chem_data3.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=30., reverse_rate=5.) 

# Reaction A <-> S (slow shunt, excellent thermodynamical energetic advantage)
chem_data3.add_reaction(reactants=["A"], products=["S"],
                       forward_rate=3., reverse_rate=0.1)

chem_data3.describe_reactions()

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics3 = ReactionDynamics(reaction_data=chem_data3)
dynamics3.set_conc([50., 0, 0.], snapshot=True)
dynamics3.describe_state()

# %% [markdown] tags=[]
# ### Run the reaction

# %%
dynamics3.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()
#dynamics3.verbose_list = [1, 2, 3]  # Uncomment for detailed run information (meant for debugging the adaptive variable time step)

# The changes of concentrations vary very rapidly early on; 
# so, we'll be using the dynamic_substeps option to increase time resolution,
# as long as the reaction remains "fast" (based on a threshold of % change, as specified by fast_threshold)
dynamics3.single_compartment_react(time_step=0.005, reaction_duration=0.3,
                                   snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                   dynamic_substeps=5, rel_fast_threshold=10)

# %%
# Continue running the reaction at lover resolution
dynamics3.single_compartment_react(time_step=0.25, reaction_duration=6.7,
                                   snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                   dynamic_substeps=5, rel_fast_threshold=10)

# %%
df3 = dynamics3.get_history()
df3

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics3.is_in_equilibrium(tolerance=13)

# %% [markdown] tags=[]
# ## <a name="react_5_plot"></a>Plots of changes of concentration with time

# %%
dynamics3.plot_curves(colors=["blue", "green", "red"], 
                      title="Coupled reactions A <-> B and A <-> S (slow but with energetic advantage)")

# %% [markdown]
# ### Notice how the "alternate (shunt) path" of the reaction, i.e. S (red)   
# ### has a SLOW START (slow kinetics),
# ### but EVENTUALLY DOMINATES (energy advantage)

# %% [markdown]
# #### Please note the much-longer timescale from the earlier plots
# If we look at the initial [0-0.1] interval, this is what it looks like:

# %%
fig = px.line(data_frame=dynamics3.get_history().loc[:96], x="SYSTEM TIME", y=["A", "B", "S"], 
              title="Same as above, both only showing initial detail",
              color_discrete_sequence = ['blue', 'green', 'red'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
dynamics3.curve_intersection(t_start=0, t_end=0.08, var1="A", var2="B")

# %%
# Verify that the stoichiometry is respected at every reaction step/substep (NOTE: it requires earlier activation of saving diagnostic data)
dynamics.stoichiometry_checker_entire_run()

# %%
