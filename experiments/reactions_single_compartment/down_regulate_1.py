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
# LAST REVISED: June 4, 2023

# %% [markdown]
# ## Bathtub analogy:
# A is initially full, while B and S are empty.  
# If the "shunt" S is present, scenario 2 corresponds to a large pipe and a small elevation change...  
# while scenario 3 corresponds to a narrow pipe and a large elevation change.

# %% [markdown]
# ![Downregulated by shunt](../../docs/down_regulate_1.png)

# %%

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.chemicals.chem_data import ChemData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import plotly.express as px
from src.modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_1"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %%

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
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc([50., 0.], snapshot=True)
dynamics.describe_state()

# %% [markdown]
# ### Run the reaction

# %%
dynamics.set_diagnostics()          # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.5, high=1.0, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
dynamics.set_step_factors(upshift=1.5, downshift=0.5, abort=0.5)
dynamics.set_error_step_factor(0.5)

# The changes of concentrations vary very rapidly early on; automated variable timesteps will take care of that
dynamics.single_compartment_react(initial_step=0.001, reaction_duration=0.3,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True, explain_variable_steps=False)

# %%
dynamics.get_history()

# %%
dynamics.explain_time_advance()

# %%
dynamics.plot_curves(colors=["blue", "green"], title="Single reaction A <-> B (no downregulation)", 
                     show_intervals=True)

# %% [markdown]
# #### Notice the intersection at the exact midpoint of the 2 initial concentrations (50 and 0):

# %%
dynamics.curve_intersection('A', 'B', t_start=0, t_end=0.1)

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium(tolerance=2)

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
dynamics = ReactionDynamics(chem_data=chem_data)   # Notice we're over-writing the earlier "dynamics" object
dynamics.set_conc([50., 0, 0.], snapshot=True)
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Run the reaction

# %%
dynamics.set_diagnostics()         # To save diagnostic information about the call to single_compartment_react()

# All of these settings are currently close to the default values... but subject to change; set for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.5, high=1.0, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.05, high=0.5, abort=1.5)
dynamics.set_step_factors(upshift=1.4, downshift=0.5, abort=0.5)
dynamics.set_error_step_factor(0.333)

# The changes of concentrations vary very rapidly early on; automated variable timesteps will take care of that
dynamics.single_compartment_react(initial_step=0.001, reaction_duration=0.3,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True, explain_variable_steps=False)

# %%
dynamics.plot_curves(colors=["blue", "green", "red"], 
                     title="Coupled reactions A <-> B and A <-> S (fast but disadvantaged energetically)",
                     show_intervals=True)

# %% [markdown]
# ### Notice how the "alternate (shunt) path" of the reaction, i.e. S (red)   
# ### has a FAST START (fast kinetics),
# ### but EVENTUALLY PETERS OUT (energy dis-advantage)

# %%
dynamics.explain_time_advance()

# %%
dynamics.get_history()

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium()

# %%

# %% [markdown]
# # <a name="down_regulate_1_scenario_3"></a> Scenario 3: 
# ### downregulated by shunt:   
# ### kinetically slow,   
# ### but with thermodynamical advantage (i.e. energetically favored)

# %%
# Specify the chemicals  (notice that we're starting with new objects)
chem_data = chem(names=["A", "B", "S"])

# Reaction A <-> B (as before)
chem_data.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=30., reverse_rate=5.) 

# Reaction A <-> S (slow shunt, excellent thermodynamical energetic advantage)
chem_data.add_reaction(reactants=["A"], products=["S"],
                       forward_rate=3., reverse_rate=0.1)

chem_data.describe_reactions()

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc([50., 0, 0.], snapshot=True)
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Run the reaction

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# These settings can be tweaked to make the time resolution finer or coarser
dynamics.set_thresholds(norm="norm_A", low=2.0, high=5.0, abort=10.0)
dynamics.set_thresholds(norm="norm_B", low=0.008, high=0.5, abort=2.0)    # The "low" value here seems especially critical to fend off instabilities
dynamics.set_step_factors(upshift=1.5, downshift=0.25, abort=0.25)
dynamics.set_error_step_factor(0.2)

# The changes of concentrations vary very rapidly early on; automated variable timesteps will take care of that
dynamics.single_compartment_react(initial_step=0.005, reaction_duration=7.0,
                                   snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                   variable_steps=True, explain_variable_steps=False)

# %%
dynamics.plot_curves(colors=["blue", "green", "red"], 
                      title="Coupled reactions A <-> B and A <-> S (slow but with energetic advantage)")

# %% [markdown]
# ### Notice how the "alternate (shunt) path" of the reaction, i.e. S (red)   
# ### has a SLOW START (slow kinetics),
# ### but EVENTUALLY DOMINATES (energy advantage)

# %%
dynamics.explain_time_advance()

# %%
dynamics.get_history()

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium(tolerance=12)

# %% [markdown]
# #### Please note the much-longer timescale from the earlier plots
# If we look at early time interval, this is what it looks like:

# %%
fig = px.line(data_frame=dynamics.get_history().loc[:250], x="SYSTEM TIME", y=["A", "B", "S"], 
              title="Same plot as above, both only showing initial detail",
              color_discrete_sequence = ["blue", "green", "red"],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
# Look at where the curves intersect
dynamics.curve_intersection("A", "B", t_start=0, t_end=0.1)

# %%
dynamics.curve_intersection("A", "S", t_start=0.1, t_end=0.2)

# %%
