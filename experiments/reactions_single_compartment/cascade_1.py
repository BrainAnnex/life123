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
# ## 2 COUPLED reactions of different speeds, forming a "cascade":  
# ### `A <-> B` (fast) and `B <-> C` (slow)
# All 1st order. Taken to equilibrium. Both reactions are mostly forward.
# The concentration of the intermediate product B manifests 1 oscillation (transient "overshoot")
#
# (Adaptive variable time resolution is used, with extensive diagnostics.)
#
# LAST REVISED: June 4, 2023

# %% [markdown]
# ## Bathtub analogy:
# A is initially full, while B and C are empty.  
# Tubs are progressively lower (reactions are mostly forward.)  
# A BIG pipe connects A and B: fast kinetics.  A small pipe connects B and C: slow kinetics. 
#
# INTUITION: B, unable to quickly drain into C while at the same time being blasted by a hefty inflow from A,  
# will experience a transient surge, in excess of its final equilibrium level.
#
# * [Compare with the final reaction plot (the red line is B)](#cascade_1_plot)

# %% [markdown]
# ![2 Coupled Reactions](../../docs/2_coupled_reactions.png)

# %%

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.chemicals.chem_data import ChemData
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

# %% [markdown]
# # Initialize the System
# Specify the chemicals and the reactions

# %% tags=[]
# Specify the chemicals
chem_data = ChemData(names=["A", "B", "C"])

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
dynamics.get_history()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()         # To save diagnostic information about the call to single_compartment_react()

# These settings can be tweaked to make the time resolution finer or coarser
dynamics.set_thresholds(norm="norm_A", low=0.5, high=1.0, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.2, high=0.5, abort=1.5)
dynamics.set_step_factors(upshift=1.4, downshift=0.5, abort=0.5)
dynamics.set_error_step_factor(0.333)

dynamics.single_compartment_react(initial_step=0.02, reaction_duration=0.4,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True, explain_variable_steps=False)

# %% [markdown]
# ### <a name="cascade_1_plot"></a> Plots of changes of concentration with time

# %%
dynamics.plot_curves(title="Coupled reactions A <-> B and B <-> C",
                     colors=['blue', 'red', 'green'])

# %%
dynamics.curve_intersection("A", "B", t_start=0, t_end=0.05)

# %%
dynamics.curve_intersection("A", "C", t_start=0, t_end=0.05)

# %%
dynamics.curve_intersection("B", "C", t_start=0.05, t_end=0.1)

# %%
dynamics.get_history()

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium()

# %%

# %% [markdown]
# # EVERYTHING BELOW IS FOR DIAGNOSTIC INSIGHT
#
# ### Perform some verification

# %% [markdown]
# ### Take a peek at the diagnostic data saved during the earlier reaction simulation

# %%
# Concentration increments due to reaction 0 (A <-> B)
# Note that [C] is not affected
dynamics.get_diagnostic_rxn_data(rxn_index=0)

# %%
# Concentration increments due to reaction 1 (B <-> C)
# Note that [A] is not affected.  
# Also notice that the 0-th row from the A <-> B reaction isn't seen here, because that step was aborted
# early on, before even getting to THIS reaction
dynamics.get_diagnostic_rxn_data(rxn_index=1)

# %%

# %% [markdown]
# ### Provide a detailed explanation of all the steps/substeps of the reactions, from the saved diagnostic data

# %%

# %% [markdown] tags=[]
# # Re-run with very small constanst steps

# %% [markdown]
# We'll use constant steps of size 0.0005, which is 1/4 of the smallest steps (the "substep" size) previously used in the variable-step run

# %%
dynamics2 = ReactionDynamics(reaction_data=chem_data)

# %% tags=[]
dynamics2.set_conc([50., 0, 0.], snapshot=True)

# %%
dynamics2.single_compartment_react(initial_step=0.0005, reaction_duration=0.4,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  )      

# %%
fig = px.line(data_frame=dynamics2.get_history(), x="SYSTEM TIME", y=["A", "B", "C"], 
              title="Coupled reactions A <-> B and B <-> C (run with constant steps)",
              color_discrete_sequence = ['blue', 'red', 'green'],
              labels={"value":"concentration", "variable":"Chemical"})
fig.show()

# %%
dynamics2.curve_intersection(t_start=0, t_end=0.05, chem1="A", chem2="B")

# %%
dynamics2.curve_intersection(t_start=0, t_end=0.05, chem1="A", chem2="C")

# %%
dynamics2.curve_intersection(t_start=0.05, t_end=0.1, chem1="B", chem2="C")

# %%
df2 = dynamics2.get_history()
df2

# %% [markdown]
# ## Notice that we now did 801 steps - vs. the 56 of the earlier variable-resolution run!

# %% [markdown]
# ### Let's compare some entries with the coarser previous variable-time run

# %% [markdown]
# ## Let's compare the plots of [B] from the earlier (variable-step) run, and the latest (high-precision, fixed-step) one:

# %%
# Earlier run (using variable time steps)
fig1 = px.line(data_frame=dynamics.get_history(), x="SYSTEM TIME", y=["B"], 
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Variable-step run"})
fig1.show()

# %%
# Latest run (high-precision result from fine fixed-resolution run)
fig2 = px.line(data_frame=dynamics2.get_history(), x="SYSTEM TIME", y=["B"], 
              color_discrete_sequence = ['orange'],
              labels={"value":"concentration", "variable":"Fine fixed-step run"})
fig2.show()

# %%
import plotly.graph_objects as go

all_fig = go.Figure(data=fig1.data + fig2.data)    # Note that the + is concatenating lists
all_fig.update_layout(title="The 2 runs contrasted")
all_fig['data'][0]['name']="B (variable steps)"
all_fig['data'][1]['name']="B (fixed high precision)"
all_fig.show()

# %%
