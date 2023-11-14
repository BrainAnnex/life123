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
# Taken to equilibrium. Both reactions are mostly forward. All 1st order.  
#
# CONTINUATION of experiment "cascade_1"; please refer to it for background information
#
# LAST REVISED: Nov. 11, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

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

chem_data.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
graph_data = chem_data.prepare_graph_network()
GraphicLog.export_plot(graph_data, "vue_cytoscape_1")

# %% [markdown]
# # Start the simulation

# %%
dynamics = ReactionDynamics(chem_data=chem_data)

# %% [markdown]
# ### Set the initial concentrations of all the chemicals, in their index order

# %%
dynamics.set_conc([50., 0, 0.], snapshot=True)

# %%
dynamics.describe_state()

# %% [markdown] tags=[]
# ## Run the reaction

# %%
dynamics.set_diagnostics()         # To save diagnostic information about the call to single_compartment_react()

# These settings can be tweaked to make the time resolution finer or coarser.  
# The values are the current defaults, but setting them here for repeatability
dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
dynamics.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4)
dynamics.set_error_step_factor(0.25)

dynamics.single_compartment_react(initial_step=0.01, reaction_duration=0.5,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True, explain_variable_steps=False)

# %% [markdown]
# ### <a name="cascade_1_plot"> Plots of changes of concentration with time</a>
# Notice the variable time steps (vertical dashed lines)

# %%
dynamics.plot_curves(title="Coupled reactions A <-> B and B <-> C",
                     colors=['blue', 'red', 'green'], show_intervals=True)

# %% [markdown]
# # Now, let's consider the scenario that we don't know anything about the intermediate state (B), and we only see the overall reaction `A <-> C'

# %% [markdown]
# Let's look at just the plots of the varying concentrations with time of `A` and `C` :

# %%
dynamics.plot_curves(chemicals=['A', 'C'], title="Overall reactions A <-> C",
                     colors=['blue', 'green'], show_intervals=False)

# %%
dynamics.curve_intersection("A", "C", t_start=0, t_end=0.05)

# %% [markdown]
# ### Let's look at the final equilibrium

# %%
A_final = dynamics.get_chem_conc("A")
A_final

# %%
C_final = dynamics.get_chem_conc("C")
C_final

# %% [markdown]
# #### Their ratio:

# %%
C_final / A_final

# %% [markdown]
# ### As expected the equilibrium constant for the overall reaction `A <-> C` (approx. 48, assuming everything is still first-order, which we'll discuss below) is indeed the product of the equilibrium constants of the two elementary reactions (K = 8 and K = 6, respectively) that we saw earlier; to repeat:

# %%
chem_data.describe_reactions()

# %%

# %%
dynamics.get_history(columns = ["SYSTEM TIME", "A", "C"])

# %%

# %%
