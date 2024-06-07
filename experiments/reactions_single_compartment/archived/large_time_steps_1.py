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
# ## Excessively large time steps - and remediation
# ### Based on experiment `cycles_1`
# #### A cycle of reactions `A <-> B <-> C <-> A` and `C + E_High <-> A + E_Low`
#
# ### [RUN 1](#large_time_steps_1_run1) - Initially, very large time steps are taken - at the edge of pushing some concentrations into negative values.
# The system automatic detects those problems, intercepts the problematic steps and re-runs them with 1/2 the time step.   
# Negative concentrations are automatically avoided, but nonetheless the plots are ragged... and the solutions are eventually unstable.
#
# ### [RUN 2](#large_time_steps_1_run2) - Same primary steps as for run #1, but with the option of using 1/2 substeps as needed, 
# with thresholds that lead to those substeps being actually utilized a fair part of the time.   
# The raggedness and instabilities are now eliminated.
# (Note: the 1/2 substeps are on a per-reaction basis)
#
# LAST REVISED: Dec. 3, 2023     *** THIS IS AN ARCHIVED EXPERIMENT ***  (lightly updated May 6, 2024)

# %% [markdown]
# # IMPORTANT: DO NOT ATTEMPT TO RUN THIS NOTEBOOK!   
# ## This is a **"frozen run"** that depends on an old version of Life123, for demonstration purposes  
# #### (newer versions of Life123 tend to recover from instability more gracefully!!)  
# If you bypass the execution exit in the first cell, and run the other cells, you WILL NOT REPLICATE the results below!

# %%
# To stop the current and subsequent cells: USED TO PREVENT ACCIDENTAL RUNS OF THIS NOTEBOOK!

class StopExecution(Exception):
    def _render_traceback_(self):
        return []

raise StopExecution     # See: https://stackoverflow.com/a/56953105/5478830

# %%

# %%

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from src.modules.chemicals.chem_data import ChemData as chem
from src.modules.reactions.reaction_dynamics import UniformCompartment

from src.modules.visualization.graphic_log import GraphicLog

# %% tags=[]
# Initialize the HTML logging
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
                  extra_js="https://cdnjs.cloudflare.com/ajax/libs/cytoscape/3.21.2/cytoscape.umd.js")

# %% [markdown]
# ### Initialize the system

# %%
# Initialize the system
chem_data = chem(names=["A", "B", "C", "E_high", "E_low"])

# Reaction A <-> B, mostly in forward direction (favored energetically)
# Note: all reactions in this experiment have 1st-order kinetics for all species
chem_data.add_reaction(reactants="A", products="B",
                       forward_rate=9., reverse_rate=3.)

# Reaction B <-> C, also favored energetically
chem_data.add_reaction(reactants="B", products="C",
                       forward_rate=8., reverse_rate=4.)

# Reaction C + E_High <-> A + E_Low, also favored energetically, but kinetically slow
# Note that, thanks to the energy donation from E, we can go "upstream" from C, to the higher-energy level of "A"
chem_data.add_reaction(reactants=["C" , "E_high"], products=["A", "E_low"],
                       forward_rate=1., reverse_rate=0.2)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}  # Note the abundant energy source "E_high"
initial_conc

# %% [markdown] tags=[]
# ### We'll split each simulation in three segments (apparent from the graphs of the runs, further down):  
# Time [0-0.03] fast changes   
# Time [0.03-5.] medium changes   
# Time [5.-8.] slow changes, as we approach equilibrium  

# %% [markdown]
# ## <a name="large_time_steps_1_run1"></a> RUN 1 - Initially, very large time steps are taken - at the edge of pushing some concentrations into negative values.
# The system automatic detects those problems, intercepts the problematic steps and re-runs them with 1/2 the time step.   
# Negative concentrations are automatically avoided, but nonetheless the plots are ragged... and the solutions are eventually unstable.

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc=initial_conc, snapshot=True)
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# %%
#dynamics.verbose_list = ["substeps"]     # Uncomment for debugging information

# %%
dynamics.single_compartment_react(time_step=0.0012, stop_time=0.03)
#dynamics.get_history()

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Notice how the system automatically intervened and reduced some steps in 1/2, to prevent negative concentrations

# %%
dynamics.single_compartment_react(time_step=0.0025, stop_time=5.)
#dynamics.get_history()

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Notice how the system automatically periodically intervened and reduced some steps in 1/2, to prevent negative concentrations.
# The time step that was originally requested hovers near the cusp of being so large as to lead to negative concentrations

# %%
dynamics.single_compartment_react(time_step=0.008, stop_time=8.)

# %%
df = dynamics.get_history()
df

# %% [markdown]
# ### Notice that with our aggressively-large time steps, we ended up with _far fewer_ data points than in any of the runs in the experiment `cycles_1`, which this analysis is based on

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Yet again, the system automatically keeps periodically throttling down the time steps, to prevent negative concentrations.
# At times, especially near the end, the throttling gets very frequent ("frantic"!)

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_history(chemicals=["E_high", "E_low"], colors=["red", "grey"])

# %%
dynamics.plot_history(chemicals=["A", "B", "C"])

# %% [markdown]
# # The negative concentrations were automatically avoided, but nonetheless the plots are ragged... and the solutions display instabilities

# %%

# %% [markdown]
# ## <a name="large_time_steps_1_run2"></a> RUN 2 - Same primary steps as for run #1, but with the option of using 1/2 substeps as needed, 
# with thresholds that lead to those substeps being actually utilized a fair part of the time.   
# The raggedness and instabilities are now eliminated.
# (Note: the 1/2 substeps are on a per-reaction basis)

# %%
dynamics = UniformCompartment(chem_data=chem_data)   # Note: OVER-WRITING the "dynamics" object
dynamics.set_conc(conc=initial_conc, snapshot=True) 
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# %%
#dynamics.verbose_list = ["substeps"]     # Uncomment for debugging information

# %%
dynamics.single_compartment_react(time_step=0.0012, stop_time=0.03,
                                  dynamic_substeps=2, abs_fast_threshold=750.)
#dynamics.get_history()

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Compared to the earlier run, the system is now at times reaching for smaller 1/4 steps 
# (1/2 factors arise from either the negative concentrations or the optional substeps triggered by large changes; sometime both occur)

# %%
dynamics.single_compartment_react(time_step=0.0025, stop_time=5.,
                                  dynamic_substeps=2, abs_fast_threshold=250.)
#dynamics.get_history()

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Again, we see many of those smaller 1/4 steps

# %%
dynamics.single_compartment_react(time_step=0.008, stop_time=8.,
                                  dynamic_substeps=2, abs_fast_threshold=2.)

# %%
df = dynamics.get_history()
df

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### In the later times (the 5.-8. range) there isn't that "frantic throttling up and down" that we had before

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_history(chemicals=["E_high", "E_low"], colors=["red", "grey"])

# %%
dynamics.plot_history(chemicals=["A", "B", "C"])

# %% [markdown]
# # The raggedness, and instability in the solutions, is now gone :)

# %%
