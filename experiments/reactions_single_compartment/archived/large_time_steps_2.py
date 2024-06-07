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
# # From finer to coarser resolution in advancing the two coupled reactions: 
# ### `2 S <-> U` and `S <-> X` (both mostly forward)
#
# 1st-order kinetics throughout.  
#
# Notes:  
# * for an exploration of instabilities, see the experiment `negative_concentrations_1`
# * for an accurate longer run of the same reactions (with the same initial conditions), see the experiment "up_regulate_3"
#
# LAST REVISED: Dec. 3, 2023  *** THIS IS AN ARCHIVED EXPERIMENT ***  (lightly updated May 6, 2024)

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
chem_data = chem(names=["U", "X", "S"])

# Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants=[(2, "S")], products="U",
                       forward_rate=8., reverse_rate=2.)

# Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants="S", products="X",
                       forward_rate=6., reverse_rate=3.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %% [markdown]
# # Run 1 : extremely small fixed time steps (no substeps)

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
#dynamics.describe_state()

dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.001, stop_time=0.8)

df = dynamics.get_history()
#df
dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=['green', 'orange', 'blue'])

# %% [markdown]
# # Run 2 : very small time steps, with dynamic substeps

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
#dynamics.describe_state()

dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.01, stop_time=0.8,
                                  dynamic_substeps=4, rel_fast_threshold=50)

df = dynamics.get_history()
#df
dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=['green', 'orange', 'blue'])

# %% [markdown]
# # Run 3 : small-ish time steps, with dynamic substeps

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
#dynamics.describe_state()

dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.08, stop_time=0.8,
                                  dynamic_substeps=4, rel_fast_threshold=250)

df = dynamics.get_history()
#df
dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=['green', 'orange', 'blue'])

# %% [markdown]
# # Run 4 : same as previous run, but fewer dynamic substeps

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
#dynamics.describe_state()

dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.08, stop_time=0.8,
                                  dynamic_substeps=2, rel_fast_threshold=150)

df = dynamics.get_history()
#df
dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=['green', 'orange', 'blue'])

# %% [markdown]
# # Run 5 : same as previous run, but slightly larger primary step

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
#dynamics.describe_state()

dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.1, stop_time=0.8,
                                  dynamic_substeps=2, abs_fast_threshold=80.0)

df = dynamics.get_history()
#df
dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ### Note: the above run is identical to the last run in the experiment `negative_concentrations_1`

# %% [markdown]
# # Run 6 : same as previous run, but no substeps

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
#dynamics.describe_state()

dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.1, stop_time=0.8)

df = dynamics.get_history()
#df
dynamics.explain_time_advance()

# %% [markdown]
# ### Notice the automated detection - and correction - of negative concentrations arising from the excessively-large time steps

# %%
dynamics.plot_history(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ### Note: for an exploration of instabilities, see the experiment `negative_concentrations_1`

# %%
