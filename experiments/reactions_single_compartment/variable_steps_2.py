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
# ## Accurate results in the simulation of the 2 coupled reactions `2 S <-> U` and `S <-> X`   
# Both mostly forward.  1st-order kinetics throughout.   
#
# Same as `variable_steps_1`, but with FIXED TIME STEPS: a lot of TINY steps - as a proxy for the "exact value"
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta34.1)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData as chem
from life123 import UniformCompartment

from life123 import GraphicLog

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
chem_data.add_reaction(reactants=[(2, "S", 1)], products="U",
                       forward_rate=8., reverse_rate=2.)

# Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
chem_data.add_reaction(reactants="S", products="X",
                       forward_rate=6., reverse_rate=3.)

chem_data.describe_reactions()

# Send the plot of the reaction network to the HTML log file
chem_data.plot_reaction_network("vue_cytoscape_2")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100.})
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(initial_step=0.0001, target_end_time=0.35, 
                                  variable_steps=False)

df = dynamics.get_history()
df

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Look up some values

# %%
df.iloc[400]

# %%
df.iloc[1850]

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_history(colors=['green', 'orange', 'darkturquoise'])

# %%
dynamics.curve_intersect("U", "X", t_start=0.3, t_end=0.35)

# %%
