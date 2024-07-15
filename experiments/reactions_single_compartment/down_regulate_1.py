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
# All reactions are 1st order, mostly forward.  Taken to equilibrium.
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta36)

# %% [markdown]
# ## Bathtub analogy:
# A is initially full, while B and S are empty.  
# If the "shunt" S is present, scenario 2 corresponds to a large pipe and a small elevation change...  
# while scenario 3 corresponds to a narrow pipe and a large elevation change.

# %% [markdown]
# ![Downregulated by shunt](../../docs/down_regulate_1.png)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData as chem
from life123 import UniformCompartment

import plotly.express as px
from life123 import GraphicLog

# %% tags=[]
# Initialize the HTML logging (for the graphics)
log_file = get_notebook_basename() + ".log.htm"    # Use the notebook base filename for the log file

# Set up the use of some specified graphic (Vue) components
GraphicLog.config(filename=log_file,
                  components=["vue_cytoscape_2"],
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
# ### Set the initial concentrations of all the chemicals

# %%
dynamics = UniformCompartment(chem_data=chem_data, preset="fast")
dynamics.set_conc(conc={"A": 50.}, snapshot=True)
dynamics.describe_state()

# %% [markdown]
# ### Run the reaction

# %%
dynamics.enable_diagnostics()          # To save diagnostic information about the call to single_compartment_react()

# The changes of concentrations vary very rapidly early on; automated variable timesteps will take care of that
dynamics.single_compartment_react(initial_step=0.001, duration=0.3,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True)

# %%
dynamics.get_history()

# %%
dynamics.explain_time_advance()

# %%
dynamics.plot_history(colors=["darkturquoise", "green"], title="Single reaction A <-> B (no downregulation)",
                      show_intervals=True)

# %% [markdown]
# #### Notice the intersection at the exact midpoint of the 2 initial concentrations (50 and 0):

# %%
dynamics.curve_intersect('A', 'B', t_start=0, t_end=0.1)

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium()

# %%

# %% [markdown]
# # <a name="down_regulate_1_scenario_2"></a>  Scenario 2 - downregulated by shunt: 
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
chem_data.plot_reaction_network("vue_cytoscape_2")

# %%
dynamics = UniformCompartment(chem_data=chem_data, preset="mid")   # Notice we're over-writing the earlier "dynamics" object
dynamics.set_conc(conc={"A": 50.}, snapshot=True)
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Run the reaction

# %%
dynamics.enable_diagnostics()         # To save diagnostic information about the call to single_compartment_react()

# The changes of concentrations vary very rapidly early on; automated variable timesteps will take care of that
dynamics.single_compartment_react(initial_step=0.001, duration=0.3,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True)

# %%
dynamics.plot_history(colors=["darkturquoise", "green", "red"],
                      title="Coupled reactions A <-> B and A <-> S (fast but disadvantaged energetically)",
                      show_intervals=True)

# %% [markdown]
# ### Notice how the "alternate (shunt) path" of the reaction, i.e. S (red)   
# ### has a FAST START (fast kinetics),
# ### but EVENTUALLY PETERS OUT (energy dis-advantage)

# %%
# Verify that all the reactions have reached equilibrium
dynamics.is_in_equilibrium()

# %%

# %% [markdown]
# # <a name="down_regulate_1_scenario_3"></a> Scenario 3 - downregulated by shunt:   
# ### kinetically slow,   
# ### but with thermodynamical advantage (i.e. energetically favored)

# %%
# Specify the chemicals  (notice that we're starting with new objects)
new_chem_data = chem(names=["A", "B", "S"])

# Reaction A <-> B (as before)
new_chem_data.add_reaction(reactants=["A"], products=["B"],
                       forward_rate=30., reverse_rate=5.) 

# Reaction A <-> S (slow shunt, excellent thermodynamical energetic advantage)
new_chem_data.add_reaction(reactants=["A"], products=["S"],
                       forward_rate=3., reverse_rate=0.1)

new_chem_data.describe_reactions()

# %%
dynamics = UniformCompartment(chem_data=new_chem_data, preset="small_rel_change")
dynamics.set_conc(conc={"A": 50.}, snapshot=True)
dynamics.describe_state()

# %% [markdown] tags=[]
# ### Run the reaction

# %%
dynamics.enable_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# The changes of concentrations vary very rapidly early on; automated variable timesteps will take care of that
dynamics.single_compartment_react(initial_step=0.005, duration=7.0,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"},
                                  variable_steps=True)

# %%
dynamics.plot_history(colors=["darkturquoise", "green", "red"],
                      title="Coupled reactions A <-> B and A <-> S (slow but with energetic advantage)")

# %% [markdown]
# ### Notice how the "alternate (shunt) path" of the reaction, i.e. S (red)   
# ### has a SLOW START (slow kinetics),
# ### but EVENTUALLY DOMINATES (energy advantage)

# %%
dynamics.explain_time_advance()

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions are close to equilibrium
dynamics.is_in_equilibrium(tolerance=12)

# %% [markdown]
# ### Please note the **much-longer** timescale from the earlier plots
# If we look at early time interval, this is what it looks like:

# %%
dynamics.plot_history(colors=["darkturquoise", "green", "red"],
                      title="Same plot as above, both only showing initial detail", xrange=[0, 0.3])

# %%
# Look at where the curves intersect
dynamics.curve_intersect("A", "B", t_start=0, t_end=0.1)

# %%
dynamics.curve_intersect("A", "S", t_start=0.1, t_end=0.2)

# %%
