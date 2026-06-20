# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.2
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

# %% [markdown]
# ### TAGS :  "uniform compartment"

# %%

# %% [markdown]
# ## Bathtub analogy:
# A is initially full, while B and S are empty.  
# If the "shunt" S is present, scenario 2 corresponds to a large pipe and a small elevation change...  
# while scenario 3 corresponds to a narrow pipe and a large elevation change.

# %% [markdown]
# ![Downregulated by shunt](../../docs/down_regulate_1.png)

# %%
LAST_REVISED = "June 20, 2026"
LIFE123_VERSION = "1.0.0rc8"     # Library version this experiment is based on

# %%
#import set_path                 # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import ipynbname
from IPython.display import IFrame

from life123 import check_version, SpeciesRegistry, UniformCompartment

import plotly.express as px

# %%
check_version(LIFE123_VERSION)

# %%
# Initialize the HTML logging (for the graphics)
log_file = ipynbname.name() + ".log.htm"    # Use the notebook base filename for the log file
                                            # IN CASE OF PROBLEMS, set manually to any desired name
log_file

# %%

# %%

# %% [markdown]
# # <a name="down_regulate_1_scenario_1"></a>  Scenario 1: A <-> B in the absence of the 2nd reaction

# %% [markdown]
# ### Initialize the System

# %%
# Instantiate the simulator and specify the chemicals.
# The diagnostics will be give insight into the inner workings of the simulation
uc = UniformCompartment(names=["A", "B"], preset="fast", enable_diagnostics=True)

# Reaction A <-> B , with 1st-order kinetics in both directions
uc.add_reaction(reactants="A", products="B", kF=30., kR=5.)

uc.describe_reactions()

# %%

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
uc.set_conc(conc={"A": 50.})
uc.describe_state()

# %%

# %% [markdown]
# ### Run the reaction

# %%
# The changes of concentrations vary very rapidly early on; automated variable timesteps will take care of that
uc.single_compartment_react(initial_step=0.001, duration=0.3)

# %%
uc.get_history()

# %%
uc.get_diagnostics().explain_time_advance()

# %%
uc.plot_history(colors=["darkturquoise", "green"],
                title="Single reaction A <-> B (no downregulation)",
                show_intervals=True)

# %% [markdown]
# #### Notice the intersection at the exact midpoint of the 2 initial concentrations (50 and 0):

# %%
uc.curve_intersect('A', 'B', t_start=0, t_end=0.1)

# %%
# Verify that all the reactions have reached equilibrium
uc.is_in_equilibrium()

# %%
# Notice that the plot colors that we assigned got saved with the species data
uc.get_species_data().as_dataframe()

# %%

# %%

# %% [markdown]
# # <a name="down_regulate_1_scenario_2"></a>  Scenario 2 - downregulated by shunt: 
# ### kinetically fast,   
# ### but with pooor thermodynamical advantage (i.e. not very favored energetically)

# %%
# Add the reaction A <-> S (fast shunt, poor thermodynical energetic advantage)
uc.add_reaction(reactants="A", products="S",
                kF=150., kR=100.) 

uc.describe_reactions()

# %%
# Send a plot of the network of reactions to the HTML log file
uc.plot_reaction_network(log_file=log_file)

# %%
IFrame(log_file, width=1000, height=600)         # You may also open the log file in a browser

# %%

# %%
reaction_data = uc.get_reactions()

# %%
uc = UniformCompartment(reactions=reaction_data, preset="mid", enable_diagnostics=True)   
# Notice we're over-writing the earlier "uc" object, and re-using the reaction data

uc.set_conc(conc={"A": 50.})
uc.describe_state()

# %% [markdown]
# ### Run the reaction

# %%
# The changes of concentrations vary very rapidly early on; automated variable timesteps will take care of that
uc.single_compartment_react(initial_step=0.001, duration=0.3)

# %%
#species_data.set_value(species_id="S", field_name="plot_color", value="red")

# %%
uc.plot_history(add_colors="red", title="Coupled reactions A <-> B and A <-> S (fast but disadvantaged energetically)",
                show_intervals=True)

# %% [markdown]
# ### Notice how the "alternate (shunt) path" of the reaction, i.e. S (red)   
# ### has a FAST START (fast kinetics),
# ### but EVENTUALLY PETERS OUT (energy dis-advantage)

# %%
# Verify that all the reactions have reached equilibrium
uc.is_in_equilibrium()

# %%

# %%

# %% [markdown]
# # <a name="down_regulate_1_scenario_3"></a> Scenario 3 - downregulated by shunt:   
# ### kinetically slow,   
# ### but with thermodynamical advantage (i.e. energetically favored)

# %%
# Start with new UniformCompartment object, but re-use the earlier species data
species_data = uc.get_species_data()

species_data.as_dataframe()
# Notice that the new color that we assigned for the species `S`, in the last plot, got saved with its species data

# %%
uc_new = UniformCompartment(species_data=species_data, preset="small_rel_change",
                            enable_diagnostics=True)

# Reaction A <-> B (just as before)
uc_new.add_reaction(reactants="A", products="B",
                    kF=30., kR=5.) 

# Reaction A <-> S (slow shunt, excellent thermodynamical energetic advantage)
uc_new.add_reaction(reactants="A", products="S",
                    kF=3., kR=0.1)

uc_new.describe_reactions()


# %%
uc_new.set_conc(conc={"A": 50.})
uc_new.describe_state()

# %% [markdown]
# ### Run the reaction

# %%
# The changes of concentrations vary very rapidly early on; automated variable timesteps will take care of that
uc_new.single_compartment_react(initial_step=0.005, duration=7.0)

# %%
# The plot now simply uses the previously-assigned colors
uc_new.plot_history(title="Coupled reactions A <-> B and A <-> S (slow but with energetic advantage)")

# %% [markdown]
# ### Notice how the "alternate (shunt) path" of the reaction, i.e. S (red)   
# ### has a SLOW START (slow kinetics),
# ### but EVENTUALLY DOMINATES (energy advantage)

# %%
uc_new.get_diagnostics().explain_time_advance()

# %% [markdown]
# ### Check the final equilibrium

# %%
# Verify that all the reactions are close to equilibrium
uc_new.is_in_equilibrium(tolerance=12)

# %% [markdown]
# ### Please note the **much-longer** timescale from the earlier plots
# If we look at early time interval, this is what it looks like:

# %%
uc_new.plot_history(title="Same plot as above, both only showing initial detail", range_x=[0, 0.3])

# %%
# Look at where the curves intersect
uc_new.curve_intersect("A", "B", t_start=0, t_end=0.1)

# %%
uc_new.curve_intersect("A", "S", t_start=0.1, t_end=0.2)

# %%
