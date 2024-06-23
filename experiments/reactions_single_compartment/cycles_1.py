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
# ## A cycle of reactions between `A`, `B` and `C` :  (1) `A <-> B <-> C `
# #### the "closing" of the above cycle (the "return" path from `C` to `A`) is coupled with an "energy donor" reaction:
# ## (2) `C + E_High <-> A + E_Low`
# #### where `E_High` and `E_Low` are, respectively, the high- and low- energy molecules that drive the cycle (for example, think of ATP/ADP).   
# Comparisons are made between results obtained with 3 different time resolutions.
#
# All 1st-order kinetics.    
#
# LAST REVISED: June 23, 2024 (using v. 1.0 beta36)

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %% tags=[]
from experiments.get_notebook_info import get_notebook_basename

from life123 import ChemData
from life123 import UniformCompartment
from life123 import Numerical as num

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
chem_data = ChemData(names=["A", "B", "C", "E_high", "E_low"])

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
initial_conc = {"A": 100., "E_high": 1000.}  # Note the abundant energy source "E_high"; anything not specified will default to zero
initial_conc

# %%

# %%

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true jp-MarkdownHeadingCollapsed=true tags=[] jp-MarkdownHeadingCollapsed=true tags=[]
# # Run # 1 : FIXED time resolution, with COARSE time steps - broken up in 3 time intervals    
# (trial and error, not shown, reveals that increasing any of the time steps below, leads to "excessive time step" errors;  
# note: fixed time resolution is generelly NOT recommended, except for double-checks and error analysis)

# %%
dynamics = UniformCompartment(chem_data=chem_data)
dynamics.set_conc(conc=initial_conc, snapshot=True)
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# %% [markdown] tags=[]
# ### We'll split the simulation in three segments (apparent from the graphs of the runs, further down):  
# Time [0-0.03] fast changes   
# Time [0.03-5.] medium changes   
# Time [5.-8.] slow changes, as we approach equilibrium  

# %%
dynamics.single_compartment_react(initial_step=0.0008, target_end_time=0.03, variable_steps=False)
#dynamics.get_history()

# %%
dynamics.single_compartment_react(initial_step=0.001, target_end_time=5., variable_steps=False)
#dynamics.get_history()

# %%
dynamics.single_compartment_react(initial_step=0.005, target_end_time=8., variable_steps=False)

# %%
dynamics.get_history()

# %% [markdown]
# ### Notice we created 5,609 data points, in 3 hand-picked segments at different time resolutions

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_history(chemicals=["E_high", "E_low"], colors=["red", "grey"])

# %%
dynamics.plot_history(chemicals=["A", "B", "C"])

# %% [markdown]
# ### The plots has 4 distinctive intersections; locate them and save them for later comparisons across repeated runs:

# %%
run1 = []

# %%
run1.append(dynamics.curve_intersect("E_high", "E_low", t_start=1., t_end=2.))

# %%
run1.append(dynamics.curve_intersect("A", "B", t_start=2.31, t_end=2.33))

# %%
run1.append(dynamics.curve_intersect("A", "C", t_start=3., t_end=4.))

# %%
run1.append(dynamics.curve_intersect("B", "C", t_start=3., t_end=4.))

# %%
run1

# %% [markdown]
# ### Verify the final equilibrium state:

# %%
dynamics.is_in_equilibrium()

# %%

# %%

# %%

# %% [markdown] tags=[]
# # _NOTE: Everything below, to the end, is JUST A REPEAT of the same experiment, with different time steps, for accuracy comparisons_

# %%

# %% [markdown]
# # Run # 2. VARIABLE time resolution

# %%
dynamics = UniformCompartment(chem_data=chem_data, preset="mid_inclusive")   # Note: OVER-WRITING the "dynamics" object  # heavy_brakes
dynamics.set_conc(conc=initial_conc, snapshot=True) 
dynamics.describe_state()

# %%
# These will affect the dynamic automated choices of time steps, 
# and were specified by the preset used in instantiating the "dynamics" object
dynamics.show_adaptive_parameters()  

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# %%
dynamics.single_compartment_react(initial_step=0.0001, target_end_time=8.0,
                                  variable_steps=True, explain_variable_steps=None)

# %% [markdown] tags=[]
# ### Notice we created 1,911 data points, a good deal less than in run #1

# %%
dynamics.get_history()

# %%
dynamics.plot_history(chemicals=["E_high", "E_low"], colors=["red", "grey"])

# %%
dynamics.plot_history(chemicals=["A", "B", "C"])

# %%
# Show the timestepe taken (vertical dashed lines) in a small section of the plot
dynamics.plot_history(chemicals=["A", "B", "C"], show_intervals=True, xrange=[2.5, 3.5])

# %%
dynamics.explain_time_advance()  # Notice a lot of timestep adjustments, as needed

# %% [markdown]
# ### Just like we saw in the earlier run, the two plots have 4 distinctive intersections among them; locate and save them:

# %%
run2 = []

# %%
run2.append(dynamics.curve_intersect(t_start=1., t_end=2., chem1="E_high", chem2="E_low"))  # This can be seen in the 1st plot

# %%
run2.append(dynamics.curve_intersect(t_start=2.31, t_end=2.33, chem1="A", chem2="B"))       # This, and later, can be seen in the 2nd plot

# %%
run2.append(dynamics.curve_intersect(t_start=3., t_end=4., chem1="A", chem2="C"))

# %%
run2.append(dynamics.curve_intersect(t_start=3., t_end=4., chem1="B", chem2="C"))

# %%
run2

# %%

# %%

# %% [markdown]
# # Run # 3. FIXED time resolution, using FINER fixed resolution than in run #1 

# %%
dynamics = UniformCompartment(chem_data=chem_data)   # Note: OVER-WRITING the "dynamics" object
dynamics.set_conc(conc=initial_conc, snapshot=True) 
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# %% tags=[]
dynamics.single_compartment_react(initial_step=0.0004, target_end_time=0.03, variable_steps=False)

# %%
dynamics.single_compartment_react(initial_step=0.0005, target_end_time=5., variable_steps=False)

# %%
dynamics.single_compartment_react(initial_step=0.0025, target_end_time=8., variable_steps=False)

# %%
dynamics.get_history()

# %% [markdown]
# ### Notice we created 11,215 data points, A LOT MORE than in runs #1 and #2

# %%
dynamics.explain_time_advance()

# %% [markdown] tags=[]
# ## Plots of changes of concentration with time

# %%
dynamics.plot_history(chemicals=["E_high", "E_low"], colors=["red", "grey"])

# %%
dynamics.plot_history(chemicals=["A", "B", "C"])

# %% [markdown]
# ### Yet again, the two plots have 4 distinctive intersections among them; locate and save them:

# %%
run3 = []

# %%
run3.append(dynamics.curve_intersect(t_start=1., t_end=2., chem1="E_high", chem2="E_low"))

# %%
run3.append(dynamics.curve_intersect(t_start=2.31, t_end=2.33, chem1="A", chem2="B"))

# %%
run3.append(dynamics.curve_intersect(t_start=3., t_end=4., chem1="A", chem2="C"))

# %%
run3.append(dynamics.curve_intersect(t_start=3., t_end=4., chem1="B", chem2="C"))

# %%
run3

# %%

# %%

# %% [markdown]
# # Finally, compare (using a Euclidean metric) the discrepancy between the runs
# Run #3, at high resolution, could be thought of as "the closest to the actual values"

# %%
# Discrepancy of run1 from run3
num.compare_results(run1, run3)

# %%
# Discrepancy of run2 from run3
num.compare_results(run2, run3)

# %% [markdown]
# Run 2 (with a lot fewer points than run 1), not surprisingly deviates a bit more from the (presumably) highly-accurate run 3

# %% [markdown]
# #### The coordinates of the 4 critical points, from the 3 different runs, are pretty similar to one another - as can be easily seen:

# %%
run1

# %%
run2

# %%
run3

# %%
