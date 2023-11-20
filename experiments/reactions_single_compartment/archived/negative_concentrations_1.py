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
# ## Tour of Life123's process of automatic detection - and stepsize reduction - in case any concentration goes negative
# ### Case study of two coupled reactions: `2 S <-> U` and `S <-> X` (both mostly forward)
#
# 1st-order kinetics throughout.   
#
# LAST REVISED: June 4, 2023   **THIS IS AN ARCHIVED EXPERIMENT**
#
# * [RUN 1 : negative concentrations are detected but otherwise ignored](#negative_concentrations_1_run_1)  
# (_"DEMONIC"_ in figure below)
#
# * [RUN 2 : we restored some, but not all, of the code that had been disabled in Run #1](#negative_concentrations_1_run_2)   
# Negative concentrations from *individual* reactions are now automatically corrected - but negative concentrations from *combined* (synergistic) reactions can still slip thru   
# (_"POSSESSED"_ in figure below)
#
# * [RUN 3 : we restored ALL the code that had been disabled in the earlier runs](#negative_concentrations_1_run_3)   
# Negative concentrations from *individual* reactions are now automatically corrected - and so are negative concentrations from *combined* (synergistic) reactions   
# (_"DISTURBED"_ in figure below)
#
# * [RUN 4 : same as the previous run, but with slightly finer time resolution](#negative_concentrations_1_run_4)   
# (_"HEALING"_ in figure below)
#
# For even more accurate solutions - _"HEALTHY"_ in figure below - see the experiment `large_time_steps_2`
#

# %% [markdown]
# ![Exorcising Instabilities](../../../docs/negative_concentrations_1.png)

# %% [markdown]
# # IMPORTANT: DO NOT ATTEMPT TO RUN THIS NOTEBOOK!   
# ## This is a **"frozen run"** that depends on Life123 _code changes_ for demonstration purposes  
# (a code change that turned off the automatic correction of time steps that lead to negative concentrations.)  
# If you bypass the execution exit in the first cell, and run the other cells, you WILL NOT REPLICATE the results below!

# %%
# To stop the current and subsequent cells: USED TO PREVENT ACCIDENTAL RUNS OF THIS NOTEBOOK!

class StopExecution(Exception):
    def _render_traceback_(self):
        return []

raise StopExecution     # See: https://stackoverflow.com/a/56953105/5478830

# %% [markdown]
# ## The `StopExecution` above is on purpose, TO PREVENT ACCIDENTAL RUNS OF THIS NOTEBOOK!

# %%

# %% tags=[]

from src.modules.chemicals.chem_data import ChemData as chem
from src.modules.reactions.reaction_dynamics import ReactionDynamics

import plotly.express as px
import plotly.graph_objects as go

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

# %% [markdown]
#

# %% [markdown]
# # <a name="negative_concentrations_1_run_1"></a>RUN 1 : negative concentrations are detected but otherwise ignored
# This run required the disabling of multiple software features that detect, and automatically correct, such issues.   
# It cannot be replicated with current versions of Life123.  DO NOT ATTEMPT TO RE-RUN!

# %%
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# ************  DO NOT RUN : Life123 automatically catches and corrects these errors;
#                            this run required a code change for demonstration purposes! 
#                            If you run it, you WILL NOT REPLICATE the results below.
#                            This is meant to be a "frozen run"
dynamics.single_compartment_react(time_step=0.1, stop_time=0.8)

dynamics.explain_time_advance()

# %%
df = dynamics.get_history()
df

# %% [markdown]
# ## Notice how all concentrations dips into negative values various times!

# %%
fig0 = dynamics.plot_history(colors=['green', 'orange', 'blue'], suppress=True)   # Prepare, but don't show, the main plot

# Add a second plot, with a horizontal red line at concentration = 0
fig1 = px.line(x=[0,0.9], y=[0,0], color_discrete_sequence = ['red'])

# Combine the plots, and display them
all_fig = go.Figure(data=fig0.data + fig1.data, layout = fig0.layout)    # Note that the + is concatenating lists
all_fig.update_layout(title="Changes in concentrations; notice the several dips into negative concentrations (red line)")
all_fig.show()

# %% [markdown]
# # There are 3 separate scenarios that lead to negative concentrations
#
# [More information](https://life123.science/reactions)
#
# #### Scenario 1 : A reaction causes a dip into the negative, and the combined other reactions fail to remedy it
#
# #### Scenario 2 : A reaction causes a dip into the negative, and the other reactions counterbalance it enough to remedy it
# (though "counterbalanced", this is still regarded as a sign of instability)
#
# #### Scenario 3 : No single reaction causes any negative concentration, but - combined - they do

# %% [markdown]
# **Example of scenario 1** (all the necessary diagnostic data in the next several rows)     
# The step from t=0.1 to 0.2 :   
# [S] is initially 50   
# * Reaction 0 causes a Delta_S of -64.000000 , which would send it into the negative    
# (leading to message _"*** DETECTING NEGATIVE CONCENTRATION in chemical `S` from reaction 2 S <-> U"_)      
# * Reaction 1 causes a Delta_S of -9.000000 , which is fine by itself (no error message) but fails to remedy the negative (in fact, makes it worse)  
#
# The net result is [S] = 50 -64 -9 = -23 at the final t=0.2   
# That's what had led to the message _"+++ SYSTEM STATE ERROR: FAILED TO CATCH negative concentration upon advancing reactions from system time t=0.1"_

# %% [markdown]
# **Example of scenario 2** (all the necessary diagnostic data in the next several rows)     
# The step from t=0.4 to 0.5 :  
# [S] is initially -67.99  
# * Reaction 0 causes a Delta_S of 146.96 , which would remedy the negative when combined with the initial value      
# * Reaction 1 causes a Delta_S of 63.927 , sending [S] into the negative (or, in this case, keeping it into the negative);   
# (leading to message _"*** DETECTING NEGATIVE CONCENTRATION in chemical `S` from reaction S <-> X"_)
#
# The net result is [S] = -67.99 +146.96 +63.927 = 142.897 at the final t=0.5    
# The value is now positive, so there's no SYSTEM STATE ERROR

# %% [markdown]
# **Example of scenario 3**    
# This scenario doesn't occur in this run, but scroll down to RUN 2 - where it occurs.   
# An example would be if the initial concentration were 20, and each of the two reactions caused a delta_concentration of -15 ;   
# combined, they'll bring to concentration to  20 -15 -15 = -10

# %%
df

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=1)

# %%

# %% [markdown]
# ##################################################################################

# %% [markdown]
# # <a name="negative_concentrations_1_run_2"></a>RUN 2 : we restored some, but not all, of the code that had been disabled in Run #1.  
# ## Negative concentrations from *individual* reactions are now automatically corrected - but negative concentrations from *combined* (synergistic) reactions can still slip thru
# ### (With the restored code) Life123 automatically detects whenever an individual reaction would cause any chemical concentration to go negative.   
# At that point, it raises an `ExcessiveTimeStep` exception.
# That exception gets caught - which results in the following:   
# * the (partial) execution of the last step gets discarded
# * the step size gets halved
# * the run automatically resumes from the previous time
#
# Such an automatic detection and remediation eliminates "Scenarios 1 and 2" (see earlier in notebook for definitions) BUT NOT "scenario 3"  
# IMPORTANT: scenario 3 normally gets caught and remedied, too - but that feature got disabled by a code change in this run, _FOR DEMONSTRATION PURPOSES_
#
# This run required the disabling of multiple software features that detect, and automatically correct, such issues.   
# It cannot be replicated with current versions of Life123.  **DO NOT ATTEMPT TO RE-RUN!**

# %%
# Same as for Run #1
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

# ************  DO NOT RUN : Life123 automatically catches and corrects these errors;
#                            this run required a code change for demonstration purposes! 
#                            If you run it, you WILL NOT REPLICATE the results below.
#                            This is meant to be a "frozen run"
dynamics.single_compartment_react(time_step=0.1, stop_time=0.8)

dynamics.explain_time_advance()

# %% [markdown]
# ### Notice how the system automatically temporarily slowed down from the requested time step of 0.1 at t=0.1 and again at t=0.65  
# Those actions intercepted, and automatically remedied, some - but NOT all - of the negative concentrations (explained below in the notebook.)  
# We now took a total of 10 steps, instead of the 9 ones of Run #1

# %%
df = dynamics.get_history()
df

# %% [markdown]
# ## Notice how [S] dips into negative values at t=0.55!
# But all the multitude of other negative dips we had in Run #1 are now gone :)

# %%
fig0 = dynamics.plot_history(colors=['green', 'orange', 'blue'], suppress=True)   # Prepare, but don't show, the main plot

# Add a second plot, with a horizontal red line at concentration = 0
fig1 = px.line(x=[0,0.9], y=[0,0], color_discrete_sequence = ['red'])

# Combine the plots, and display them
all_fig = go.Figure(data=fig0.data + fig1.data, layout = fig0.layout)    # Note that the + is concatenating lists
all_fig.update_layout(title="Changes in concentrations; notice how S dips into negative concentrations (red line) at t=0.55")
all_fig.show()

# %% [markdown]
# ## With this partially-restored functionality, we're averting what we called "scenarios 1 and 2" earlier in the notebook - but are still afflicted by "scenario 3" (no single reaction causes any negative concentration, but - combined - they do)

# %% [markdown]
# **That's exactly what happens to [S] at t=0.55**   (all the necessary diagnostic data in the next few rows)       
# The step from t=0.45 to 0.55 :   
# [S] is initially 37.453500    
# * Reaction 0 causes a Delta_S of -36.445600 , which by itself would cause no problem    
# * Reaction 1 causes a Delta_S of  -8.928150 , by itself would cause no problem - but leads to a negative value when combined with reaction 0!    
#
# The net result is [S] = 37.453500 -36.445600 -8.928150 = -7.92025 at the final t=0.55   
# That's what had led to the message _"+++ SYSTEM STATE ERROR: FAILED TO CATCH negative concentration upon advancing reactions from system time t=0.45"_

# %%
df

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=0)

# %%
dynamics.get_diagnostic_rxn_data(rxn_index=1)

# %%

# %% [markdown]
# ##################################################################################

# %% [markdown]
# # <a name="negative_concentrations_1_run_3"></a>RUN 3 : we restored ALL the code that had been disabled in the earlier runs.
# ## Negative concentrations from *individual* reactions are now automatically corrected - and so are negative concentrations from *combined* (synergistic) reactions
# ### (With the restored code) Life123 automatically detects whenever any chemical concentration goes negative from any cause.   
# At that point, it raises an `ExcessiveTimeStep` exception.
# That exception gets caught - which results in the following:   
# * the (partial) execution of the last step gets discarded
# * the step size gets halved
# * the run automatically resumes from the previous time

# %%
# Same as for the previous runs
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
dynamics.describe_state()

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.1, stop_time=0.8)

dynamics.explain_time_advance()

# %% [markdown]
# ### Notice how the system automatically temporarily slowed down from the requested time step of 0.1 at t=0.1 and again at t=0.45  
# Those actions intercepted, and automatically remedied, ALL the negative concentrations, whether caused by any single reaction, or by the cumulative effect of multiple ones.  
# We now took a total of 10 steps, instead of the 9 ones of Run #1

# %%
df = dynamics.get_history()
df

# %% [markdown]
# ## Notice how negative concentrations are no longer seen anywhere :)

# %%
fig0 = dynamics.plot_history(colors=['green', 'orange', 'blue'], suppress=True)   # Prepare, but don't show, the main plot

# Add a second plot, with a horizontal red line at concentration = 0
fig1 = px.line(x=[0,0.9], y=[0,0], color_discrete_sequence = ['red'])

# Combine the plots, and display them
all_fig = go.Figure(data=fig0.data + fig1.data, layout = fig0.layout)    # Note that the + is concatenating lists
all_fig.update_layout(title="Changes in concentrations; notice negative concentrations are totally absent now")
all_fig.show()

# %% [markdown]
# ## With the completely-restored functionality, we're averting what we called "scenarios 1, 2 and 3" earlier in the notebook :)

# %%

# %% [markdown]
# ##################################################################################

# %% [markdown]
# # <a name="negative_concentrations_1_run_4"></a>RUN 4 : same as the previous run, but with slightly finer time resolution
#
# In run 3, we demonstrated catching - and resolving - all cases of negative concentrations that would slip in if it weren't for Life123's detection and automatic remediation.  But we still have a lot of instability, especially in [S].
#
# A slight improvement in the time steps will make a profound difference in this run

# %%
# Same as for the previous runs
dynamics = ReactionDynamics(chem_data=chem_data)
dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})
dynamics.describe_state()

# %% [markdown]
# ## Rather than decreasing (for the whole run) the time step of 0.1, we'll simply opt to AUTOMATICALLY take 2 substeps whenever any of the reactions are "fast changing" based on a threshold we provide

# %%
dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

dynamics.single_compartment_react(time_step=0.1, stop_time=0.8, 
                                  dynamic_substeps=2, abs_fast_threshold=80.0)

dynamics.explain_time_advance()

# %% [markdown]
# ### Note how the system automatically splits the time steps into substeps for up to t=0.3 
# that's when most of the change is occurring

# %%
df = dynamics.get_history()
df

# %%
dynamics.plot_history(colors=['green', 'orange', 'blue'])

# %% [markdown]
# ## The resolution is still coarse, especially up to t=0.1, but all the instability is gone :)

# %% [markdown]
# ### Note: For even more accurate solutions, see the experiment `large_time_steps_2`

# %%
