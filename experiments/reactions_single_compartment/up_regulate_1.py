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
# ## `A` up-regulates `B` , 
# ### by being *the limiting reagent* in the reaction `A + X <-> 2B` (mostly forward), where `X` is plentiful
# 1st-order kinetics.   
# If [A] is low, [B] remains low, too.  Then, if [A] goes high, then so does [B].  However, at that point, A can no longer bring B down to any substantial extent.
#
# See also the experiment `1D/reactions/up_regulation_1`

# %%
LAST_REVISED = "Jan. 12, 2025"
LIFE123_VERSION = "1.0.0rc2"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

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
chem_data = chem(names=["A", "X", "B"], plot_colors=['red', 'darkorange', 'green'])

# %%
uc = UniformCompartment(chem_data=chem_data, preset="fast", enable_diagnostics=True)     # We'll use diagnostics to explain the adaptive steps taken during the simulation

# %%
# Reaction A + X <-> 2B , with 1st-order kinetics for all species
uc.add_reaction(reactants=["A" , "X"], products=(2, "B", 1),
                forward_rate=8., reverse_rate=2.)

uc.describe_reactions()

# Send the plot of the reaction network to the HTML log file
uc.plot_reaction_network("vue_cytoscape_2")

# %% [markdown]
# ### Set the initial concentrations of all the chemicals

# %%
uc.set_conc(conc={"A": 5., "X": 100.})      # A is scarce, X is plentiful, B is absent
uc.describe_state()

# %% [markdown] tags=[]
# ### Take the initial system to equilibrium

# %%
uc.single_compartment_react(initial_step=0.0005, duration=0.015)

# %%
uc.plot_history(show_intervals=True)

# %% [markdown]
# **A, as the scarse limiting reagent, stops the reaction.  
# As long as A is low, B also remains low.**

# %%
uc.get_history()

# %%
uc.get_diagnostics().explain_time_advance()

# %% [markdown]
# ### <a name="sec_equilibrium"></a>Equilibrium

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium()

# %%

# %%

# %% [markdown] tags=[]
# # STEP 2. Now, let's suddenly increase [A]

# %%
uc.set_single_conc(species_name="A", conc=50.)
uc.describe_state()

# %%
uc.get_history(tail=5)

# %% [markdown]
# Notice the discontinuity in [A] in the last 2 rows

# %% [markdown]
# ### Again, take the system to equilibrium

# %%
uc.single_compartment_react(initial_step=0.0005, target_end_time=0.035)

# %%
uc.plot_history(show_intervals=True)

# %% [markdown]
# `A`, still the **limiting reagent**, is again eventually stopping the reaction...    
# The (transiently) high value of [A] led to a high value of [B]  
# Notice how the variable time steps automatically adapt.

# %%
#uc.get_history()

#uc.get_diagnostics().explain_time_advance()

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium()

# %%

# %%

# %% [markdown] tags=[]
# # STEP 3.  Let's again suddenly increase [A]

# %%
uc.set_single_conc(species_name="A", conc=30.)
uc.describe_state()

# %%
uc.get_history(tail=5)

# %% [markdown]
# ### Yet again, take the system to equilibrium

# %%
uc.single_compartment_react(initial_step=0.0005, target_end_time=0.070)

# %%
uc.plot_history(show_intervals=True)

# %% [markdown]
# `A`, again the scarce limiting reagent, stops the reaction yet again.  
# And, again, the (transiently) high value of [A] up-regulated [B]
#
# Notes:   
# `A` can up-regulate `B`, but it cannot bring it down.  
# `X` will soon need to be replenished, if `A` is to continue being the limiting reagent.**

# %%
# Zoom in on the early part of the last plot
uc.plot_history(range_x=[0, 0.02])

# %%
# Look up the first intersection of the [A] and [B] curves
uc.curve_intersect("A", "B", t_start=0, t_end=0.01)

# %%
# The next intersection isn't particulary meaning; just an artifact of the sudden jump in [A]
uc.curve_intersect("A", "B", t_start=0.01, t_end=0.02)

# %%
# The final intersection
uc.curve_intersect("A", "B", t_start=0.016, t_end=0.02)

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium()

# %%
