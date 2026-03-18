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
# ### `A <-> 2B` elementary reversible decomposition reaction, taken to equilibrium.  
# ### Examine State Space trajectory, using [A] and [B] as state variables

# %% [markdown]
# ### TAGS :  "uniform compartment"

# %%
LAST_REVISED = "Mar. 17, 2026"
LIFE123_VERSION = "1.0.0rc7"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import check_version, UniformCompartment, ReactionKinetics, PlotlyHelper

# %%

# %% [markdown]
# ### Initialize the system

# %%
uc = UniformCompartment(names=["A", "B"])
uc.set_conc(conc={"A": 40., "B": 0.})

# %%
# Elementary reaction A <-> 2B
uc.add_reaction(reactants="A", products=[(2,"B")], kF=3., kR=2.)

uc.describe_reactions()

# %%

# %%

# %% [markdown]
# ### Simulate the reaction to equilibrium

# %%
uc.single_compartment_react(initial_step=0.005, duration=0.08, variable_steps=False)

# %%
df = uc.get_history()
df

# %%
# Verify that the reaction has reached equilibrium
uc.is_in_equilibrium(tolerance=2)

# %%
uc.plot_history(colors=['navy', 'orange'])

# %%

# %% [markdown]
# ## Same data, but shown differently

# %%
PlotlyHelper.plot_pandas(df=uc.get_history(), x_var="B", fields="A", 
                         title="State space of reaction A <-> 2B : [A] vs. [B]", 
                         colors="#C83778", 
                         annotation_field="SYSTEM TIME", show_points=True, 
                         annotate={"dx": 5, "dy": -13, "x_offset": -50, "y_offset": -8, "frequency": 2})

# %% [markdown]
# ### **Note how the trajectory is progressively slowing down towards the dynamical system's "attractor" (equilibrium state of the reaction)**

# %%
