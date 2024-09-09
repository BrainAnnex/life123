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
# ## The reversible Synthesis/Decomposition reaction `A + B <-> C`
# #### with 1st-order kinetics for each species, taken to equilibrium.
# #### Comparison of 2 approximate solutions and exact solutions 

# %% [markdown]
# ### TAGS :  "uniform compartment", "numerical"

# %%
LAST_REVISED = "Sep. 8, 2024"
LIFE123_VERSION = "1.0.0.beta.38"      # Version this experiment is based on

# %%
#import set_path            # Using MyBinder?  Uncomment this before running the next cell!
                            # Importing this local file will add the project's home directory to sys.path

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

import pandas as pd
from life123 import check_version, UniformCompartment, ReactionDynamics, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# # 1. Run a simulation with low accuracy

# %% [markdown]
# ### Initialize the System
# Specify the chemicals, the reactions, and the initial concentrations

# %%
# Instantiate the simulator and specify the chemicals
# Here we use the "fast" preset for the variable steps, which leads to fewer steps, but generally less-accurate results
uc_fast = UniformCompartment(preset="fast")    

# %%
# Reaction A + B <-> C , with 1st-order kinetics for each species
uc_fast.add_reaction(reactants=["A" , "B"], products="C",
                     forward_rate=5., reverse_rate=2.)

uc_fast.describe_reactions()

# %%
# Set the initial concentrations of all the chemicals
uc_fast.set_conc({"A": 10., "B": 50., "C": 20.}, snapshot=True)

# %%

# %% [markdown] tags=[]
# ### Run the reaction

# %%
uc_fast.single_compartment_react(initial_step=0.004, duration=0.06, variable_steps=True)

# %%
df_fast = uc_fast.get_history()
df_fast

# %% [markdown] tags=[]
# ### Plots changes of concentrations with time

# %%
uc_fast.plot_history(colors=['red', 'darkorange', 'green'], show_intervals=True)

# %%

# %%

# %% [markdown]
# # 2. Let's now repeat the simulation with the "slow" preset, which yields more data points
# ### and, generally, more accuracy

# %%
# Instantiate the simulator and specify the chemicals
uc_slow = UniformCompartment(preset="slow")

# %%
# Reaction A + B <-> C , with 1st-order kinetics for each species
uc_slow.add_reaction(reactants=["A" , "B"], products="C",
                     forward_rate=5., reverse_rate=2.)

# Set the initial concentrations of all the chemicals
uc_slow.set_conc({"A": 10., "B": 50., "C": 20.}, snapshot=True)

# %%
uc_slow.single_compartment_react(initial_step=0.004, duration=0.06, variable_steps=True)

# %% [markdown]
# ### Note that **80 steps** were now used, instead of the earlier 23

# %%
df_slow = uc_slow.get_history()
df_slow

# %%

# %%

# %% [markdown]
# # 3. And, finally, get the EXACT analytical solution

# %%
rxn = uc_slow.get_chem_data().get_reaction(0)

# %%
reactants, products, kF, kR = rxn.unpack_for_dynamics()

# %%
# We'll use the same, larger, number of time points as in the "slow" simulation of step 2
t_arr = uc_slow.get_history(columns="SYSTEM TIME").to_numpy()
t_arr

# %%
# The EXACT, analytical solution
A_exact, B_exact, C_exact = ReactionDynamics.exact_solution_combination_rxn(kF, kR, A0=10., B0=50., C0=20., t_arr=t_arr)

# %%
C_exact

# %%
df_exact = pd.DataFrame({
    'SYSTEM TIME': t_arr,
    'C_exact': C_exact
})

df_exact

# %%

# %%

# %% [markdown]
# # 4. Let's compare plots for the concentration of `C`, as a function of time, from the earlier 3 simulations

# %%
p1 = uc_fast.plot_history(chemicals="C", colors=['#F5B914'], title="fast (less precise)", show=False)
p2 = uc_slow.plot_history(chemicals="C", colors=['#DBF514'], title="slow (more precise)",  show=False)
p3 = PlotlyHelper.plot_pandas(df=df_exact, fields=["SYSTEM TIME", "C_exact"], colors=['forestgreen'], title="exact", show=False)

# %%
PlotlyHelper.combine_plots([p1, p2, p3], title="Comparison of simulation accuracies", xrange=[0, 0.02])

# %% [markdown]
# ## The comparison reveals a gradient "less precise" -> "more precise" -> "exact" solution

# %%
