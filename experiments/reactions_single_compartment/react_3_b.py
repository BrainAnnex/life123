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
# ## Synthesis/Decomposition  reaction `A + B <-> C`
# #### with 1st-order kinetics for each species, taken to equilibrium.
# #### Compatison of approximate and exact solutions 

# %% [markdown]
# ### TAGS :  "uniform compartment", "numerical"

# %%
LAST_REVISED = "Sep. 5, 2024"
LIFE123_VERSION = "1.0.0.beta.38"      # Version this experiment is based on

# %%
#import set_path            # Using MyBinder?  Uncomment this before running the next cell!
                            # Importing this local file will add the project's home directory to sys.path

# %% tags=[]
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

from life123 import check_version, UniformCompartment, ReactionDynamics

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# # Initialize the System
# Specify the chemicals, the reactions, and the initial concentrations

# %%
# Instantiate the simulator and specify the chemicals
uc = UniformCompartment(preset="fast")

# %%
# Reaction A + B <-> C , with 1st-order kinetics for each species
uc.add_reaction(reactants=["A" , "B"], products="C",
                forward_rate=5., reverse_rate=2.)

uc.describe_reactions()

# %%
# Set the initial concentrations of all the chemicals
uc.set_conc({"A": 10., "B": 50., "C": 20.}, snapshot=True)

# %%

# %%

# %% [markdown] tags=[]
# # Run the reaction

# %%
uc.single_compartment_react(initial_step=0.004, duration=0.06,
                                  variable_steps=True,
                                  snapshots={"initial_caption": "1st reaction step",
                                             "final_caption": "last reaction step"})

# %%
uc.get_history()

# %% [markdown] tags=[]
# ## Plots changes of concentration with time

# %%
uc.plot_history(colors=['red', 'darkorange', 'green'], show_intervals=True)

# %%

# %%

# %% [markdown]
# ## Compare to EXACT, analytic solution

# %%
rxn = uc.get_chem_data().get_reaction(0)

# %%
reactants, products, kF, kR = rxn.unpack_for_dynamics()

# %%
t_arr = uc.get_history(columns="SYSTEM TIME").to_numpy()
t_arr

# %%
A_exact, B_exact, C_exact = ReactionDynamics.approx_solution_combination_rxn(kF, kR, A0=10., B0=50., C0=20., t_arr=t_arr)

# %%
df = uc.get_history()

# %%
df["A_exact"] = A_exact
df

# %%
uc.plot_data(df=df, fields=["A", "A_exact"], xrange=[0, 0.02], colors=["red", "purple"], legend_header="Simulation") 

# %% [markdown]
# ### The match is perhaps passable, but not spectacular.
# That's a result of our simulation run done with the "fast" preset, which yields fewer data point.  
# Let's now repeat the simulation with the "slow" preset

# %%
# Start a new simulator, with the same chemicals and same reactions
uc_slow = UniformCompartment(preset="slow", chem_data=uc.get_chem_data())

# %%
# Again, set the initial concentrations of all the chemicals
uc_slow.set_conc({"A": 10., "B": 50., "C": 20.}, snapshot=True)

# %%
# Re-run the simulation, with the new preset
uc_slow.single_compartment_react(initial_step=0.004, duration=0.06, variable_steps=True)

# %% [markdown]
# ### Notice that we now took 80 steps, vs. the 23 with the "fast" preset

# %%
uc_slow.get_history()

# %%
t_arr = uc_slow.get_history(columns="SYSTEM TIME").to_numpy()
t_arr

# %%
A_exact, B_exact, C_exact = ReactionDynamics.approx_solution_combination_rxn(kF, kR, A0=10., B0=50., C0=20., t_arr=t_arr)

# %%
df = uc_slow.get_history()
df["A_exact"] = A_exact
df

# %%
uc.plot_data(df=df, fields=["A", "A_exact"], xrange=[0, 0.02], colors=["pink", "purple"], legend_header="Simulation") 

# %%

# %%
# Start a new simulator, with the same chemicals and same reactions
uc_fixed = UniformCompartment(chem_data=uc.get_chem_data())

# %%
# Again, set the initial concentrations of all the chemicals
uc_fixed.set_conc({"A": 10., "B": 50., "C": 20.}, snapshot=True)

# %%
# Re-run the simulation, with the new preset
uc_fixed.single_compartment_react(initial_step=0.00001, duration=0.06, variable_steps=False)

# %%
df = uc_fixed.get_history()
df

# %%
t_arr = uc_fixed.get_history(columns="SYSTEM TIME").to_numpy()

# %%

# %%
df.iloc[1234]

# %%

# %%

# %%
A_exact, B_exact, C_exact = ReactionDynamics.approx_solution_combination_rxn(kF, kR, A0=10., B0=50., C0=20., t_arr=t_arr)

# %%
df = uc_fixed.get_history()
df["A_exact"] = A_exact
df

# %%
uc.plot_data(df=df, fields=["A", "A_exact"], xrange=[0, 0.02], colors=["orange", "purple"], legend_header="Simulation") 

# %%
uc.plot_data(df=df, fields=["A", "A_exact"], xrange=[0, 0.02], colors=["orange", "purple"], legend_header="Simulation") 

# %%
