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
# ## Carry out a few single indivual steps of diffusion, 
# ### and directly verify that the values satisfy the diffusion equation
#
# In this "PART 3", we perform all the steps done in part2,
# with an even finer resolution, and more complex initial concentrations,
# repeated for 2 different diffusion algorithms.
#
# LAST REVISED: Nov. 28, 2022

# %%
# Extend the sys.path variable, to contain the project's root directory
import set_path
set_path.add_ancestor_dir_to_syspath(3)  # The number of levels to go up 
                                         # to reach the project's home, from the folder containing this notebook

# %%
from experiments.get_notebook_info import get_notebook_basename

from life_1D.bio_sim_1d import BioSim1D
from modules.reactions.reaction_data import ReactionData as chem
from modules.movies.movies import MovieArray
from modules.numerical.numerical import Numerical as num

import numpy as np

import plotly.express as px
import plotly.graph_objects as go

# %%
# Parameters of the simulation run.  We'll be considering just 1 chemical species, "A"
diffusion_rate = 10.
delta_t = 0.01
n_bins = 5000
delta_x = 2       # Note that the number of bins also define the fraction of the sine wave cycle in each bin

# %%
chem_data = chem(diffusion_rates=[diffusion_rate], names=["A"])

# %% [markdown]
# # ALGORITHM 1

# %%
algorithm = None   # "Explicit, with 3+1 stencil"

# %%
# Initialize the system
bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

# Initialize the concentrations to 2 superposed sine waves
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=12, bias=40)
bio.inject_sine_conc(species_name="A", frequency=2, amplitude=10)
bio.inject_sine_conc(species_name="A", frequency=16, amplitude=5)

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= "Initial System State (for the tiny system)",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# #### Now do 4 rounds of single-step diffusion, to collect the system state at a total of 5 time points: t0 (the initial state), plus t1, t2, t3 and t4

# %%
history = MovieArray()   # All the system state will get collected in this object
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %% tags=[]
# Do the 4 rounds of single-step diffusion; accumulate all data in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history.shape   

# %%
# Compute time derivatives (for each bin), using 5-point stencils
df_dt_all_bins = np.apply_along_axis(num.gradient_order4_1d, 0, all_history, delta_t)

# Let's consider the state at the midpoint in time (t2)
f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
f_at_t2.shape

# %%
# Computer the second spacial derivative, using 5-point stencils
gradient_x_at_t2 = num.gradient_order4_1d(arr=f_at_t2, dx=delta_x)
second_gradient_x_at_t2 = num.gradient_order4_1d(arr=gradient_x_at_t2, dx=delta_x)
second_gradient_x_at_t2.shape

# %% tags=[]
# Compare the left and right hand sides of the diffusion equation
lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
rhs = diffusion_rate*second_gradient_x_at_t2

num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end

# %% [markdown]
# The above number is a measure of the discrepancy from the perfect match (zero distance) that an ideal solution would provide. 

# %% [markdown]
# # ALGORITHM 2

# %%
algorithm = "5_1_explicit"   # "Explicit, with 5+1 stencil"

# %%
# Initialize the system
bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

# Initialize the concentrations to 2 superposed sine waves
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=12, bias=40)
bio.inject_sine_conc(species_name="A", frequency=2, amplitude=10)
bio.inject_sine_conc(species_name="A", frequency=16, amplitude=5)

# %% [markdown]
# #### Now do 4 rounds of single-step diffusion, to collect the system state at a total of 5 time points: t0 (the initial state), plus t1, t2, t3 and t4

# %%
history = MovieArray()   # All the system state will get collected in this object
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %% tags=[]
# Do the 4 rounds of single-step diffusion; accumulate all data in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history.shape   

# %%
# Compute time derivatives (for each bin), using 5-point stencils
df_dt_all_bins = np.apply_along_axis(num.gradient_order4_1d, 0, all_history, delta_t)

# Let's consider the state at the midpoint in time (t2)
f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
f_at_t2.shape

# %%
# Computer the second spacial derivative, using 5-point stencils
gradient_x_at_t2 = num.gradient_order4_1d(arr=f_at_t2, dx=delta_x)
second_gradient_x_at_t2 = num.gradient_order4_1d(arr=gradient_x_at_t2, dx=delta_x)
second_gradient_x_at_t2.shape

# %% tags=[]
# Compare the left and right hand sides of the diffusion equation
lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
rhs = diffusion_rate*second_gradient_x_at_t2

num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end

# %% [markdown]
# # Both algorithms show good measures of accuracy

# %%
