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
# In this "PART 2", we quickly repeat all the steps discussed at length in part1,
# but using a much-finer resolution.
# This time, we'll start with slightly more complex initial concentrations built from 2 superposed sine waves.
#
# **We'll also explore the effects of:**  
# -Spacial resolution ("delta x")  
# -Temporal resolution ("delta t")    
# -Alternate methods of estimating numerical derivatives
#
# LAST REVISED: Oct. 6, 2022

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
# We'll be considering just 1 chemical species, "A"
diffusion_rate = 10.

chem_data = chem(diffusion_rates=[diffusion_rate], names=["A"])

# %% [markdown]
# # BASELINE
# This will be our initial system, whose adherence to the diffusion equation we'll test.
# Afterwards, we'll tweak individual parameters - and observed their effect on the closeness of the approximation

# %%
# Parameters of the simulation run (the diffusion rate got set earlier, and will never vary)
delta_t = 0.01
n_bins = 300
delta_x = 2       # Note that the number of bins also define the fraction of the sine wave cycle in each bin
algorithm = None  # "Explicit, with 3+1 stencil"

# %%
# Initialize the system
bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

# Initialize the concentrations to 2 superposed sine waves
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= "Initial System State",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %% [markdown]
# #### Now do 4 rounds of single-step diffusion, to collect the system state at a total of 5 time points: t0 (the initial state), plus t1, t2, t3 and t4

# %%
# All the system state will get collected in this object
history = MovieArray()

# %%
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Do the 4 rounds of single-step diffusion; accumulate all data in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history.shape   

# %%
# Compute time derivatives (for each bin)
df_dt_all_bins = np.apply_along_axis(np.gradient, 0, all_history, delta_t)
df_dt_all_bins.shape

# %%
# Let's consider the state at the midpoint in time (t2)
f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
f_at_t2.shape

# %%
# Computer the second spacial derivative (simple method)
gradient_x_at_t2 = np.gradient(f_at_t2, delta_x)
second_gradient_x_at_t2 = np.gradient(gradient_x_at_t2, delta_x)
second_gradient_x_at_t2.shape

# %% [markdown]
# #### Now, let's use the computed values to see how the left- and right-hand side of the diffusion equation compare
# at all bins (except 2 at each edge), at the middle simulation time t2

# %%
lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
lhs.shape

# %%
rhs = diffusion_rate*second_gradient_x_at_t2
rhs.shape

# %%
num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end

# %% [markdown]
# The above number is a measure of the discrepancy from the perfect match (zero distance) that an ideal solution would provide. 
#
# Notice how vastly closer to zero it is than the counterpart value we got with the very coarse simulation in experiment "validate_diffusion_1"

# %% [markdown]
# # VARIATIONS on the BASELINE
# We'll now tweak some parameters, and observe the effect on the estimate of the discrepancy from the exact diffusion equation
#
# The baseline distance was **0.023163289760024783**

# %% [markdown]
# ## Variation 1 : reduce spacial resolution
# Expectation: greater discrepancy

# %%
n_bins = 100          # Reducing the spacial resolution

# %%
# Initialize the system
bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

# Initialize the concentrations to 2 superposed sine waves
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)

# %%
history = MovieArray()   # All the system state will get collected in this object
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Do the 4 rounds of single-step diffusion; accumulate all data in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history.shape   

# %%
# Compute time derivatives (for each bin) : simple method, using central differences
df_dt_all_bins = np.apply_along_axis(np.gradient, 0, all_history, delta_t)

# Let's consider the state at the midpoint in time (t2)
f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
f_at_t2.shape

# %%
# Computer the second spacial derivative (simple method, using central differences)
gradient_x_at_t2 = np.gradient(f_at_t2, delta_x)
second_gradient_x_at_t2 = np.gradient(gradient_x_at_t2, delta_x)
second_gradient_x_at_t2.shape

# %%
# Compare the left and right hand sides of the diffusion equation
lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
rhs = diffusion_rate*second_gradient_x_at_t2

num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end

# %% [markdown]
# #### The discrepancy value has indeed gone up from the baseline 0.023163289760024783

# %% [markdown]
# ## Variation 2 : increase spacial resolution
# Expectation: smaller discrepancy

# %%
n_bins = 900          # Increasing the spacial resolution (from the baseline of 300)

# %%
# Initialize the system
bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

# Initialize the concentrations to 2 superposed sine waves
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)

# %%
history = MovieArray()   # All the system state will get collected in this object
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Do the 4 rounds of single-step diffusion; accumulate all data in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history.shape   

# %%
# Compute time derivatives (for each bin) : simple method, using central differences
df_dt_all_bins = np.apply_along_axis(np.gradient, 0, all_history, delta_t)

# Let's consider the state at the midpoint in time (t2)
f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
f_at_t2.shape

# %%
# Computer the second spacial derivative (simple method, using central differences)
gradient_x_at_t2 = np.gradient(f_at_t2, delta_x)
second_gradient_x_at_t2 = np.gradient(gradient_x_at_t2, delta_x)
second_gradient_x_at_t2.shape

# %%
# Compare the left and right hand sides of the diffusion equation
lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
rhs = diffusion_rate*second_gradient_x_at_t2

num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end

# %% [markdown]
# #### The discrepancy value has indeed gone down from the baseline 0.023163289760024783

# %% [markdown]
# ## Variation 3 : reduce time resolution
# Expectation: greater discrepancy

# %%
n_bins = 300          # Restoring the baseline number of bins
delta_t = 0.02        # Reducing the time resolution (from baselin 0.01)

# %%
# Initialize the system
bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

# Initialize the concentrations to 2 superposed sine waves
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)

# %%
history = MovieArray()   # All the system state will get collected in this object
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Do the 4 rounds of single-step diffusion; accumulate all data in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history.shape   

# %%
# Compute time derivatives (for each bin) : simple method, using central differences
df_dt_all_bins = np.apply_along_axis(np.gradient, 0, all_history, delta_t)

# Let's consider the state at the midpoint in time (t2)
f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
f_at_t2.shape

# %%
# Computer the second spacial derivative (simple method, using central differences)
gradient_x_at_t2 = np.gradient(f_at_t2, delta_x)
second_gradient_x_at_t2 = np.gradient(gradient_x_at_t2, delta_x)
second_gradient_x_at_t2.shape

# %%
# Compare the left and right hand sides of the diffusion equation
lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
rhs = diffusion_rate*second_gradient_x_at_t2

num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end

# %% [markdown]
# #### The discrepancy value has indeed gone up from the baseline 0.023163289760024783

# %% [markdown]
# ## Variation 4 : increase time resolution
# Expectation: smaller discrepancy

# %%
delta_t = 0.005        # Increasing the time resolution (from baselin 0.01)

# %%
# Initialize the system
bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

# Initialize the concentrations to 2 superposed sine waves
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)

# %%
history = MovieArray()   # All the system state will get collected in this object
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Do the 4 rounds of single-step diffusion; accumulate all data in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history.shape   

# %%
# Compute time derivatives (for each bin) : simple method, using central differences
df_dt_all_bins = np.apply_along_axis(np.gradient, 0, all_history, delta_t)

# Let's consider the state at the midpoint in time (t2)
f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
f_at_t2.shape

# %%
# Computer the second spacial derivative (simple method, using central differences)
gradient_x_at_t2 = np.gradient(f_at_t2, delta_x)
second_gradient_x_at_t2 = np.gradient(gradient_x_at_t2, delta_x)
second_gradient_x_at_t2.shape

# %%
# Compare the left and right hand sides of the diffusion equation
lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
rhs = diffusion_rate*second_gradient_x_at_t2

num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end

# %% [markdown]
# #### The discrepancy value has indeed gone down from the baseline 0.023163289760024783

# %% [markdown]
# ## Variation 5 : Use a better (higher-order) method to numerically estimate the SPACIAL derivatives
# Expectation: presumably smaller discrepancy (unless dwarfed by errors from approximations in the simulation)

# %%
n_bins = 300          # Restoring the baseline number of bins
delta_t = 0.01        # Reducing the baseline time resolution

# %%
# Initialize the system
bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

# Initialize the concentrations to 2 superposed sine waves
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)

# %%
history = MovieArray()   # All the system state will get collected in this object
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Do the 4 rounds of single-step diffusion; accumulate all data in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history.shape   

# %%
# Compute time derivatives (for each bin) : simple method, using central differences
df_dt_all_bins = np.apply_along_axis(np.gradient, 0, all_history, delta_t)

# Let's consider the state at the midpoint in time (t2)
f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
f_at_t2.shape

# %%
# Computer the second spacial derivative (MORE SOPHISTICATED METHOD, using 5-point stencils)
gradient_x_at_t2 = num.gradient_order4_1d(arr=f_at_t2, dx=delta_x)
second_gradient_x_at_t2 = num.gradient_order4_1d(arr=gradient_x_at_t2, dx=delta_x)
second_gradient_x_at_t2.shape

# %%
# Compare the left and right hand sides of the diffusion equation
lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
rhs = diffusion_rate*second_gradient_x_at_t2

num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end

# %% [markdown]
# #### Not surprisingly, the discrepancy value (in part caused by poor numeric estimation of spacial derivatives) has indeed gone down from the baseline 0.023163289760024783

# %% [markdown]
# ## Variation 6 : Use a better (higher-order) method to numerically estimate the TIME derivatives
# Expectation: presumably smaller discrepancy (unless dwarfed by errors from approximations in the simulation)
#
# (Note: we revert to the original method for the *spacial* derivatives)

# %%
# Initialize the system
bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

# Initialize the concentrations to 2 superposed sine waves
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)

# %%
history = MovieArray()   # All the system state will get collected in this object
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Do the 4 rounds of single-step diffusion; accumulate all data in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history.shape   

# %%
# Compute time derivatives (for each bin) : MORE SOPHISTICATED METHOD, using 5-point stencils
df_dt_all_bins = np.apply_along_axis(num.gradient_order4_1d, 0, all_history, delta_t)

# Let's consider the state at the midpoint in time (t2)
f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
f_at_t2.shape

# %%
# Computer the second spacial derivative (BACK TO SIMPLE METHOD, using central differences)
gradient_x_at_t2 = np.gradient(f_at_t2, delta_x)
second_gradient_x_at_t2 = np.gradient(gradient_x_at_t2, delta_x)
second_gradient_x_at_t2.shape

# %%
# Compare the left and right hand sides of the diffusion equation
lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
rhs = diffusion_rate*second_gradient_x_at_t2

num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end

# %% [markdown]
# #### Here, the discrepancy value has remained largely the same from the baseline 0.023163289760024783

# %%
