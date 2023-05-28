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
# In this "PART 1", we'll be looking at a tiny system with just 10 bins, to easily inspect the numbers directly.
# We'll use concentrations initialized to an upward-shifted sine wave with 1 cycle over the length system
#
# LAST REVISED: May 27, 2023

# %%
import set_path      # Importing this module will add the project's home directory to sys.path

# %%
from experiments.get_notebook_info import get_notebook_basename

from src.life_1D.bio_sim_1d import BioSim1D
from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.movies.movies import MovieArray
from src.modules.numerical.numerical import Numerical as num

import numpy as np

import plotly.express as px
import plotly.graph_objects as go

# %%
# We'll be considering just 1 chemical species, "A"
diffusion_rate = 10.

chem_data = chem(diffusion_rates=[diffusion_rate], names=["A"])

# %%
# Initialize the system with just a few bins
bio = BioSim1D(n_bins=10, chem_data=chem_data)

# %%
# Initialize the concentrations to a sine wave with 1 cycle over the system
# (with a bias to always keep it > 0)
bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)

# %%
bio.show_system_snapshot()

# %%
fig = px.line(data_frame=bio.system_snapshot(), y=["A"], 
              title= "Initial System State (for the tiny system)",
              color_discrete_sequence = ['red'],
              labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})
fig.show()

# %%
# Show as heatmap
fig = px.imshow(bio.system_snapshot().T, 
                title= "Initial System State (for the tiny system)", 
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto=False, color_continuous_scale="gray_r") 

fig.data[0].xgap=1
fig.data[0].ygap=1

fig.show()

# %% [markdown]
# ### Now do 4 rounds of single-step diffusion, to collect the system state at a total of 5 time points: 
# #### t0 (the initial state), plus t1, t2, t3 and t4

# %%
# All the system states will get collected in this object
history = MovieArray()

# %%
# Store the initial state
arr = bio.lookup_species(species_index=0, copy=True)
history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Take a look at what got stored so far (a matrix whose only row is the initial state)
history.get_array()

# %%
# Additional parameters of the simulation run (the diffusion rate got set earlier)
delta_t = 0.01
delta_x = 2       # Note that the number of bins also define the fraction of the sine wave cycle in each bin
algorithm = None  # "Explicit, with 3+1 stencil"

# %%
# Do the 4 rounds of single-step diffusion; show the system state after each step, and accumulate all data
# in the history object
for _ in range(4):
    bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)
    bio.describe_state(concise=True)

    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

# %%
# Now, let's examine the data collected at the 5 time points
all_history = history.get_array()
all_history

# %% [markdown]
# #### Each row in the above matrix is a state snapshot at a different time:
# first row is the initial state at time 0; the successive rows are at t1, t2, t3, t4

# %%
# Let's consider the state, i.e. the concentration across the bins, at the midpoint in time (t2)
f_at_t2 = all_history[2]
f_at_t2

# %% [markdown]
# If one compares the above state (at t2) with the initial state (the top row in the matrix), 
# one can see, for example, that the leftmost bin's concentration is increasing ("pulled up" by its neighbor to the right):
# 50. has now become 50.28881576
#
# The rightmost bin's concentration is decreasing ("pulled down" by its neighbor to the left):
# 44.12214748 has now become 43.94505274

# %% [markdown]
# ## The diffusion equation states that the partial derivative of the concentration values with respect to time  must equal the (diffusion rate) times (the 2nd partial derivative  with respect to space).
# Let's see if that is the case for the values at time t2!  We are picking t2 because we need at least a value at the earlier time, and a value at the later time

# %% [markdown]
# ## A. The 2nd partial derivative  with respect to space

# %%
# A simple-minded way of computing the 2nd spacial derivative is to use the Numpy gradient function TWICE across the x value
# (that function approximates the derivative using differences)
gradient_x_at_t2 = np.gradient(f_at_t2, delta_x)   # Start with taking the first derivative
gradient_x_at_t2    # This will be the partial derivative of the concentration with respect to x, at time t2

# %% [markdown]
# For example, the 2nd entry in the above array of estimated derivatives, can be manually checked from the 1st and 3rd value in the function f_at_t2(x), as follows:

# %%
(59.32979676 - 50.28881576) / (2*delta_x)    # This way of numerically estimating derivatives is called "Central Differences"

# %% [markdown]
# #### Now take the derivative again, with the Numpy gradient function, to arrive at a coarse estimate of the 2nd derivative with respect to x:

# %%
second_gradient_x_at_t2 = np.gradient(gradient_x_at_t2, delta_x)
second_gradient_x_at_t2

# %% [markdown]
# Note how the 2nd derivative is 0 at bin 5 (bins are numbered 0 thru 9): if you look at the earlier sine plot, x5 is the inflection point.

# %% [markdown]
# ### B. The partial derivative with respect to time
# Now, let's look at how concentrations change with time.  Let's first revisit the full history:

# %%
all_history

# %% [markdown]
# #### For simplicity, let's start by just inspecting how the values change over time at the *3rd bin* from the left, 
# i.e. the 3rd *column* (index 2 because counting starts at 0) of the above matrix:

# %%
f_of_t_at_x2 = all_history[ : , 2]  # The column with index 2
f_of_t_at_x2

# %% [markdown]
# ### The above a function of time 
# (we took an entry from each row - the rows being the system states at t0, t1, t2, t3, t4); let's look at its time derivative:

# %%
gradient_t_at_x2 = np.gradient(f_of_t_at_x2, delta_t)
gradient_t_at_x2

# %% [markdown]
# ### The above is the rate of change of the concentration, as the diffusion proceeds, at the position x2
# At time t2, the midpoint in the simulation, the value is:

# %%
gradient_t_at_x2[2]

# %% [markdown]
# ### C. All said and done, we have collected the time derivative and the 2nd spacial derivative, at the point (x2, t2).
# Do those values satisfy the diffusion equation??  Does the time derivative indeed equal the (diffusion rate) x (the 
# 2nd spacial derivative)?  Let's see:

# %%
gradient_t_at_x2[2]

# %%
diffusion_rate * second_gradient_x_at_t2[2]

# %% [markdown]
# ## D. The 2 value indeed roughly match - considering the coarseness of the large spacial grid, and the coarseness of estimating the derivatives numerically. 
# We have evidence that the diffusion equation may indeed be satisfied by the values we obtained from our diffusion algorithm, in the proximity of the point (x2, t2)

# %% [markdown]
# ### E. Finally, instead of just scrutining the match at the point x2, let's do that for all the points in space at time t2,
# WITH THE EXCEPTION of the outmost points (because the numeric estimation of the derivatives gets very crummy at the boundary).
#  
# Let's first re-visit all the data once again:

# %%
all_history

# %%
# The following is just an expanded version of what we did before for the time derivative; 
# instead of just considering 1 column,
# like we did before, we're now repeating the computations along all columns
# (the computations are applied vertically, "along axis 0" in Numpy-speak)
gradient_t = np.apply_along_axis(np.gradient, 0, all_history, delta_t)
gradient_t

# %%
# Again, we focus on time t2 (the 3rd row), to stay away from the edges
gradient_t_at_t2 = gradient_t[2]
gradient_t_at_t2

# %% [markdown]
# Note the value -8.94751865, 3rd from left : that's the single value we looked at before (at point x2)

# %% [markdown]
# #### Time to again check the match of the two sides of the diffusion equation

# %%
lhs = gradient_t_at_t2   # The LEFT-hand side of the diffusion equation, as a vector for all the spacial points at time t2
lhs

# %%
rhs = diffusion_rate * second_gradient_x_at_t2  # The RIGHT-hand side of the diffusion equation, again as a vector
rhs

# %% [markdown]
# ## The left-hand side and the right-hand side of the diffusion equation appear to generally agree, except at the boundary points, where our approximations are just too crummy

# %%
lhs - rhs

# %%
# Here we use a handy function to compare two equal-sized vectors,
# while opting to disregarding a specified number of entries at each edge.
# It returns the Euclidean distance ("L2 norm") of the shortened vectors
num.compare_vectors(lhs, rhs, trim_edges=1)

# %% [markdown]
# #### IMPORTANT: all values in this experiment are VERY coarse, because of the large effective delta_x (the tiny number of bins.)
# In part2 (notebook "validate_diffusion_2"), much-better approximations will get looked at!

# %%
