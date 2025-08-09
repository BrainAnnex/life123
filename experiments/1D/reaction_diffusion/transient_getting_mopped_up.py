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
# ## TBA
#
# Diffusion of a narrow, bell-shaped, initial concentration of a single chemical.   
# With increasing distance from the location of the transient signal, hardly any change with time is detected, as the system goes to equilibrium

# %% [markdown]
# ### TAGS : "reactions 1D", "diffusion 1D"

# %%
LAST_REVISED = "Aug. 7, 2025"
LIFE123_VERSION = "1.0.0rc5"       # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim1D, ChemData, Reactions, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %%
# Initialize the system
chem_data = ChemData(names=["A", "B", "C"], diffusion_rates=[10., 2., 30.],  
                     plot_colors=["red", "turquoise", "green"]) 

rxns = Reactions(chem_data=chem_data)

# Reaction A + B <-> C , with 1st-order kinetics for each species; note that it's mostly in the forward direction
rxns.add_reaction(reactants=["A", "B"], products="C", forward_rate=20., reverse_rate=1.)
rxns.describe_reactions()

# %%
bio = BioSim1D(n_bins=300, chem_data=chem_data, reactions=rxns)

# %%
# Set up the initial bell-shape concentration of `A`, with the very narrow peak close to one end of the system,
# centered at 1/10 of the width of the system, i.e. at bin 30
bio.inject_bell_curve(chem_label="A", center=0.1, sd=0.005, amplitude=0.1, bias=10.)

# %%
# Chemical `B`, by contrast, is uniformly distributed
bio.set_uniform_concentration(chem_label="B", conc=20.)

# %%

# %%
df = bio.describe_state()
df

# %%
df[df.columns[22:44:2]]  # Zoom in where the action is

# %%

# %%
# Visualize the system state so far
bio.visualize_system(title_prefix="Initial strong, localized transient of chemical `A`")

# %%
# Show as heatmap
bio.system_heatmaps(title_prefix="Initial strong, localized transient")

# %%

# %% [markdown]
# ## Request history-keeping for some bins

# %%
# Request to save the concentration history at the bins with the initial concentration injection, 
# and the bins at, or near, both ends of the system
bio.enable_history(bins=[0, 10, 20, 30, 100, 200, 299], frequency=25, take_snapshot=True)    

# %%

# %%
bio.get_reactions().describe_reactions()

# %% [markdown]
# ## Start the reaction-diffusion

# %%
bio.react_diffuse(time_step=0.0000001, n_steps=1)

# %%

# %%

# %%

# %%
bio.react_diffuse(total_duration=5, time_step=0.000000001)

# %%
bio.visualize_system()

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ### Scrutinize the changes with time, at bins increasingly further away from the transient peak of bin 50

# %%
bio.plot_history_single_bin(title_prefix="Bin with the concentration injection.", bin_address=50)

# %%
bio.plot_history_single_bin(title_prefix="Bin very close to the location of the concentration injection.", bin_address=48)

# %%
bio.plot_history_single_bin(title_prefix="Bin relatively far from the location of the concentration injection.", bin_address=25)

# %%

# %% [markdown]
# ### Continue the diffusion, to equilibrium

# %%
# Do several rounds of diffusion, over relatively small time steps
for _ in range(5):
    bio.diffuse(total_duration=25, time_step=0.02)
    bio.visualize_system(show=True)

# %% [markdown]
# #### Notice how the wave of diffusion hits the left edge of the system

# %%

# %%
# Do more rounds of diffusion, over larger time steps
for _ in range(5):
    bio.diffuse(total_duration=400, time_step=0.025)
    bio.visualize_system(show=True)

# %%
bio.system_heatmaps()

# %% [markdown]
# #### We're now close to equilibrium

# %%

# %% [markdown]
# ### Scrutinize the changes with time, at bins increasingly further away from the transient peak of bin 50

# %%
bio.plot_history_single_bin(title_prefix="Bin with the concentration injection.", bin_address=50)

# %%
bio.plot_history_single_bin(title_prefix="Bin very close to the location of the concentration injection.", bin_address=48)

# %%
bio.plot_history_single_bin(title_prefix="Bin relatively far from the location of the concentration injection.", bin_address=25)

# %%
bio.plot_history_single_bin(title_prefix="Bin also relatively far from the location of the concentration injection (but on opposite side).", 
                            bin_address=0)

# %%
bio.plot_history_single_bin(title_prefix="Bin quite far from the location of the concentration injection.", bin_address=100)

# %%
bio.plot_history_single_bin(title_prefix="Bin very far from the location of the concentration injection.", bin_address=200)

# %%
bio.plot_history_single_bin(title_prefix="Bin hugely far from the location of the concentration injection.", bin_address=400)

# %% [markdown]
# # Notice how faraway locations barely register that the distant transient ever happened!

# %%
