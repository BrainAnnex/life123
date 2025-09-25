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
# ## A chemical that can cross membranes, if not diffusing fast enough thru the length of a compartment, will get converted by a compartment enzyme into a form that cannot cross membranes, and thus gets trapped inside.  
#
# As in experiment `membranes_racing_condition_1`, **a racing condition between diffusion and reaction.**   
# Like before, the "window of opportunity to escape" closes fast!  
# This time, however, the initial concentration of `A` is outside the compartment, and the reaction is an enzyme-catalyed unimolecular `A + E -> B + E` .  The enzyme `E` cannot cross membranes.   
#
# If, metaphorically speaking, we think of `A` as a "traveler", then:
#
# #### SCENARIO 1 ("The Hotel California") - `A` diffuses so **slowly** that, in its journey across the compartment, it "falls under the spell of the enzyme" (gets converted into a form `B` that cannot cross membranes), and thus **"can never leave"**...
#
# #### SCENARIO 2 ("Ulysses' escape") - `A` diffuses so **fast** that it **"stays a step ahead of trouble"**, and is largely out of the compartment  prior to getting entangled into it by the "magic spell" (reaction) that would have "trapped it into the island"
#
# **Recommended background:**  
#
# * experiment `1D/reaction_diffusion/membranes_racing_condition_1`

# %% [markdown]
# ### TAGS : "reactions 1D", "diffusion 1D", "membranes 1D"

# %%
LAST_REVISED = "Sep. 25, 2025"
LIFE123_VERSION = "1.0.0rc7"       # Library version this experiment is based on

# %%
#import set_path              # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim1D, ChemData, ReactionRegistry, ReactionEnzyme, PlotlyHelper, check_version

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## Initialize the Chemical Data and the Reactions.  They will be RE-USED (except for a change in the diffusion rate of `A`) in both scenarios

# %%
# Initialize the chemical data.  
# The enzyme `E`, typically a large molecule, is given a relatively sall diffusion rate (not too important,
# because of our initial condition later, of `E` already uniformly diffused)
chem_data = ChemData(names=["A", "B", "E"], 
                     diffusion_rates=[200., 250., 80.],        # The diffusion rate of `A` will later be increased in scenario 2
                     plot_colors=["red", "green", "violet"]) 

rxns = ReactionRegistry(chem_data=chem_data)

# Enzymatic reaction A + E -> B + E
r = ReactionEnzyme(substrate="A", product="B", enzyme="E", 
                   k1_F=10., k1_R=2., k2_F=50.)

rxns.register_reaction(r)

rxns.describe_reactions()

# %%
chem_data.all_chemicals()   # We'll let a color for `EA` get automatically assigned later

# %%

# %%

# %% [markdown]
# # SCENARIO 1 - `A` diffuses slowly, relatively to the Enzymatic reaction `A + E -> B + E`

# %% [markdown]
# ### Initialize the 1D System, including Membranes

# %%
bio = BioSim1D(n_bins=30, chem_data=chem_data, reactions=rxns)

# %%
bio.membranes().set_membranes(membranes=[ (11, 21) ])
bio.membranes().membrane_list

# %%
# We'll use 1/2 of the diffusion rate of `A` as its membrane permeability (by passive transport)
# The conversion product `B` and the enzyne `E`, by constrast, keep their default 0 permeability (i.e., can't cross membranes)
bio.membranes().change_permeability("A", 50.)

# %%
bio.membranes().show_permeability()   # Values not shown mean 0

# %%

# %%

# %% [markdown]
# ### Initialize the initial concentrations

# %%
# Set up the initial bell-shape concentration of `A`, with a narrow peak close to one end of the system,
# centered at 20% of the width of the system, i.e. at bin 6
bio.inject_bell_curve(chem_label="A", center=0.166666, sd=0.05, max_amplitude=200., bias=0., clip=(0,9))

# %%
# The enzyme `E`, by contrast, is uniformly distributed within the membranes of the 1st (and only) compartment
bio.set_compartment_uniform_concentration(compartment_id=0, chem_label="E", conc=15.)

# %%
# Show as heatmap (including the membranes, shown in brown)
bio.system_heatmaps(title_prefix="Initial strong, localized transient of chemical `A` (membranes shown in brown)")

# %%

# %% [markdown]
# ### The initial transient of `A` is localized to the left of the compartment that spans the space between bins 10 and 20  
# The enzyme `E` is uniformly localized in that compartment.  
# `B` and the reaction intermerdiate `EA` are not present anywhere.

# %%
# Visualize the system state so far
bio.visualize_system(title_prefix=["Initial strong, localized transient of chemical `A` (membranes shown in brown).", 
                                   "Notice that the enzyme `E` is found inside the compartment"])

# %%
df = bio.describe_state()
df

# %%
df[df.columns[4:24]]  # Zoom in where the action is

# %%

# %% [markdown]
# ### Request history-keeping for some bins

# %%
# Request to save the concentration history 
# at the bin with the initial concentration injection, on the left of the system, 
# plus 2 bins inside the compartment, and one on the far side on the right
bio.enable_history(bins=[5,11,15,25], frequency=15, take_snapshot=True)     

# %%

# %%

# %% [markdown]
# ## Start the simulation of the reaction-diffusion

# %%
# The first round of reaction-diffusion, over a small time duration
bio.react_diffuse(total_duration=0.05, fraction_max_step=0.5, show_status=True)
bio.visualize_system(title_prefix=["The transient `A` starts diffusing into the compartment, and turning into `B` by `A + E -> B + E` ",
                                   "Notice the production of `B`, which can't cross the membranes"])

# %% [markdown]
# At this point, some of the enzyme `E` is in the form of the intermediate `EA`, especially in the left of the compartment (near bin 11), where `A` is diffusing in from the left, and initiating the reaction

# %%
# SAME IN HEATMAP VIEW
bio.system_heatmaps(title_prefix=["The transient `A` starts diffusing into the compartment, and turning into `B` by `A + E -> B + E`"])

# %%

# %% [markdown]
# ### Let's continue the reaction-diffusion

# %%
bio.react_diffuse(total_duration=0.05, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.1, fraction_max_step=0.5, show_status=True)
bio.visualize_system(title_prefix="Here we turned off the curve smoothing in the plot!",smoothed=False)

# %%
bio.react_diffuse(total_duration=0.2, fraction_max_step=0.5, show_status=True)
bio.visualize_system(smoothed=False)

# %%
bio.react_diffuse(total_duration=0.2, fraction_max_step=0.5, show_status=True)
bio.visualize_system(smoothed=False)

# %%
bio.react_diffuse(total_duration=1, fraction_max_step=0.5, show_status=True)
bio.visualize_system(smoothed=False)

# %%
bio.react_diffuse(total_duration=1.4, fraction_max_step=0.5, show_status=True)
bio.visualize_system(smoothed=False)

# %% [markdown]
# ### `A` is crossing to some extent the nearby left membrane, but not making it in time to reach the right membrane, before getting consumed   

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ### Now, let's look at a few individual bins, and their concentration change with time

# %%
bio.plot_history_single_bin(title_prefix=["Diffusion, membrane passive transport, and reaction `A + E -> B + E`",
                                          "Time evolution at bin where the transient originates."], 
                             bin_address=5)

# %%
bio.plot_history_single_bin(title_prefix=["Time evolution inside the compartment, at the very left of it"], 
                            bin_address=11)

# %% [markdown]
# ### Notice the arrival of `A`, albeit short-lasting.  Also notice the transient buildup of the intermediary `EA` (and its corresponding transient dip in `E`)

# %%
bio.plot_history_single_bin(title_prefix=["Time evolution at a bin in the middle of the compartment"], 
                            bin_address=15)

# %% [markdown]
# ### `A` enters the compartment from the left (bin 11), thru passive transport across the membrane), but barely has a chance to diffuse even to the middle of the compartment (bin 15) because of its quick conversion to `B`

# %%
bio.plot_history_single_bin(title_prefix=["Time evolution several bins to the right of the compartment"], 
                            bin_address=25)

# %% [markdown]
# ## Only PUNY amounts of `A` ever reach the faraway bins to the right of the compartment!

# %%

# %%

# %%

# %% [markdown]
# # SCENARIO 2 - `A` diffuses quickly
# Let's repeat, starting like before, but with a much higher diffusion rate of `A`

# %%
# Initial conditioned just like before, with same reaction, membrane layoyt and concentrations

bio = BioSim1D(n_bins=30, chem_data=chem_data, reactions=rxns)

bio.membranes().set_membranes(membranes=[ (11, 21) ])

bio.inject_bell_curve(chem_label="A", center=0.166666, sd=0.05, max_amplitude=200., bias=0., clip=(0,9))
bio.set_compartment_uniform_concentration(compartment_id=0, chem_label="E", conc=15.)

# %%
# Re-using all the chemical data from before, EXCEPT for substantially increasing the diffusion rate of `A`
# and its membrane permeability, compared to the earlier scenario

# %%
# 15 times faster diffusion rate of 'A', and of its membrane permeability, than before

chem_data.set_diffusion_rate(chem_label="A", diff_rate = 3000)  # **** x15 from scenario 1

bio.membranes().change_permeability("A", 750.)                 # **** x15 from scenario 1 

# %%
chem_data.all_chemicals()   # Notice that the plot_color for `EA` was automatically-assigned earlier

# %%
bio.membranes().membrane_list

# %%
bio.membranes().show_permeability()   # Values not shown mean 0

# %%

# %%
# Request to save the concentration history 
# at the bin with the initial concentration injection, on the left of the system, 
# plus 2 bins inside the compartment, and one on the far side on the right
bio.enable_history(bins=[5,11,15,25], frequency=15, take_snapshot=True)    

# %%

# %%

# %%
# Visualize the system state so far
bio.visualize_system(title_prefix=["Initial strong, localized transient of chemical `A` (membranes shown in brown).", 
                                   "Notice that the enzyme `E` is found inside the compartment"])

# %%

# %% [markdown]
# ### Repeat the simulation

# %%
# The first round of reaction-diffusion, over a small time duration
bio.react_diffuse(total_duration=0.05, fraction_max_step=0.5, show_status=True)
bio.visualize_system(title_prefix=["The transient `A` starts diffusing into the compartment, and turning into `B` by `A + E -> B + E` ",
                                   "Notice the production of `B`, which can't cross the membranes"])

# %%
# SAME IN HEATMAP VIEW
bio.system_heatmaps(title_prefix=["The transient `A` starts diffusing into the compartment, and turning into `B` by `A + E -> B + E`"])

# %%

# %% [markdown]
# ### Let's continue the reaction-diffusion

# %%
bio.react_diffuse(total_duration=0.05, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.1, fraction_max_step=0.5, show_status=True)
bio.visualize_system(title_prefix="Here we turned off the curve smoothing in the plot!",smoothed=False)

# %%
bio.react_diffuse(total_duration=2.8, fraction_max_step=0.5, show_status=True)
bio.visualize_system(smoothed=False)

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ### Now, let's look at a few individual bins, and their concentration change with time

# %%
bio.plot_history_single_bin(title_prefix=["Diffusion, membrane passive transport, and reaction `A + E -> B + E`",
                                          "Time evolution at bin where the transient originates."], 
                             bin_address=5)

# %%
bio.plot_history_single_bin(title_prefix=["Time evolution inside the compartment, at the very left of it"], 
                            bin_address=11)

# %%
bio.plot_history_single_bin(title_prefix=["Time evolution at a bin in the middle of the compartment"], 
                            bin_address=15)

# %%
bio.plot_history_single_bin(title_prefix=["Time evolution several bins to the right of the compartment"], 
                            bin_address=25)

# %%
