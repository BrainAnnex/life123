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
# ## A chemical that can cross membranes, if not diffusing fast enough thru the length of an organelle, will get converted by an organelle enzyme into a form that cannot cross membranes, and thus gets trapped by the organelle.  
#
# As in experiment `membranes_racing_condition_1`, **a racing condition between diffusion and reaction.**   
# Like before, the "window of opportunity to escape" closes fast!  
# This time, however, the initial concentration of `A` is outside the organelle, and the reaction is an enzyme-catalyed unimolecular `A + E -> B + E` .  The enzyme `E` cannot cross membranes.   
#
# If, metaphorically speaking, we think of `A` as a "traveler":
#
# #### SCENARIO 1 ("The Hotel California") - `A` diffuses so **slowly** that, in its journey across the organelle, it "falls under the spell of the enzyme" (gets converted into a form `B` that cannot cross membranes), and thus "can never leave"...
#
# #### SCENARIO 2 ("Ulysses' escape") - `A` diffuses so **fast** that it stays "a step ahead of trouble", and is largely out of the organelle  prior to getting entangled into it by the "magic spell" (reaction) that would have "trapped it into the island"
#
# **Recommended background:**  
#
# * experiment `1D/reaction_diffusion/membranes_racing_condition_1`

# %% [markdown]
# ### TAGS : "reactions 1D", "diffusion 1D", "membranes 1D"

# %%
LAST_REVISED = "Sep. 10, 2025"
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
# ## Initialize the Chemical Data and the Reactions.  They will be re-used in both scenarios

# %%
# Initialize the chemical data
chem_data = ChemData(names=["A", "B", "E"], 
                     diffusion_rates=[100., 0, 80.],        # The diffusion rate of `A` will later be increased in scenario 2
                     plot_colors=["red", "green", "purple"]) 

rxns = ReactionRegistry(chem_data=chem_data)

# Enzymatic reaction A + E <-> B + E
r = ReactionEnzyme(substrate="A", product="B", enzyme="E", 
                  k1_F=10., k1_R=2., k2_F=5.)

rxns.register_reaction(r)

rxns.describe_reactions()

# %%

# %%
rxns.labels_of_active_chemicals()      # TODO: FIX!!!

# %%

# %% [markdown]
# # SCENARIO 1 - `A` diffuses slowly, relatively to the Enzymatic reaction A + E <-> B + E

# %% [markdown]
# ### Initialize the 1D System, including Membranes

# %%
bio = BioSim1D(n_bins=30, chem_data=chem_data, reactions=rxns)

# %%
bio.membranes().set_membranes(membranes=[ (11, 21) ])
bio.membranes().membrane_list

# %%
# We'll use 1/2 of the diffusion rate of `A` as its membrane permeability (by passive transport)
# The conversion product `B` and the enzyne `E`, by constrast, keeps their default 0 permeability (i.e., can't cross membranes)
bio.membranes().change_permeability("A", 50.)

# %%

# %%

# %% [markdown]
# ### Initialize the initial concentrations

# %%
# Set up the initial bell-shape concentration of `A`, with a narrow peak close to one end of the system,
# centered at 20% of the width of the system, i.e. at bin 6
bio.inject_bell_curve(chem_label="A", center=0.166666, sd=0.05, max_amplitude=200., bias=0., clip=(0,9))

# %%
# The enzyme `E`, by contrast, is uniformly distributed within the organelle
for addr in range (11, 21):
    bio.set_bin_conc(bin_address=addr, conc=15., chem_label="E")

# %%
# Show as heatmap (including the membranes, shown in brown)
bio.system_heatmaps(title_prefix="Initial strong, localized transient of chemical `A` (membranes shown in brown)")

# %% [markdown]
# ### The initial transient of `A` is localized to the left of the compartment (organelle) that spans the space between bins 10 and 20

# %%
# Visualize the system state so far
bio.visualize_system(title_prefix=["Initial strong, localized transient of chemical `A` (membranes shown in brown).", 
                                   "Notice that the enzyme `E` is found inside the organelle"])

# %%
df = bio.describe_state()
df

# %%
df[df.columns[4:24]]  # Zoom in where the action is

# %%

# %% [markdown]
# ### Request history-keeping for some bins

# %%
# Request to save the concentration history at the bins with the initial concentration injection, 
# and the bins at, or near, both ends of the system
bio.enable_history(bins=[5,15,25], frequency=15, take_snapshot=True)    

# %%

# %%

# %% [markdown]
# ## Start the simulation of the reaction-diffusion

# %%
# The first round of reaction-diffusion, over a small time duration
bio.react_diffuse(total_duration=0.025, fraction_max_step=0.5, show_status=True)
bio.visualize_system(title_prefix=["The localized transient `A` starts turning into `C` by `A + B <-> C`, ",
                                   "before it can diffuse away much.  Notice the production of `C`, which can't cross the membrane"])

# %%
# SAME IN HEATMAP VIEW
bio.system_heatmaps(title_prefix=["The localized transient `A` starts turning into `C` by `A + B <-> C`, ",
                                  "before it can diffuse away much.  Notice the production of `C`, which can't cross the membrane"])

# %%

# %% [markdown]
# ### Let's continue the reaction-diffusion

# %%
bio.react_diffuse(total_duration=0.025, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %% [markdown]
# ### `A` is crossing to some extent the nearby left membrane, but not making it in time to reach the right membrane, before getting consumed   

# %%
bio.react_diffuse(total_duration=0.025, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.025, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.05, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.15, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=1.2, fraction_max_step=0.9, show_status=True)
bio.visualize_system()

# %%
bio.system_heatmaps()

# %%

# %% [markdown]
# ### Now, let's look at a few individual bins, and their concentration change with time

# %%
bio.plot_history_single_bin(title_prefix=["Diffusion, membrane passive transport, and reaction `A + B <-> C`",
                                          "Bin where the transient originates."], 
                             bin_address=6)

# %%
bio.plot_history_single_bin(title_prefix=["Diffusion, membrane passive transport, and reaction `A + B <-> C`",
                                          "At left-most edge of system."], 
                            bin_address=0)

# %% [markdown]
# #### Only a small amount of `A` reaches bin 0 early on, and later converts to `C`

# %%
# Save this plot, for later comparison with the counterpart from scenario 2
scenario_1 = bio.plot_history_single_bin(title_prefix=["Diffusion, membrane passive transport, and reaction `A + B <-> C`",
                                                       "Faraway bin."], 
                                         bin_address=29)
scenario_1 

# %% [markdown]
# ## Virtually no `A` ever reaches the faraway bins!

# %%

# %%

# %%

# %% [markdown]
# # SCENARIO 2 - `A` diffuses quickly
# Let's repeat, starting like before, but with a much higher diffusion rate of `A`

# %%
# Initial conditioned just like before

bio = BioSim1D(n_bins=30, chem_data=chem_data, reactions=rxns)

bio.membranes().set_membranes(membranes=[ (2, 18) ])

bio.inject_bell_curve(chem_label="A", center=0.2, sd=0.05, max_amplitude=200., bias=0., clip=(2,17))
bio.set_uniform_concentration(chem_label="B", conc=80.)

# %%
# 15 times faster diffusion rate of 'A', and of its membrane permeability, than before

chem_data.set_diffusion_rate(chem_label="A", diff_rate = 1500)  # **** x15 from scenario 1

bio.membranes().change_permeability("A", 750.)  # **** x15 from scenario 1 
bio.membranes().change_permeability("B", 400.)  # Same as before
# `C` keeps the default 0 permeability (i.e., can't cross membranes)

# %%
# Request history-keeping for some bins
bio.enable_history(bins=[0, 6, 29], frequency=15, take_snapshot=True)   

# %%
# Visualize the system state so far
bio.visualize_system(title_prefix=["Initial strong, localized transient of chemical `A` (membranes shown in brown)", 
                                   "Same as in the earlier scenario."])

# %%

# %% [markdown]
# ### Repeat the simulation

# %%
# The first round of reaction-diffusion, over a small time duration
bio.react_diffuse(total_duration=0.025, fraction_max_step=0.5, show_status=True)
bio.visualize_system(title_prefix=["This time, the localized transient `A` diffuses more substantially across the membrane,",
                                   "before being converted into `C`, which can't cross the membrane"])

# %%
bio.react_diffuse(total_duration=0.025, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %% [markdown]
# ### In this scenario, `A` is managing to cross both membranes, to some extent, before getting consumed   

# %%
bio.react_diffuse(total_duration=0.025, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.025, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.05, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %%
bio.react_diffuse(total_duration=0.15, fraction_max_step=0.5, show_status=True)
bio.visualize_system()

# %% [markdown]
# #### The increasingly prominent `C` (green) on the right is from `A` that managed to cross the membrane, and then converted to `C`

# %%
bio.react_diffuse(total_duration=1.2, fraction_max_step=0.9, show_status=True)
bio.visualize_system()

# %%

# %%

# %% [markdown]
# ### Now, let's look at the faraway bin 29, and it concentration change with time

# %%
scenario_2 = bio.plot_history_single_bin(title_prefix=["Diffusion, membrane passive transport, and reaction `A + B <-> C`",
                                                       "Faraway bin."], 
                                         bin_address=29)
scenario_2

# %% [markdown]
# ## A lot more `A` (eventually converted into `C`) reaches the faraway bins, compared to the earlier scenario!

# %%
PlotlyHelper.combine_plots(fig_list = [scenario_1, scenario_2], 
                           title="Diffusion, membrane passive transport, and reaction `A + B <-> C`<br>Concentrations in faraway bin 29<br>Comparison of scenario 1 (dotted) and 2 (solid)", 
                           layout_index=0, modify = {0: "dot"})

# %% [markdown]
# Magnify plot to separate the red and green dots (very close to each other)!  
#
# Side note: `B` eventually equilibrates to the same value under either scenario - since `B` diffuses freely across membranes, and is the excess reactant: it diffuses throughout to bind to `A`, wherever `A` might be...

# %%
