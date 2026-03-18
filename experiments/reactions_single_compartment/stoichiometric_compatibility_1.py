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
# ## Stoichiometric Compatibility in reaction network `A -> 2 B` and `2 B -> A` 
#
# (i.e. reversible synthesis/decomposition reaction)
#
# From different starting conditions, we explore the stoichiometric compatibility class for this reaction, and observe that the trajectories in phase space follow 1-dimensional paths, as expected.  
#
# **Background**  
# We'll be following "Foundations of Chemical Reaction Network Theory" (2019), by Martin Feinberg, section 3.4
#
# See also: experiment _"react_4"_ , for an intro to this reaction

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

from life123 import check_version, ReactionRegistry, UniformCompartment, PlotlyHelper, Colors

# %%
check_version(LIFE123_VERSION)

# %%

# %%

# %% [markdown]
# # PART 1 - simulation of the reaction

# %%
# Set up the reactions and their chemicals (common for all the simulations below)
rxns = ReactionRegistry()

# Elementary reaction A <-> 2 B
rxns.add_reaction(reactants="A", products=[(2, "B")], kF=3., kR=2.)

rxns.describe_reactions()

# %%
# Instantiate the simulator with the predefined chemical species and set of reactions
uc = UniformCompartment(reactions=rxns, preset="mid")

# Initial concentrations of all the chemicals
uc.set_conc({"A": 40., "B": 0.})
uc.describe_state()

# %%
uc.get_history()

# %%

# %%

# %% [markdown]
# ### Let's look ahead at the final equilibrium point

# %%
uc.find_equilibrium_conc(rxn_index=0)    # This is an EXACT equilibrium solution, 
                                         # for 1 reaction (the only reaction)

# %%

# %%

# %% [markdown]
# ### Run the reaction

# %%
uc.react_to_equilibrium(initial_step=0.01)

# %%
uc.is_in_equilibrium()

# %%
df = uc.get_history()
df

# %%

# %% [markdown]
# ### Visualize

# %%
color_palette = Colors.assign_default_colors(n=8)

# %%
PlotlyHelper.plot_pandas(df=df, x_var="B", fields="A", colors=color_palette[0], 
                         show_points=True, annotation_field="SYSTEM TIME",
                         annotate={"frequency": 3})

# %% [markdown]
# #### Because of the stoichiometry of A <-> 2B , for every molecule of `A` consumed, 2 of `B` are produced.  Hence, the linearity of `A` vs. `B` over time (in our scenario, the reaction starts in the upper left point, and attains equilibrium in the bottom right)

# %% [markdown]
# ## More formally, the **stoichiometric compatibility class** of our reaction network (`A->2B` and `2B->A`) is uni-dimensional.
# And that's because our reaction vectors are {A−2B , 2B−A}, which are linearly dependent.

# %% [markdown]
# For more detailed explanations, see *"Foundations of Chemical Reaction Network Theory" (2019),* by Martin Feinberg, section 3.4

# %%

# %%
# Let's save up a simpler version of the above plot, to later combine with others
plot = PlotlyHelper.plot_pandas(df=df, x_var="B", fields="A", title=f"A0 = 40",
                                show_points=True, annotation_field="SYSTEM TIME",)

all_plots = [plot]

# %%

# %%

# %% [markdown]
# ## Repeat, with a different start condition for A(0)
#

# %%
for i, A0 in enumerate([10, 20, 30, 50, 60, 70, 80]):    # Skipping 40, which we already did  
    # Instantiate the simulator with the predefined chemical species and set of reactions
    uc = UniformCompartment(reactions=rxns, preset="mid")

    # Initial concentrations of all the chemicals
    uc.set_conc({"A": A0, "B": 0.})
 
    # Run the reaction
    uc.react_to_equilibrium(initial_step=0.01)    # Variable steps is the default
    
    df = uc.get_history()
    
    plot = PlotlyHelper.plot_pandas(df=df, x_var="B", fields="A", title=f"A0 = {A0}",
                                    show_points=True, annotation_field="SYSTEM TIME",
                                    colors=color_palette[i+1])  # Avoiding color_palette[0], already used
    
    all_plots.append(plot)

# %%

# %% [markdown]
# ### Finally, combine all plots

# %%
PlotlyHelper.combine_plots(fig_list=all_plots,
                           title="Some stoichiometric compatibility classes<br>for reaction network A->2B and 2B->A",
                           legend_title="Initial starting point")

# %% [markdown]
# ### Always linearity, but with different starting and ending points

# %% [markdown]
# ### Compare with the phase portrait in fig. 3.1 on page 32 of "Foundations of Chemical Reaction Network Theory" (2019), by Martin Feinberg 
# (Here we're not approaching the equilibria points from the other side, i.e. with initial concentrations of just `B`)

# %% [markdown]
# #### In a simple reaction like our `A <-> 2B`, the linearity of the "A vs. B" plot is obvious ("for every `A` consumed, 2 of `B` are produced"), but for a larger network of reactions, it generally won't be 1-dimensional trajectories

# %%
