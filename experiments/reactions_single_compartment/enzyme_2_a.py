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
# ## **Enzyme Kinetics** : 
#
# #### Our model: `E + S <-> ES` (with kinetic parameters _k1_forward_ and _k1_reverse_), and  `ES -> E + P`  (_k2_forward_)  
#
# #### In experiment `enzyme_1_a`, we were given `k1_forward`, `k1_reverse` and `k2_forward`...  But what to do if we're **just given `kM` and `kcat`** ?  
#
# Background: please see experiment `enzyme_1_a`

# %% [markdown]
# #### THE REACTION:  
# the enzyme `Adenosinedeaminase`,  
# with the substrate `2,6-Diamino-9-β-D-deoxyribofuranosyl-9-H-purine`.
#
# Source of kinetic parameters:  *page 16 of "Analysis of Enzyme Reaction Kinetics, Vol. 1", by F. Xavier Malcata, Wiley, 2023*

# %% [markdown]
# ### TAGS :  "uniform compartment", "chemistry", "numerical", "enzymes"

# %%
LAST_REVISED = "Sep. 2, 2025"
LIFE123_VERSION = "1.0.0rc6"         # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

import numpy as np
import plotly.express as px

from life123 import check_version, ReactionEnzyme

# %%
check_version(LIFE123_VERSION)    # To check compatibility

# %%

# %%

# %% [markdown]
# ## Assume we're only given values for `kM` and `kcat`
# ### What values of `k1_forward`, `k1_reverse` and `k2_forward` are compatible with them?

# %%
# We'll use the following values, taken from experiment `enzyme_1_a`
kM = 8.27777777777777
kcat = 49

# %%
# k2_forward equals kcat ; not much to say here!
k2_forward = kcat
k2_forward

# %% [markdown]
# By definition: `kM = (k2_forward + k1_reverse) / k1_forward`   
#
# We are given `kM` and `k2_forward` (same as `kcat`), as those are typical quantities measured experimentally, i.e.:  
#
# `kM = (kcat + k1_reverse) / k1_forward`  
#
# But how to solve for `k1_forward` and `k1_reverse`??  We have just 1 equation and 2 variables!  **The system of equations is "underdetermined"** : what can we do?   
#
# We'll explore fixing a guess for `k1_forward`, and then computing the corresponding `k1_reverse` - or vice versa.  
#
# The Life123 class `ReactionEnzyme` conveniently provides the necessary transformations.

# %%
enz = ReactionEnzyme()

# %%
# Example, using the k1_forward=18. from experiment `enzyme_1_a`, to determine k1_reverse

k1_reverse = enz.compute_k1_reverse(kM=kM, kcat=kcat, k1_forward = 18.)
k1_reverse

# %%
# Conversely, using the k1_reverse=100. from experiment `enzyme_1`, to determine k1_forward

k1_forward = enz.compute_k1_forward(kM=kM, kcat=kcat, k1_reverse = 100.)
k1_forward

# %% [markdown]
# #### Naturally, we're getting the same values we had in experiment `enzyme_1_a`,  
# namely `k1_forward = 18.` and `k1_reverse = 100.`  
# #### But what if neither `k1_forward` nor `k1_reverse` are known?

# %%

# %% [markdown]
# ### PART 1. Let's try a variety of values for `k1_reverse`, and determine the corresponding values for `k1_forward`

# %% [markdown]
# #### `k1_reverse` must be non-negative because it's a reaction rate constant, but in other respects there's no conceptual restriction on its value, as plugged into our equation:
# `kM = (kcat + k1_reverse) / k1_forward`

# %%
k1_reverse_choices = np.linspace(0., 300., 3)    # Even grid of values
k1_reverse_choices

# %%
k1_forward_choices = enz.compute_k1_forward(kM=kM, kcat=kcat, k1_reverse = k1_reverse_choices)
k1_forward_choices

# %%
fig = px.line(x=k1_reverse_choices, y=k1_forward_choices, title="k1_forward for given k1_reverse values")

fig.update_layout(xaxis_title='k1_reverse',
                  yaxis_title='k1_forward')

fig.add_vline(x=100., line_width=1, line_dash="dot", line_color="gray")
fig.add_hline(y=18.,  line_width=1, line_dash="dot", line_color="gray")
fig.add_scatter(x=[100.], y=[18.], 
                mode="markers", marker={"color": "red"}, name="actual values")

# %%

# %%

# %% [markdown]
# ### PART 2. Conversely, let's try a variety of values for `k1_forward`, and determine the corresponding values for `k1_reverse`

# %% [markdown]
# #### There's an extra consideration in pursuing the converse approach, namely we cannot simply give `k1_forward` any non-negative value as we please, because exessively small values would cause `k1_reverse` to become negative, as may be observed from our equation:
# `kM = (kcat + k1_reverse) / k1_forward`   
# which can be re-written as:  
# `kM * k1_forward = kcat + k1_reverse`, i.e.  
# `k1_reverse = kM * k1_forward - kcat`  
#
# To insure that `k1_reverse >= 0`  we must enforce `kM * k1_forward >= kcat`, i.e.  `k1_forward >= kcat / kM`  
#
# Hence, the required minimun value for `k1_forward` is `kcat / kM`.   This may be conveniently computed with the following function call:

# %%
enz.min_k1_forward(kM=kM, kcat=kcat)

# %% [markdown]
# We'll use the above min value as the start of the range of values that we'll explore for `k1_forward`

# %%
k1_forward_choices = np.linspace(5.92, 50., 3)    # Even grid of values; notice the start value just above the min required value
k1_forward_choices

# %%
k1_reverse_choices = enz.compute_k1_reverse(kM=kM, kcat=kcat, k1_forward = k1_forward_choices)
k1_reverse_choices

# %%
fig = px.line(x=k1_forward_choices, y=k1_reverse_choices, 
              title="k1_reverse for given k1_forward values<br>(dashed lines show actual values)")

fig.update_layout(xaxis_title='k1_forward',
                  yaxis_title='k1_reverse')

fig.add_vline(x=18., line_width=1, line_dash="dot", line_color="gray")
fig.add_hline(y=100., line_width=1, line_dash="dot", line_color="gray")
fig.add_scatter(x=[18.], y=[100.], 
                mode="markers", marker={"color": "red"}, name="actual values")

# %% [markdown]
# Note that smallest value of `k1_forward` is NOT zero, but rather about 5.92, as discussed earlier

# %% [markdown]
# ### In the continuation experiment, `enzyme_2_b`, we'll explore how variations of `k1_forward` and `k1_reverse` (guesses consistent with `kM` and `kcat`) affect the kinetics of our enzymatic reaction...

# %%
