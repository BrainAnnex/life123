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
# # Tests and showcase of `PlotlyHelper.heatmap_stack_1D()`

# %%
LAST_REVISED = "Apr. 28, 2025"
LIFE123_VERSION = "1.0.0rc3"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

import numpy as np
from life123 import check_version, PlotlyHelper

# %%
check_version(LIFE123_VERSION)

# %%

# %% [markdown]
# ## A single heatmap

# %%
data_matrix = np.array([ [0.,  0., 100.,  10.,  0.,  0.,  50.,  0.,  25.,  0.] 
                       ])
data_matrix

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A"])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, 
                              labels=["A"], title="Test of 1-D heatmap stacking (with just 1 row)", 
                              data_name="Conc.", entity_name="CHEM")

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, 
                              labels=["A"], title="Test of 1-D heatmap stacking (with just 1 row)", 
                              data_name="Conc.", entity_name="CHEM", height=100)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, 
                              labels=["A"], title="Test of 1-D heatmap stacking (with just 1 row)", 
                              data_name="Conc.", entity_name="CHEM", height=500)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, 
                              labels=["A"], title="Test of 1-D heatmap stacking, color-coded", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise"])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, 
                              labels=["A"], title="Test of 1-D heatmap stacking, with color borders", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise"], color_borders=True)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix,
                              labels=["A"], title="Test of 1-D heatmap stacking, color-coded but monochromatic heatmap",
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise"], monochromatic=True)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix,
                              labels=["A"], title="Test of 1-D heatmap stacking, color-coded with borders, but monochromatic heatmap",
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise"], color_borders=True, monochromatic=True)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, 
                              labels=["A"], title="Test of 1-D heatmap stacking with barriers", 
                              data_name="Conc.", entity_name="CHEM",
                              barriers=[1, 2, 4, 9])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, 
                              labels=["A"], title="Test of 1-D heatmap stacking with barriers and color-coding", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise"], barriers=[1, 2, 4, 9])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, 
                              labels=["A"], title="1-D heatmap stacking with barriers and color-coding, incl. borders", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise"], color_borders=True, barriers=[1, 2, 4, 9])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix,
                              labels=["A"], title="1-D heatmap stacking with barriers and color-coding (gray scale)",
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise"], monochromatic=True, barriers=[1, 2, 4, 9])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix,
                              labels=["A"], title="1-D heatmap stacking with barriers and color-coding and color borders (gray scale)",
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise"], color_borders=True, monochromatic=True, barriers=[1, 2, 4, 9])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, 
                              labels=["A"], title="1-D heatmap stacking with barriers at the far edges, using a different color for barriers", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise"], 
                              barriers=[0, 1, 2, 4, 9, 10], barrier_color="pink")

# %%

# %%

# %% [markdown]
# ## Two heatmaps

# %%
data_matrix = np.array([ [0.,  0., 100.,  10.,  0.,  0.,  50.,  0.,  25.,  0.], 
                         [-5.,  3., 6.,  10.,  7.,  4.,  1.,  8.,  5., -2.]])

data_matrix

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"], 
                              title="Specifying vertical grid spacing", grid_vert_spacing=0.4)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"], 
                              title="Specifying vertical grid spacing and colorbar extension", 
                              grid_vert_spacing=0.4, colorbar_extend=2)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"], 
                              title="Specifying vertical grid spacing", grid_vert_spacing=0)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"], 
                              title="Specifying vertical grid spacing and colorbar extension", 
                              grid_vert_spacing=0, colorbar_extend=1.1)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"], 
                              title="Specifying the total height", height=550)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"], 
                              title="Specifying the total height", height=1000)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"], 
                              colors=["turquoise", "yellow"])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"], 
                              colors=["turquoise", "yellow"], color_borders=True)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"],
                              colors=["turquoise", "yellow"], monochromatic=True)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"],
                              colors=["turquoise", "yellow"], color_borders=True, monochromatic=True)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B"], 
                              colors=["turquoise", "yellow"],
                              barriers=[0, 2, 4, 9])

# %%

# %%

# %% [markdown]
# ## Three heatmaps

# %%
data_matrix = np.array([ [0.,  0., 100.,  10.,  0.,  0.,  50.,  0.,  25.,  0.], 
                         [-5.,  3., 6.,  10.,  7.,  4.,  1.,  8.,  5., -2.],
                         [8.,  23., 16.,  8.,  4.,  9.,  18.,  28.,  15., 8.]
                       ])

data_matrix

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C"])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C"], height=900)

# %%

# %%

# %% [markdown]
# ## Four heatmaps

# %%
data_matrix = np.array([ [1.234,  0., 100.,  10.,  0.,  0.,  50.,  0.,  25.,  0.], 
                         [-5.,  3., 6.,  10.,  7.,  4.,  1.,  8.,  5., -2.],
                         [8.23,  23., 16.,  8.,  4.,  9.,  18.,  28.,  15., 8.],
                         [66.,  3., 6.,  10.,  8.,  6.,  -4.,  3.,  2.,  0.]
                       ])

data_matrix

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="4 stacked heatmaps", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"],
                              barriers=[1, 2, 4, 10])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="4 stacked heatmaps", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"],
                              barriers=[1, 2, 4, 10], color_borders=True)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="4 stacked heatmaps", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"],
                              barriers=[1, 2, 4, 10], monochromatic=True)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="4 stacked heatmaps", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"],
                              barriers=[1, 2, 4, 10], monochromatic=True, color_borders=True)

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="4 stacked heatmaps", 
                              data_name="Conc.", entity_name="CHEM", height=50,
                              colors=["turquoise", "yellow", "green", "pink"],
                              barriers=[1, 2, 4, 10])

# %%
# Small number of bins

data_matrix = np.array([ [1.234,  0., 100.,  10.,  0.], 
                         [-5.987,  3., 6.,  10.,  7.],
                         [8.23,  23., 16.,  8.,  4.],
                         [66.66,  3., 6.,  10.,  8.]
                       ])

data_matrix

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="Because of the small number of bins, more text is shown on them", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"],
                              barriers=[1, 2, 4, 5])

# %%
# Larger number of bins

data_matrix = np.array([ [1.234,  0., 100.,  10.,  0.,  0.,  50.,  0.,  25.,  0., 11.56], 
                         [-5.987,  3., 6.,  10.,  7.,  4.,  1.,  8.,  5., -2., 8.34],
                         [8.23,  23., 16.,  8.,  4.,  9.,  18.,  28.,  15., 8., 15.3],
                         [66.66,  3., 6.,  10.,  8.,  6.,  -4.,  3.,  2.,  0., 33.7]
                       ])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="Because of the larger number of bins, less text is shown on them", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"],
                              barriers=[1, 2, 4, 10])

# %%
# Even larger number of bins

data_matrix = np.array([ [1.234,  0., 100.,  10.,  0.,  0.,  50.,  0.,  25.,  0., 11.56, 13.2, 0.234, 9.34,  100.,  10.,  0., 7.,  4.,  1.,  8.], 
                         [-5.987,  3., 6.,  10.,  7.,  4.,  1.,  8.,  5., -2., 8.34, 66.,  3., 6.,  10.,  8.,  6.,  -4.,  3.,  2.,  0.],
                         [8.23,  23., 16.,  8.,  4.,  9.,  18.,  28.,  15., 8., 15.3, -5.,  3., 6.,  10.,  7.,  4.,  1.,  8.,  5., -2.],
                         [66.66,  3., 6.,  10.,  8.,  6.,  -4.,  3.,  2.,  0., 33.7, 66.,  3., 6.,  10.,  8.,  6.,  -4.,  3.,  2.,  0.]
                       ])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="Because of the large number of bins, no text is shown on them", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"],
                              barriers=[1, 2, 4, 10])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="Because of the large number of bins, no text is shown on them", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"], color_borders=True, 
                              barriers=[1, 2, 4, 10])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="Because of the large number of bins, no text is shown on them", 
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"], color_borders=True, monochromatic=True, 
                              barriers=[1, 2, 4, 10])

# %%
PlotlyHelper.heatmap_stack_1D(data_matrix=data_matrix, labels=["A", "B", "C", "D"],
                              title="Because of the large number of bins, no text is shown on them",
                              data_name="Conc.", entity_name="CHEM",
                              colors=["turquoise", "yellow", "green", "pink"], monochromatic=True,
                              barriers=[1, 2, 4, 10])

# %%
