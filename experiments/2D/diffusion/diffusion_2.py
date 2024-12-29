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
# # IN-PROGRESS

# %% [markdown]
# ## An initial concentration pulse of 2 chemicals, starting near diametrically-opposite ends of the system, diffusing towards equilibrium with different rates.
#
# No reaction takes place; the system is left undisturbed, and followed to equilibrium.

# %% [markdown]
# ### TAGS :  "diffusion 2D"

# %%
LAST_REVISED = "Dec. 27, 2024"
LIFE123_VERSION = "1.0.0rc2"        # Library version this experiment is based on

# %%
#import set_path                    # Using MyBinder?  Uncomment this before running the next cell!

# %%
#import sys
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path   

from life123 import BioSim2D, ChemData, check_version

# %%
check_version(LIFE123_VERSION)

# %%
# Prepare the initial system, with a single non-zero bin, near the left edge of the system, positioned halfway vertically
chem_data = ChemData(names=["A", "B"], diffusion_rates=[0.02, 0.01])
bio = BioSim2D(n_bins=(6, 10), chem_data=chem_data)

bio.set_bin_conc(bin_x = 1, bin_y = 1, chem_label="A", conc=10.)
bio.set_bin_conc(bin_x = 4, bin_y = 8, chem_label="B", conc=20.)

bio.describe_state()

# %%
bio.system_snapshot(chem_label="A")

# %%
bio.system_snapshot(chem_label="B")

# %%
bio.heatmap_single_chem(chem_label="A", title_prefix="Diffusion of 2 chemicals", height=400)

# %%
bio.heatmap_single_chem(chem_label="B", title_prefix="Diffusion of 2 chemicals", height=400)

# %%
bio.heatmap(chem_label="A")

# %%
bio.heatmap(chem_label="A", color_name="blue")

# %%
bio.heatmap(chem_label="A", color_name="red")

# %%
bio.heatmap(chem_label="A", color_name="green")

# %%
bio.heatmap(chem_label="A", color_name="black")

# %%
heatmap_NEW(chem_label="A", color_name="green")

# %%
heatmap_NEW(chem_labels=["A", "B"], colors=["yellow", "green"])

# %%
heatmap_NEWER(chem_labels=["A", "B"], colors=["yellow", "green"])


# %%
def heatmap_NEWER(chem_labels :[str], title_prefix = "", height=500, colors=None) -> pgo.Figure:

    title = f"System state at time t={bio.system_time:.5g}"
    if title_prefix:
        title = f"{title_prefix}.  {title}"

    # Create subplots for all the
    nrows = 1  # Number of rows in subplot grid
    ncols = 2  # Number of columns in subplot grid

    fig = sp.make_subplots(rows=nrows, cols=ncols, subplot_titles=[f'{c}' for c in chem_labels],
                           horizontal_spacing=0.2, vertical_spacing=0.2)
    
    row = 1
    col = 1
    for i, chem in enumerate(chem_labels):
        color_name = colors[i]
        if color_name is None:
            color_scale = "gray_r"
        else:
            lighter_color = PlotlyHelper.lighten_color(color_name, factor=.96)
            color_scale = [
                [0.0, lighter_color],   # Light tint
                [1.0, color_name],      # Full color
            ]

        # Create the Heatmap object
        hm = pgo.Heatmap(z=bio.system_snapshot(chem_label=chem),
                         xgap=2, ygap=2,
                         hovertemplate='Conc.: %{z}<br>x bin: %{x}<br>y bin: %{y}<extra>A</extra>',
                         texttemplate = '%{z:.2f}',
                         colorscale=color_scale,
                         colorbar=dict(
                                    title="Conc",
                                    len=1.1,  # Length of the colorbar
                                    y=0.55 - (row - 1) * 0.45,  # Adjust vertical position
                                    x=0.4 + 0.6 *(col-1),         # Adjust horizontal position
                                )
                         )
        
        fig.add_trace(hm, row=row, col=col)   # row and col values start at 1
        
        # Add x-axis and y-axis titles for this subplot
        fig.update_xaxes(
                            title_text=f"x bin",
                            row=row,
                            col=col
                        )
        
        fig.update_yaxes(
                            title_text=f"y bin",
                            row=row,
                            col=col
                        )
        
        # Invert the y-axis for this subplot
        fig.update_yaxes(
                            autorange="reversed",
                            row=row,
                            col=col
                        )
            
        col += 1


    # Update layout
    fig.update_layout(
        title=title,
        height=height,
        showlegend=False
    )

    return fig


# %%
def heatmap_NEW(chem_labels :[str], title_prefix = "", height=550, width=None, colors=None) -> pgo.Figure:

    title = f"System state at time t={bio.system_time:.5g}"
    if title_prefix:
        title = f"{title_prefix}.  {title}"

    # Create subplots for all the
    rows = 1  # Number of rows in subplot grid
    cols = 2  # Number of columns in subplot grid

    fig = sp.make_subplots(rows=rows, cols=cols, subplot_titles=[f'{c}' for c in chem_labels])
    
    row = 1
    col = 1
    for i, chem in enumerate(chem_labels):
        color_name = colors[i]
        if color_name is None:
            color_scale = "gray_r"
        else:
            lighter_color = PlotlyHelper.lighten_color(color_name, factor=.96)
            color_scale = [
                [0.0, lighter_color],   # Light tint
                [1.0, color_name],      # Full color
            ]

        # Create the Heatmap object
        hm = pgo.Heatmap(z=bio.system_snapshot(chem_label=chem),
                         colorscale=color_scale,
                         xgap=2, ygap=2,
                         hovertemplate='Conc.: %{z}<br>x bin: %{x}<br>y bin: %{y}<extra>A</extra>',
                         texttemplate = '%{z:.2f}',
                         colorbar_title="Concentration"
                         )

        fig.add_trace(hm, row=row, col=col)   # row and col values start at 1
        col += 1
        
        
        
    # Create the Figure object
    #fig = pgo.Figure(data=hm)

    # Update layout
    fig.update_layout(
        title=title,
        height=height,
        width=width,

        yaxis_title='y bin',
        yaxis_autorange="reversed"    
    )

    return fig


# %%
import plotly.subplots as sp
import plotly.graph_objects as pgo
import numpy as np
from life123.visualization.plotly_helper import PlotlyHelper

# %%
import plotly.subplots as sp
import plotly.graph_objects as go
import numpy as np

# Example 3D numpy array: 10 chemicals, each with a 2D heatmap of size 20x20
data = np.random.rand(10, 20, 20)

# Create subplots for all slices
rows = 2  # Number of rows in subplot grid
cols = 5  # Number of columns in subplot grid

fig = sp.make_subplots(rows=rows, cols=cols, subplot_titles=[f'Chemical {i+1}' for i in range(data.shape[0])])

for i in range(data.shape[0]):
    row = i // cols + 1
    col = i % cols + 1
    fig.add_trace(go.Heatmap(z=data[i], colorscale='Viridis', showscale=False), row=row, col=col)

# Update layout
fig.update_layout(height=600, width=1000, title='Chemical Concentration Heatmaps')

# %%

# %% [markdown]
# # Initial Diffusion Step

# %%
delta_time = 10.

status = bio.diffuse(total_duration=delta_time, time_step=0.1)
print("\n", status)

bio.describe_state()

# %%
bio.heatmap_single_chem(chem_label="A", title_prefix="Diffusion")

# %% [markdown]
# ## This is still an early stage in the diffusion process; let's advance it more... (Visualization from results shown at selected times)

# %% tags=[]
for i in range(200):
    status = bio.diffuse(total_duration=delta_time, time_step=0.1)

    if i<2 or i==6 or i>=199:
        bio.describe_state()
        fig = bio.heatmap_single_chem(chem_label="A", title_prefix="Diffusion", height=400)
        fig.show()


# %% [markdown]
# # All bins now have essentially uniform concentration. The diffusion has reached equilibrium
#
# Notice, throughout the simulation, the continued symmetry across the mid-row (ybin 2).
#
# **Mass conservations**: the initial "10. units of concentration" are now uniformly spread across the 40 (5x8) bins, leading to a near-constant concentration of 10./40

# %%
10./40

# %%
