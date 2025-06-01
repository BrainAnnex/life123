# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# # TEST OF Plotly
# with a peek at its data structure

# %%
import plotly.express as px

# %%

# %% [markdown]
# ### Test of a line figure

# %%
fig = px.line(x=["a","b","c"], y=[1,3,2], title="sample figure")
fig.show()

# %%
print(fig)

# %%
print(fig.data)

# %%
print(fig.layout)

# %%

# %%

# %% [markdown]
# ### Test of a barplot

# %%
fig = px.bar(x=["a", "b", "c"], y=[1, 3, 2])
fig.show()

# %%
import plotly.graph_objects as go

# %%
fig_widget = go.FigureWidget(fig)
fig_widget

# %%
