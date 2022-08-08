# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
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
# **TEST OF plotly**

# %%
import plotly.express as px

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
