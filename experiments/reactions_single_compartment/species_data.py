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
# ### A demo/tutorial on the data structure for Species and SpeciesRegistry
#
# **"Species"** is very broadly defined to be a molecule or larger entity that the simulator acts on.   
# For now, small molecules, enzymes and macromolecules.   
# NOT to be confused with ecological species (homo sapiens, etc)!

# %% [markdown]
# ### TAGS :   "quick-start", "basic"

# %%
LAST_REVISED = "July 13, 2026"
LIFE123_VERSION = "1.0.0rc8"     # Library version this experiment is based on

# %%
#import set_path            # Using MyBinder?  Uncomment this before running the next cell!
                            # Importing this module will add the project's home directory to sys.path

# %%
#import sys, os
#os.getcwd()
#sys.path.append("C:/some_path/my_env_or_install")   # CHANGE to the folder containing your venv or libraries installation!
# NOTE: If any of the imports below can't find a module, uncomment the lines above, or try:  import set_path

from life123 import check_version, Species, SpeciesRegistry

# %%
check_version(LIFE123_VERSION)    # To check compatibility

# %%

# %% [markdown]
# ## 1. Individual species

# %%
# Minimal!
s1 = Species("A")   # The first, and only required, argument is a unique "id"

s1    # Notice the various attributes (fields); by default, "name" and "label" are made equal to the supplied "id"

# %%

# %%
# Additional properties may be set at initialization...

s2 = Species(
            id="metK",
            name="methionine adenosyltransferase",
            categories=["protein", "enzyme"],
            molecular_weight=43000,
            ec_number="2.5.1.6",
            metadata={"compartment": "cytosol"},
            plot_color="darkturquoise"
        )

print(s2)

# %%
# ... or added later

s2.diffusion_rate=123.   # s1 is a "Species" object (a python dataclass), and its attributes may be read/set the usual way

print(s2)

# %%
print(s2.plot_color)

# %%
# A number of validations are automatically performed

try:
    s2.diffusion_rate = -1       # Negative diffusion!
except Exception as ex:
    print("What were you thinking?!? ", ex)

# %%
print(s2)       # Unchanged

# %%

# %%

# %%

# %% [markdown]
# ## 2. Registry of species
# Keeping together all your species data

# %%
# Good for testing; auto-assigns id's, names and labels: "A", "B", "C", ....
test_registry = SpeciesRegistry(n_species=3)
test_registry

# %%
test_registry.as_dataframe()

# %%

# %%

# %%
# If you already have existing "Species" objects, you can initialize a SpeciesRegistry with them
sr = SpeciesRegistry(species=[s1,s2])
sr.as_recordset()

# %%
sr.as_dataframe()

# %%
sr.number_of_species()

# %%

# %%

# %%
# Alternatively, you can start with an empty SpeciesRegistry, and then add Species later
sr_2 = SpeciesRegistry()
sr_2.as_recordset()

# %%
sr_2.add_species("X")       # Minimal!

# %%
sr_2.as_recordset()     # `name` and `label`, not directly supplied, were by default made equal to `id`

# %%
sr_2.add_species("Y", annotation="need to re-test")

# %%

# %%

# %%
# A SpeciesRegistry can also be started with lists/tuples of data
rome_cocktail = SpeciesRegistry(id=["S", "P", "Q", "R"], molecular_weight=[1200, 600, 2900, 1500], annotation=["Senatus",  "Populusque", None, "Romanus"])
rome_cocktail.as_dataframe()

# %%
# Several functions help deal with setting and reading values
rome_cocktail.update(species_id="Q", diffusion_rate=5, plot_color="red")

rome_cocktail.set_value(species_id="R", field="diffusion_rate", value=8)   # Alternate way to set a single value

rome_cocktail.as_dataframe()

# %%
rome_cocktail.get_value(species_id="Q", field="molecular_weight")

# %%
rome_cocktail.get_all_values(field="molecular_weight")

# %%

# %%
# Some search/aggregation functions
rome_cocktail.max_value(field="molecular_weight")

# %%
rome_cocktail.has_missing_values(field="diffusion_rate")   # Some of the diffusion rates are indeed missing

# %%
