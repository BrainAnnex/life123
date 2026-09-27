from life123 import SpeciesRegistry, UniformCompartment
from life123.react_diffuse import ReactDiffuse
from life123.index_species import IndexSpecies
from life123 import BioSim1D, BioSim2D, BioSim3D



def test_constructor():
    rd0 = ReactDiffuse(type="uniform")

    assert type(rd0.module) == UniformCompartment
    assert type(rd0.species_registry) == SpeciesRegistry
    assert type(rd0.index_species) == IndexSpecies

    assert rd0.species_registry.get_all_species_ids() == []
    assert rd0.index_species.number_of_system_species() == 0


    sr = SpeciesRegistry(n_species=3)
    rd1 = ReactDiffuse(type="1d", species_registry=sr, n_bins=2)

    assert type(rd1.module) == BioSim1D
    assert type(rd1.species_registry) == SpeciesRegistry
    assert type(rd1.index_species) == IndexSpecies

    assert rd1.species_registry.get_all_species_ids() == ["A", "B", "C"]
    assert rd1.index_species.number_of_system_species() == 3
