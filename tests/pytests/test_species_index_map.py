import pytest
from life123.species_index_map import SpeciesIndexMap
from life123.species_registry import SpeciesRegistry


def test_CONSTRUCTOR():
    ind = SpeciesIndexMap()

    assert ind.index_to_species == []
    assert ind.species_to_index == {}


    sp = SpeciesRegistry()
    with pytest.raises(Exception):
        SpeciesIndexMap(species_ids={"A"}, species_registry=sp)     # Unknown species id

    sp.add_single_species("A")
    ind = SpeciesIndexMap(species_ids={"A"}, species_registry=sp)
    assert ind.index_to_species == ["A"]
    assert ind.species_to_index == {"A": 0}



def test_number_of_mapped_species():
    ind = SpeciesIndexMap()

    assert ind.number_of_mapped_species() == 0

    ind.add_species({"A"})
    assert ind.number_of_mapped_species() == 1
    assert len(ind.species_to_index) == 1

    ind.add_species(["X", "Y", "Z"])
    assert ind.number_of_mapped_species() == 4
    assert len(ind.species_to_index) == 4
    assert len(ind.species_to_index) == 4



def test_index_of():
    ind = SpeciesIndexMap()

    with pytest.raises(Exception):
        ind.index_of("A")

    ind.add_species({"A"})
    assert ind.index_of("A") == 0

    ind.add_species( ("X", "Y") )
    assert ind.index_of("X") == 1
    assert ind.index_of("Y") == 2

    with pytest.raises(Exception):
        ind.index_of("ZZZ")



def test_species_at():
    ind = SpeciesIndexMap()

    with pytest.raises(Exception):
        ind.species_at(0)

    ind.add_species({"A"})
    assert ind.species_at(0) == "A"

    ind.add_species( ["X", "Y"] )
    assert ind.species_at(1) == "X"
    assert ind.species_at(2) == "Y"

    with pytest.raises(Exception):
        ind.species_at(3)    # Out of range

    with pytest.raises(Exception):
        ind.species_at(-1)    # Out of range



def test_add_species():
    ind = SpeciesIndexMap()
    assert ind.index_to_species == []
    assert ind.species_to_index == {}

    ind.add_species({"A"})
    assert ind.index_to_species == ["A"]
    assert ind.species_to_index == {"A": 0}

    ind.add_species({"Z", "H"})
    assert ind.index_to_species == ["A", "H", "Z"]  # Notice the automatic sorting of the added set
    assert ind.species_to_index == {"A": 0, "H": 1, "Z": 2}

    ind.add_species( ("F", "E") )
    assert ind.index_to_species == ["A", "H", "Z", "F", "E"]  # Notice that the specified tuple order is respected
    assert ind.species_to_index == {"A": 0, "H": 1, "Z": 2, "F": 3, "E": 4}

    ind.add_species( ("W", "C") )
    assert ind.index_to_species == ["A", "H", "Z", "F", "E", "W", "C"]  # Notice that the specified tuple order is respected
    assert ind.species_to_index == {"A": 0, "H": 1, "Z": 2, "F": 3, "E": 4, "W": 5, "C": 6}



def test_clear_index():
    ind = SpeciesIndexMap()

    ind.add_species(["Z", "H"])
    ind.clear_index()
    assert ind.index_to_species == []
    assert ind.species_to_index == {}
