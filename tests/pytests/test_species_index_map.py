import pytest
from life123.species_index_map import SpeciesIndexMap


def test_CONSTRUCTOR():
    ind = SpeciesIndexMap()

    assert ind.index_to_species == []
    assert ind.species_to_index == {}



def test_number_of_system_species():
    ind = SpeciesIndexMap()

    assert ind.number_of_system_species() == 0

    ind.add_species({"A"})
    assert ind.number_of_system_species() == 1
    assert len(ind.species_to_index) == 1

    ind.add_species(["X", "Y", "Z"])
    assert ind.number_of_system_species() == 4
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
        assert ind.species_at(3)    # Out of range



def test_add_species():
    ind = SpeciesIndexMap()
    assert ind.index_to_species == []
    assert ind.species_to_index == {}

    ind.add_species({"A"})
    assert ind.index_to_species == ["A"]
    assert ind.species_to_index == {"A": 0}

    ind.add_species(["Z", "H"])
    assert ind.index_to_species == ["A", "H", "Z"]  # Notice the automatic sorting
    assert ind.species_to_index == {"A": 0, "H": 1, "Z": 2}

    ind.add_species( ("F", "E") )
    assert ind.index_to_species == ["A", "H", "Z", "E", "F"]  # Notice the automatic group-wise sorting
    assert ind.species_to_index == {"A": 0, "H": 1, "Z": 2, "E": 3, "F": 4}



def test_clear_index():
    ind = SpeciesIndexMap()

    ind.add_species(["Z", "H"])
    ind.clear_index()
    assert ind.index_to_species == []
    assert ind.species_to_index == {}
