import pytest
import numpy as np
from life123.species_index_map import SpeciesIndexMap
from life123.species_registry import SpeciesRegistry
from life123.uniform_compartment import UniformCompartment
from life123.species_registry import MacroMolecules



def test_CONSTRUCTOR():
    species_list = ["A", "B", "C"]
    sp = SpeciesRegistry(ids=species_list)
    ind = SpeciesIndexMap(species_ids=species_list, species_registry=sp)

    uc = UniformCompartment(species_index_map=ind)
    assert np.allclose(uc.system, np.array([0, 0, 0]))
    assert uc.species_registry == sp

    uc = UniformCompartment(species_index_map=ind, species_registry=sp)
    assert np.allclose(uc.system, np.array([0, 0, 0]))
    assert uc.species_registry == sp

    # Incompatible args `species_registry` and `species_index_map`
    with pytest.raises(Exception):
        UniformCompartment(species_index_map=ind, species_registry=SpeciesRegistry())




###############   TO READ SYSTEM DATA  ###############


def test_number_of_system_species():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])

    uc = UniformCompartment(species_index_map=ind)
    assert uc.number_of_system_species() == 3

    uc.add_species("X")
    assert uc.number_of_system_species() == 4

    uc.add_species("B")
    assert uc.number_of_system_species() == 4   # Already in place; nothing added

    ind.add_species("ZZZ")
    assert uc.species_index_map.index_to_species == ["A", "B", "C", "X", "ZZZ"]
    assert uc.number_of_system_species() == 4   # Unaffected by the larger index

    uc.system = np.array([1, 2, 3, 4, 5, 6])    # Artificially exceed the species index size
    with pytest.raises(Exception):
        uc.number_of_system_species()



def test_get_system_conc():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    uc = UniformCompartment(species_index_map=ind)
    assert np.allclose(uc.get_system_conc(), np.array([0, 0, 0]))

    uc.set_conc(conc={"B": 4, "C": 88})
    assert np.allclose(uc.get_system_conc(), np.array([0, 4, 88]))



def test_get_species_conc():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    uc = UniformCompartment(species_index_map=ind)

    uc.set_conc(conc={"B": 4, "C": 88})
    assert uc.get_species_conc("A") == 0
    assert uc.get_species_conc("B") == 4
    assert uc.get_species_conc("C") == 88

    with pytest.raises(Exception):
        uc.get_species_conc("UNKNOWN")



def test_get_conc_dict():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    uc = UniformCompartment(species_index_map=ind)

    uc.set_conc(conc={"A": 7, "B": 4, "C": 88})

    assert uc.get_conc_dict() == {"A": 7, "B": 4, "C": 88}

    assert uc.get_conc_dict(species_ids=["B"]) == {"B": 4}
    assert uc.get_conc_dict(species_ids=("B", "C")) == {"B": 4, "C": 88}
    assert uc.get_conc_dict(species_ids={"B", "A", "C"}) == {"A": 7, "B": 4, "C": 88}

    assert uc.get_conc_dict(system_data=np.array([1, 2, 3])) == {"A": 1, "B": 2, "C": 3}

    with pytest.raises(Exception):
        uc.get_conc_dict(system_data=np.array([1, 2, 3, 4]))

    assert uc.get_conc_dict(species_ids={"B", "A", "C"}, system_data=np.array([1, 2, 3, 4])) \
                    == {"A": 1, "B": 2, "C": 3}

    with pytest.raises(Exception):
        uc.get_conc_dict(species_ids={"B", "A", "C"}, system_data=np.array([1, 2]))




###############   TO SET SYSTEM DATA  ###############

def test_add_species():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    uc = UniformCompartment(species_index_map=ind)
    assert np.allclose(uc.system, np.array([0, 0, 0]))

    result = uc.add_species("X")
    assert result == 1
    assert np.allclose(uc.system, np.array([0, 0, 0, 0]))   # Added one
    assert uc.species_index_map.index_to_species == ["A", "B", "C", "X"]

    result = uc.add_species("B")
    assert result == 0
    assert np.allclose(uc.system, np.array([0, 0, 0, 0]))   # Already in place; nothing added
    assert uc.species_index_map.index_to_species == ["A", "B", "C", "X"]

    result = uc.add_species( {"X", "Z", "H"} )
    assert result == 2
    assert np.allclose(uc.system, np.array([0, 0, 0, 0, 0, 0]))   # Added two
    assert uc.species_index_map.index_to_species == ["A", "B", "C", "X", "H", "Z"]  # Notice the sorting of the newly-added set



def test_set_single_conc():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    uc = UniformCompartment(species_index_map=ind)
    assert np.allclose(uc.system, np.array([0, 0, 0]))

    with pytest.raises(Exception):
        uc.set_single_conc(conc=123)

    with pytest.raises(Exception):
        uc.set_single_conc(conc=123, species_id="B", species_index=1) # Can't pass both args


    uc.set_single_conc(conc=123, species_id="B")
    assert np.allclose(uc.system, np.array([0, 123, 0]))

    with pytest.raises(Exception):
        uc.set_single_conc(conc=3.1415, species_id="UNKNOWN")

    with pytest.raises(Exception):
        uc.set_single_conc(conc=-1, species_id="B")     # Negative!


    uc.set_single_conc(conc=88, species_index=0)
    assert np.allclose(uc.system, np.array([88, 123, 0]))

    assert uc.species_index_map.index_to_species == ["A", "B", "C"]
    with pytest.raises(Exception):
        uc.set_single_conc(conc=666, species_index=4, allow_growth=True)  # No such index in uc.species_index_map

    uc.species_index_map.add_species({"P", "Q"})

    with pytest.raises(Exception):
        # The index is now fine, but needs allow_growth=True
        uc.set_single_conc(conc=666, species_index=4, allow_growth=False)

    uc.set_single_conc(conc=666, species_index=4, allow_growth=True)
    assert np.allclose(uc.system, np.array([88, 123, 0, 0, 666]))   # The array got expanded accordingly

    with pytest.raises(Exception):
        uc.set_single_conc(conc=666, species_index=-1, allow_growth=True)


    l = []
    def foo(caption :str, arg=l):   # Use a mutable argument, to affect the passed variable l
        l.append(caption)

    uc.set_single_conc(conc=1000, species_id="C", snapshot=foo)
    assert l[0] == "Set concentration of `C`"
    assert np.allclose(uc.system, np.array([88, 123, 1000, 0, 666]))



def test_set_conc():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    uc = UniformCompartment(species_index_map=ind)
    assert np.allclose(uc.system, np.array([0, 0, 0]))

    with pytest.raises(Exception):
        uc.set_conc(conc="What's a list doing here?")

    with pytest.raises(Exception):
        uc.set_conc(conc={"UNKNOWN": 666})

    uc.set_conc(conc={"B": 123})
    assert np.allclose(uc.system, np.array([0, 123, 0]))

    uc.set_conc(conc={"B": 4, "C": 88})
    assert np.allclose(uc.system, np.array([0, 4, 88]))

    uc.species_index_map.add_species({"P", "Q"})
    uc.set_conc(conc={"Q": 999})
    assert np.allclose(uc.system, np.array([0, 4, 88, 0, 999]) )     # The system state got expanded


    l = []
    def foo_2(caption :str, arg=l):   # Use a mutable argument, to affect the passed variable l
        l.append(caption)

    uc.set_conc(conc={"P": 1, "A": 55}, snapshot=foo_2)
    assert l[0] == "Set concentration"
    assert np.allclose(uc.system, np.array([55, 4, 88, 1, 999]) )

    #TODO: test the diagnostic option
