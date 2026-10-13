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
    assert np.allclose(uc.system_state, np.array([0, 0, 0]))
    assert uc.species_registry == sp

    uc = UniformCompartment(species_index_map=ind, species_registry=sp)
    assert np.allclose(uc.system_state, np.array([0, 0, 0]))
    assert uc.species_registry == sp

    # Incompatible args `species_registry` and `species_index_map`
    with pytest.raises(Exception):
        UniformCompartment(species_index_map=ind, species_registry=SpeciesRegistry())




###############   TO READ SYSTEM DATA  ###############

def test_get_species_registry():
    pass    # TODO



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

    uc.system_state = np.array([1, 2, 3, 4, 5, 6])    # Artificially exceed the species index size
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

    with pytest.raises(Exception):
        uc.get_conc_dict(species_ids=["UNKNOWN"])





###############   TO SET SYSTEM DATA  ###############

def test_add_species():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    uc = UniformCompartment(species_index_map=ind)
    assert np.allclose(uc.system_state, np.array([0, 0, 0]))

    result = uc.add_species("X")
    assert result == 1
    assert np.allclose(uc.system_state, np.array([0, 0, 0, 0]))   # Added one
    assert uc.species_index_map.index_to_species == ["A", "B", "C", "X"]

    result = uc.add_species("B")
    assert result == 0
    assert np.allclose(uc.system_state, np.array([0, 0, 0, 0]))   # Already in place; nothing added
    assert uc.species_index_map.index_to_species == ["A", "B", "C", "X"]

    result = uc.add_species( {"X", "Z", "H"} )
    assert result == 2
    assert np.allclose(uc.system_state, np.array([0, 0, 0, 0, 0, 0]))   # Added two
    assert uc.species_index_map.index_to_species == ["A", "B", "C", "X", "H", "Z"]  # Notice the sorting of the newly-added set



def test_set_system_state():
    pass    # TODO



def test_set_single_conc():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    uc = UniformCompartment(species_index_map=ind)
    assert np.allclose(uc.system_state, np.array([0, 0, 0]))

    with pytest.raises(Exception):
        uc.set_single_conc(conc=123)

    with pytest.raises(Exception):
        uc.set_single_conc(conc=123, species_id="B", species_index=1) # Can't pass both args


    uc.set_single_conc(conc=123, species_id="B")
    assert np.allclose(uc.system_state, np.array([0, 123, 0]))

    with pytest.raises(Exception):
        uc.set_single_conc(conc=3.1415, species_id="UNKNOWN")

    with pytest.raises(Exception):
        uc.set_single_conc(conc=-1, species_id="B")     # Negative!


    uc.set_single_conc(conc=88, species_index=0)
    assert np.allclose(uc.system_state, np.array([88, 123, 0]))

    assert uc.species_index_map.index_to_species == ["A", "B", "C"]
    with pytest.raises(Exception):
        uc.set_single_conc(conc=666, species_index=4, allow_growth=True)  # No such index in uc.species_index_map

    uc.species_index_map.add_species({"P", "Q"})

    with pytest.raises(Exception):
        # The index is now fine, but needs allow_growth=True
        uc.set_single_conc(conc=666, species_index=4, allow_growth=False)

    uc.set_single_conc(conc=666, species_index=4, allow_growth=True)
    assert np.allclose(uc.system_state, np.array([88, 123, 0, 0, 666]))   # The array got expanded accordingly

    with pytest.raises(Exception):
        uc.set_single_conc(conc=666, species_index=-1, allow_growth=True)


    l = []
    def foo(caption :str, arg=l):   # Use a mutable argument, to affect the passed variable l
        l.append(caption)

    uc.set_single_conc(conc=1000, species_id="C", snapshot=foo)
    assert l[0] == "Set concentration of `C`"
    assert np.allclose(uc.system_state, np.array([88, 123, 1000, 0, 666]))



def test_set_conc():
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    uc = UniformCompartment(species_index_map=ind)
    assert np.allclose(uc.system_state, np.array([0, 0, 0]))

    with pytest.raises(Exception):
        uc.set_conc(conc="What's a list doing here?")

    with pytest.raises(Exception):
        uc.set_conc(conc={"UNKNOWN": 666})

    uc.set_conc(conc={"B": 123})
    assert np.allclose(uc.system_state, np.array([0, 123, 0]))

    uc.set_conc(conc={"B": 4, "C": 88})
    assert np.allclose(uc.system_state, np.array([0, 4, 88]))

    uc.species_index_map.add_species({"P", "Q"})
    uc.set_conc(conc={"Q": 999})
    assert np.allclose(uc.system_state, np.array([0, 4, 88, 0, 999]))     # The system state got expanded


    l = []
    def foo_2(caption :str, arg=l):   # Use a mutable argument, to affect the passed variable l
        l.append(caption)

    uc.set_conc(conc={"P": 1, "A": 55}, snapshot_funct=foo_2)
    assert l[0] == "Set concentration"
    assert np.allclose(uc.system_state, np.array([55, 4, 88, 1, 999]))

    #TODO: test the diagnostic option




####################   MACROMOLECULES   ####################

def test_set_macromolecules():
    sr = SpeciesRegistry(ids=["A", "B", "M1", "M2"])
    mm = MacroMolecules(sr)

    mm.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=3)
    mm.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=5)
    mm.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="B", Kd=11)
    mm.set_binding_site_affinity(macromolecule="M2", site_number=3, ligand="B", Kd=102)

    uc = UniformCompartment(species_registry=sr, macromolecules=mm)
    uc.set_macromolecules()    # By default, set counts to 1 for all the registered macromolecules

    assert uc.macro_system == {"M1": 1, "M2": 1}
    assert len(uc.macro_system_state) == 2
    assert uc.macro_system_state["M1"] == {1: ("A", 0.), 2: ("B", 0.)}
    assert uc.macro_system_state["M2"] == {1: ("B", 0.), 3: ("B", 0.)}


    # Over-write the previous settings
    uc.set_macromolecules({"M2": 4})

    assert uc.macro_system == {"M2": 4}
    assert len(uc.macro_system_state) == 1
    assert uc.macro_system_state["M2"] == {1: ("B", 0.), 3: ("B", 0.)}


    with pytest.raises(Exception):
        uc.set_macromolecules({"M999": 2})        # Unknown macromolecule

    sr.add_single_species("M999")

    # Over-write the previous settings
    uc.set_macromolecules({"M999": 2})
    assert uc.macro_system == {"M999": 2}
    assert len(uc.macro_system_state) == 1
    assert uc.macro_system_state["M999"] == {}     # No binding sites were registered

    mm.set_binding_site_affinity(macromolecule="M999", site_number=12, ligand="A", Kd=4.5)
    assert uc.macro_system_state["M999"] == {}     # The system state has not been updated yet
    # Over-write the previous settings
    uc.set_macromolecules({"M999": 2})
    assert uc.macro_system_state["M999"] == {12: ("A", 0.)}



def test_get_occupancy():
    sr = SpeciesRegistry(ids=["A", "B", "C", "M1", "M2"])
    mm = MacroMolecules(sr)

    uc = UniformCompartment(species_registry=sr, macromolecules=mm)

    with pytest.raises(Exception):
        # The system state for macromolecules has not been set yet
        uc.get_occupancy(macromolecule="M1", site_number=1)

    with pytest.raises(Exception):
        # No occupancy data yet set for macromolecule `Unknown`
        uc.get_occupancy(macromolecule="Unknown", site_number=1)

    with pytest.raises(Exception):
        # No occupancy data yet set for site number 1 of macromolecule `M1`
        uc.get_occupancy(macromolecule="M1", site_number=1)

    mm.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="C", Kd=23.5)
    uc.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=0.5)

    assert uc.macro_system == {"M1": 1}
    assert uc.macro_system_state == {"M1": {1: ("C", 0.5)} }

    assert uc.get_occupancy(macromolecule="M1", site_number=1) == 0.5



def test_set_occupancy():
    sr = SpeciesRegistry(ids=["A", "B", "M1", "M2"])
    mm = MacroMolecules(sr)

    rxn = UniformCompartment(species_registry=sr, macromolecules=mm)

    with pytest.raises(Exception):
        # Occupancy out of range
        rxn.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=1.1)

    with pytest.raises(Exception):
        # Occupancy out of range
        rxn.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=-0.3)

    with pytest.raises(Exception):
        # Unknown macromolecule
        rxn.set_occupancy(macromolecule="Unknown", site_number=1, fractional_occupancy=0.5)

    with pytest.raises(Exception):
        # No binding sites are defined on macromolecule M1
        rxn.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=0.5)

    mm.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="B", Kd=6.)
    rxn.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=0.5)

    assert rxn.macro_system == {"M1": 1}
    assert rxn.macro_system_state == {"M1": {1: ("B", 0.5)}}

    assert rxn.get_occupancy(macromolecule="M1", site_number=1) == 0.5
