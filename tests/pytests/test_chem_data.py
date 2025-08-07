import pytest
import numpy as np
import pandas as pd
from life123 import ChemData
from life123.chem_data import Diffusion
from tests.utilities.comparisons import *


#############  ChemCore  #############

def test_number_of_chemicals():
    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.number_of_chemicals() == 3

    chem_data = ChemData(names=[])
    assert chem_data.number_of_chemicals() == 0



def test_assert_valid_species_index():
    chem_data = ChemData(names=['A', 'B', 'C'])
    chem_data.assert_valid_chem_index(0)
    chem_data.assert_valid_chem_index(1)
    chem_data.assert_valid_chem_index(2)

    with pytest.raises(Exception):
        chem_data.assert_valid_chem_index(3)     # Too large

    with pytest.raises(Exception):
        chem_data.assert_valid_chem_index(-1)    # Too small

    with pytest.raises(Exception):
        chem_data.assert_valid_chem_index("2")   # Not an int



def test_get_index():
    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.get_index('A') == 0
    assert chem_data.get_index('B') == 1
    assert chem_data.get_index('C') == 2
    with pytest.raises(Exception):
        assert chem_data.get_index('X')     # Not found



def test_name_exists():
    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.label_exists("A")
    assert chem_data.label_exists("C")
    assert not chem_data.label_exists("X")

    chem_data.add_chemical(name="X")
    assert chem_data.label_exists("X")

    assert not chem_data.label_exists("Z")



def test_get_name():
    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.get_label(0) == 'A'
    assert chem_data.get_label(1) == 'B'
    assert chem_data.get_label(2) == 'C'

    with pytest.raises(Exception):
        chem_data.get_label(3)               # Out of bounds
    with pytest.raises(Exception):
        chem_data.get_label(-1)              # Invalid argument value
    with pytest.raises(Exception):
        chem_data.get_label("some string")   # Invalid argument type
    with pytest.raises(Exception):
        chem_data.get_label(3.14)            # Invalid argument type



def test_get_all_names():
    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.get_all_labels() == ['A', 'B', 'C']

    # Reach into the internal data structure, to mess up some names
    chem_data.chemical_data = [{'name': 'A'}, {'name': ''}, {'name': 'C'}]

    with pytest.raises(Exception):
        chem_data.get_all_labels()   # The former name 'B' is now a blank

    chem_data.chemical_data = [{'name': 'A'}, {'name': 'B'}, {}]

    with pytest.raises(Exception):
        chem_data.get_all_labels()   # The former name 'B' is now missing



def test_all_chemicals():
    chem_data = ChemData(names=['A', 'B'])
    chem_data.add_chemical(label="C", note="this is C", name="CH4", cost="200")
    chem_data.add_chemical_with_diffusion("D", diff_rate=10, note="this is D")

    result = chem_data.all_chemicals()
    expected_recordset = [{'name': 'A', 'label': 'A'}, {'name': 'B', 'label': 'B'},
                          {'name': 'CH4', 'label': 'C', 'note': 'this is C', 'cost': '200'},
                          {'name': 'D', 'label': 'D', 'note': 'this is D'}]
    expected_df = pd.DataFrame(expected_recordset)

    assert expected_df.equals(result)



def test_add_chemical():
    chem_data = ChemData()
    assert chem_data.number_of_chemicals() == 0
    assert chem_data.chemical_data == []

    result = chem_data.add_chemical(name="A")
    assert result == 0
    assert chem_data.number_of_chemicals() == 1
    assert chem_data.chemical_data == [{"name": "A", "label": "A"}]
    assert chem_data.label_dict == {"A": 0}

    with pytest.raises(Exception):
        chem_data.add_chemical(name="A")     # Duplicate!

    result = chem_data.add_chemical(name="B", note="some note")
    assert result == 1
    assert chem_data.number_of_chemicals() == 2
    assert chem_data.chemical_data == [ {"name": "A", "label": "A"},
                                        {"name": "B", "label": "B", "note": "some note"}]
    assert chem_data.label_dict == {"A": 0, "B": 1}

    result = chem_data.add_chemical(name="C")
    assert result == 2
    assert chem_data.number_of_chemicals() == 3
    assert chem_data.chemical_data == [{"name": "A", "label": "A"},
                                       {"name": "B", "label": "B", "note": "some note"},
                                       {"name": "C", "label": "C"}]
    assert chem_data.label_dict == {"A": 0, "B": 1, "C": 2}

    result = chem_data.add_chemical(name="Some long name", label="D")
    assert result == 3
    assert chem_data.number_of_chemicals() == 4
    assert chem_data.chemical_data == [{"name": "A", "label": "A"},
                                       {"name": "B", "label": "B", "note": "some note"},
                                       {"name": "C", "label": "C"},
                                       {"name": "Some long name", "label": "D"}]
    assert chem_data.label_dict == {"A": 0, "B": 1, "C": 2, "D": 3}

    with pytest.raises(Exception):
        chem_data.add_chemical(name="Some long name")   # Duplicate name!

    with pytest.raises(Exception):
        chem_data.add_chemical(name="D")                # Name cannot be same as an existing label!

    with pytest.raises(Exception):
        chem_data.add_chemical(name="Z", label="A")     # Duplicate label!

    with pytest.raises(Exception):
        chem_data.add_chemical(name="Z", label="Some long name")    # Label cannot be same as an existing name!

    with pytest.raises(Exception):
        chem_data.add_chemical(name=123)    # Name is not a string

    with pytest.raises(Exception):
        chem_data.add_chemical(name="")     # Missing name


    # Re-start
    chem_data = ChemData(names="X")

    result = chem_data.add_chemical("Y", note="test")
    assert result == 1
    assert chem_data.number_of_chemicals() == 2
    assert chem_data.chemical_data == [{"name": "X", "label": "X"}, {"name": "Y", "label": "Y", "note": "test"}]
    assert chem_data.label_dict == {"X": 0, "Y": 1}

    result = chem_data.add_chemical(label="Z", name="CH3OH")
    assert result == 2
    assert chem_data.number_of_chemicals() == 3
    assert chem_data.chemical_data == [{"name": "X", "label": "X"},
                                       {"name": "Y", "label": "Y", "note": "test"},
                                       {"name": "CH3OH", "label": "Z"}]
    assert chem_data.label_dict == {"X": 0, "Y": 1, "Z": 2}





#############  Diffusion  #############


def test_add_chemical_with_diffusion():
    chem_data = ChemData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2, 3.14])
    assert chem_data.number_of_chemicals() == 3
    assert chem_data.get_all_labels() == ['A', 'B', 'C']
    assert chem_data.label_dict == {'A': 0, 'B': 1, 'C': 2}
    assert chem_data.chemical_data == [ {'name': 'A', 'label': 'A'}, {'name': 'B', 'label': 'B'}, {'name': 'C', 'label': 'C'} ]
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2, 3.14])

    chem_data.add_chemical_with_diffusion(name="D", diff_rate=8)
    assert chem_data.number_of_chemicals() == 4
    assert chem_data.get_all_labels() == ['A', 'B', 'C', 'D']
    assert chem_data.label_dict == {'A': 0, 'B': 1, 'C': 2, 'D': 3}
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2, 3.14, 8])
    assert chem_data.chemical_data == [ {'name': 'A', 'label': 'A'}, {'name': 'B', 'label': 'B'}, {'name': 'C', 'label': 'C'},
                                        {'name': 'D', 'label': 'D'}]

    with pytest.raises(Exception):
        chem_data.add_chemical_with_diffusion(name="E", diff_rate="I'm not a number")

    with pytest.raises(Exception):
        chem_data.add_chemical_with_diffusion(name="E", diff_rate=-666.)

    with pytest.raises(Exception):
        chem_data.add_chemical_with_diffusion(name=666, diff_rate=25.)    # Wrong type for name

    chem_data.add_chemical_with_diffusion(name="E", diff_rate=25., note="my note", some_field="some value")
    assert chem_data.number_of_chemicals() == 5
    assert chem_data.get_all_labels() == ['A', 'B', 'C', 'D', 'E']
    assert chem_data.label_dict == {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2, 3.14, 8, 25])
    assert chem_data.chemical_data == [ {'name': 'A', 'label': 'A'}, {'name': 'B', 'label': 'B'}, {'name': 'C', 'label': 'C'},
                                        {'name': 'D', 'label': 'D'},
                                        {'name': 'E',  'label': 'E', 'note': 'my note', 'some_field': 'some value'}]

    with pytest.raises(Exception):
        chem_data.add_chemical_with_diffusion(name="", diff_rate=25.)



def test_set_diffusion_rate():
    pass



def test_assert_valid_diffusion():
    chem_data = ChemData()

    chem_data.assert_valid_diffusion(3.14)
    with pytest.raises(Exception):
        chem_data.assert_valid_diffusion("Not a number!")

    with pytest.raises(Exception):
        chem_data.assert_valid_diffusion(-0.001)



def test_get_diffusion_rate():
    chem_data = ChemData(names=["A", "B", "C"], diffusion_rates=[10, 11, 12])

    assert chem_data.get_diffusion_rate(name="A") == 10
    assert chem_data.get_diffusion_rate(chem_index=0) == 10
    with pytest.raises(Exception):
        assert chem_data.get_diffusion_rate(name="A", chem_index=0)

    assert chem_data.get_diffusion_rate(name="C") == 12
    assert chem_data.get_diffusion_rate(chem_index=2) == 12

    with pytest.raises(Exception):
        chem_data.get_diffusion_rate(name="X")          # X doesn't exist

    with pytest.raises(Exception):
        chem_data.get_diffusion_rate(chem_index=3)   # Index 3 doesn't exist

    chem_data.add_chemical(name="Z")                    # This will be given index 3
    assert chem_data.get_diffusion_rate(name="Z") is None       # No diffusion value assigned
    assert chem_data.get_diffusion_rate(chem_index=3) is None

    chem_data.set_diffusion_rate(label="Z", diff_rate=8)
    assert chem_data.get_diffusion_rate(name="Z") == 8
    assert chem_data.get_diffusion_rate(chem_index=3) == 8



def test_get_all_diffusion_rates():
    diff = Diffusion()

    result = diff.get_all_diffusion_rates()
    assert result == []

    diff.add_chemical_with_diffusion(name="A", diff_rate=1)
    result = diff.get_all_diffusion_rates()
    assert result == [1]

    diff.set_diffusion_rate(label="A", diff_rate=10)
    result = diff.get_all_diffusion_rates()
    assert result == [10]

    diff.add_chemical(name="B")
    result = diff.get_all_diffusion_rates()
    assert result == [10, None]

    diff.add_chemical_with_diffusion(name="C", diff_rate=8)
    result = diff.get_all_diffusion_rates()
    assert result == [10, None, 8]

    diff.set_diffusion_rate(label="B", diff_rate=3)
    result = diff.get_all_diffusion_rates()
    assert result == [10, 3, 8]



def test_missing_diffusion_rate():
    diff = Diffusion()

    result = diff.missing_diffusion_rate()
    assert result == False

    diff.add_chemical_with_diffusion(name="A", diff_rate=1)
    result = diff.missing_diffusion_rate()
    assert result == False

    diff.set_diffusion_rate(label="A", diff_rate=10)
    result = diff.missing_diffusion_rate()
    assert result == False

    diff.add_chemical(name="B")
    result = diff.missing_diffusion_rate()
    assert result == True       # B's diffusion rate is missing

    diff.add_chemical_with_diffusion(name="C", diff_rate=8)
    result = diff.missing_diffusion_rate()
    assert result == True       # B's diffusion rate is still missing

    diff.set_diffusion_rate(label="B", diff_rate=3)
    result = diff.missing_diffusion_rate()
    assert result == False






#############  ChemData  #############

def test_constructor():
    chem_data = ChemData()
    assert chem_data.chemical_data == []
    assert chem_data.label_dict == {}
    assert chem_data.color_dict == {}
    assert chem_data.diffusion_rates == {}


    with pytest.raises(Exception):
        ChemData(names=123)         # Not a list/tuple/str


    with pytest.raises(Exception):
        ChemData(names=["A", 2])    # Some of the names aren't strings


    chem_data = ChemData(names='A')
    assert chem_data.chemical_data == [{"label": "A", "name": "A"}]
    assert chem_data.label_dict == {'A': 0}
    assert chem_data.color_dict == {}
    assert chem_data.diffusion_rates == {}

    chem_data = ChemData(names=['A', 'B'])
    assert chem_data.chemical_data == [{"label": "A", "name": "A"}, {"label": "B", "name": "B"}]
    assert chem_data.label_dict == {'A': 0, 'B': 1}
    assert chem_data.color_dict == {}
    assert chem_data.diffusion_rates == {}


    with pytest.raises(Exception):
        ChemData(labels=True)     # Not a list/tuple/str

    with pytest.raises(Exception):
        ChemData(labels=["L1", 2])    # Some of the labels aren't strings

    with pytest.raises(Exception):
        ChemData(names=["A", "B"], labels=["L1", 2])    # Some of the labels aren't strings

    chem_data = ChemData(labels='A')
    assert chem_data.chemical_data == [{"label": "A", "name": "A"}]
    assert chem_data.label_dict == {'A': 0}
    assert chem_data.color_dict == {}
    assert chem_data.diffusion_rates == {}

    chem_data = ChemData(labels=['A', 'B'])
    assert chem_data.chemical_data == [{"label": "A", "name": "A"}, {"label": "B", "name": "B"}]
    assert chem_data.label_dict == {'A': 0, 'B': 1}
    assert chem_data.color_dict == {}
    assert chem_data.diffusion_rates == {}


    chem_data = ChemData(names=['name1', 'name2'], labels=['A', 'B'])
    assert chem_data.chemical_data == [{"label": "A", "name": "name1"}, {"label": "B", "name": "name2"}]
    assert chem_data.label_dict == {'A': 0, 'B': 1}
    assert chem_data.color_dict == {}
    assert chem_data.diffusion_rates == {}


    with pytest.raises(Exception):
        ChemData(diffusion_rates=False)     # Not a list/tuple/numpy array/number

    with pytest.raises(Exception):
        ChemData(diffusion_rates=[3.1, "I'm not a number"])     # One of the values isn't a number

    with pytest.raises(Exception):
        ChemData(diffusion_rates=[3.1, -6.66])                  # Values cannot be negative

    chem_data = ChemData(diffusion_rates=33)    # No names, nor labels, specified (defaults are used)
    assert chem_data.chemical_data == [{"label": "A", "name": "A"}]
    assert chem_data.label_dict == {"A": 0}
    assert chem_data.color_dict == {}
    assert chem_data.diffusion_rates == {"A": 33}

    chem_data = ChemData(diffusion_rates=[0.15, 3.2])   # No names, nor labels, specified (defaults are used)
    assert chem_data.chemical_data == [{"label": "A", "name": "A"},
                                       {"label": "B", "name": "B"}]
    assert chem_data.label_dict == {"A": 0, "B": 1}
    assert chem_data.color_dict == {}
    assert len(chem_data.diffusion_rates) == 2
    assert np.allclose(0.15, chem_data.diffusion_rates["A"])
    assert np.allclose(3.2, chem_data.diffusion_rates["B"])


    with pytest.raises(Exception):
        ChemData(plot_colors=123)     # Not a list/tuple/str

    with pytest.raises(Exception):
        ChemData(plot_colors=["red", 22])    # Some of the colors aren't strings

    chem_data = ChemData(plot_colors="red")
    assert chem_data.chemical_data == [{"label": "A", "name": "A"}]
    assert chem_data.label_dict == {"A": 0}
    assert chem_data.color_dict == {"A": "red"}
    assert chem_data.diffusion_rates == {}

    chem_data = ChemData(plot_colors=["red", "blue"])
    assert chem_data.chemical_data == [{"label": "A", "name": "A"},
                                       {"label": "B", "name": "B"}]
    assert chem_data.label_dict == {"A": 0, "B": 1}
    assert chem_data.color_dict == {"A": "red", "B": "blue"}
    assert chem_data.diffusion_rates == {}


    # Various mismatch in counts
    with pytest.raises(Exception):
        ChemData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2])

    with pytest.raises(Exception):
        ChemData(names=['A', 'B', 'C'], labels="L")

    with pytest.raises(Exception):
        ChemData(diffusion_rates=[0.15, 1.2], labels="L")

    with pytest.raises(Exception):
        ChemData(diffusion_rates=[0.15, 1.2], plot_colors="red")

    with pytest.raises(Exception):
        ChemData(names=["Name1"], plot_colors=["red", "blue"])


    chem_data = ChemData(diffusion_rates=[10, 20], plot_colors=["red", "blue"])
    assert chem_data.chemical_data == [{"label": "A", "name": "A"},
                                       {"label": "B", "name": "B"}]
    assert chem_data.label_dict == {"A": 0, "B": 1}
    assert chem_data.color_dict == {"A": "red", "B": "blue"}
    assert chem_data.diffusion_rates == {"A": 10, "B": 20}

    chem_data = ChemData(names=["Name1", "Name2"], diffusion_rates=[10, 20], plot_colors=["red", "blue"])
    assert chem_data.chemical_data == [{"label": "Name1", "name": "Name1"},
                                       {"label": "Name2", "name": "Name2"}]
    assert chem_data.label_dict == {"Name1": 0, "Name2": 1}
    assert chem_data.color_dict == {"Name1": "red", "Name2": "blue"}
    assert chem_data.diffusion_rates == {"Name1": 10, "Name2": 20}

    chem_data = ChemData(names=["Name1", "Name2"], labels=["L1", "L2"],
                         diffusion_rates=[10, 20], plot_colors=["red", "blue"])
    assert chem_data.chemical_data == [{"label": "L1", "name": "Name1"},
                                       {"label": "L2", "name": "Name2"}]
    assert chem_data.label_dict == {"L1": 0, "L2": 1}
    assert chem_data.color_dict == {"L1": "red", "L2": "blue"}
    assert chem_data.diffusion_rates == {"L1": 10, "L2": 20}


    chem_data = ChemData(names=['A', 'B', 'C'], diffusion_rates=np.array([0.15, 2, 3.14]))
    assert len(chem_data.chemical_data) == 3
    assert chem_data.get_all_labels() == ['A', 'B', 'C']
    assert chem_data.label_dict == {'A': 0, 'B': 1, 'C': 2}
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 2, 3.14])





##########  MACRO-MOLECULES  ##########

def test_add_macromolecules():
    chem_data = ChemData()

    chem_data.add_macromolecules(["M1", "M2"])
    assert chem_data.get_macromolecules() == ["M1", "M2"]

    chem_data.add_macromolecules("M3")
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3"]

    chem_data.add_macromolecules(["M4", "M5"])
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3", "M4", "M5"]

    chem_data.add_macromolecules("M2")          # Redundant
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3", "M4", "M5"]

    chem_data.add_macromolecules(["M1", "M6"])  # Partially redundant
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3", "M4", "M5", "M6"]



def test_get_macromolecules():
    chem_data = ChemData()

    assert chem_data.get_macromolecules() == []
    chem_data.add_macromolecules(["M1", "M2"])
    assert chem_data.get_macromolecules() == ["M1", "M2"]
    chem_data.add_macromolecules(["M3", "M1"])
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3"]



def test_set_binding_site_affinity():
    chem_data = ChemData(["A", "B", "ZZZ"])
    chem_data.add_macromolecules(["M1", "M2"])

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=3)
    result = chem_data.get_binding_site_affinity(macromolecule="M1", site_number=1)
    assert result == ("A", 3)
    assert result.chemical == "A"
    assert result.Kd == 3

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=5)
    assert chem_data.get_binding_site_affinity("M1", site_number=1) == ("A", 3)
    assert chem_data.get_binding_site_affinity("M1", site_number=2) == ("B", 5)

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="B", Kd=11)
    assert chem_data.get_binding_site_affinity("M1", site_number=1) == ("A", 3)
    assert chem_data.get_binding_site_affinity("M1", site_number=2) == ("B", 5)
    assert chem_data.get_binding_site_affinity("M2", site_number=1) == ("B",11)

    # Over-write previous value of affinity of "A" to site 1 of macromolecule "M1"
    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=999)
    assert chem_data.get_binding_site_affinity("M1", site_number=1) == ("A", 999)
    assert chem_data.get_binding_site_affinity("M1", site_number=2) == ("B", 5)
    assert chem_data.get_binding_site_affinity("M2", site_number=1) == ("B",11)

    # Attempting to associate another chemical to site 1 of macromolecule "M1", will result in error
    with pytest.raises(Exception):
        chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="ZZZ", Kd=999)

    chem_data.set_binding_site_affinity(macromolecule="M3", site_number=10, ligand="A", Kd= 4)
    assert chem_data.get_binding_site_affinity("M3", site_number=10) == ("A", 4)
    assert chem_data.get_macromolecules() == ["M1", "M2", "M3"]     # "M3" got automatically added
    assert chem_data.get_binding_site_affinity("M1", site_number=1) == ("A", 999)
    assert chem_data.get_binding_site_affinity("M1", site_number=2) == ("B", 5)
    assert chem_data.get_binding_site_affinity("M2", site_number=1) == ("B",11)

    # Unknown chemical "X"
    with pytest.raises(Exception):
        chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="X", Kd=100)



def test_get_binding_site_affinity():
    pass  # TODO



def test_get_binding_sites():
    chem_data = ChemData(["A", "B"])
    chem_data.add_macromolecules(["M1", "M2"])

    with pytest.raises(Exception):
        chem_data.get_binding_sites("M9999")    # Unknown macromolecule

    assert chem_data.get_binding_sites("M1") == []
    assert chem_data.get_binding_sites("M2") == []

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=3)
    assert chem_data.get_binding_sites("M1") == [1]
    assert chem_data.get_binding_sites("M2") == []

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=5)
    assert chem_data.get_binding_sites("M1") == [1, 2]
    assert chem_data.get_binding_sites("M2") == []

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="B", Kd=11)
    assert chem_data.get_binding_sites("M1") == [1, 2]
    assert chem_data.get_binding_sites("M2") == [1]

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=3, ligand="B", Kd=102)
    assert chem_data.get_binding_sites("M1") == [1, 2]
    assert chem_data.get_binding_sites("M2") == [1, 3]



def test_get_binding_sites_and_ligands():
    chem_data = ChemData(["A", "B"])
    chem_data.add_macromolecules(["M1", "M2"])

    with pytest.raises(Exception):
        chem_data.get_binding_sites_and_ligands("M9999")    # Unknown macromolecule

    assert chem_data.get_binding_sites_and_ligands("M1") == {}
    assert chem_data.get_binding_sites_and_ligands("M2") == {}

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=3)
    assert chem_data.get_binding_sites_and_ligands("M1") == {1: "A"}
    assert chem_data.get_binding_sites_and_ligands("M2") == {}

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=5)
    assert chem_data.get_binding_sites_and_ligands("M1") == {1: "A", 2: "B"}
    assert chem_data.get_binding_sites_and_ligands("M2") == {}

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="B", Kd=11)
    assert chem_data.get_binding_sites_and_ligands("M1") == {1: "A", 2: "B"}
    assert chem_data.get_binding_sites_and_ligands("M2") == {1: "B"}

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=3, ligand="B", Kd=102)
    assert chem_data.get_binding_sites_and_ligands("M1") == {1: "A", 2: "B"}
    assert chem_data.get_binding_sites_and_ligands("M2") == {1: "B", 3: "B"}



def test_get_ligand_name():
    chem_data = ChemData(["A", "B"])
    chem_data.add_macromolecules(["M1", "M2"])

    with pytest.raises(Exception):
        chem_data.get_ligand_name(macromolecule="M9999", site_number=123)    # Unknown macromolecule

    with pytest.raises(Exception):
        chem_data.get_ligand_name(macromolecule="M1", site_number=1)    # M1 has no sites associated to it

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="A", Kd=3.1)

    with pytest.raises(Exception):
        chem_data.get_ligand_name(macromolecule="M1", site_number=1)    # No site number 1 on M1

    assert chem_data.get_ligand_name(macromolecule="M1", site_number=2) == "A"

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="B", Kd=6.)
    assert chem_data.get_ligand_name(macromolecule="M1", site_number=1) == "B"

    chem_data.add_chemical("C")
    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=5, ligand="C", Kd=11.4)
    assert chem_data.get_ligand_name(macromolecule="M2", site_number=5) == "C"



def test_reset_macromolecule():
    pass   # TODO


def test_clear_macromolecules():
    chem_data = ChemData()
    chem_data.add_macromolecules(["M1", "M2"])
    chem_data.clear_macromolecules()
    assert chem_data.macromolecules == []
    assert chem_data.binding_sites == {}





#############  PRIVATE METHODS  #############

def test__internal_reactions_data():
    pass   # TODO



def test__generate_generic_names():
    chem_data = ChemData()

    assert chem_data._generate_generic_names(1) == ["A"]
    assert chem_data._generate_generic_names(2) == ["A", "B"]
    assert chem_data._generate_generic_names(26) == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                                                     'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
    assert chem_data._generate_generic_names(27) == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                                                     'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
                                                     'Z2']
    assert chem_data._generate_generic_names(28) == ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
                                                     'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
                                                     'Z2', 'Z3']