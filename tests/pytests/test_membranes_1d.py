import pytest
import numpy as np
from life123.bio_sim_1d_experimental import Membranes1D
from life123 import ChemData



#########   TESTS OF INITIALIZATION, SETTING AND VIEWING    #########

def test_constructor():
    with pytest.raises(Exception):
        Membranes1D()                  # Missing required arguments

    with pytest.raises(Exception):
        Membranes1D(n_bins=3)          # Missing required arguments


    chem_data = ChemData(names="A")

    with pytest.raises(Exception):
        Membranes1D(n_bins="I'm not an integer", chem_data=chem_data)

    with pytest.raises(Exception):
        Membranes1D(n_bins=0, chem_data=chem_data)      # Must be at least 1



    bio = Membranes1D(n_bins=5, chem_data=chem_data)

    assert bio.n_bins == 5
    assert bio.n_species == 1
    assert bio.chem_data == chem_data
    assert bio.global_Dx == 1
    assert np.allclose(bio.system_time, 0)
    expected = np.zeros((1, 5), dtype=float)
    assert np.allclose(bio.system, expected)
    assert bio.membranes == []
    assert bio.permeability == {}


    # New test
    chem_data = ChemData(names=["A", "B", "C"])
    bio = Membranes1D(n_bins=15, chem_data=chem_data)
    assert bio.n_bins == 15
    assert bio.n_species == 3
    assert bio.chem_data == chem_data
    expected = np.zeros((3,15), dtype=float)
    assert np.allclose(bio.system, expected)
    assert bio.membranes == []
    assert bio.permeability == {}



def test_describe_state():
    pass



def test_uses_membranes():
    chem_data = ChemData(names=["A", "B"])
    bio = Membranes1D(n_bins=5, chem_data=chem_data)

    assert not bio.uses_membranes()
    bio.set_membranes(membranes=[(2, 4)])
    assert bio.uses_membranes()
    
    
    bio = Membranes1D(n_bins=123, chem_data=ChemData(labels="A"))
    assert not bio.uses_membranes()

    bio.set_membranes(membranes=[(0, 4)])
    assert bio.uses_membranes()

    bio.set_membranes(membranes=[])
    assert not bio.uses_membranes()



def test_set_membranes():
    bio = Membranes1D(n_bins=10, chem_data=ChemData("A")) # 10 bins, numbered 0 thru 9

    with pytest.raises(Exception):
        bio.set_membranes()     # Missing required arg

    with pytest.raises(Exception):
        bio.set_membranes(membranes=123)        # Wrong arg type

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[8])        # Element isn't a pair

    with pytest.raises(Exception):
        bio.set_membranes([1, 2, 3])
        
    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(1,2,3)])  # Element isn't a pair

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(1,"x")])  # Element isn't a pair of integers

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(8,3)])    # Value not in sorted order

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(-1,3)])   # Bad left value

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(2,11)])   # Bad right value

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(2,4) , (4,9)])    # Touching

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(2,7) , (3,9)])    # Overlapping

    bio.set_membranes(membranes=[(2,4) , (5,9)])
    assert bio.membranes == [(2,4) , (5,9)]

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(5,9), (2,4)])     # Not in sorted order

    bio.set_membranes(membranes=[(0,3) , (4,6), (8,10)])
    assert bio.membranes == [(0,3) , (4,6), (8,10)]
    assert bio.permeability == {}



def test_set_permeability():
    pass

def test_change_permeability():
    pass



def test_membranes_list():
    bio = Membranes1D(n_bins=40, chem_data=ChemData(labels="A"))

    assert bio.membranes_list() == []

    bio.set_membranes([ (8,10) ])
    assert bio.membranes_list() == [8,10]

    bio.set_membranes([ (1,2), (3,4) ])
    assert bio.membranes_list() == [1,2,3,4]



def test_membrane_on_left():
    bio = Membranes1D(n_bins=40, chem_data=ChemData(labels="A"))

    assert not bio.membrane_on_left(10)

    bio.set_membranes([ (6,10) ])

    assert not bio.membrane_on_left(5)
    assert bio.membrane_on_left(6)
    assert not bio.membrane_on_left(7)

    assert not bio.membrane_on_left(9)
    assert bio.membrane_on_left(10)
    assert not bio.membrane_on_left(11)

    bio.set_membranes([ (6,10), (11,18) ])

    assert not bio.membrane_on_left(9)
    assert bio.membrane_on_left(10)
    assert bio.membrane_on_left(11)
    assert not bio.membrane_on_left(12)
    assert not bio.membrane_on_left(17)
    assert bio.membrane_on_left(18)
    assert not bio.membrane_on_left(19)



def test_membrane_on_right():
    bio = Membranes1D(n_bins=40, chem_data=ChemData(labels="A"))

    assert not bio.membrane_on_right(10)

    bio.set_membranes([ (6,10) ])

    assert not bio.membrane_on_right(4)
    assert bio.membrane_on_right(5)
    assert not bio.membrane_on_right(6)

    assert not bio.membrane_on_right(8)
    assert bio.membrane_on_right(9)
    assert not bio.membrane_on_right(10)

    bio.set_membranes([ (6,10), (11,18) ])

    assert not bio.membrane_on_right(8)
    assert bio.membrane_on_right(9)
    assert bio.membrane_on_right(10)
    assert not bio.membrane_on_right(11)
    assert not bio.membrane_on_right(16)
    assert bio.membrane_on_right(17)
    assert not bio.membrane_on_right(18)
