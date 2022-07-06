import pytest
import numpy as np
from modules.chemicals.chemicals import Chemicals as chem



def test_initialize():

    chem_data = chem(names=['A', 'B', 'C'])
    assert chem_data.n_species == 3
    assert chem_data.names == ['A', 'B', 'C']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2}

    with pytest.raises(Exception):
        chem(names='this is not a list')   # Not a list/tuple

    chem_data = chem(diffusion_rates=[0.15, 1.2])
    assert chem_data.n_species == 2
    assert np.allclose(chem_data.diffusion_rates, [0.15, 1.2])

    with pytest.raises(Exception):
        chem(diffusion_rates=[0.15, 1.2, "I'm not a number"])   # Bad value

    with pytest.raises(Exception):
        chem(diffusion_rates=123.456)   # Not a list/tuple


    with pytest.raises(Exception):
        chem(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2])  # mismatch in count

    chem_data = chem(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2, 3.14])
    assert chem_data.n_species == 3
    assert chem_data.names == ['A', 'B', 'C']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2}
    assert np.allclose(chem_data.diffusion_rates, [0.15, 1.2, 3.14])



def test_set_diffusion_rates():
    pass


def test_set_names():
    pass



def test_get_name():
    chem_data = chem(names=['A', 'B', 'C'])
    assert chem_data.get_name(0) == 'A'
    assert chem_data.get_name(1) == 'B'
    assert chem_data.get_name(2) == 'C'
    assert chem_data.get_name(3) is None

    with pytest.raises(Exception):
        chem_data.get_name(-1)              # Invalid argument value
    with pytest.raises(Exception):
        chem_data.get_name("some string")   # Invalid argument type
    with pytest.raises(Exception):
        chem_data.get_name(3.14)            # Invalid argument type


def test_get_index():
    chem_data = chem(names=['A', 'B', 'C'])
    assert chem_data.get_index('A') == 0
    assert chem_data.get_index('B') == 1
    assert chem_data.get_index('C') == 2
    assert chem_data.get_index('X') is None



def test_add_chemical():
    chem_data = chem(names=['A', 'B', 'C'])

    with pytest.raises(Exception):
        chem_data.add_chemical(name="D", diffusion_rate=0.8)    # The diffusion rates of the earlier elements weren't set

    chem_data = chem(diffusion_rates=[0.15, 1.2])
    with pytest.raises(Exception):
        chem_data.add_chemical(name="X", diffusion_rate=666.)    # The names of the earlier elements weren't set

    chem_data = chem(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2, 3.14])
    chem_data.add_chemical(name="D", diffusion_rate=8)
    assert chem_data.n_species == 4
    assert chem_data.names == ['A', 'B', 'C', 'D']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2, 'D': 3}
    assert np.allclose(chem_data.diffusion_rates, [0.15, 1.2, 3.14, 8])
    with pytest.raises(Exception):
        chem_data.add_chemical(name="E", diffusion_rate="I'm not a number")

    with pytest.raises(Exception):
        chem_data.add_chemical(name="E", diffusion_rate=-666.)

    with pytest.raises(Exception):
        chem_data.add_chemical(name=666, diffusion_rate=25.)    # Wrong type for name

    chem_data.add_chemical(name="E", diffusion_rate=25.)
    assert chem_data.n_species == 5
    assert chem_data.names == ['A', 'B', 'C', 'D', 'E']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}
    assert np.allclose(chem_data.diffusion_rates, [0.15, 1.2, 3.14, 8, 25])

