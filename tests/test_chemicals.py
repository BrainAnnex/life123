import pytest
import numpy as np
from modules.chemicals.chemicals import Chemicals as chem



def test_initialize():
    chem_data = chem(n_species=4)
    assert chem_data.n_species == 4

    chem_data = chem(n_species=3, names=['A', 'B', 'C'])
    assert chem_data.n_species == 3
    assert chem_data.names == ['A', 'B', 'C']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2}

    chem_data = chem(names=['A', 'B', 'C'])
    assert chem_data.n_species == 3     # This number was inferred
    assert chem_data.names == ['A', 'B', 'C']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2}

    with pytest.raises(Exception):
        chem(n_species=100, names=['A', 'B', 'C'])   # mismatch

    chem_data = chem(n_species=2, diffusion_rates=[0.15, 1.2])
    assert chem_data.n_species == 2
    assert np.allclose(chem_data.diffusion_rates, [0.15, 1.2])

    chem_data = chem(diffusion_rates=[0.15, 1.2])
    assert chem_data.n_species == 2    # This number was inferred
    assert np.allclose(chem_data.diffusion_rates, [0.15, 1.2])

    with pytest.raises(Exception):
        chem(n_species=100, diffusion_rates=[0.15, 1.2])   # mismatch



def test_set_diffusion_rates():
    pass


def test_set_names():
    pass



def test_get_name():
    chem_data = chem(n_species=3, names=['A', 'B', 'C'])
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
    chem_data = chem(n_species=3, names=['A', 'B', 'C'])
    assert chem_data.get_index('A') == 0
    assert chem_data.get_index('B') == 1
    assert chem_data.get_index('C') == 2
    assert chem_data.get_index('X') is None
