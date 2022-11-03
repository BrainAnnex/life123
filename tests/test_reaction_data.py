import pytest
import numpy as np
from modules.reactions.reaction_data import ReactionData



def test_initialize():
    chem_data = ReactionData(names=['A', 'B', 'C'])
    assert chem_data.n_species == 3
    assert chem_data.get_all_names() == ['A', 'B', 'C']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2}

    with pytest.raises(Exception):
        ReactionData(names='this is not a list')   # Not a list/tuple

    chem_data = ReactionData(diffusion_rates=[0.15, 1.2])
    assert chem_data.n_species == 2
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2])

    with pytest.raises(Exception):
        ReactionData(diffusion_rates=[0.15, 1.2, "I'm not a number"])   # Bad value

    with pytest.raises(Exception):
        ReactionData(diffusion_rates=123.456)   # Not a list/tuple


    with pytest.raises(Exception):
        ReactionData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2])  # mismatch in count

    chem_data = ReactionData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2, 3.14])
    assert chem_data.n_species == 3
    assert chem_data.get_all_names() == ['A', 'B', 'C']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2}
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2, 3.14])



def test_number_of_chemicals():
    chem_data = ReactionData(names=['A', 'B', 'C'])
    assert chem_data.number_of_chemicals() == 3

    chem_data = ReactionData(names=[])
    assert chem_data.number_of_chemicals() == 0



def test_init_chemical_data():
    chem_data = ReactionData()
    chem_data.init_chemical_data(names=["A", "B", "C", "D", "E", "F"])

    assert chem_data.get_all_names() == ["A", "B", "C", "D", "E", "F"]
    assert chem_data.name_dict == {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5}

    with pytest.raises(Exception):
        chem_data = ReactionData()
        chem_data.init_chemical_data(names="Do I look like a list??")

    with pytest.raises(Exception):
        chem_data = ReactionData(names=['A', 'B', 'C'])
        chem_data.init_chemical_data(names=['X', 'Y', 'Z'])

    with pytest.raises(Exception):
        chem_data = ReactionData()
        chem_data.init_chemical_data(names=['A', 'B'], diffusion_rates=[1, 2, 3])

    chem_data = ReactionData()
    chem_data.init_chemical_data(names=['A', 'B', 'C'], diffusion_rates=[1, 2, 3])
    assert chem_data.number_of_chemicals() == 3
    assert chem_data.get_all_names() == ['A', 'B', 'C']
    # Verify that the name index also got created successfully
    assert chem_data.name_dict['A'] == 0
    assert chem_data.name_dict['B'] == 1
    assert chem_data.name_dict['C'] == 2



def test_get_name():
    chem_data = ReactionData(names=['A', 'B', 'C'])
    assert chem_data.get_name(0) == 'A'
    assert chem_data.get_name(1) == 'B'
    assert chem_data.get_name(2) == 'C'

    with pytest.raises(Exception):
        chem_data.get_name(3)               # Out of bounds
    with pytest.raises(Exception):
        chem_data.get_name(-1)              # Invalid argument value
    with pytest.raises(Exception):
        chem_data.get_name("some string")   # Invalid argument type
    with pytest.raises(Exception):
        chem_data.get_name(3.14)            # Invalid argument type



def test_get_index():
    chem_data = ReactionData(names=['A', 'B', 'C'])
    assert chem_data.get_index('A') == 0
    assert chem_data.get_index('B') == 1
    assert chem_data.get_index('C') == 2
    with pytest.raises(Exception):
        assert chem_data.get_index('X')     # Not found



def test_add_chemical():
    chem_data = ReactionData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2, 3.14])
    assert chem_data.n_species == 3
    assert chem_data.get_all_names() == ['A', 'B', 'C']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2}
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2, 3.14])

    chem_data.add_chemical(name="D", diffusion_rate=8)
    assert chem_data.n_species == 4
    assert chem_data.get_all_names() == ['A', 'B', 'C', 'D']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2, 'D': 3}
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2, 3.14, 8])
    with pytest.raises(Exception):
        chem_data.add_chemical(name="E", diffusion_rate="I'm not a number")

    with pytest.raises(Exception):
        chem_data.add_chemical(name="E", diffusion_rate=-666.)

    with pytest.raises(Exception):
        chem_data.add_chemical(name=666, diffusion_rate=25.)    # Wrong type for name

    chem_data.add_chemical(name="E", diffusion_rate=25.)
    assert chem_data.n_species == 5
    assert chem_data.get_all_names() == ['A', 'B', 'C', 'D', 'E']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4}
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2, 3.14, 8, 25])
