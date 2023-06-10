import pytest
import numpy as np
from src.modules.chemicals.chem_data import ChemData



def test_initialize():
    chem_data = ChemData()
    assert chem_data.n_species == 0

    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.n_species == 3
    assert chem_data.get_all_names() == ['A', 'B', 'C']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2}
    assert chem_data.get_all_diffusion_rates() == [None, None, None]

    with pytest.raises(Exception):
        ChemData(names='this is not a list')   # Not a list/tuple

    with pytest.raises(Exception):
        ChemData(names=[1, 2])   # The names aren't strings

    chem_data = ChemData(diffusion_rates=[0.15, 1.2])
    assert chem_data.n_species == 2
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2])
    assert chem_data.get_all_names() == ['Chemical 1', 'Chemical 2']

    with pytest.raises(Exception):
        ChemData(diffusion_rates=123.456)   # Not a list/tuple

    with pytest.raises(Exception):
        ChemData(diffusion_rates=[0.15, 1.2, "I'm not a number"])   # Bad value

    with pytest.raises(Exception):
        ChemData(diffusion_rates=[-6.66])   # Values cannot be negative

    with pytest.raises(Exception):
        ChemData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2])  # mismatch in count

    chem_data = ChemData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2, 3.14])
    assert chem_data.n_species == 3
    assert chem_data.get_all_names() == ['A', 'B', 'C']
    assert chem_data.name_dict == {'A': 0, 'B': 1, 'C': 2}
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2, 3.14])



def test_init_chemical_data():
    chem_data = ChemData()
    chem_data.init_chemical_data(names=["A", "B", "C", "D", "E", "F"])

    assert chem_data.get_all_names() == ["A", "B", "C", "D", "E", "F"]
    assert chem_data.name_dict == {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5}

    with pytest.raises(Exception):
        chem_data = ChemData()
        chem_data.init_chemical_data(names="Do I look like a list??")

    with pytest.raises(Exception):
        chem_data = ChemData(names=['A', 'B', 'C'])
        chem_data.init_chemical_data(names=['X', 'Y', 'Z'])

    with pytest.raises(Exception):
        chem_data = ChemData()
        chem_data.init_chemical_data(names=['A', 'B'], diffusion_rates=[1, 2, 3])

    chem_data = ChemData()
    chem_data.init_chemical_data(names=['A', 'B', 'C'], diffusion_rates=[1, 2, 3])
    assert chem_data.number_of_chemicals() == 3
    assert chem_data.get_all_names() == ['A', 'B', 'C']
    # Verify that the name index also got created successfully
    assert chem_data.name_dict['A'] == 0
    assert chem_data.name_dict['B'] == 1
    assert chem_data.name_dict['C'] == 2



def test_number_of_chemicals():
    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.number_of_chemicals() == 3

    chem_data = ChemData(names=[])
    assert chem_data.number_of_chemicals() == 0



def test_get_name():
    chem_data = ChemData(names=['A', 'B', 'C'])
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
    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.get_index('A') == 0
    assert chem_data.get_index('B') == 1
    assert chem_data.get_index('C') == 2
    with pytest.raises(Exception):
        assert chem_data.get_index('X')     # Not found



def test_add_chemical():
    chem_data = ChemData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2, 3.14])
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



def test_add_reaction():
    chem = ChemData(names=["A", "B", "C", "D", "E", "F"])
    chem.set_temp(None)

    # Reactants and the products can't be the same
    with pytest.raises(Exception):
        chem.add_reaction(reactants=["A"], products=["A"])
    with pytest.raises(Exception):
        chem.add_reaction(reactants=["A"], products=[("A")])
    with pytest.raises(Exception):
        chem.add_reaction(reactants=["A"], products=[(1, "A")])
    with pytest.raises(Exception):
        chem.add_reaction(reactants=["A"], products=[(1, "A", 1)])
    with pytest.raises(Exception):
        chem.add_reaction(reactants=[(2, "B")], products=[(2, "B")])
    with pytest.raises(Exception):
        chem.add_reaction(reactants=[(2, "B")], products=[(2, "B", 1)])
    with pytest.raises(Exception):
        chem.add_reaction(reactants=["A", "B"], products=["A", "B"])
    with pytest.raises(Exception):
        chem.add_reaction(reactants=["A", (3, "B")], products=["A", (3, "B")])
    with pytest.raises(Exception):
        chem.add_reaction(reactants=["A", "B"], products=["B", "A"])


    # Add the first (0-th) reaction
    chem.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert chem.number_of_reactions() == 1
    r = chem.get_reaction(0)
    assert np.allclose(r.kF , 3.)
    assert np.allclose(r.kR , 2.)
    assert np.allclose(r.K , 3./2.)
    assert r.reactants == [(1, 0, 1)]
    assert r.products == [(1, 1, 1)]
    assert r.Delta_H is None
    assert r.Delta_S is None
    assert r.Delta_G is None

    assert chem.get_reactants(0) == [(1, 0, 1)]
    assert chem.get_reactants_formula(0) == "A"

    assert chem.get_products(0) == [(1, 1, 1)]
    assert chem.get_products_formula(0) == "B"

    assert np.allclose(chem.get_forward_rate(0), 3.)
    assert np.allclose(chem.get_reverse_rate(0),  2.)


    # Another reaction (reaction 1)
    chem.add_reaction(reactants=[(2, "B")], products=[(5, "C")], forward_rate=9., reverse_rate=7.)

    assert chem.number_of_reactions() == 2

    r = chem.get_reaction(0)
    assert np.allclose(r.kF , 3.)
    assert np.allclose(r.kR , 2.)
    assert np.allclose(r.K , 3./2.)
    assert r.reactants == [(1, 0, 1)]
    assert r.products == [(1, 1, 1)]
    assert r.Delta_H is None
    assert r.Delta_S is None
    assert r.Delta_G is None

    r = chem.get_reaction(1)
    assert np.allclose(r.kF , 9.)
    assert np.allclose(r.kR , 7.)
    assert np.allclose(r.K , 9./7.)
    assert r.reactants == [(2, 1, 1)]
    assert r.products == [(5, 2, 1)]
    assert r.Delta_H is None
    assert r.Delta_S is None
    assert r.Delta_G is None


    # Add another reaction (reaction index 2).  This time, first set the temperature
    chem.temp = 200

    chem.add_reaction(reactants=[(2, "D", 3)], products=[(1, "C", 2)], forward_rate=11., reverse_rate=13.)
    assert chem.number_of_reactions() == 3

    r = chem.get_reaction(2)
    assert np.allclose(r.kF , 11.)
    assert np.allclose(r.kR , 13.)
    assert np.allclose(r.K , 11./13.)
    assert r.reactants == [(2, 3, 3)]
    assert r.products == [(1, 2, 2)]
    assert r.Delta_H is None
    assert r.Delta_S is None
    assert np.allclose(r.Delta_G , 277.7928942715384)   # - RT log(K)


    # Add a multi-term reaction (reaction index 3)
    chem.add_reaction(reactants=["A", (2, "B")], products=[(3, "C", 2), "D"], forward_rate=5., reverse_rate=1.)
    assert chem.number_of_reactions() == 4

    r = chem.get_reaction(3)
    assert np.allclose(r.kF , 5.)
    assert np.allclose(r.kR , 1.)
    assert np.allclose(r.K , 5./1.)
    assert r.reactants == [(1, 0, 1), (2, 1, 1)]
    assert r.products == [(3, 2, 2), (1, 3, 1)]
    assert r.Delta_H is None
    assert r.Delta_S is None
    assert np.allclose(r.Delta_G , -2676.321364705849)   # - RT log(K)

    # Check the descriptions we has so far
    rxn_info = chem.multiple_reactions_describe()
    assert rxn_info[0] == '0: A <-> B  (kF = 3 / kR = 2 / K = 1.5) | 1st order in all reactants & products'
    assert rxn_info[1] == '1: 2 B <-> 5 C  (kF = 9 / kR = 7 / K = 1.28571) | 1st order in all reactants & products'
    assert rxn_info[2] == '2: 2 D <-> C  (kF = 11 / kR = 13 / Delta_G = 277.793 / K = 0.846154) | 3-th order in reactant D | 2-th order in product C'
    assert rxn_info[3] == '3: A + 2 B <-> 3 C + D  (kF = 5 / kR = 1 / Delta_G = -2,676.32 / K = 5) | 2-th order in product C'


    # Add another reaction (reaction index 4), this time with thermodynamic data;
    # the reverse reaction rate will get computed from the thermodynamic data
    chem.add_reaction(reactants=["A"], products=[(2, "B")], forward_rate=10.,
                      delta_H= 5., delta_S= 0.4)
    assert chem.number_of_reactions() == 5

    r = chem.get_reaction(4)
    assert r.reactants == [(1, 0, 1)]
    assert r.products == [(2, 1, 1)]
    assert r.Delta_H == 5.
    assert r.Delta_S == 0.4
    assert np.allclose(r.Delta_G , -75.0)         # 5 - 200 * 0.4
    assert np.allclose(r.K , 1.0461347154679432)  # exp(75/(8.3144598 * 200))
    assert np.allclose(r.kF , 10.)
    assert np.allclose(r.kR , 9.558998331803693) # 10. / 1.0461347154679432
    

    
def test_clear_reactions():
    pass



# TO DESCRIBE THE DATA

def test_describe_reactions():
    pass



def test_get_chemicals_in_reaction():
    chem = ChemData(names=["A", "B"])

    with pytest.raises(Exception):
        chem.get_chemicals_in_reaction(0)   # There are no reactions defined yet

    chem.add_reaction(reactants="A", products="B")  # Reaction 0 : A <-> B
    assert chem.get_chemicals_in_reaction(0) == {0, 1}

    with pytest.raises(Exception):
        chem.get_chemicals_in_reaction(1)   # There is no reaction 1

    chem.add_chemical("C")

    chem.add_reaction(reactants=["B"], products=[(2, "C")])  # Reaction 1 : B <-> 2C
    assert chem.get_chemicals_in_reaction(0) == {0, 1}
    assert chem.get_chemicals_in_reaction(1) == {1, 2}

    chem.add_reaction(reactants=["A"], products=["C"])      # Reaction 2 : A <-> C
    assert chem.get_chemicals_in_reaction(0) == {0, 1}
    assert chem.get_chemicals_in_reaction(1) == {1, 2}
    assert chem.get_chemicals_in_reaction(2) == {0, 2}

    chem.add_chemical("D")
    chem.add_reaction(reactants=["A", "B"], products="D")    # Reaction 3 : A + B <-> D
    assert chem.get_chemicals_in_reaction(0) == {0, 1}
    assert chem.get_chemicals_in_reaction(1) == {1, 2}
    assert chem.get_chemicals_in_reaction(2) == {0, 2}
    assert chem.get_chemicals_in_reaction(3) == {0, 1, 3}




###  SUPPORT FOR CREATION OF NETWORK DIAGRAMS  ###

def test_create_graph_network_data():
    chem = ChemData(names=["A", "B"])
    chem.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    network_data = chem.create_graph_network_data()
    #print(network_data)
    expected = [{'id': 2, 'label': 'Reaction', 'name': 'RXN', 'kF': '3', 'kR': '2', 'Delta_G': '-1,005.13','K': '1.5'},
                {'id': 1, 'label': 'Chemical', 'name': 'B', 'diff_rate': None, 'stoich': 1, 'rxn_order': 1},
                {'id': 3, 'name': 'produces', 'source': 2, 'target': 1},
                {'id': 0, 'label': 'Chemical', 'name': 'A', 'diff_rate': None, 'stoich': 1, 'rxn_order': 1},
                {'id': 4, 'name': 'reacts', 'source': 0, 'target': 2}
                ]
    assert network_data == expected
