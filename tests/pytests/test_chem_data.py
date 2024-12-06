import pytest
import numpy as np
import pandas as pd
from life123 import ChemData
from tests.utilities.comparisons import *


#############  ChemCore  #############

def test_number_of_chemicals():
    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.number_of_chemicals() == 3

    chem_data = ChemData(names=[])
    assert chem_data.number_of_chemicals() == 0



def test_assert_valid_species_index():
    chem_data = ChemData(names=['A', 'B', 'C'])
    chem_data.assert_valid_species_index(0)
    chem_data.assert_valid_species_index(1)
    chem_data.assert_valid_species_index(2)

    with pytest.raises(Exception):
        chem_data.assert_valid_species_index(3)     # Too large

    with pytest.raises(Exception):
        chem_data.assert_valid_species_index(-1)    # Too small

    with pytest.raises(Exception):
        chem_data.assert_valid_species_index("2")   # Not an int



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

    result = chem_data.add_chemical("A")
    assert result == 0
    assert chem_data.number_of_chemicals() == 1
    assert chem_data.chemical_data == [{"name": "A", "label": "A"}]
    assert chem_data.label_dict == {"A": 0}

    with pytest.raises(Exception):
        chem_data.add_chemical("A") # Duplicate!

    result = chem_data.add_chemical("B", note="some note")
    assert result == 1
    assert chem_data.number_of_chemicals() == 2
    assert chem_data.chemical_data == [ {"name": "A", "label": "A"},
                                        {"name": "B", "label": "B", "note": "some note"}]
    assert chem_data.label_dict == {"A": 0, "B": 1}

    result = chem_data.add_chemical("C")
    assert result == 2
    assert chem_data.number_of_chemicals() == 3
    assert chem_data.chemical_data == [{"name": "A", "label": "A"},
                                       {"name": "B", "label": "B", "note": "some note"},
                                       {"name": "C", "label": "C"}]
    assert chem_data.label_dict == {"A": 0, "B": 1, "C": 2}


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

    with pytest.raises(Exception):
        chem_data.add_chemical(name=123)    # Name is not a string

    with pytest.raises(Exception):
        chem_data.add_chemical(name="")




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
    assert chem_data.get_diffusion_rate(species_index=0) == 10
    with pytest.raises(Exception):
        assert chem_data.get_diffusion_rate(name="A", species_index=0)

    assert chem_data.get_diffusion_rate(name="C") == 12
    assert chem_data.get_diffusion_rate(species_index=2) == 12

    with pytest.raises(Exception):
        chem_data.get_diffusion_rate(name="X")          # X doesn't exist

    with pytest.raises(Exception):
        chem_data.get_diffusion_rate(species_index=3)   # Index 3 doesn't exist

    chem_data.add_chemical(name="Z")                    # This will be given index 3
    assert chem_data.get_diffusion_rate(name="Z") is None       # No diffusion value assigned
    assert chem_data.get_diffusion_rate(species_index=3) is None

    chem_data.set_diffusion_rate(label="Z", diff_rate=8)
    assert chem_data.get_diffusion_rate(name="Z") == 8
    assert chem_data.get_diffusion_rate(species_index=3) == 8



def test_get_all_diffusion_rates():
    pass

def test_missing_diffusion_rate():
    pass





#############  AllReactions  #############

def test_number_of_reactions():
    chem = ChemData(names=["A", "B", "C"])

    assert chem.number_of_reactions() == 0

    chem.add_reaction(reactants=["A"], products=["B"])
    assert chem.number_of_reactions() == 1

    chem.add_reaction(reactants=["A", (2, "B")], products=["C"])
    assert chem.number_of_reactions() == 2



def test_assert_valid_rxn_index():
    pass   # TODO

def test_get_reaction():
    pass   # TODO

def test_get_reactants():
    pass   # TODO

def test_get_reactants_formula():
    pass   # TODO

def test_get_products():
    pass   # TODO

def test_get_products_formula():
    pass   # TODO

def test_get_forward_rate():
    pass   # TODO

def test_get_reverse_rate():
    pass   # TODO



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



def test_get_chemicals_indexes_in_reaction():
    chem = ChemData(names=["A", "B"])

    with pytest.raises(Exception):
        chem.get_chemicals_indexes_in_reaction(0)   # There are no reactions defined yet

    chem.add_reaction(reactants="A", products="B")  # Reaction 0 : A <-> B
    assert chem.get_chemicals_indexes_in_reaction(0) == [0, 1]

    with pytest.raises(Exception):
        chem.get_chemicals_indexes_in_reaction(1)   # There is no reaction 1

    chem.add_chemical("C")

    chem.add_reaction(reactants=["B"], products=[(2, "C")])  # Reaction 1 : B <-> 2C
    assert chem.get_chemicals_indexes_in_reaction(0) == [0, 1]
    assert chem.get_chemicals_indexes_in_reaction(1) == [1, 2]

    chem.add_reaction(reactants=["A"], products=["C"])      # Reaction 2 : A <-> C
    assert chem.get_chemicals_indexes_in_reaction(0) == [0, 1]
    assert chem.get_chemicals_indexes_in_reaction(1) == [1, 2]
    assert chem.get_chemicals_indexes_in_reaction(2) == [0, 2]

    chem.add_chemical("D")
    chem.add_reaction(reactants=["A", "B"], products="D")    # Reaction 3 : A + B <-> D
    assert chem.get_chemicals_indexes_in_reaction(0) == [0, 1]
    assert chem.get_chemicals_indexes_in_reaction(1) == [1, 2]
    assert chem.get_chemicals_indexes_in_reaction(2) == [0, 2]
    assert chem.get_chemicals_indexes_in_reaction(3) == [0, 1, 3]



def test_get_reactions_participating_in():
    pass   # TODO



def test_set_temp():
    chem_data = ChemData()
    chem_data.set_temp(123)
    assert chem_data.temp == 123



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
        chem.add_reaction(reactants=[(2, "B")], products=[(2, "B", 2)])
    with pytest.raises(Exception):
        chem.add_reaction(reactants=[(2, "B")], products=[(2, "B", 2)])
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
    assert r.reactants == [(1,"A", 1)]
    assert r.products == [(1, "B", 1)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert r.delta_G is None

    assert chem.get_reactants(0) == [(1, "A", 1)]
    assert chem.get_reactants_formula(0) == "A"

    assert chem.get_products(0) == [(1, "B", 1)]
    assert chem.get_products_formula(0) == "B"

    assert np.allclose(chem.get_forward_rate(0), 3.)
    assert np.allclose(chem.get_reverse_rate(0),  2.)

    assert chem.active_chemicals == {"A", "B"}
    assert chem.active_enzymes == set()


    # Another reaction (reaction 1)
    chem.add_reaction(reactants=[(2, "B", 1)], products=[(5, "C", 1)], forward_rate=9., reverse_rate=7.)

    assert chem.number_of_reactions() == 2

    r = chem.get_reaction(0)
    assert np.allclose(r.kF , 3.)
    assert np.allclose(r.kR , 2.)
    assert np.allclose(r.K , 3./2.)
    assert r.reactants == [(1, "A", 1)]
    assert r.products == [(1, "B", 1)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert r.delta_G is None

    r = chem.get_reaction(1)
    assert np.allclose(r.kF , 9.)
    assert np.allclose(r.kR , 7.)
    assert np.allclose(r.K , 9./7.)
    assert r.reactants == [(2, "B", 1)]
    assert r.products == [(5, "C", 1)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert r.delta_G is None

    assert chem.active_chemicals == {"A", "B", "C"}
    assert chem.active_enzymes == set()


    # Add another reaction (reaction index 2).  This time, first set the temperature
    chem.temp = 200

    chem.add_reaction(reactants=[(2, "D", 3)], products=[(1, "C", 2)],
                      forward_rate=11., reverse_rate=13.)
    assert chem.number_of_reactions() == 3

    r = chem.get_reaction(2)
    assert np.allclose(r.kF , 11.)
    assert np.allclose(r.kR , 13.)
    assert np.allclose(r.K , 11./13.)
    assert r.reactants == [(2, "D", 3)]
    assert r.products == [(1, "C", 2)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert np.allclose(r.delta_G, 277.7928942715384)   # - RT log(K)

    assert chem.active_chemicals == {"A", "B", "C", "D"}
    assert chem.active_enzymes == set()


    # Add a multi-term reaction (reaction index 3)
    chem.add_reaction(reactants=["A", (2, "B", 1)], products=[(3, "C", 2), "D"],
                      forward_rate=5., reverse_rate=1.)
    assert chem.number_of_reactions() == 4

    r = chem.get_reaction(3)
    assert np.allclose(r.kF , 5.)
    assert np.allclose(r.kR , 1.)
    assert np.allclose(r.K , 5./1.)
    assert r.reactants == [(1, "A", 1), (2, "B", 1)]
    assert r.products == [(3, "C", 2), (1, "D", 1)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert np.allclose(r.delta_G, -2676.321364705849)   # - RT log(K)

    assert chem.active_chemicals == {"A", "B", "C", "D"}
    assert chem.active_enzymes == set()


    # Check the descriptions we has so far
    rxn_info = chem.multiple_reactions_describe()
    assert rxn_info[0] == '0: A <-> B  (kF = 3 / kR = 2 / K = 1.5) | 1st order in all reactants & products'
    assert rxn_info[1] == '1: 2 B <-> 5 C  (kF = 9 / kR = 7 / K = 1.2857) | 1st order in all reactants & products'
    assert rxn_info[2] == '2: 2 D <-> C  (kF = 11 / kR = 13 / delta_G = 277.79 / K = 0.84615) | 3-th order in reactant D | 2-th order in product C'
    assert rxn_info[3] == '3: A + 2 B <-> 3 C + D  (kF = 5 / kR = 1 / delta_G = -2,676.3 / K = 5) | 2-th order in product C'


    # Add another reaction (reaction index 4), this time with thermodynamic data;
    # the reverse reaction rate will get computed from the thermodynamic data
    chem.add_reaction(reactants=["A"], products=[(2, "B", 1)], forward_rate=10.,
                      delta_H= 5., delta_S= 0.4)
    assert chem.number_of_reactions() == 5

    r = chem.get_reaction(4)
    assert r.reactants == [(1, "A", 1)]
    assert r.products == [(2, "B", 1)]
    assert r.delta_H == 5.
    assert r.delta_S == 0.4
    assert np.allclose(r.delta_G, -75.0)         # 5 - 200 * 0.4
    assert np.allclose(r.K , 1.0461347154679432)  # exp(75/(8.3144598 * 200))
    assert np.allclose(r.kF , 10.)
    assert np.allclose(r.kR , 9.558998331803693) # 10. / 1.0461347154679432

    assert chem.active_chemicals == {"A", "B", "C", "D"}
    assert chem.active_enzymes == set()


    # Add another reaction (reaction index 5), this time involving a catalytic term
    chem.add_reaction(reactants=[(3, "A", 1), "D"], products=["D", (2, "B", 1)])
    assert chem.number_of_reactions() == 6
    assert chem.active_chemicals == {"A", "B", "C", "D"}
    assert chem.active_enzymes == {"D"}         # Notice that chemical "D is an enzyme in this reaction,
                                                # but an active participant in others

    # Add another reaction (reaction index 6), again involving a catalytic term
    chem.add_reaction(reactants=["F", (2, "C", 1), (3, "A", 1)], products=["D", "F", (2, "D", 1)])
    assert chem.number_of_reactions() == 7
    assert chem.active_chemicals == {"A", "B", "C", "D"}
    assert chem.active_enzymes == {"D", "F"}

    # Add another reaction (reaction index 7), where "F" (reaction index 5) is NOT a catalyst
    chem.add_reaction(reactants=[(2, "B", 1), "F"], products=[(3, "A", 1)])
    assert chem.number_of_reactions() == 8
    assert chem.active_chemicals == {"A", "B", "C", "D", "F"}
    assert chem.active_enzymes == {"D", "F"}



def test_clear_reactions_data():
    pass   # TODO




#################  TO DESCRIBE THE DATA  #################

def test_multiple_reactions_describe():
    pass   # TODO

def test_single_reaction_describe():
    pass   # TODO



def test_names_of_active_chemicals():
    chem = ChemData(names=['A', 'B', 'C', 'X', 'Y'])
    assert chem.names_of_active_chemicals() == set()    # No reactions yet

    chem.add_reaction(reactants="A", products="B")
    assert chem.names_of_active_chemicals() == {"A", "B"}

    chem.add_reaction(reactants=["B", "X"], products=["C", "X"])
    assert chem.names_of_active_chemicals() == {"A", "B", "C"}  # "X" is an enzyme

    chem.add_reaction(reactants="X", products="Y")
    assert chem.names_of_active_chemicals() == {"A", "B", "C", "X", "Y"}
    # "X" is now involved in some reactions in a non-enzymatic role

    chem.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"])
    assert chem.names_of_active_chemicals() == {"A", "B", "C", "X", "Y"}  # "Z" is an enzyme



def test_indexes_of_active_chemicals():
    chem = ChemData(names=['Y', 'X', 'C', 'B', 'A'])
    assert chem.indexes_of_active_chemicals() == []                 # No reactions yet

    chem.add_reaction(reactants="A", products="B")
    assert chem.indexes_of_active_chemicals() == [3, 4]             # ["A", "B"]

    chem.add_reaction(reactants=["B", "X"], products=["C", "X"])
    assert chem.indexes_of_active_chemicals() == [2, 3, 4]          # ["A", "B", "C"]

    chem.add_reaction(reactants="X", products="Y")
    assert chem.indexes_of_active_chemicals() == [0, 1, 2, 3, 4]    # All

    chem.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"])
    assert chem.indexes_of_active_chemicals() == [0, 1, 2, 3, 4]    # All



def test_names_of_enzymes():
    chem = ChemData(names=['A', 'B', 'C', 'X', 'Y'])
    assert chem.names_of_enzymes() == set()    # No enzymes yet

    chem.add_reaction(reactants="A", products="B")
    assert chem.names_of_enzymes() == set()     # No enzymes yet

    chem.add_reaction(reactants=["B", "X"], products=["C", "X"])
    assert chem.names_of_enzymes() == {"X"}     # "X" is an enzyme

    chem.add_reaction(reactants="X", products="Y")
    assert chem.names_of_enzymes() == {"X"}      # "X" is an enzyme in SOME reaction

    chem.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"])
    assert chem.names_of_enzymes() == {"X", "Z"}





#############  ChemData  #############

def test_initialize():
    chem_data = ChemData()
    assert chem_data.number_of_chemicals() == 0

    chem_data = ChemData(names=['A', 'B', 'C'])
    assert chem_data.number_of_chemicals() == 3
    assert chem_data.get_all_labels() == ['A', 'B', 'C']
    assert chem_data.label_dict == {'A': 0, 'B': 1, 'C': 2}
    assert chem_data.diffusion_rates == {}

    with pytest.raises(Exception):
        ChemData(names=123)     # Not a list/tuple/str

    with pytest.raises(Exception):
        ChemData(names=[1, 2])  # The names aren't strings

    chem_data = ChemData(diffusion_rates=[0.15, 1.2])
    assert chem_data.number_of_chemicals() == 2
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2])
    assert chem_data.get_all_labels() == ['Chemical 1', 'Chemical 2']

    with pytest.raises(Exception):
        ChemData(diffusion_rates=123.456)   # Not a list/tuple

    with pytest.raises(Exception):
        ChemData(diffusion_rates=[0.15, 1.2, "I'm not a number"])   # Bad value

    with pytest.raises(Exception):
        ChemData(diffusion_rates=[-6.66])   # Values cannot be negative

    with pytest.raises(Exception):
        ChemData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2])  # mismatch in count

    chem_data = ChemData(names=['A', 'B', 'C'], diffusion_rates=[0.15, 1.2, 3.14])
    assert chem_data.number_of_chemicals() == 3
    assert chem_data.get_all_labels() == ['A', 'B', 'C']
    assert chem_data.label_dict == {'A': 0, 'B': 1, 'C': 2}
    assert np.allclose(chem_data.get_all_diffusion_rates(), [0.15, 1.2, 3.14])

    assert np.allclose(chem_data.temp, 298.15)      # The default temperature



def test_init_chemical_data():
    chem_data = ChemData()
    chem_data.init_chemical_data(names=["A", "B", "C", "D", "E", "F"])

    assert chem_data.get_all_labels() == ["A", "B", "C", "D", "E", "F"]
    assert chem_data.label_dict == {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5}

    with pytest.raises(Exception):
        chem_data = ChemData()
        chem_data.init_chemical_data(names=123) # Not a list/tuple/str

    with pytest.raises(Exception):
        chem_data = ChemData(names=['A', 'B', 'C'])
        chem_data.init_chemical_data(names=['X', 'Y', 'Z'])

    with pytest.raises(Exception):
        chem_data = ChemData()
        chem_data.init_chemical_data(names=['A', 'B'], diffusion_rates=[1, 2, 3])

    chem_data = ChemData()
    chem_data.init_chemical_data(names=['A', 'B', 'C'], diffusion_rates=[1, 2, 3])
    assert chem_data.number_of_chemicals() == 3
    assert chem_data.get_all_labels() == ['A', 'B', 'C']
    # Verify that the name index also got created successfully
    assert chem_data.label_dict['A'] == 0
    assert chem_data.label_dict['B'] == 1
    assert chem_data.label_dict['C'] == 2




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





#############  FOR CREATION OF NETWORK DIAGRAMS  #############

def test_prepare_graph_network():
    # Set up an A <-> B reaction
    chem = ChemData(names=["A", "B"])
    chem.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    graph_data = chem.prepare_graph_network()

    expected_structure = [{'name': 'A', 'diff_rate': None, 'id': 'C-0', 'labels': ['Chemical']},
                          {'name': 'B', 'diff_rate': None, 'id': 'C-1', 'labels': ['Chemical']},
                          {'name': 'RXN', 'kF': '3', 'kR': '2', 'delta_G': '-1,005.13', 'K': '1.5', 'id': 'R-0', 'labels': ['Reaction']},

                          {'name': 'produces', 'source': 'R-0', 'target': 'C-1', 'id': 'edge-1', 'stoich': 1, 'rxn_order': 1},
                          {'name': 'reacts', 'source': 'C-0', 'target': 'R-0', 'id': 'edge-2', 'stoich': 1, 'rxn_order': 1}
                          ]

    assert compare_recordsets(graph_data["structure"], expected_structure)
    assert graph_data["color_mapping"] == {'Chemical': '#8DCC92', 'Reaction': '#D9C8AD'}
    assert graph_data["caption_mapping"] == {'Chemical': 'name', 'Reaction': 'name'}







#############  PRIVATE METHODS  #############

def test__internal_reactions_data():
    pass   # TODO