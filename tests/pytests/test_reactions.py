import pytest
import numpy as np
from life123 import ChemData
from life123.reactions import Reactions
from tests.utilities.comparisons import *



def test_constructor():
    chem_data = ChemData()
    rxns = Reactions(chem_data=chem_data)

    assert rxns.reaction_list == []
    assert np.allclose(rxns.temp, 298.15)       # The default temperature

    assert rxns.active_chemicals == set()       # Empty set
    assert rxns.active_enzymes == set()         # Empty set



def test_number_of_reactions():
    chem_data = ChemData(names=["A", "B", "C"])
    rxns = Reactions(chem_data=chem_data)

    assert rxns.number_of_reactions() == 0

    rxns.add_reaction(reactants="A", products="B")
    assert rxns.number_of_reactions() == 1

    rxns.add_reaction(reactants=["A", (2, "B")], products="C")
    assert rxns.number_of_reactions() == 2



def test_active_reaction_indices():
    chem_data = ChemData(names=["A", "B", "C"])
    rxns = Reactions(chem_data)

    assert rxns.active_reaction_indices() == []

    rxns.add_reaction(reactants="A", products="B")
    assert rxns.active_reaction_indices() == [0]

    rxns.add_reaction(reactants=["A", (2, "B")], products="C")
    assert rxns.active_reaction_indices() == [0, 1]



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
    chem_data = ChemData(names=["A", "B"])
    rxns = Reactions(chem_data)

    with pytest.raises(Exception):
        rxns.get_chemicals_in_reaction(0)   # There are no reactions defined yet

    rxns.add_reaction(reactants="A", products="B")  # Reaction 0 : A <-> B
    assert rxns.get_chemicals_in_reaction(0) == {0, 1}

    with pytest.raises(Exception):
        rxns.get_chemicals_in_reaction(1)   # There is no reaction 1

    chem_data.add_chemical("C")

    rxns.add_reaction(reactants=["B"], products=[(2, "C")])  # Reaction 1 : B <-> 2C
    assert rxns.get_chemicals_in_reaction(0) == {0, 1}
    assert rxns.get_chemicals_in_reaction(1) == {1, 2}

    rxns.add_reaction(reactants=["A"], products=["C"])      # Reaction 2 : A <-> C
    assert rxns.get_chemicals_in_reaction(0) == {0, 1}
    assert rxns.get_chemicals_in_reaction(1) == {1, 2}
    assert rxns.get_chemicals_in_reaction(2) == {0, 2}

    chem_data.add_chemical("D")
    rxns.add_reaction(reactants=["A", "B"], products="D")    # Reaction 3 : A + B <-> D
    assert rxns.get_chemicals_in_reaction(0) == {0, 1}
    assert rxns.get_chemicals_in_reaction(1) == {1, 2}
    assert rxns.get_chemicals_in_reaction(2) == {0, 2}
    assert rxns.get_chemicals_in_reaction(3) == {0, 1, 3}



def test_get_chemicals_indexes_in_reaction():
    chem_data = ChemData(names=["A", "B"])
    rxns = Reactions(chem_data)

    with pytest.raises(Exception):
        rxns.get_chemicals_indexes_in_reaction(0)   # There are no reactions defined yet

    rxns.add_reaction(reactants="A", products="B")  # Reaction 0 : A <-> B
    assert rxns.get_chemicals_indexes_in_reaction(0) == [0, 1]

    with pytest.raises(Exception):
        rxns.get_chemicals_indexes_in_reaction(1)   # There is no reaction 1

    chem_data.add_chemical("C")

    rxns.add_reaction(reactants=["B"], products=[(2, "C")])  # Reaction 1 : B <-> 2C
    assert rxns.get_chemicals_indexes_in_reaction(0) == [0, 1]
    assert rxns.get_chemicals_indexes_in_reaction(1) == [1, 2]

    rxns.add_reaction(reactants=["A"], products=["C"])      # Reaction 2 : A <-> C
    assert rxns.get_chemicals_indexes_in_reaction(0) == [0, 1]
    assert rxns.get_chemicals_indexes_in_reaction(1) == [1, 2]
    assert rxns.get_chemicals_indexes_in_reaction(2) == [0, 2]

    chem_data.add_chemical("D")
    rxns.add_reaction(reactants=["A", "B"], products="D")    # Reaction 3 : A + B <-> D
    assert rxns.get_chemicals_indexes_in_reaction(0) == [0, 1]
    assert rxns.get_chemicals_indexes_in_reaction(1) == [1, 2]
    assert rxns.get_chemicals_indexes_in_reaction(2) == [0, 2]
    assert rxns.get_chemicals_indexes_in_reaction(3) == [0, 1, 3]



def test_get_reactions_participating_in():
    pass   # TODO



def test_set_temp():
    rxns = Reactions(ChemData())
    rxns.set_temp(123)
    assert rxns.temp == 123



def test_add_reaction():
    chem_data = ChemData(names=["A", "B", "C", "D", "E", "F"])
    rxns = Reactions(chem_data)
    rxns.temp = None


    # Reactants and the products can't be the same
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A"], products=["A"])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A"], products=[("A")])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A"], products=[(1, "A")])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A"], products=[(1, "A", 1)])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=[(2, "B")], products=[(2, "B", 2)])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=[(2, "B")], products=[(2, "B", 2)])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A", "B"], products=["A", "B"])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A", (3, "B")], products=["A", (3, "B")])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A", "B"], products=["B", "A"])


    # Add the first (0-th) reaction
    rxns.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert rxns.number_of_reactions() == 1
    r = rxns.get_reaction(0)
    assert np.allclose(r.kF , 3.)
    assert np.allclose(r.kR , 2.)
    assert np.allclose(r.K , 3./2.)
    assert r.reactants == [(1,"A", 1)]
    assert r.products == [(1, "B", 1)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert r.delta_G is None

    assert rxns.get_reactants(0) == [(1, "A", 1)]
    assert rxns.get_reactants_formula(0) == "A"

    assert rxns.get_products(0) == [(1, "B", 1)]
    assert rxns.get_products_formula(0) == "B"

    assert np.allclose(rxns.get_forward_rate(0), 3.)
    assert np.allclose(rxns.get_reverse_rate(0),  2.)

    assert rxns.active_chemicals == {"A", "B"}
    assert rxns.active_enzymes == set()


    # Another reaction (reaction 1)
    rxns.add_reaction(reactants=[(2, "B", 1)], products=[(5, "C", 1)], forward_rate=9., reverse_rate=7.)

    assert rxns.number_of_reactions() == 2

    r = rxns.get_reaction(0)
    assert np.allclose(r.kF , 3.)
    assert np.allclose(r.kR , 2.)
    assert np.allclose(r.K , 3./2.)
    assert r.reactants == [(1, "A", 1)]
    assert r.products == [(1, "B", 1)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert r.delta_G is None

    r = rxns.get_reaction(1)
    assert np.allclose(r.kF , 9.)
    assert np.allclose(r.kR , 7.)
    assert np.allclose(r.K , 9./7.)
    assert r.reactants == [(2, "B", 1)]
    assert r.products == [(5, "C", 1)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert r.delta_G is None

    assert rxns.active_chemicals == {"A", "B", "C"}
    assert rxns.active_enzymes == set()


    # Add another reaction (reaction index 2).  This time, first set the temperature
    rxns.temp = 200

    rxns.add_reaction(reactants=[(2, "D", 3)], products=[(1, "C", 2)],
                      forward_rate=11., reverse_rate=13.)
    assert rxns.number_of_reactions() == 3

    r = rxns.get_reaction(2)
    assert np.allclose(r.kF , 11.)
    assert np.allclose(r.kR , 13.)
    assert np.allclose(r.K , 11./13.)
    assert r.reactants == [(2, "D", 3)]
    assert r.products == [(1, "C", 2)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert np.allclose(r.delta_G, 277.7928942715384)   # - RT log(K)

    assert rxns.active_chemicals == {"A", "B", "C", "D"}
    assert rxns.active_enzymes == set()


    # Add a multi-term reaction (reaction index 3)
    rxns.add_reaction(reactants=["A", (2, "B", 1)], products=[(3, "C", 2), "D"],
                      forward_rate=5., reverse_rate=1.)
    assert rxns.number_of_reactions() == 4

    r = rxns.get_reaction(3)
    assert np.allclose(r.kF , 5.)
    assert np.allclose(r.kR , 1.)
    assert np.allclose(r.K , 5./1.)
    assert r.reactants == [(1, "A", 1), (2, "B", 1)]
    assert r.products == [(3, "C", 2), (1, "D", 1)]
    assert r.delta_H is None
    assert r.delta_S is None
    assert np.allclose(r.delta_G, -2676.321364705849)   # - RT log(K)

    assert rxns.active_chemicals == {"A", "B", "C", "D"}
    assert rxns.active_enzymes == set()


    # Check the descriptions we has so far
    rxn_info = rxns.multiple_reactions_describe()
    assert rxn_info[0] == '0: A <-> B  (kF = 3 / kR = 2 / K = 1.5) | 1st order in all reactants & products'
    assert rxn_info[1] == '1: 2 B <-> 5 C  (kF = 9 / kR = 7 / K = 1.2857) | 1st order in all reactants & products'
    assert rxn_info[2] == '2: 2 D <-> C  (kF = 11 / kR = 13 / delta_G = 277.79 / K = 0.84615) | 3-th order in reactant D | 2-th order in product C'
    assert rxn_info[3] == '3: A + 2 B <-> 3 C + D  (kF = 5 / kR = 1 / delta_G = -2,676.3 / K = 5) | 2-th order in product C'


    # Add another reaction (reaction index 4), this time with thermodynamic data;
    # the reverse reaction rate will get computed from the thermodynamic data
    rxns.add_reaction(reactants=["A"], products=[(2, "B", 1)], forward_rate=10.,
                      delta_H= 5., delta_S= 0.4)
    assert rxns.number_of_reactions() == 5

    r = rxns.get_reaction(4)
    assert r.reactants == [(1, "A", 1)]
    assert r.products == [(2, "B", 1)]
    assert r.delta_H == 5.
    assert r.delta_S == 0.4
    assert np.allclose(r.delta_G, -75.0)         # 5 - 200 * 0.4
    assert np.allclose(r.K , 1.0461347154679432)  # exp(75/(8.3144598 * 200))
    assert np.allclose(r.kF , 10.)
    assert np.allclose(r.kR , 9.558998331803693) # 10. / 1.0461347154679432

    assert rxns.active_chemicals == {"A", "B", "C", "D"}
    assert rxns.active_enzymes == set()


    # Add another reaction (reaction index 5), this time involving a catalytic term
    rxns.add_reaction(reactants=[(3, "A", 1), "D"], products=["D", (2, "B", 1)])
    assert rxns.number_of_reactions() == 6
    assert rxns.active_chemicals == {"A", "B", "C", "D"}
    assert rxns.active_enzymes == {"D"}         # Notice that chemical "D is an enzyme in this reaction,
                                                # but an active participant in others

    # Add another reaction (reaction index 6), again involving a catalytic term
    rxns.add_reaction(reactants=["F", (2, "C", 1), (3, "A", 1)], products=["D", "F", (2, "D", 1)])
    assert rxns.number_of_reactions() == 7
    assert rxns.active_chemicals == {"A", "B", "C", "D"}
    assert rxns.active_enzymes == {"D", "F"}

    # Add another reaction (reaction index 7), where "F" (reaction index 5) is NOT a catalyst
    rxns.add_reaction(reactants=[(2, "B", 1), "F"], products=[(3, "A", 1)])
    assert rxns.number_of_reactions() == 8
    assert rxns.active_chemicals == {"A", "B", "C", "D", "F"}
    assert rxns.active_enzymes == {"D", "F"}



def test_clear_reactions_data():
    pass   # TODO




#################  TO DESCRIBE THE DATA  #################

def test_multiple_reactions_describe():
    pass   # TODO

def test_single_reaction_describe():
    pass   # TODO



def test_names_of_active_chemicals():
    chem_data = ChemData(names=['A', 'B', 'C', 'X', 'Y'])
    rxns = Reactions(chem_data)

    assert rxns.labels_of_active_chemicals() == set()    # No reactions yet

    rxns.add_reaction(reactants="A", products="B")
    assert rxns.labels_of_active_chemicals() == {"A", "B"}

    rxns.add_reaction(reactants=["B", "X"], products=["C", "X"])
    assert rxns.labels_of_active_chemicals() == {"A", "B", "C"}  # "X" is an enzyme

    rxns.add_reaction(reactants="X", products="Y")
    assert rxns.labels_of_active_chemicals() == {"A", "B", "C", "X", "Y"}
    # "X" is now involved in some reactions in a non-enzymatic role

    rxns.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"])
    assert rxns.labels_of_active_chemicals() == {"A", "B", "C", "X", "Y"}  # "Z" is an enzyme



def test_indexes_of_active_chemicals():
    chem_data = ChemData(names=['Y', 'X', 'C', 'B', 'A'])
    rxns = Reactions(chem_data)

    assert rxns.indexes_of_active_chemicals() == []                 # No reactions yet

    rxns.add_reaction(reactants="A", products="B")
    assert rxns.indexes_of_active_chemicals() == [3, 4]             # ["A", "B"]

    rxns.add_reaction(reactants=["B", "X"], products=["C", "X"])
    assert rxns.indexes_of_active_chemicals() == [2, 3, 4]          # ["A", "B", "C"]

    rxns.add_reaction(reactants="X", products="Y")
    assert rxns.indexes_of_active_chemicals() == [0, 1, 2, 3, 4]    # All

    rxns.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"])
    assert rxns.indexes_of_active_chemicals() == [0, 1, 2, 3, 4]    # All



def test_names_of_enzymes():
    chem_data = ChemData(names=['A', 'B', 'C', 'X', 'Y'])
    rxns = Reactions(chem_data)

    assert rxns.names_of_enzymes() == set()    # No enzymes yet

    rxns.add_reaction(reactants="A", products="B")
    assert rxns.names_of_enzymes() == set()     # No enzymes yet

    rxns.add_reaction(reactants=["B", "X"], products=["C", "X"])
    assert rxns.names_of_enzymes() == {"X"}     # "X" is an enzyme

    rxns.add_reaction(reactants="X", products="Y")
    assert rxns.names_of_enzymes() == {"X"}      # "X" is an enzyme in SOME reaction

    rxns.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"])
    assert rxns.names_of_enzymes() == {"X", "Z"}





#############  FOR CREATION OF NETWORK DIAGRAMS  #############

def test_prepare_graph_network():
    # Set up an A <-> B reaction
    chem_data = ChemData(names=["A", "B"])
    rxns = Reactions(chem_data)

    rxns.add_reaction(reactants="A", products="B", forward_rate=3., reverse_rate=2.)

    graph_data = rxns.prepare_graph_network()

    expected_structure = [{'name': 'A', 'diff_rate': None, 'id': 'C-0', 'labels': ['Chemical']},
                          {'name': 'B', 'diff_rate': None, 'id': 'C-1', 'labels': ['Chemical']},
                          {'name': 'RXN', 'kF': '3', 'kR': '2', 'delta_G': '-1,005.13', 'K': '1.5', 'id': 'R-0', 'labels': ['Reaction']},

                          {'name': 'produces', 'source': 'R-0', 'target': 'C-1', 'id': 'edge-1', 'stoich': 1, 'rxn_order': 1},
                          {'name': 'reacts', 'source': 'C-0', 'target': 'R-0', 'id': 'edge-2', 'stoich': 1, 'rxn_order': 1}
                          ]

    assert compare_recordsets(graph_data["structure"], expected_structure)
    assert graph_data["color_mapping"] == {'Chemical': '#8DCC92', 'Reaction': '#D9C8AD'}
    assert graph_data["caption_mapping"] == {'Chemical': 'name', 'Reaction': 'name'}
