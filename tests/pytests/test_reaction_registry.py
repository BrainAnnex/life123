import pytest
import numpy as np
from life123 import ChemData, ReactionRegistry, ReactionKinetics, \
                ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition, ReactionGeneric
from tests.utilities.comparisons import *



def test_constructor_ReactionRegistry():
    chem_data = ChemData()
    rxns = ReactionRegistry(chem_data=chem_data)

    assert rxns.reaction_list == []

    assert rxns.active_chemicals == set()       # Empty set



def test_number_of_reactions():
    chem_data = ChemData(names=["A", "B", "C"])
    rxns = ReactionRegistry(chem_data=chem_data)

    assert rxns.number_of_reactions() == 0

    rxns.add_reaction(reactants="A", products="B")
    assert rxns.number_of_reactions() == 1

    rxns.add_reaction(reactants=["A", (2, "B")], products="C")
    assert rxns.number_of_reactions() == 2



def test_active_reaction_indices():
    chem_data = ChemData(names=["A", "B", "C"])
    rxns = ReactionRegistry(chem_data)

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



def test_get_reactant_species():
    rxns = ReactionRegistry()

    rxns.add_reaction(reactants="A", products="B")
    assert rxns.get_reactant_species(0) == ["A"]

    rxns.add_reaction(reactants=[(2, "F")], products="X")
    assert rxns.get_reactant_species(0) == ["A"]
    assert rxns.get_reactant_species(1) == ["F"]

    rxns.add_reaction(reactants=["B", (2, "H")], products=[(2, "C"), "L"])
    assert rxns.get_reactant_species(0) == ["A"]
    assert rxns.get_reactant_species(1) == ["F"]
    assert rxns.get_reactant_species(2) == ["B", "H"]


def test_get_reactants_formula():
    pass   # TODO


def test_get_products():
    pass   # TODO



def test_get_product_species():
    rxns = ReactionRegistry()

    rxns.add_reaction(reactants="A", products="B")
    assert rxns.get_product_species(0) == ["B"]

    rxns.add_reaction(reactants=[(2, "F")], products="X")
    assert rxns.get_product_species(0) == ["B"]
    assert rxns.get_product_species(1) == ["X"]

    rxns.add_reaction(reactants=["B", (2, "H")], products=[(2, "C"), "L"])
    assert rxns.get_product_species(0) == ["B"]
    assert rxns.get_product_species(1) == ["X"]
    assert rxns.get_product_species(2) == ["C", "L"]



def test_get_products_formula():
    pass   # TODO

def test_get_forward_rate():
    pass   # TODO

def test_get_reverse_rate():
    pass   # TODO



def test_get_chemicals_in_reaction():
    chem_data = ChemData(names=["A", "B"])
    rxns = ReactionRegistry(chem_data)

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
    rxns = ReactionRegistry(chem_data)

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



def test_add_reaction():
    chem_data = ChemData(names=["A", "B", "C", "D", "E", "F"])
    rxns = ReactionRegistry(chem_data)

    # Reactants and the products can't be the same
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A"], products=["A"])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=[1, "A"], products=[("A")])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A"], products=[(1, "A")])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A", "B"], products=["A", "B"])
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A", "B"], products=["B", "A"])


    # Add the first (0-th) reaction : A <-> B  ('ReactionUnimolecular')
    rxns.add_reaction(reactants=["A"], products=["B"], kF=3., kR=2.)

    assert rxns.number_of_reactions() == 1
    r = rxns.get_reaction(0)
    assert type(r) == ReactionUnimolecular
    assert np.allclose(r.kF , 3.)
    assert np.allclose(r.kR , 2.)
    assert np.allclose(r.K , 3./2.)
    assert r.reactant == "A"
    assert r.product ==  "B"
    assert r.delta_H is None
    assert r.delta_S is None
    assert r.delta_G is None

    assert rxns.get_reactants(0) == [(1, "A")]
    assert rxns.get_reactants_formula(0) == "A"

    assert rxns.get_products(0) == [(1, "B")]
    assert rxns.get_products_formula(0) == "B"

    assert np.allclose(rxns.get_forward_rate(0), 3.)
    assert np.allclose(rxns.get_reverse_rate(0),  2.)

    assert rxns.active_chemicals == {"A", "B"}
    #assert rxns.active_enzymes == set()


    # Another reaction (reaction 1):  2 B <-> 5 C     ('ReactionGeneric')
    rxns.add_reaction(reactants=[(2, "B")], products=[(5, "C")], kF=9., kR=7.)

    assert rxns.number_of_reactions() == 2

    r = rxns.get_reaction(0)
    assert type(r) == ReactionUnimolecular
    assert np.allclose(r.kF , 3.)
    assert np.allclose(r.kR , 2.)
    assert np.allclose(r.K , 3./2.)
    assert r.reactant == "A"
    assert r.product ==  "B"
    assert r.delta_H is None
    assert r.delta_S is None
    assert r.delta_G is None

    r = rxns.get_reaction(1)
    assert type(r) == ReactionGeneric
    assert np.allclose(r.kF , 9.)
    assert np.allclose(r.kR , 7.)
    assert np.allclose(r.K , 9./7.)
    assert r.reactants == [(2, "B")]
    assert r.products == [(5, "C")]
    assert r.delta_H is None
    assert r.delta_S is None
    assert r.delta_G is None

    assert rxns.active_chemicals == {"A", "B", "C"}
    #assert rxns.active_enzymes == set()


    # Add another reaction (reaction index 2) : 2 D <-> C  ('ReactionSynthesis')
    # This time, specify the temperature
    rxns.add_reaction(reactants=[(2, "D")], products=[(1, "C")],
                      kF=11., kR=13., temp=200)
    assert rxns.number_of_reactions() == 3

    r = rxns.get_reaction(2)
    assert type(r) == ReactionSynthesis
    assert np.allclose(r.kF , 11.)
    assert np.allclose(r.kR , 13.)
    assert np.allclose(r.K , 11./13.)
    assert r.reactant_1 == "D"
    assert r.reactant_2 == "D"
    assert r.product == "C"
    assert r.delta_H is None
    assert r.delta_S is None
    assert np.allclose(r.delta_G, 277.7928942715384)   # - RT log(K)

    assert rxns.active_chemicals == {"A", "B", "C", "D"}
    #assert rxns.active_enzymes == set()


    # Add a multi-term reaction (reaction index 3):  A + 2B <-> 3C + D  ('ReactionGeneric')
    rxns.add_reaction(reactants=["A", (2, "B")], products=[(3, "C"), "D"],
                      kF=5., kR=1., temp=200)
    assert rxns.number_of_reactions() == 4

    r = rxns.get_reaction(3)
    assert type(r) == ReactionGeneric
    assert np.allclose(r.kF , 5.)
    assert np.allclose(r.kR , 1.)
    assert np.allclose(r.K , 5./1.)
    assert r.reactants == [(1, "A"), (2, "B")]
    assert r.products == [(3, "C"), (1, "D")]
    assert r.delta_H is None
    assert r.delta_S is None
    assert np.allclose(r.delta_G, -2676.321364705849)   # - RT log(K)

    assert rxns.active_chemicals == {"A", "B", "C", "D"}
    #assert rxns.active_enzymes == set()


    # Check the descriptions we has so far
    rxn_info = rxns.multiple_reactions_describe()
    assert rxn_info[0] == '0: A <-> B  Elementary Unimolecular reaction  (kF = 3 / kR = 2 / K = 1.5)'
    assert rxn_info[1] == '1: 2 B <-> 5 C  (kF = 9 / kR = 7 / K = 1.2857) | Kinetic rate function : `compute_rate_mass_action_kinetics()`'
    assert rxn_info[2] == '2: 2 D  <-> C  Elementary Synthesis reaction  (kF = 11 / kR = 13 / delta_G = 277.79 / K = 0.84615 / Temp = -73.15 C)'
    assert rxn_info[3] == '3: A + 2 B <-> 3 C + D  (kF = 5 / kR = 1 / delta_G = -2,676.3 / K = 5 / Temp = -73.15 C) | Kinetic rate function : `compute_rate_mass_action_kinetics()`'


    # Add another reaction (reaction index 4), this time with thermodynamic data;
    # the reverse reaction rate will get computed from the thermodynamic data
    rxns.add_reaction(reactants=["A"], products=[(2, "B")], kF=10.,
                      delta_H= 5., delta_S= 0.4, temp=200)      # A <-> 2 B     ('ReactionDecomposition')
    assert rxns.number_of_reactions() == 5

    r = rxns.get_reaction(4)
    assert type(r) == ReactionDecomposition
    assert r.reactant == "A"
    assert r.product_1 == "B"
    assert r.product_2 == "B"
    assert r.delta_H == 5.
    assert r.delta_S == 0.4
    assert np.allclose(r.delta_G, -75.0)            # 5 - 200 * 0.4
    assert np.allclose(r.K , 1.0461347154679432)    # exp(75/(8.3144598 * 200))
    assert np.allclose(r.kF , 10.)
    assert np.allclose(r.kR , 9.558998331803693)    # 10. / 1.0461347154679432

    assert rxns.active_chemicals == {"A", "B", "C", "D"}


    # Add another reaction (reaction index 5)
    i = rxns.add_reaction(reactants=[(2, "B"), "F"], products=[(3, "A")])   # 2 B + F = 3 A     ('ReactionGeneric')
    assert rxns.number_of_reactions() == 6
    assert rxns.active_chemicals == {"A", "B", "C", "D", "F"}

    r = rxns.get_reaction(i)
    assert type(r) == ReactionGeneric



def test_add_elementary_reaction():
    chem_data = ChemData(names=["A", "B", "C"])
    rxns = ReactionRegistry(chem_data)

    i = rxns.add_elementary_reaction(reactants="A", products="B")
    assert type(rxns.get_reaction(i)) == ReactionUnimolecular

    i = rxns.add_elementary_reaction(reactants=["A", "B"], products="C")
    assert type(rxns.get_reaction(i)) == ReactionSynthesis

    i = rxns.add_elementary_reaction(reactants=["A", "A"], products=["B"])
    assert type(rxns.get_reaction(i)) == ReactionSynthesis

    i = rxns.add_elementary_reaction(reactants="A", products=["B", "C"])
    assert type(rxns.get_reaction(i)) == ReactionDecomposition

    with pytest.raises(Exception):
        rxns.add_elementary_reaction(reactants=["A", "A"], products=["B", "C"])



def test_register_reaction():
    chem_data = ChemData()
    rxns = ReactionRegistry(chem_data)

    r_uni_AB = ReactionUnimolecular(reactant="A", product="B",
                                    kF=11., kR=13., temp=200)
    rxns.register_reaction(r_uni_AB)

    assert rxns.number_of_reactions() == 1

    r = rxns.get_reaction(0)        # Get the 0-th (and only so far) reaction
    assert r == r_uni_AB
    assert r.reactant == "A"
    assert r.product == "B"
    assert np.allclose(r.kF, 11.)
    assert np.allclose(r.kR, 13.)
    assert np.allclose(r.K , 11./13.)
    assert r.delta_H is None
    assert r.delta_S is None
    assert np.allclose(r.delta_G, 277.7928942715384)   # - RT log(K)
    assert rxns.active_chemicals == {"A", "B"}


    r_uni_CD = ReactionUnimolecular(reactant="C", product="D",
                                    kF=10., delta_H= 5., delta_S= 0.4, temp=200)
    rxns.register_reaction(r_uni_CD)

    assert rxns.number_of_reactions() == 2

    r = rxns.get_reaction(1)
    assert r == r_uni_CD

    assert r.reactant == "C"
    assert r.product == "D"
    assert np.allclose(r.kF, 10.)                   # The given value
    assert r.delta_H == 5.                          # The given value
    assert r.delta_S == 0.4                         # The given value
    assert np.allclose(r.delta_G, -75.0)            # 5 - 200 * 0.4
    assert np.allclose(r.K , 1.0461347154679432)    # exp(75/(8.3144598 * 200))
    assert np.allclose(r.kR, 9.558998331803693)     # 10. / 1.0461347154679432
    assert rxns.active_chemicals == {"A", "B", "C", "D"}



def test_clear_reactions_data():
    pass   # TODO



def test_determine_reaction_rate_ReactionGeneric():
    # Reaction A <-> B
    rxn = ReactionGeneric(reactants="A", products="B",
                          kF=20., kR=2., reversible=True)
    rxn.set_rate_function(ReactionKinetics.compute_rate_mass_action_kinetics)  # "standard rate law"
    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5. - 2. * 8.)  # 84.0

    rxn.kR = 0

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5.)            # 100.0


    # Switch to using the following test function for the reaction rate
    def my_rate_law_1(reactant_terms :[(int, str)], product_terms :[(int, str)],
                      kF :float, kR :float,
                      conc_dict :dict) -> float:
        return 1234.


    rxn.set_rate_function(my_rate_law_1)

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 1234.)



    # Reaction 5A <-> 2B , hypothetically with 1st-order kinetics in both directions.
    rxn = ReactionGeneric(reactants=[(5, "A")], products=[(2, "B")],
                          kF=20., kR=2.)

    rxn.set_rate_function(ReactionKinetics.compute_rate_first_order)

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})

    assert np.allclose(result, 20. * 5. - 2. * 8.)      # 84.0



    # Reaction  2A <-> B , with 2nd-ORDER kinetics in the forward direction
    rxn = ReactionGeneric(reactants=[(2, "A")], products="B",
                          kF=5., kR=2., reversible=True)
    rxn.set_rate_function(ReactionKinetics.compute_rate_mass_action_kinetics)  # "standard rate law"
    result = rxn.determine_reaction_rate(conc_dict={"A": 4.5, "B": 6.})
    assert np.allclose(result, 5. * 4.5 **2 - 2. * 6.)  # 89.25

    rxn.kR = 0
    result = rxn.determine_reaction_rate(conc_dict={"A": 4.5, "B": 6.})
    assert np.allclose(result, 5. * 4.5 **2)            # 101.25

    rxn.kR = 2.
    rxn.kF = 3.
    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 59.)


    # Reaction  B <-> 2C , with 2nd-ORDER kinetics in the reverse direction
    rxn = ReactionGeneric(reactants="B", products=[(2, "C")], kF=4., kR=2.)
    rxn.set_rate_function(ReactionKinetics.compute_rate_mass_action_kinetics)  # "standard rate law"
    result = rxn.determine_reaction_rate(conc_dict={"B": 5., "C": 4})
    assert np.allclose(result, 4. * 5. - 2. * 4. **2)           # -12.0


    # Reaction 2B <-> 3C , hypothetically with 1st-order kinetics in both directions
    rxn = ReactionGeneric(reactants=[(2, "B")], products=[(3, "C")],
                          kF=10., kR=25.)
    rxn.set_rate_function(ReactionKinetics.compute_rate_first_order)
    result = rxn.determine_reaction_rate(conc_dict={"B": 8., "C": 15.})
    assert np.allclose(result,  10. * 8. - 25. * 15.)   # -295.0


    # Reaction 2A + 5B <-> 4C + 3D , hypothetically with 1st-order kinetics for each species
    rxn = ReactionGeneric(reactants=[(2, "A") , (5, "B")],
                          products=[(4, "C") , (3, "D")],
                          kF=5., kR=2.)
    rxn.set_rate_function(ReactionKinetics.compute_rate_first_order)
    result = rxn.determine_reaction_rate(conc_dict={"A": 3.5, "B": 9., "C": 11., "D": 7.})
    assert np.allclose(result,  5. * 3.5 * 9. - 2. * 11. * 7.)  # 3.5

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 15., "D": 7.})
    assert np.allclose(result,  -10.)




#################  TO DESCRIBE THE DATA  #################

def test_multiple_reactions_describe():
    pass   # TODO

def test_single_reaction_describe():
    pass   # TODO



def test_labels_of_active_chemicals():
    chem_data = ChemData(names=['A', 'B', 'C', 'X', 'Y'])
    rxns = ReactionRegistry(chem_data)

    assert rxns.labels_of_active_chemicals() == []   # No reactions yet

    rxns.add_reaction(reactants="A", products="B")
    assert set(rxns.labels_of_active_chemicals()) == {"A", "B"}
    assert rxns.labels_of_active_chemicals(sort_by_index=True) == ["A", "B"]

    rxns.add_reaction(reactants=["B", "X"], products=["C", "X"])
    assert set(rxns.labels_of_active_chemicals()) == {"A", "B", "C", "X"}
    assert rxns.labels_of_active_chemicals(sort_by_index=True) == ["A", "B", "C", "X"]

    rxns.add_reaction(reactants="X", products="Y")
    assert set(rxns.labels_of_active_chemicals()) == {"A", "B", "C", "X", "Y"}
    assert rxns.labels_of_active_chemicals(sort_by_index=True) == ["A", "B", "C", "X", "Y"]

    rxns.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"])
    assert set(rxns.labels_of_active_chemicals()) == {"A", "B", "C", "X", "Y", "Z"}
    assert rxns.labels_of_active_chemicals(sort_by_index=True) == ["A", "B", "C", "X", "Y", "Z"]



def test_indexes_of_active_chemicals():
    chem_data = ChemData(names=['Y', 'X', 'C', 'B', 'A'])
    rxns = ReactionRegistry(chem_data)

    assert rxns.indexes_of_active_chemicals() == []                 # No reactions yet

    rxns.add_reaction(reactants="A", products="B")
    assert rxns.indexes_of_active_chemicals() == [3, 4]             # ["A", "B"]

    rxns.add_reaction(reactants=["B", "X"], products=["C", "X"])
    assert rxns.indexes_of_active_chemicals() == [1, 2, 3, 4]          # ["X", "A", "B", "C"]

    rxns.add_reaction(reactants="X", products="Y")
    assert rxns.indexes_of_active_chemicals() == [0, 1, 2, 3, 4]    # All

    rxns.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"])
    assert rxns.indexes_of_active_chemicals() == [0, 1, 2, 3, 4, 5]    # All



def test__parse_reaction_term():
    rxn = ReactionRegistry(chem_data=ChemData())     # Won't actually use the reactants/products

    with pytest.raises(Exception):
        rxn._parse_reaction_term(5)    # The argument is not a string nor a tuple nor a list

    assert rxn._parse_reaction_term("F") == (1, "F")


    with pytest.raises(Exception):
        rxn._parse_reaction_term( (2, 5) )   # The last item in the pair is not a string

    assert rxn._parse_reaction_term( (2, "F") ) == (2, "F")    # order defaults to stoichiometry
    assert rxn._parse_reaction_term( [2, "F"] ) == (2, "F")

    with pytest.raises(Exception):
        rxn._parse_reaction_term( (2, 5) )   # The mid-item in the triplet is not a string

    assert rxn._parse_reaction_term( (2, "F") ) == (2, "F")
    assert rxn._parse_reaction_term( [2, "F"] ) == (2, "F")

    with pytest.raises(Exception):
        rxn._parse_reaction_term( (3, "F", 2, 123) )     # Extra element in tuple





#############  FOR CREATION OF NETWORK DIAGRAMS  #############


def test_prepare_graph_network():
    # Set up an A <-> B reaction
    chem_data = ChemData(names=["A", "B"])
    rxns = ReactionRegistry(chem_data)

    rxns.add_reaction(reactants="A", products="B", kF=3., kR=2., temp=298.15)

    graph_data = rxns.prepare_graph_network()

    expected_nodes = [{'name': 'A', 'id': 'C-0', '_node_labels': ['Chemical']},
                      {'name': 'B', 'id': 'C-1', '_node_labels': ['Chemical']},
                      {'name': 'RXN', 'kF': '3', 'kR': '2', 'delta_G': '-1,005.13', 'K': '1.5', 'id': 'RXN-0',
                       '_node_labels': ['Reaction'], 'formula': 'A <-> B'}
                      ]     # Note: 'diff_rate': None   is NOT included

    expected_edges = [
                        {'name': 'produces', 'source': 'RXN-0', 'target': 'C-1', 'id': 'edge-1', 'stoich': 1},
                        {'name': 'reacts', 'source': 'C-0', 'target': 'RXN-0', 'id': 'edge-2', 'stoich': 1}
                      ]

    assert compare_recordsets(graph_data["nodes"], expected_nodes)
    assert compare_recordsets(graph_data["edges"], expected_edges)
    assert graph_data["color_mapping"] == {'Chemical': '#8DCC92', 'Reaction': '#D9C8AD'}
    assert graph_data["caption_mapping"] == {'Chemical': 'name', 'Reaction': 'id'}
