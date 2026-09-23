import pytest
import numpy as np
from life123 import SpeciesRegistry, ReactionRegistry, ReactionKinetics
from life123.reactions import Stoichiometry, ReactionDefinition, SimulationReaction
from tests.utilities.comparisons import *



def test_constructor_ReactionRegistry():
    chem_data = SpeciesRegistry()
    rxns = ReactionRegistry(species_data=chem_data)

    assert rxns.reaction_list == []

    assert rxns.active_chemicals == set()       # Empty set



def test_number_of_reactions():
    chem_data = SpeciesRegistry(ids=["A", "B", "C"])
    rxns = ReactionRegistry(species_data=chem_data)

    assert rxns.number_of_reactions() == 0

    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action")
    assert rxns.number_of_reactions() == 1

    rxns.add_reaction(reactants=["A", (2, "B")], products="C", reaction_model="mass action")
    assert rxns.number_of_reactions() == 2



def test_active_reaction_indices():
    chem_data = SpeciesRegistry(ids=["A", "B", "C"])
    rxns = ReactionRegistry(chem_data)

    assert rxns.active_reaction_indices() == []

    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action")
    assert rxns.active_reaction_indices() == [0]

    rxns.add_reaction(reactants=["A", (2, "B")], products="C", reaction_model="mass action")
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
    chem_data = SpeciesRegistry(ids=["A", "B"])
    rxns = ReactionRegistry(chem_data)

    with pytest.raises(Exception):
        rxns.get_chemicals_in_reaction(0)   # There are no reactions defined yet

    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action")  # Reaction 0 : A <-> B
    assert rxns.get_chemicals_in_reaction(0) == {0, 1}

    with pytest.raises(Exception):
        rxns.get_chemicals_in_reaction(1)   # There is no reaction 1

    chem_data.add_species("C")

    rxns.add_reaction(reactants=["B"], products=[(2, "C")], reaction_model="mass action")  # Reaction 1 : B <-> 2C
    assert rxns.get_chemicals_in_reaction(0) == {0, 1}
    assert rxns.get_chemicals_in_reaction(1) == {1, 2}

    rxns.add_reaction(reactants=["A"], products=["C"], reaction_model="mass action")      # Reaction 2 : A <-> C
    assert rxns.get_chemicals_in_reaction(0) == {0, 1}
    assert rxns.get_chemicals_in_reaction(1) == {1, 2}
    assert rxns.get_chemicals_in_reaction(2) == {0, 2}

    chem_data.add_species("D")
    rxns.add_reaction(reactants=["A", "B"], products="D", reaction_model="mass action")    # Reaction 3 : A + B <-> D
    assert rxns.get_chemicals_in_reaction(0) == {0, 1}
    assert rxns.get_chemicals_in_reaction(1) == {1, 2}
    assert rxns.get_chemicals_in_reaction(2) == {0, 2}
    assert rxns.get_chemicals_in_reaction(3) == {0, 1, 3}



def test_get_chemicals_indexes_in_reaction():
    chem_data = SpeciesRegistry(ids=["A", "B"])
    rxns = ReactionRegistry(chem_data)

    with pytest.raises(Exception):
        rxns.get_chemicals_indexes_in_reaction(0)   # There are no reactions defined yet

    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action")  # Reaction 0 : A <-> B
    assert rxns.get_chemicals_indexes_in_reaction(0) == [0, 1]

    with pytest.raises(Exception):
        rxns.get_chemicals_indexes_in_reaction(1)   # There is no reaction 1

    chem_data.add_species("C")

    rxns.add_reaction(reactants=["B"], products=[(2, "C")], reaction_model="mass action")  # Reaction 1 : B <-> 2C
    assert rxns.get_chemicals_indexes_in_reaction(0) == [0, 1]
    assert rxns.get_chemicals_indexes_in_reaction(1) == [1, 2]

    rxns.add_reaction(reactants=["A"], products=["C"], reaction_model="mass action")      # Reaction 2 : A <-> C
    assert rxns.get_chemicals_indexes_in_reaction(0) == [0, 1]
    assert rxns.get_chemicals_indexes_in_reaction(1) == [1, 2]
    assert rxns.get_chemicals_indexes_in_reaction(2) == [0, 2]

    chem_data.add_species("D")
    rxns.add_reaction(reactants=["A", "B"], products="D", reaction_model="mass action")    # Reaction 3 : A + B <-> D
    assert rxns.get_chemicals_indexes_in_reaction(0) == [0, 1]
    assert rxns.get_chemicals_indexes_in_reaction(1) == [1, 2]
    assert rxns.get_chemicals_indexes_in_reaction(2) == [0, 2]
    assert rxns.get_chemicals_indexes_in_reaction(3) == [0, 1, 3]



def test_get_reactions_participating_in():
    rxns = ReactionRegistry()

    result = rxns.get_reactions_participating_in(species_id="A", side="reagent")
    assert result == []

    # A -> B  (reaction 0)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action")

    result = rxns.get_reactions_participating_in(species_id="A", side="reagent")
    assert len(result) == 1
    rnx_defn = result[0]
    assert type(rnx_defn) == ReactionDefinition
    assert (rnx_defn.id) == 0

    result = rxns.get_reactions_participating_in(species_id="A", side="product")
    assert result == []

    result = rxns.get_reactions_participating_in(species_id="B", side="product")
    rnx_defn = result[0]
    assert type(rnx_defn) == ReactionDefinition
    assert (rnx_defn.id) == 0


    # A + C -> P  (reaction 1)
    rxns.add_reaction(reactants=["A", "C"], products="P", reaction_model="mass action")

    result = rxns.get_reactions_participating_in(species_id="A", side="reagent")
    assert len(result) == 2
    assert type(result[0]) == ReactionDefinition
    assert result[0].id == 0
    assert type(result[1]) == ReactionDefinition
    assert result[1].id == 1

    result = rxns.get_reactions_participating_in(species_id="C", side="reagent")
    assert len(result) == 1
    assert type(result[0]) == ReactionDefinition
    assert result[0].id == 1

    result = rxns.get_reactions_participating_in(species_id="C", side="product")
    assert result == []

    result = rxns.get_reactions_participating_in(species_id="P", side="product")
    assert len(result) == 1
    assert type(result[0]) == ReactionDefinition
    assert result[0].id == 1



def test_add_reaction():
    chem_data = SpeciesRegistry(ids=["A", "B", "C", "D", "E", "F"])
    rxns = ReactionRegistry(chem_data)

    assert rxns.number_of_reactions() == 0
    assert len(rxns.reaction_defn_list) == 0

    # Reactants and the products can't be the same
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A"], products=["A"], reaction_model="mass action")
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=[1, "A"], products=[("A")], reaction_model="mass action")
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A"], products=[(1, "A")], reaction_model="mass action")
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A", "B"], products=["A", "B"], reaction_model="mass action")
    with pytest.raises(Exception):
        rxns.add_reaction(reactants=["A", "B"], products=["B", "A"], reaction_model="mass action")


    assert rxns.number_of_reactions() == 0
    assert len(rxns.reaction_defn_list) == 0

    # Add the first (0-th) reaction : A <-> B
    result = rxns.add_reaction(reactants="A", products="B",
                              reaction_model="mass action",
                              kinetic_parameters={"kF": 3., "kR": 2.})
    assert result == 0
    assert rxns.number_of_reactions() == 1
    assert len(rxns.reaction_defn_list) == 1

    r = rxns.get_reaction(0)        # Get the 0-th (and only so far) reaction
    assert type(r) == SimulationReaction
    assert r.stoichiometry == Stoichiometry({"A": -1, "B": 1})

    r_defn = rxns.reaction_defn_list[0]
    assert type(r_defn) == ReactionDefinition
    assert r_defn.id == 0
    assert r_defn.stoichiometry == Stoichiometry({"A": -1, "B": 1})

    assert rxns.active_chemicals == {"A", "B"}

    expected_kinetics = {'kF': 3, 'kR': 2, 'K': 1.5, 'reversible': True}
    # Note:  K derived as:  3/2
    assert compare_dicts(r.get_parameters(), expected_kinetics)
    assert r.model.derived_pars == {'K', 'reversible'}

    expected_thermo = {'K_eq':  1.5}
    # Note:  K_eq derived from kinetics data
    assert compare_dicts(r_defn.thermodynamics.to_dict(), expected_thermo)
    assert r_defn.thermodynamics.derived_pars == {"K_eq"}


    # Another reaction (reaction 1):  2 B <-> 5 C
    result = rxns.add_reaction(reactants=[(2, "B")], products=[(5, "C")],
                              reaction_model="custom",
                              kinetic_parameters={"kF": 9., "kR": 7.})
    assert result == 1
    assert rxns.number_of_reactions() == 2
    assert len(rxns.reaction_defn_list) == 2

    assert rxns.active_chemicals == {"A", "B", "C"}

    assert rxns.get_reaction(0) == r
    assert rxns.reaction_defn_list[0] == r_defn

    r_1 = rxns.get_reaction(1)
    r_defn_1 = rxns.reaction_defn_list[1]

    assert r_1.stoichiometry == Stoichiometry({"B": -2, "C": 5})
    assert r_defn_1.stoichiometry == Stoichiometry({"B": -2, "C": 5})

    expected_kinetics = {'kF': 9, 'kR': 7, 'K': 1.285714285714, 'reversible': True, 'rate_function': None}
    # Note:  K derived as:  9/7
    assert compare_dicts(r_1.get_parameters(), expected_kinetics)
    assert r.model.derived_pars == {'K', 'reversible'}

    expected_thermo = {'K_eq':  1.285714285714}
    # Note:  K_eq derived from kinetics data
    assert compare_dicts(r_defn_1.thermodynamics.to_dict(), expected_thermo)
    assert r_defn_1.thermodynamics.derived_pars == {"K_eq"}


    # Add another reaction (reaction index 2) : 2 D <-> C  ('ReactionSynthesis')
    # This time, specify the temperature
    result = rxns.add_reaction(reactants=[(2, "D")], products=[(1, "C")],
                               reaction_model="mass action",
                               kinetic_parameters={"kF": 11., "kR": 13.},
                               thermodynamic_parameters={"temp": 200})
    assert result == 2
    assert rxns.number_of_reactions() == 3

    assert rxns.active_chemicals == {"A", "B", "C", "D"}

    r_2 = rxns.get_reaction(2)
    r_defn_2 = rxns.reaction_defn_list[2]

    assert r_2.stoichiometry == Stoichiometry({"D": -2, "C": 1})
    assert r_defn_2.stoichiometry == Stoichiometry({"D": -2, "C": 1})

    expected_kinetics = {'kF': 11, 'kR': 13, 'K': 0.84615384615, 'reversible': True}
    # Note:  K derived as:  11/13
    assert compare_dicts(r_2.get_parameters(), expected_kinetics)
    assert r_2.model.derived_pars == {'K', 'reversible'}

    expected_thermo = {'K_eq':  0.84615384615, 'delta_G': 0.2777929884283404, 'temp': 200}
    # Note:  K_eq derived from kinetics data ;  delta_G derived from K_eq and temp : - RT log(K)
    assert compare_dicts(r_defn_2.thermodynamics.to_dict(), expected_thermo)
    assert r_defn_2.thermodynamics.derived_pars == {"K_eq", "delta_G"}


    # Add a multi-term reaction (reaction index 3):  A + 2B <-> 3C + D
    result = rxns.add_reaction(reactants=["A", (2, "B")], products=[(3, "C"), "D"],
                               reaction_model="custom",
                               kinetic_parameters={"kF": 5., "kR": 1.},
                               thermodynamic_parameters={"temp": 200})

    assert result == 3
    assert rxns.number_of_reactions() == 4

    assert rxns.active_chemicals == {"A", "B", "C", "D"}

    r_3 = rxns.get_reaction(3)
    assert r_3.stoichiometry == Stoichiometry({"A": -1, "B": -2, "C": 3, "D": 1})


    # Check the descriptions we has so far
    rxn_info = rxns.multiple_reactions_describe(concise=True)
    assert rxn_info[0] == '0: A <-> B'
    assert rxn_info[1] == '1: 2 B <-> 5 C'
    assert rxn_info[2] == '2: 2 D <-> C'
    assert rxn_info[3] == '3: A + 2 B <-> 3 C + D'



def test_add_elementary_reaction():
    chem_data = SpeciesRegistry(ids=["A", "B", "C"])
    rxns = ReactionRegistry(chem_data)

    assert rxns.number_of_reactions() == 0

    i = rxns.add_elementary_reaction(reactants="A", products="B")
    assert i == 0
    assert rxns.number_of_reactions() == 1
    assert type(rxns.reaction_defn_list[i] is ReactionDefinition)
    assert rxns.reaction_defn_list[i].id == 0
    assert rxns.reaction_defn_list[i].reaction_category == "Unimolecular rearrangement/isomerization"
    assert type(rxns.reaction_list[i] is SimulationReaction)

    i = rxns.add_elementary_reaction(reactants=["A", "B"], products="C")
    assert i == 1
    assert rxns.number_of_reactions() == 2
    assert type(rxns.reaction_defn_list[i] is ReactionDefinition)
    assert rxns.reaction_defn_list[i].id == 1
    assert rxns.reaction_defn_list[i].reaction_category == "Bimolecular synthesis"
    assert type(rxns.reaction_list[i] is SimulationReaction)

    i = rxns.add_elementary_reaction(reactants=["A", "A"], products=["B"])
    assert i == 2
    assert rxns.reaction_defn_list[i].id == 2
    assert rxns.reaction_defn_list[i].reaction_category == "Bimolecular synthesis"

    i = rxns.add_elementary_reaction(reactants="A", products=["B", "C"])
    assert i == 3
    assert rxns.reaction_defn_list[i].id == 3
    assert rxns.reaction_defn_list[i].reaction_category == "Unimolecular decomposition"

    with pytest.raises(Exception):
        rxns.add_elementary_reaction(reactants=["A", "A"], products=["B", "C"])



def test_register_reaction():
    chem_data = SpeciesRegistry()
    rxns = ReactionRegistry(chem_data)

    # A -> B
    r_uni_AB = ReactionDefinition(reactants="A", products="B",
                                  species_registry=chem_data, autoregister_species=True,
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF": 11., "kR":13.},
                                  thermodynamic_parameters={"temp": 200})
    rxns.register_reaction(r_uni_AB)

    assert rxns.number_of_reactions() == 1

    r = rxns.get_reaction(0)        # Get the 0-th (and only so far) reaction
    assert type(r) == SimulationReaction
    assert r.stoichiometry == Stoichiometry({"A": -1, "B": 1})

    r_defn = rxns.reaction_defn_list[0]
    assert type(r_defn) == ReactionDefinition
    assert r_defn == r_uni_AB
    assert r_defn.id == 0

    assert rxns.active_chemicals == {"A", "B"}

    expected_thermo = {'K_eq': 0.8461538461538461, 'delta_G': 0.2777929884283404, 'temp': 200}
    # Note:  K_eq derived from kinetics data  ;  delta_G derived as:  - RT log(K)   , in kJ/mol
    assert compare_dicts(r_defn.thermodynamics.to_dict(), expected_thermo)
    assert r_defn.thermodynamics.derived_pars == {"K_eq", "delta_G"}

    expected_kinetics = {'kF': 11.0, 'kR': 13.0, 'K': 0.8461538461538461, 'reversible': True}
    # Note:  K derived as:  11. / 13.
    assert compare_dicts(r.get_parameters(), expected_kinetics)
    assert r.model.derived_pars == {'K', 'reversible'}


    # C -> D
    r_uni_CD = ReactionDefinition(reactants="C", products="D",
                                  species_registry=chem_data, autoregister_species=True,
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF": 10.},
                                  thermodynamic_parameters={"delta_H": 0.005, "delta_S": 0.4, "temp": 200})
    rxns.register_reaction(r_uni_CD)

    assert rxns.number_of_reactions() == 2

    r = rxns.get_reaction(1)        # Get the 0-th (and only so far) reaction
    assert type(r) == SimulationReaction
    assert r.stoichiometry == Stoichiometry({"C": -1, "D": 1})

    r_defn = rxns.reaction_defn_list[1]
    assert type(r_defn) == ReactionDefinition
    assert r_defn == r_uni_CD
    assert r_defn.id == 1

    assert rxns.active_chemicals == {"A", "B", "C", "D"}

    expected_thermo = {"K_eq": 1.0461346994754837, "delta_G": -0.075, "temp": 200, "delta_H": 0.005, "delta_S": 0.4}
    # Note -  K_eq derived as: exp(75/(8.3144598 * 200))  ;  delta_G derived as: 0.005 - 200 * 0.4/1000 , in kJ/mol
    assert compare_dicts(r_defn.thermodynamics.to_dict(), expected_thermo)
    assert r_defn.thermodynamics.derived_pars == {"K_eq", "delta_G"}

    expected_kinetics = {'kF': 10, 'kR': 9.558998477933914, 'K': 1.0461346994754837, 'reversible': True}
    # Note -  K derived from thermodynamics;   kR derived from: 10. / 1.0461347154679432
    assert compare_dicts(r.get_parameters(), expected_kinetics)
    assert r.model.derived_pars == {'K', 'kR', 'reversible'}



def test_clear_reactions_data():
    pass   # TODO




#######################  TO DESCRIBE THE DATA  #######################

def test_multiple_reactions_describe():
    pass   # TODO

def test_single_reaction_describe():
    pass   # TODO



def test_labels_of_active_chemicals():
    chem_data = SpeciesRegistry(ids=['A', 'B', 'C', 'X', 'Y'])
    rxns = ReactionRegistry(chem_data)

    assert rxns.labels_of_active_chemicals() == []   # No reactions yet

    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action")
    assert set(rxns.labels_of_active_chemicals()) == {"A", "B"}
    assert rxns.labels_of_active_chemicals(sort_by_index=True) == ["A", "B"]

    rxns.add_reaction(reactants=["B", "X"], products=["C", "X"], reaction_model="mass action")
    assert set(rxns.labels_of_active_chemicals()) == {"A", "B", "C", "X"}
    assert rxns.labels_of_active_chemicals(sort_by_index=True) == ["A", "B", "C", "X"]

    rxns.add_reaction(reactants="X", products="Y", reaction_model="mass action")
    assert set(rxns.labels_of_active_chemicals()) == {"A", "B", "C", "X", "Y"}
    assert rxns.labels_of_active_chemicals(sort_by_index=True) == ["A", "B", "C", "X", "Y"]

    rxns.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"], reaction_model="mass action")
    assert set(rxns.labels_of_active_chemicals()) == {"A", "B", "C", "X", "Y", "Z"}
    assert rxns.labels_of_active_chemicals(sort_by_index=True) == ["A", "B", "C", "X", "Y", "Z"]



def test_indexes_of_active_chemicals():
    chem_data = SpeciesRegistry(ids=['Y', 'X', 'C', 'B', 'A'])
    rxns = ReactionRegistry(chem_data)

    assert rxns.indexes_of_active_chemicals() == []                 # No reactions yet

    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action")
    assert rxns.indexes_of_active_chemicals() == [3, 4]             # ["A", "B"]

    rxns.add_reaction(reactants=["B", "X"], products=["C", "X"], reaction_model="mass action")
    assert rxns.indexes_of_active_chemicals() == [1, 2, 3, 4]          # ["X", "A", "B", "C"]

    rxns.add_reaction(reactants="X", products="Y", reaction_model="mass action")
    assert rxns.indexes_of_active_chemicals() == [0, 1, 2, 3, 4]    # All

    rxns.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"], reaction_model="mass action")
    assert rxns.indexes_of_active_chemicals() == [0, 1, 2, 3, 4, 5]    # All



def test__parse_reaction_term():
    rxn = ReactionRegistry(species_data=SpeciesRegistry())     # Won't actually use the reactants/products

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
    chem_data = SpeciesRegistry(ids=["A", "B"])
    rxns = ReactionRegistry(chem_data)

    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.},
                      thermodynamic_parameters={"temp": 298.15})

    graph_data = rxns.prepare_graph_network()

    expected_nodes = [{'name': 'A', 'id': 'C-0', '_node_labels': ['Chemical']},
                      {'name': 'B', 'id': 'C-1', '_node_labels': ['Chemical']},
                      {'name': 'RXN', 'kF': '3', 'kR': '2', 'K': '1.5', 'delta_G': '-1.00513', 'K_eq': '1.5', 'temp': '298.15', 'id': 'RXN-0',
                       'reversible': 'True', 'reaction_model': 'mass action', '_node_labels': ['Reaction'], 'formula': 'A <-> B'}
                      ]     # Note: 'diff_rate': None   is NOT included

    expected_edges = [
                        {'name': 'produces', 'source': 'RXN-0', 'target': 'C-1', 'id': 'edge-1', 'stoich': 1},
                        {'name': 'reacts', 'source': 'C-0', 'target': 'RXN-0', 'id': 'edge-2', 'stoich': 1}
                      ]

    assert compare_recordsets(graph_data["nodes"], expected_nodes)
    assert compare_recordsets(graph_data["edges"], expected_edges)
    assert graph_data["color_mapping"] == {'Chemical': '#8DCC92', 'Reaction': '#D9C8AD'}
    assert graph_data["caption_mapping"] == {'Chemical': 'name', 'Reaction': 'id'}
