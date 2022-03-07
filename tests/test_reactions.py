import pytest
from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio



# Provide an instantiated object that can be used by the various tests that need it
@pytest.fixture(scope="module")
def rxn():
    chem_data = chem()
    chem_data.set_names(["A", "B", "C", "D", "E", "F"])

    bio.initialize_universe(n_bins=10, chem_data=chem_data)

    rnx_obj = Reactions(chem_data)
    yield rnx_obj



def test_chem_data(rxn):
    assert rxn.chem_data.names == ["A", "B", "C", "D", "E", "F"]
    assert rxn.chem_data.name_dict == {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5}



def test_parse_reaction_term(rxn):
    assert rxn._parse_reaction_term(5) == (1, 5, 1)
    assert rxn._parse_reaction_term("F") == (1, 5, 1)

    assert rxn._parse_reaction_term( (3, 5) ) == (3, 5, 1)
    assert rxn._parse_reaction_term( (3, "F") ) == (3, 5, 1)

    assert rxn._parse_reaction_term( (3, 5, 2) ) == (3, 5, 2)
    assert rxn._parse_reaction_term( (3, "F", 2) ) == (3, 5, 2)



def test_add_reaction(rxn):
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert rxn.number_of_reactions() == 1
    assert rxn.get_reaction(0) == {'reactants': [(1, 0, 1)], 'products': [(1, 1, 1)], 'Rf': 3.0, 'Rb': 2.0}

    assert rxn.get_reactants(0) == [(1, 0, 1)]
    assert rxn.get_reactants_formula(0) == "A"

    assert rxn.get_products(0) == [(1, 1, 1)]
    assert rxn.get_products_formula(0) == "B"

    assert rxn.get_forward_rate(0) == 3.
    assert rxn.get_reverse_rate(0) == 2.

    # A 2nd reaction
    rxn.add_reaction(reactants=[(2, "B")], products=[(5, "C")], forward_rate=9., reverse_rate=7.)

    assert rxn.number_of_reactions() == 2
    assert rxn.get_reaction(0) == {'reactants': [(1, 0, 1)], 'products': [(1, 1, 1)], 'Rf': 3.0, 'Rb': 2.0}
    assert rxn.get_reaction(1) == {'reactants': [(2, 1, 1)], 'products': [(5, 2, 1)], 'Rf': 9.0, 'Rb': 7.0}

    # A 3rd reaction
    rxn.add_reaction(reactants=[(2, "D", 3)], products=[(1, "C", 2)], forward_rate=11., reverse_rate=13.)
    assert rxn.number_of_reactions() == 3
    assert rxn.get_reaction(2) == {'reactants': [(2, 3, 3)], 'products': [(1, 2, 2)], 'Rf': 11.0, 'Rb': 13.0}

    # A multi-term reaction
    rxn.add_reaction(reactants=["A", (2, "B")], products=[(3, "C", 2), "D"], forward_rate=5., reverse_rate=1.)
    assert rxn.number_of_reactions() == 4
    assert rxn.get_reaction(3) == {'reactants': [(1, 0, 1), (2, 1, 1)], 'products': [(3, 2, 2), (1, 3, 1)], 'Rf': 5.0, 'Rb': 1.0}

    rxn_list = rxn.describe_reactions()
    assert rxn_list[0] == '0: A <-> B  (Rf = 3.0 / Rb = 2.0)'
    assert rxn_list[1] == '1: 2 B <-> 5 C  (Rf = 9.0 / Rb = 7.0)'
    assert rxn_list[2] == '2: 2 D <-> C  (Rf = 11.0 / Rb = 13.0)'
    assert rxn_list[3] == '3: A + 2 B <-> 3 C + D  (Rf = 5.0 / Rb = 1.0)'
