import pytest
import numpy as np
from life123.random_reactions import RandomReactionNetwork, ReactionTerm
from life123.reactions import SimulationReaction
from life123.stoichiometry import Stoichiometry


def test_ReactionTerm():
    t = ReactionTerm(2, "R", 20)

    assert t[0] == 2
    assert t.coeff == 2

    assert t[1] == "R"
    assert t.species_id == "R"

    assert t[2] == 20
    assert t.atom_count == 20


def test_RandomReactionNetwork_CONSTRUCTOR():
    net = RandomReactionNetwork(n_species=4, n_rxns=0)
    assert net.species_registry.get_all_species_ids(sort=True) == ["A", "B", "C", "D"]

    # TODO: continue



def test_set_number_atoms():
    net = RandomReactionNetwork(n_species=3, n_rxns=0)
    assert net.species_registry.get_all_species_ids(sort=True) == ["A", "B", "C"]

    assert net.set_number_atoms(species_id="A", atom_count=15) == 15
    sp = net.species_registry.get_species("A")
    assert sp.annotation == "15"

    assert net.set_number_atoms(species_id="C", atom_count=7) == 7
    sp = net.species_registry.get_species("C")
    assert sp.annotation == "7"

    assert net.set_number_atoms(species_id="A", atom_count=3) == 3
    sp = net.species_registry.get_species("A")
    assert sp.annotation == "3"



def test_get_number_atoms():
    net = RandomReactionNetwork(n_species=4, n_rxns=0)

    net.set_number_atoms(species_id="A", atom_count=15)
    assert net.get_number_atoms("A") == 15

    net.set_number_atoms(species_id="C", atom_count=7)
    assert net.get_number_atoms("C") == 7

    net.set_number_atoms(species_id="A", atom_count=3)
    assert net.get_number_atoms("A") == 3

    assert net.get_number_atoms("B") == 0



def test__assign_species_to_rxn():
    net = RandomReactionNetwork(n_species=4, n_rxns=0, seed=123456)

    rxn = net._assign_species_to_rxn()
    assert type(rxn) is SimulationReaction
    assert rxn.stoichiometry.vector == {'A': -1, 'C': -1, 'D': 1}   # A + C -> D
    assert rxn.stoichiometry.catalysts == []

    rxn = net._assign_species_to_rxn()
    assert rxn.stoichiometry.vector == {'B': -2, 'A': 1}            # 2 B -> A

    rxn = net._assign_species_to_rxn()
    assert rxn.stoichiometry.vector == {'A': -1, 'D': -1, 'B': 1}   # A + D -> B

    for i in range(20):
        rxn = net._assign_species_to_rxn()
        stoic = rxn.stoichiometry
        pattern = stoic.reaction_pattern()
        # 1 or 2 reactant species, 1 product, 0 catalysts
        assert pattern == (2, 1, 0)  or pattern == (1, 1, 0)

        r_comp, p_comp = stoic.get_reaction_complexes()
        # EXAMPLE: r_comp = {"A": 1, "B": 1}   , p_comp = {"P": 1}
        assert sum(c for sp, c in r_comp.items()) == 2  # Exactly 2 reactants total
        assert sum(c for sp, c in p_comp.items()) == 1  # Exactly 1 product total



def test_solve_lhs():
    net = RandomReactionNetwork(n_species=3, n_rxns=0)

    with pytest.raises(Exception):
        net.solve_lhs(a_list=[], alpha=10, beta=20)


    # Length 1
    result = net.solve_lhs(a_list=[1], alpha=10, beta=10)
    assert result is None

    result = net.solve_lhs(a_list=[1], alpha=10, beta=11)
    assert result == (1,)

    result = net.solve_lhs(a_list=[2], alpha=10, beta=11)
    assert result is None

    result = net.solve_lhs(a_list=[2], alpha=10, beta=12)
    assert result == (1,)

    result = net.solve_lhs(a_list=[2], alpha=10, beta=13)
    assert result is None

    result = net.solve_lhs(a_list=[2], alpha=10, beta=20)
    assert result == (5,)

    result = net.solve_lhs(a_list=[3], alpha=10, beta=22)
    assert result == (4,)

    result = net.solve_lhs(a_list=[3], alpha=10, beta=23)
    assert result is None


    # Length 2
    result = net.solve_lhs(a_list=[1,1], alpha=10, beta=10)
    assert result is None

    with pytest.raises(Exception):
        net.solve_lhs(a_list=[3,1], alpha=10, beta=20)

    with pytest.raises(Exception):
        net.solve_lhs(a_list=[1,2], alpha=10, beta=20)

    result = net.solve_lhs(a_list=[1,1], alpha=10, beta=11)
    assert result is None

    result = net.solve_lhs(a_list=[1,1], alpha=10, beta=12)
    assert result == (1,1)

    result = net.solve_lhs(a_list=[1,1], alpha=10, beta=20)
    assert result == (5,5)

    result = net.solve_lhs(a_list=[1,1], alpha=10, beta=25)
    assert result == (7,8)



def test_solve_both_sides():
    net = RandomReactionNetwork(n_species=3, n_rxns=0)

    with pytest.raises(Exception):
        net.solve_both_sides(a_list=[], alpha=10, b_list=[1], beta=20)

    with pytest.raises(Exception):
        net.solve_both_sides(a_list=[1,1], alpha=10, b_list=[], beta=20)

    with pytest.raises(Exception):
        net.solve_both_sides(a_list=[1], alpha=10, b_list=[1], beta=20)

    with pytest.raises(Exception):
        net.solve_both_sides(a_list=[1,1], alpha=10, b_list=[2], beta=20)


    # `a_list` of size 1 : only allowed value is [2]
    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=20)
    assert result == ((11,), (12,))     # Both sides add up to 32

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=11)
    assert result == ((2,), (3,))       # Both sides add up to 14

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=10)
    assert result == ((1,), (2,))       # Both sides add up to 12

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=9)
    assert result == ((2,), (5,))       # Both sides add up to 12

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=8)
    assert result == ((3,), (8,))       # Both sides add up to 16

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=0)
    assert result == ((11,), (32,))       # Both sides add up to 32


    # `a_list` of size 1 : only allowed value is [1,1]
    result = net.solve_both_sides(a_list=[1,1], alpha=10, b_list=[1], beta=20)
    assert result == ((10, 10), (10,))       # Both sides add up to 30

    result = net.solve_both_sides(a_list=[1,1], alpha=10, b_list=[1], beta=11)
    assert result == ((1, 1), (1,))       # Both sides add up to 12

    result = net.solve_both_sides(a_list=[1,1], alpha=10, b_list=[1], beta=10)
    assert result == ((16, 16), (32,))       # Both sides add up to 42

    result = net.solve_both_sides(a_list=[1,1], alpha=10, b_list=[1], beta=9)
    assert result == ((1, 1), (3,))       # Both sides add up to 12

    result = net.solve_both_sides(a_list=[1,1], alpha=10, b_list=[1], beta=0)
    assert result == ((10, 10), (30,))       # Both sides add up to 30

    result = net.solve_both_sides(a_list=[1,1], alpha=0, b_list=[1], beta=0)
    assert result == ((16, 16), (32,))       # Both sides add up to 32
    #print(result)



def test___fulfill_reaction():
    net = RandomReactionNetwork(n_species=3, n_rxns=0)
    # Will initialize 3 species: A, B, C

    stoic = Stoichiometry(vector={"A": -1, "B": -1, "C": 1})    # A + B -> C , all values unset
    result = net._fulfill_reaction(stoic)
    assert result
    check = sum(coeff * int(net.species_registry.get_species(sp).annotation)
                    for sp, coeff in stoic.vector.items())
    assert check == 0   # Atom accounting is balanced

    result = net._fulfill_reaction(stoic)       # Same reaction, but with the previously set atom counts
    assert result
    check = sum(coeff * int(net.species_registry.get_species(sp).annotation)
                    for sp, coeff in stoic.vector.items())
    assert check == 0   # Atom accounting is balanced

    stoic = Stoichiometry(vector={"A": 1, "B": -1, "C": -1})    # B + C -> A, with the previously set atom counts
    result = net._fulfill_reaction(stoic)
    assert not result   # Incompatible counts
