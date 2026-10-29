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



def test_random_reaction_enthalpy():
    net = RandomReactionNetwork(n_species=3, n_rxns=0, seed=1492042)

    result_1 = net.random_reaction_enthalpy(reaction_vector={"A": -1, "B": -1, "C": 1})
    #print(net.standard_species_enthalpy)

    assert np.allclose(net.standard_species_enthalpy["A"], -12.674757661590291)
    assert np.allclose(net.standard_species_enthalpy["B"], 6.06243962532564)
    assert np.allclose(net.standard_species_enthalpy["C"], -6.702926824946904)

    assert np.allclose(result_1, -0.09060878868225242)
    assert np.allclose(result_1,
                       net.get_species_enthalpy("C")
                       - net.get_species_enthalpy("A") - net.get_species_enthalpy("B"))
    # ΔH = HC − HA − HB


    result_2 = net.random_reaction_enthalpy(reaction_vector={"C": -1, "A": 1, "B": 1})
    assert np.allclose(result_2, -result_1)     # The reverse reaction has opposite delta_H


    result_1 = net.random_reaction_enthalpy(reaction_vector={"A": -2, "C": 1})
    assert np.allclose(result_1, 18.64658849823368)
    assert np.allclose(result_1,
                       net.get_species_enthalpy("C")
                       - 2 * net.get_species_enthalpy("A"))
    # ΔH = HC - 2 HA

    result_2 = net.random_reaction_enthalpy(reaction_vector={"C": -1, "A": 2})
    assert np.allclose(result_2, -result_1)     # The reverse reaction has opposite delta_H



def test_assign_thermodynamics():
    net = RandomReactionNetwork(n_species=300, n_rxns=0, seed=8492342, verbose=False)

    net.generate_multiple_random_reactions(n_rxns=150)

    net.assign_thermodynamics()

    # Extract the delta_H values for all the reactions
    H_list = [rxn.get_thermodynamics().delta_H
              for rxn in net.get_reaction_data().get_all_reactions()]

    assert len(H_list) == 150
    #print(H_list)

    #print(np.mean(H_list))
    #print(np.std(H_list))
    # Note: the values are subject to change, if changes are made in RandomReactionNetwork's usage of random calls
    assert np.allclose(np.mean(H_list), -4.217131270343261)  # This will approach 0
    assert np.allclose(np.std(H_list),  76.38392872589338)   # Roughly around 70 - 90 (much closer to 70)



def test_random_species_enthalpy():
    net = RandomReactionNetwork(n_species=3, n_rxns=0, seed=1492042)
    H = net.random_species_enthalpy(sigma=40.4145, n=3)
    assert np.allclose(H, [-12.67475766, 6.06243963, -6.70292682])

    # Produce a linear combination of the form X1 + X2 - X3 from 3 of our normal distributions
    linear_combo = [np.sum(net.random_species_enthalpy(sigma=40.4145, n=2)) -
                           net.random_species_enthalpy(sigma=40.4145, n=1)
                        for _ in range(100)]
    np.allclose(np.mean(linear_combo), 0.8995432)   # This will approach 0
    np.allclose(np.std(linear_combo), 69.654486)    # This will approach 70

    H = net.random_species_enthalpy(sigma=35., n=1)
    assert type(H) == np.ndarray
    assert np.allclose(H, [41.60853089])

    H = net.random_species_enthalpy(sigma=35., n=None)
    assert type(H) == float
    assert np.allclose(H, -52.64502861779662)



def test_generate_multiple_random_reactions():

    net = RandomReactionNetwork(n_species=4, n_rxns=0, seed=1234)
    net.generate_multiple_random_reactions(0)
    assert net.reaction_registry.number_of_reactions() == 0


    net = RandomReactionNetwork(n_species=4, n_rxns=0, seed=1234)
    net.generate_multiple_random_reactions(1)
    assert net.reaction_registry.number_of_reactions() == 1
    rxns = net.reaction_registry.get_all_reactions()
    r = rxns[0]
    assert r.stoichiometry == Stoichiometry(vector={'D': -2, 'C': 1}, species_registry=net.species_registry)   # 2 D -> C
    assert net.species_registry.get_species('C').annotation == "32"
    assert net.species_registry.get_species('D').annotation == "16"


    net = RandomReactionNetwork(n_species=4, n_rxns=0, seed=1234)
    net.generate_multiple_random_reactions(2)
    assert net.reaction_registry.number_of_reactions() == 2
    rxns = net.reaction_registry.get_all_reactions()
    r = rxns[0]
    assert r.stoichiometry == Stoichiometry(vector={'D': -2, 'C': 1}, species_registry=net.species_registry)           # 2 D -> C
    r = rxns[1]
    assert r.stoichiometry == Stoichiometry(vector={'B': -1, 'A': -1, 'D': 1}, species_registry=net.species_registry)  # B + A -> D
    assert net.species_registry.get_species('A').annotation == "8"
    assert net.species_registry.get_species('B').annotation == "8"
    assert net.species_registry.get_species('C').annotation == "32"
    assert net.species_registry.get_species('D').annotation == "16"


    net = RandomReactionNetwork(n_species=4, n_rxns=0, seed=1234)
    net.generate_multiple_random_reactions(3)
    assert net.reaction_registry.number_of_reactions() == 3
    rxns = net.reaction_registry.get_all_reactions()
    r = rxns[0]
    assert r.stoichiometry == Stoichiometry(vector={'D': -2, 'C': 1}, species_registry=net.species_registry)           # 2 D -> C
    r = rxns[1]
    assert r.stoichiometry == Stoichiometry(vector={'B': -1, 'A': -1, 'D': 1}, species_registry=net.species_registry)  # B + A -> D
    r = rxns[2]
    assert r.stoichiometry == Stoichiometry(vector={'B': -2, 'D': 1}, species_registry=net.species_registry)           # 2 B -> D
    # All atom counts same as before
    assert net.species_registry.get_species('A').annotation == "8"
    assert net.species_registry.get_species('B').annotation == "8"
    assert net.species_registry.get_species('C').annotation == "32"
    assert net.species_registry.get_species('D').annotation == "16"


    net = RandomReactionNetwork(n_species=6, n_rxns=0, seed=987654321)
    net.generate_multiple_random_reactions(4)
    assert net.reaction_registry.number_of_reactions() == 4
    rxns = net.reaction_registry.get_all_reactions()
    assert rxns[0].stoichiometry.vector == {'E': -1, 'B': -1, 'F': 1}  # E + B -> F
    assert rxns[1].stoichiometry.vector == {'B': -1, 'F': -1, 'A': 1}  # B + F -> A
    assert rxns[2].stoichiometry.vector == {'B': -1, 'C': -1, 'F': 1}  # B + C -> F
    assert rxns[3].stoichiometry.vector == {'F': -1, 'C': -1, 'A': 1}  # F + C -> A

    assert net.species_registry.get_species('A').annotation == "48"
    assert net.species_registry.get_species('B').annotation == "16"
    assert net.species_registry.get_species('C').annotation == "16"
    assert net.species_registry.get_species('D').annotation == ""   # Not participating in any reaction
    assert net.species_registry.get_species('E').annotation == "16"
    assert net.species_registry.get_species('F').annotation == "32"

    assert net.species_registry.get_species('A').label == "A(48)"
    assert net.species_registry.get_species('D').label == "D"
    assert net.species_registry.get_species('F').label == "F(32)"

    """
    net = RandomReactionNetwork(n_species=16, n_rxns=40, seed=987654321)
    net.generate_multiple_random_reactions()
    #print(net.species_registry.as_recordset())
    net.reaction_registry.describe_reactions(detail="minimal", use_labels=True)
    #net.reaction_registry.describe_reactions(detail="summary", use_labels=True)
    """



def test_generate_new_consistent_unique_reaction():
    net = RandomReactionNetwork(n_species=4, n_rxns=0, seed=1234)

    st = net.generate_new_consistent_unique_reaction()
    assert st == Stoichiometry(vector={'D': -2, 'C': 1}, species_registry=net.species_registry) # 2 D -> C
    assert net.check_atom_count_consistency(st)
    net.reaction_registry.add_elementary_reaction(reactants=[(2, "D")], products="C")

    st = net.generate_new_consistent_unique_reaction()
    assert st.vector == {'B': -1, 'A': -1, 'D': 1}      # B + A -> D
    assert net.check_atom_count_consistency(st)
    net.reaction_registry.add_elementary_reaction(reactants=["B", "A"], products="D")

    st = net.generate_new_consistent_unique_reaction()
    assert st.vector == {'B': -2, 'D': 1}               # 2 B -> D
    assert net.check_atom_count_consistency(st)
    net.reaction_registry.add_elementary_reaction(reactants=[(2, "B")], products="D")

    st = net.generate_new_consistent_unique_reaction()
    assert st.vector == {'A': -2, 'D': 1}           # 2 A <-> D
    assert net.check_atom_count_consistency(st)
    net.reaction_registry.add_elementary_reaction(reactants=[(2, "A")], products="D")


    with pytest.raises(Exception):
        # Repeated failures to find a new reaction consistent with atom-count assignments so far
        net.generate_new_consistent_unique_reaction()



def test_generate_new_unique_reaction():
    net = RandomReactionNetwork(n_species=4, n_rxns=0, seed=1234)

    st = net.generate_new_unique_reaction()
    assert st == Stoichiometry(vector={'D': -2, 'C': 1}, catalysts=[], species_registry=net.species_registry)  # 2 D -> C

    st = net.generate_new_unique_reaction()
    assert st.vector == {'B': -1, 'A': -1, 'D': 1}      # B + A -> D

    st = net.generate_new_unique_reaction()
    assert st.vector == {'A': -1, 'B': -1, 'C': 1}      # A + B -> C



def test__assign_species_to_rxn():
    net = RandomReactionNetwork(n_species=4, n_rxns=0, seed=123456)

    reactants, products = net._assign_species_to_rxn()
    # Generates A + C -> D
    assert reactants == ['A', 'C']
    assert products == ['D']


    reactants, products = net._assign_species_to_rxn()
    # Generates 2 B -> A
    assert reactants == ['B', 'B']
    assert products == ['A']

    reactants, products = net._assign_species_to_rxn()
    # Generates  A + D -> B
    assert reactants == ['A', 'D']
    assert products == ['B']


    for i in range(20):
        reactants, products = net._assign_species_to_rxn()

        # 2 reactants (same or different species), 1 product
        assert len(reactants) == 2
        assert len(products) == 1




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
        net.solve_both_sides(a_list=[1,1], alpha=10, b_list=[2], beta=20)


    # `a_list` of size 1
    # Solving a0 * x0 + alpha = b0 * y0 + beta
    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=20)
    assert result == ((26,), (42,))     # Both sides add up to 62
    assert 2 * 26 + 10 == 1 * 42 + 20

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=11)
    assert result == ((17,), (33,))       # Both sides add up to 44

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=10)
    assert result == ((16,), (32,))

    result = net.solve_both_sides(a_list=[2], alpha=0, b_list=[1], beta=0)
    assert result == ((16,), (32,))

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=9)
    assert result == ((17,), (35,))

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=8)
    assert result == ((18,), (38,))

    result = net.solve_both_sides(a_list=[2], alpha=10, b_list=[1], beta=0)
    assert result == ((26,), (62,))


    result = net.solve_both_sides(a_list=[1], alpha=10, b_list=[1], beta=20)
    assert result == ((26,), (16,))     # Both sides add up to 36

    result = net.solve_both_sides(a_list=[1], alpha=16, b_list=[1], beta=0)
    assert result == ((32,), (48,))     # Both sides add up to 48


    # `a_list` of size 2 : only allowed value is [1,1]
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
    assert net.check_atom_count_consistency(stoic)  # Atom accounting is balanced
    """
    check = sum(coeff * int(net.species_registry.get_species(sp).annotation)
                    for sp, coeff in stoic.vector.items())
    assert check == 0   # Atom accounting is balanced
    """

    stoic = Stoichiometry(vector={"A": 1, "B": -1, "C": -1})    # B + C -> A, with the previously set atom counts
    result = net._fulfill_reaction(stoic)
    assert not result   # Incompatible counts



def test__already_used():
    net = RandomReactionNetwork(n_species=3, n_rxns=1, seed=1)
    all_chem_labels = net.species_registry.get_all_species_ids()
    assert all_chem_labels == ["A", "B", "C"]

    assert not net._already_used(reactants="A", products="B")
    assert not net._already_used(reactants="X", products="Y")

    net.reaction_registry.add_elementary_reaction(reactants="A", products="B")  # Reaction 0 : A <-> B
    assert net._already_used(reactants="A", products="B")
    assert net._already_used(reactants="B", products="A")

    assert net._already_used(reactants="A", products="B")
    assert net._already_used(reactants="B", products="A")
    assert not net._already_used(reactants=["A", "B"], products="C")

    net.reaction_registry.add_elementary_reaction(reactants=["A", "B"], products="C")  # Reaction 1 : A + B <-> C
    assert net._already_used(reactants=["A", "B"], products="C")
    assert net._already_used(reactants=["B", "A"], products="C")
    assert net._already_used(reactants="C", products=["A", "B"])

    assert not net._already_used(reactants=["A", "C"], products="B")

    net.reaction_registry.species_registry.add_single_species(species_id="D")
    assert not net._already_used(reactants=["C", "D"], products="B")

    net.reaction_registry.add_elementary_reaction(reactants=["D", "C"], products="B")  # Reaction 2 : D + C <-> B
    assert net._already_used(reactants=["D", "C"], products="B")
    assert net._already_used(reactants=["C", "D"], products="B")
    assert net._already_used(reactants="B", products=["C", "D"])

    assert not net._already_used(reactants="C", products=["D", "A"])
    net.reaction_registry.add_elementary_reaction(reactants="C", products=["D", "A"])  # Reaction 3 : C <-> D + A
    assert net._already_used(reactants="C", products=["D", "A"])
    assert net._already_used(reactants=["A", "D"], products="C")

    net.reaction_registry.add_elementary_reaction(reactants="D", products=["B", "B"])  # Reaction 3 : D <-> 2 B
    assert net._already_used(reactants="D", products=["B"])
    assert net._already_used(reactants="D", products=["B", "B"])
    assert net._already_used(reactants="B", products=["D", "D"])
