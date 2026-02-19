import pytest
import numpy as np
from life123 import RandomReactionNetwork




def test_constructor():

    with pytest.raises(Exception):
        RandomReactionNetwork()     # Missing arg

    with pytest.raises(Exception):
        RandomReactionNetwork(n_chems="not_an_int!", n_rxns=3)  # Bad arg

    with pytest.raises(Exception):
        RandomReactionNetwork(n_chems=2, n_rxns=3)  # Too few chems

    with pytest.raises(Exception):
        RandomReactionNetwork(n_chems=3, n_rxns=[1,3])

    with pytest.raises(Exception):
        RandomReactionNetwork(n_chems=3, n_rxns=0)

    for _ in range(25):
        # In some cases, the instantiation of RandomReactionNetwork will give up, because of
        # the difficulty of creating so many unique reactions from so few chemicals
        try:
            net = RandomReactionNetwork(n_chems=3, n_rxns=6, relative_rxn_prob=[0.3, 0.7])
            break
        except Exception:
            pass

    assert net.chem_data.number_of_chemicals() == 3
    assert net.chem_data.get_all_labels() == ['A', 'B', 'C']

    # Verify that all possible reactions (or their reverse ones) were indeed generated,
    # with no consideration for stoichiometry
    assert net.already_used(reactants="A", products=["B", "C"])
    assert net.already_used(reactants="B", products=["A", "C"])
    assert net.already_used(reactants="C", products=["A", "B"])
    assert net.already_used(reactants="A", products="B")    # This covers A <--> 2 B as well as 2 A <-> B
    assert net.already_used(reactants="A", products="C")
    assert net.already_used(reactants="B", products="C")



def test__assign_chems_to_rxn():
    net = RandomReactionNetwork(n_chems=3, n_rxns=1, seed=1)
    all_chem_labels = net.chem_data.get_all_labels()
    assert all_chem_labels == ["A", "B", "C"]

    for _ in range(20):
        result = net._assign_chems_to_rxn("ReactionUnimolecular", all_chem_labels=all_chem_labels)
        assert len(reactants := result[0]) == 1
        r1 = reactants[0]
        assert r1 in all_chem_labels

        assert len(products := result[1]) == 1
        p1 = products[0]
        assert p1 in all_chem_labels

        assert r1 != p1


        result = net._assign_chems_to_rxn("ReactionSynthesis", all_chem_labels=all_chem_labels)
        assert len(reactants := result[0]) == 2
        r1, r2 = reactants
        assert r1 in all_chem_labels
        assert r2 in all_chem_labels

        assert len(products := result[1]) == 1
        p1 = products[0]
        assert p1 in all_chem_labels

        assert r1 != p1
        assert r2 != p1


        result = net._assign_chems_to_rxn("ReactionDecomposition", all_chem_labels=all_chem_labels)
        assert len(reactants := result[0]) == 1
        r1 = reactants[0]
        assert r1 in all_chem_labels

        assert len(products := result[1]) == 2
        p1, p2 = products
        assert p1 in all_chem_labels
        assert p2 in all_chem_labels

        assert r1 != p1
        assert r1 != p2



def test_already_used():
    net = RandomReactionNetwork(n_chems=3, n_rxns=1, seed=1)
    all_chem_labels = net.chem_data.get_all_labels()
    assert all_chem_labels == ["A", "B", "C"]

    assert not net.already_used(reactants="A", products="B")
    assert not net.already_used(reactants="X", products="Y")

    net.registry.add_elementary_reaction(reactants="A", products="B")  # Reaction 0 : A <-> B
    assert net.already_used(reactants="A", products="B")
    assert net.already_used(reactants="B", products="A")

    assert net.already_used(reactants="A", products="B")
    assert net.already_used(reactants="B", products="A")
    assert not net.already_used(reactants=["A", "B"], products="C")

    net.registry.add_elementary_reaction(reactants=["A", "B"], products="C")  # Reaction 1 : A + B <-> C
    assert net.already_used(reactants=["A", "B"], products="C")
    assert net.already_used(reactants=["B", "A"], products="C")
    assert net.already_used(reactants="C", products=["A", "B"])

    assert not net.already_used(reactants=["A", "C"], products="B")

    net.registry.chem_data.add_chemical(name="D")
    assert not net.already_used(reactants=["C", "D"], products="B")

    net.registry.add_elementary_reaction(reactants=["D", "C"], products="B")  # Reaction 2 : D + C <-> B
    assert net.already_used(reactants=["D", "C"], products="B")
    assert net.already_used(reactants=["C", "D"], products="B")
    assert net.already_used(reactants="B", products=["C", "D"])

    assert not net.already_used(reactants="C", products=["D", "A"])
    net.registry.add_elementary_reaction(reactants="C", products=["D", "A"])  # Reaction 3 : C <-> D + A
    assert net.already_used(reactants="C", products=["D", "A"])
    assert net.already_used(reactants=["A", "D"], products="C")

    net.registry.add_elementary_reaction(reactants="D", products=["B", "B"])  # Reaction 3 : D <-> 2 B
    assert net.already_used(reactants="D", products=["B"])
    assert net.already_used(reactants="D", products=["B", "B"])



def test_random_species_enthalpy():
    net = RandomReactionNetwork(n_chems=3, n_rxns=0, seed=1492042)
    H = net.random_species_enthalpy(sigma=40.4145, n=3)
    np.allclose(H, [-12.67475766, 6.06243963, -6.70292682])

    # Produce a linear combination of the form X1 + X2 - X3 from 3 of our normal distributions
    linear_combo = [np.sum(net.random_species_enthalpy(sigma=40.4145, n=2)) -
                           net.random_species_enthalpy(sigma=40.4145, n=1)
                        for _ in range(100)]
    np.allclose(np.mean(linear_combo), 0.8995432)   # This will approach 0
    np.allclose(np.std(linear_combo), 69.654486)    # This will approach 70



def test_random_reaction_enthalpy():
    net = RandomReactionNetwork(n_chems=3, n_rxns=0, seed=1492042)

    result_1 = net.random_reaction_enthalpy(reactants=["A", "B"], products="C")
    print(net.standard_species_enthalpy)

    assert np.allclose(net.standard_species_enthalpy["A"], -12.674757661590291)
    assert np.allclose(net.standard_species_enthalpy["B"], 6.06243962532564)
    assert np.allclose(net.standard_species_enthalpy["C"], -6.702926824946904)

    assert np.allclose(result_1, -0.09060878868225242)
    assert np.allclose(result_1,
                       net.get_species_enthalpy("C")
                       - net.get_species_enthalpy("A") - net.get_species_enthalpy("B"))
    # ΔH = HC − HA − HB

    result_2 = net.random_reaction_enthalpy(reactants="C", products=["A", "B"])
    assert np.allclose(result_2, -result_1)     # The reverse reaction has opposite delta_H


    result_1 = net.random_reaction_enthalpy(reactants=["A", "A"], products="C")
    assert np.allclose(result_1, 18.64658849823368)
    assert np.allclose(result_1,
                       net.get_species_enthalpy("C")
                       - 2 * net.get_species_enthalpy("A"))
    # ΔH = HC - 2 HA

    result_2 = net.random_reaction_enthalpy(reactants="C", products=["A", "A"])
    assert np.allclose(result_2, -result_1)     # The reverse reaction has opposite delta_H



def test_random_reaction_enthalpy_2():
    net = RandomReactionNetwork(n_chems=35, n_rxns=65, seed=1492042)

    H_list = [rxn.delta_H
                for rxn in net.registry.get_all_reactions()]

    #print(np.mean(H_list))
    #print(np.std(H_list))
    np.allclose(np.mean(H_list), 3.6438823)   # This will approach 0
    np.allclose(np.std(H_list), 72.886196)    # Roughly around 70
