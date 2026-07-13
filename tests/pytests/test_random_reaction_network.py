import pytest
import numpy as np
from life123 import RandomReactionNetwork




def test_constructor():

    with pytest.raises(Exception):
        RandomReactionNetwork()     # Missing arg

    with pytest.raises(Exception):
        RandomReactionNetwork(n_species="not_an_int!", n_rxns=3)  # Bad arg

    with pytest.raises(Exception):
        RandomReactionNetwork(n_species=2, n_rxns=3)  # Too few species

    with pytest.raises(Exception):
        RandomReactionNetwork(n_species=3, n_rxns=[1, 3])


    for _ in range(25):
        # In some cases, the instantiation of RandomReactionNetwork will give up, because of
        # the difficulty of creating so many unique reactions from so few chemicals
        try:
            net = RandomReactionNetwork(n_species=3, n_rxns=6, relative_rxn_prob=[0.3, 0.7])
            break
        except Exception:
            pass

    assert net.species_data.number_of_species() == 3
    assert net.species_data.get_all_species_ids() == ['A', 'B', 'C']

    # Verify that all possible reactions (or their reverse ones) were indeed generated,
    # with no consideration for stoichiometry
    assert net._already_used(reactants="A", products=["B", "C"])
    assert net._already_used(reactants="B", products=["A", "C"])
    assert net._already_used(reactants="C", products=["A", "B"])
    assert net._already_used(reactants="A", products="B")    # This covers A <--> 2 B as well as 2 A <-> B
    assert net._already_used(reactants="A", products="C")
    assert net._already_used(reactants="B", products="C")


def test_constructor_2():
    net = RandomReactionNetwork(n_species=5, n_rxns=25, seed=888)
    for rxn in net.reaction_data.get_all_reactions():
        reactants = rxn.extract_reactant_labels()
        products = rxn.extract_product_labels()
        '''
        print()
        print(reactants)
        print(products)
        '''
        assert set(reactants) & set(products) == set()      # The intersection is an empty set



def test__assign_chems_to_rxn():
    net = RandomReactionNetwork(n_species=3, n_rxns=1, seed=1)
    all_chem_labels = net.species_data.get_all_species_ids()
    assert all_chem_labels == ["A", "B", "C"]

    for _ in range(20):
        '''
        # "ReactionUnimolecular" is currently NOT supported
        result = net._assign_chems_to_rxn("ReactionUnimolecular", all_chem_labels=all_chem_labels)
        assert len(reactants := result[0]) == 1
        r1 = reactants[0]
        assert r1 in all_chem_labels

        assert len(products := result[1]) == 1
        p1 = products[0]
        assert p1 in all_chem_labels

        assert r1 != p1
        '''

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
    net = RandomReactionNetwork(n_species=3, n_rxns=1, seed=1)
    all_chem_labels = net.species_data.get_all_species_ids()
    assert all_chem_labels == ["A", "B", "C"]

    assert not net._already_used(reactants="A", products="B")
    assert not net._already_used(reactants="X", products="Y")

    net.reaction_data.add_elementary_reaction(reactants="A", products="B")  # Reaction 0 : A <-> B
    assert net._already_used(reactants="A", products="B")
    assert net._already_used(reactants="B", products="A")

    assert net._already_used(reactants="A", products="B")
    assert net._already_used(reactants="B", products="A")
    assert not net._already_used(reactants=["A", "B"], products="C")

    net.reaction_data.add_elementary_reaction(reactants=["A", "B"], products="C")  # Reaction 1 : A + B <-> C
    assert net._already_used(reactants=["A", "B"], products="C")
    assert net._already_used(reactants=["B", "A"], products="C")
    assert net._already_used(reactants="C", products=["A", "B"])

    assert not net._already_used(reactants=["A", "C"], products="B")

    net.reaction_data.species_data.add_species(id="D")
    assert not net._already_used(reactants=["C", "D"], products="B")

    net.reaction_data.add_elementary_reaction(reactants=["D", "C"], products="B")  # Reaction 2 : D + C <-> B
    assert net._already_used(reactants=["D", "C"], products="B")
    assert net._already_used(reactants=["C", "D"], products="B")
    assert net._already_used(reactants="B", products=["C", "D"])

    assert not net._already_used(reactants="C", products=["D", "A"])
    net.reaction_data.add_elementary_reaction(reactants="C", products=["D", "A"])  # Reaction 3 : C <-> D + A
    assert net._already_used(reactants="C", products=["D", "A"])
    assert net._already_used(reactants=["A", "D"], products="C")

    net.reaction_data.add_elementary_reaction(reactants="D", products=["B", "B"])  # Reaction 3 : D <-> 2 B
    assert net._already_used(reactants="D", products=["B"])
    assert net._already_used(reactants="D", products=["B", "B"])



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



def test_random_reaction_enthalpy():
    net = RandomReactionNetwork(n_species=3, n_rxns=0, seed=1492042)

    result_1 = net.random_reaction_enthalpy(reactants=["A", "B"], products="C")
    #print(net.standard_species_enthalpy)

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
    net = RandomReactionNetwork(n_species=15, n_rxns=25, seed=1492042, verbose=False)

    H_list = [rxn.delta_H
              for rxn in net.reaction_data.get_all_reactions()]

    #print(np.mean(H_list))
    #print(np.std(H_list))

    assert np.allclose(np.mean(H_list), 7.511580706834767)   # This will approach 0
    assert np.allclose(np.std(H_list), 76.07987788990208)    # Roughly around 70



def test_random_entropy():
    net = RandomReactionNetwork(n_species=3, n_rxns=0, seed=1492042)

    with pytest.raises(Exception):
        net.random_reaction_entropy("Unknown_rxn_type")

    assert np.allclose(net.random_reaction_entropy("ReactionDecomposition"), 33.72761872021661)
    assert np.allclose(net.random_reaction_entropy("ReactionSynthesis"), -36.99986904436495)

    entropy_list = [net.random_reaction_entropy("ReactionDecomposition") for i in range(100)]
    assert np.allclose(np.mean(entropy_list), 38.876713159925636)       # This will approach 40
    assert np.allclose(np.std(entropy_list), 19.606579641379344)        # This will approach 20
    assert min(entropy_list) >= 0
    assert max(entropy_list) <= 100

    entropy_list = [net.random_reaction_entropy("ReactionSynthesis") for i in range(100)]
    assert np.allclose(np.mean(entropy_list), -39.72751726151384)       # This will approach -40
    assert np.allclose(np.std(entropy_list), 18.754802472647917)        # This will approach 20
    assert min(entropy_list) >= -100
    assert max(entropy_list) <= 0



def test_forward_activation_gibbs_energy_normal():
    net = RandomReactionNetwork(n_species=3, n_rxns=0, seed=9538625)


    assert np.allclose(net.forward_activation_gibbs_energy_normal(delta_G=40), 69.40655722019214)

    assert np.allclose(net.forward_activation_gibbs_energy_normal(delta_G=20), 71.37605795182265)

    assert np.allclose(net.forward_activation_gibbs_energy_normal(delta_G=0), 77.09026832411178)

    assert np.allclose(net.forward_activation_gibbs_energy_normal(delta_G=-20), 71.46432473759808)

    assert np.allclose(net.forward_activation_gibbs_energy_normal(delta_G=-40), 74.86038387981017)

    for _ in range(20):
         delta_G = net.rng.normal(loc=0, scale=100)
         result = net.forward_activation_gibbs_energy_normal(delta_G=delta_G)
         assert result >= 0
         assert result >= delta_G



def test_forward_activation_gibbs_energy_BEP():
    net = RandomReactionNetwork(n_species=3, n_rxns=0, seed=5447061)

    assert net.forward_activation_gibbs_energy_BEP(delta_G=40, sigma=0) == 90
    assert np.allclose(net.forward_activation_gibbs_energy_BEP(delta_G=40, sigma=12), 79.5437172573611)

    assert net.forward_activation_gibbs_energy_BEP(delta_G=20, sigma=0) == 80
    assert np.allclose(net.forward_activation_gibbs_energy_BEP(delta_G=20, sigma=12), 75.03655557989254)

    assert net.forward_activation_gibbs_energy_BEP(delta_G=0, sigma=0) == 70
    assert np.allclose(net.forward_activation_gibbs_energy_BEP(delta_G=0, sigma=12), 44.18724275817719)

    assert net.forward_activation_gibbs_energy_BEP(delta_G=-20, sigma=0) == 60
    assert np.allclose(net.forward_activation_gibbs_energy_BEP(delta_G=-20, sigma=12), 53.11334188826283)

    assert net.forward_activation_gibbs_energy_BEP(delta_G=-40, sigma=0) == 50
    assert np.allclose(net.forward_activation_gibbs_energy_BEP(delta_G=-40, sigma=12), 31.658362710826765)

    for _ in range(20):
         delta_G = net.rng.normal(loc=0, scale=100)
         result = net.forward_activation_gibbs_energy_BEP(delta_G=delta_G, sigma=12)
         assert result >= 0
         assert result >= delta_G



def test_rate_constant_from_activation_gibbs_energyy():
    net = RandomReactionNetwork(n_species=3, n_rxns=0, seed=123)

    assert np.allclose(net.rate_constant_from_activation_gibbs_energy(activation_delta_G=50, temp=300), 12312.469127155957)
    assert np.allclose(net.rate_constant_from_activation_gibbs_energy(activation_delta_G=60, temp=300), 223.472703949188)
    assert np.allclose(net.rate_constant_from_activation_gibbs_energy(activation_delta_G=65.743427304, temp=300), 22.3472703949188)
    assert np.allclose(net.rate_constant_from_activation_gibbs_energy(activation_delta_G=70, temp=300), 4.056054792471754)
    assert np.allclose(net.rate_constant_from_activation_gibbs_energy(activation_delta_G=75, temp=300), 0.5464412521361471)
    assert np.allclose(net.rate_constant_from_activation_gibbs_energy(activation_delta_G=80, temp=300), 0.07361785215286852)
    assert np.allclose(net.rate_constant_from_activation_gibbs_energy(activation_delta_G=90, temp=300), 0.0013361723233277472)
    assert np.allclose(net.rate_constant_from_activation_gibbs_energy(activation_delta_G=100, temp=300), 2.4251678436906198e-05)
