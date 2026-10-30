import pytest
import numpy as np
from life123.random_reactions import RandomReactionNetwork



def test_constructor():




    for _ in range(25):
        # In some cases, the instantiation of RandomReactionNetwork will give up, because of
        # the difficulty of creating so many unique reactions from so few chemicals
        try:
            net = RandomReactionNetwork(n_species=3, n_rxns=6, relative_rxn_prob=[0.3, 0.7])
            break
        except Exception:
            pass

    assert net.species_registry.number_of_species() == 3
    assert net.species_registry.get_all_species_ids() == ['A', 'B', 'C']

    # Verify that all possible reactions (or their reverse ones) were indeed generated,
    # with no consideration for stoichiometry
    assert net._already_used(reactants="A", products=["B", "C"])
    assert net._already_used(reactants="B", products=["A", "C"])
    assert net._already_used(reactants="C", products=["A", "B"])
    assert net._already_used(reactants="A", products="B")    # This covers A <--> 2 B as well as 2 A <-> B
    assert net._already_used(reactants="A", products="C")
    assert net._already_used(reactants="B", products="C")

