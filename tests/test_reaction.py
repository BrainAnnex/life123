import pytest
import numpy as np
from src.modules.reactions.reaction_data import ReactionData
from src.modules.reactions.reaction import Reaction



def test_initialize():
    chem_data = ReactionData(names=['A', 'B', 'C'])

    rxn = Reaction(chem_data)

    assert rxn.reactants is None
    assert rxn.products is None
    assert rxn.kF is None
    assert rxn.kR is None
    assert rxn.Delta_H is None
    assert rxn.Delta_S is None
    assert rxn.Delta_G is None
    assert rxn.K is None
    assert rxn.enzymes is None



def test_define_reaction():
    chem_data = ReactionData(names=["A", "B", "C", "D", "E", "F"])
    chem_data.set_temp(None)

    rxn = Reaction(chem_data)

    # Reactants and the products can't be the same
    with pytest.raises(Exception):
        rxn.define_reaction(reactants=["A"], products=["A"])
    with pytest.raises(Exception):
        rxn.define_reaction(reactants=["A"], products=[("A")])
    with pytest.raises(Exception):
        rxn.define_reaction(reactants=["A"], products=[(1, "A")])
    with pytest.raises(Exception):
        rxn.define_reaction(reactants=["A"], products=[(1, "A", 1)])
    with pytest.raises(Exception):
        rxn.define_reaction(reactants=[(2, "B")], products=[(2, "B")])
    with pytest.raises(Exception):
        rxn.define_reaction(reactants=[(2, "B")], products=[(2, "B", 1)])
    with pytest.raises(Exception):
        rxn.define_reaction(reactants=["A", "B"], products=["A", "B"])
    with pytest.raises(Exception):
        rxn.define_reaction(reactants=["A", (3, "B")], products=["A", (3, "B")])
    with pytest.raises(Exception):
        rxn.define_reaction(reactants=["A", "B"], products=["B", "A"])


    rxn.define_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert np.allclose(rxn.extract_forward_rate() , 3.)
    assert np.allclose(rxn.extract_reverse_rate() , 2.)
    assert np.allclose(rxn.K , 3./2.)
    assert rxn.extract_reactants() == [(1, 0, 1)]
    assert rxn.extract_products() == [(1, 1, 1)]
    assert rxn.Delta_H is None
    assert rxn.Delta_S is None
    assert rxn.Delta_G is None

    assert rxn.extract_reactants() == [(1, 0, 1)]
    assert rxn.extract_reactants_formula() == "A"

    assert rxn.extract_products() == [(1, 1, 1)]
    assert rxn.extract_products_formula() == "B"
 

    # Another reaction
    rxn.define_reaction(reactants=[(2, "B")], products=[(5, "C")], forward_rate=9., reverse_rate=7.)

    assert np.allclose(rxn.extract_forward_rate() , 9.)
    assert np.allclose(rxn.extract_reverse_rate() , 7.)
    assert np.allclose(rxn.K , 9./7.)
    assert rxn.extract_reactants() == [(2, 1, 1)]
    assert rxn.extract_products() == [(5, 2, 1)]
    assert rxn.Delta_H is None
    assert rxn.Delta_S is None
    assert rxn.Delta_G is None


    # Add another reaction.  This time, first set the temperature
    chem_data.set_temp(200)
    rxn.define_reaction(reactants=[(2, "D", 3)], products=[(1, "C", 2)], forward_rate=11., reverse_rate=13.)

    assert np.allclose(rxn.extract_forward_rate() , 11.)
    assert np.allclose(rxn.extract_reverse_rate() , 13.)
    assert np.allclose(rxn.K , 11./13.)
    assert rxn.extract_reactants() == [(2, 3, 3)]
    assert rxn.extract_products() == [(1, 2, 2)]
    assert rxn.Delta_H is None
    assert rxn.Delta_S is None
    assert np.allclose(rxn.Delta_G, 277.7928942715384)   # - RT log(K)


    # Add a multi-term reaction
    rxn.define_reaction(reactants=["A", (2, "B")], products=[(3, "C", 2), "D"], forward_rate=5., reverse_rate=1.)

    assert np.allclose(rxn.extract_forward_rate() , 5.)
    assert np.allclose(rxn.extract_reverse_rate() , 1.)
    assert np.allclose(rxn.K , 5./1.)
    assert rxn.extract_reactants() == [(1, 0, 1), (2, 1, 1)]
    assert rxn.extract_products() == [(3, 2, 2), (1, 3, 1)]
    assert rxn.Delta_H is None
    assert rxn.Delta_S is None
    assert np.allclose(rxn.Delta_G, -2676.321364705849)   # - RT log(K)


    # Add a reaction with thermodynamic data;
    # the reverse reaction rate will get computed from the thermodynamic data
    rxn = Reaction(chem_data)
    rxn.define_reaction(reactants=["A"], products=[(2, "B")], forward_rate=10.,
                        Delta_H = 5., Delta_S = 0.4)

    assert rxn.extract_reactants() == [(1, 0, 1)]
    assert rxn.extract_products() == [(2, 1, 1)]
    assert rxn.Delta_H == 5.
    assert rxn.Delta_S == 0.4
    assert np.allclose(rxn.Delta_G, -75.0)         # 5 - 200 * 0.4
    assert np.allclose(rxn.K, 1.0461347154679432)  # exp(75/(8.3144598 * 200))
    assert np.allclose(rxn.extract_forward_rate() , 10.)
    assert np.allclose(rxn.extract_reverse_rate() , 9.558998331803693) # 10. / 1.0461347154679432



def test_describe():
    chem_data = ReactionData(names=["CH4", "O2", "CO2", "H2O"])

    rxn = Reaction(chem_data)

    rxn.define_reaction(reactants=["CH4", (2, "O2")], products=["CO2", (2, "H2O")])

    assert rxn.describe(concise=True) == "CH4 + 2 O2 <-> CO2 + 2 H2O"
    assert rxn.describe(concise=False) == "0: CH4 + 2 O2 <-> CO2 + 2 H2O  () | 1st order in all reactants & products"   # TODO: eliminate the ()


    # Start over
    chem_data = ReactionData(names=["A", "B", "C", "D", "E", "F"])
    chem_data.set_temp(None)

    rxn = Reaction(chem_data)
    rxn.define_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    assert rxn.describe(concise=True) == "A <-> B"
    assert rxn.describe(concise=False) == "0: A <-> B  (kF = 3 / kR = 2 / K = 1.5) | 1st order in all reactants & products"

    rxn = Reaction(chem_data)
    rxn.define_reaction(reactants=[(2, "B")], products=[(5, "C")], forward_rate=9., reverse_rate=7.)
    assert rxn.describe(concise=True) == "2 B <-> 5 C"
    assert rxn.describe(concise=False, rxn_index=1) == "1: 2 B <-> 5 C  (kF = 9 / kR = 7 / K = 1.28571) | 1st order in all reactants & products"

    rxn = Reaction(chem_data)
    chem_data.set_temp(200)
    rxn.define_reaction(reactants=[(2, "D", 3)], products=[(1, "C", 2)], forward_rate=11., reverse_rate=13.)
    assert rxn.describe(concise=True) == "2 D <-> C"
    assert rxn.describe(concise=False, rxn_index=2) == "2: 2 D <-> C  (kF = 11 / kR = 13 / Delta_G = 277.793 / K = 0.846154) | 3-th order in reactant D | 2-th order in product C"

    rxn = Reaction(chem_data)
    rxn.define_reaction(reactants=["A", (2, "B")], products=[(3, "C", 2), "D"], forward_rate=5., reverse_rate=1.)
    assert rxn.describe(concise=True) == "A + 2 B <-> 3 C + D"
    assert rxn.describe(concise=False, rxn_index=3) == "3: A + 2 B <-> 3 C + D  (kF = 5 / kR = 1 / Delta_G = -2,676.32 / K = 5) | 2-th order in product C"



def test__standard_form_chem_eqn():
    chem_data = ReactionData(names=['Fe', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'Cl'])

    rxn = Reaction(chem_data)

    assert rxn._standard_form_chem_eqn([(1, 0, 1), (2, 8, 1)]) == "Fe + 2 Cl"

    assert rxn._standard_form_chem_eqn([(3, 0, 1), (5, 7, 2)]) == "3 Fe + 5 G"

