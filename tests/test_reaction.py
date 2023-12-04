import pytest
import numpy as np
from src.modules.chemicals.chem_data import ChemData
from src.modules.reactions.reaction import Reaction



def test_initialize():
    chem_data = ChemData(names=['A', 'B'])

    rxn = Reaction(chem_data, reactants="A", products="B")

    assert rxn.reactants == [(1, 0, 1)]
    assert rxn.products == [(1, 1, 1)]
    assert rxn.kF is None
    assert rxn.kR is None
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None
    assert rxn.K is None
    assert rxn.enzyme is None


    # More complex scenarios
    chem_data = ChemData(names=["A", "B", "C", "D", "E", "F"])
    chem_data.set_temp(None)

    # Missing products or reactants
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=["A"], products=None)
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=None, products="B")

    # Reactants and the products can't be the same
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=["A"], products=["A"])
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=["A"], products=[("A")])
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=["A"], products=[(1, "A")])
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=["A"], products=[(1, "A", 1)])
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=[(2, "B")], products=[(2, "B")])
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=[(2, "B")], products=[(2, "B", 2)])
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=["A", "B"], products=["A", "B"])
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=["A", (3, "B")], products=["A", (3, "B")])
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=["A", "B"], products=["B", "A"])
    with pytest.raises(Exception):
        Reaction(chem_data, reactants=[(2, "A"), "B", "C"], products=["B", (1, "C"), (2, "A")])


    rxn = Reaction(chem_data, reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert np.allclose(rxn.extract_forward_rate() , 3.)
    assert np.allclose(rxn.extract_reverse_rate() , 2.)
    assert np.allclose(rxn.K , 3./2.)
    assert rxn.extract_reactants() == [(1, 0, 1)]
    assert rxn.extract_products() == [(1, 1, 1)]
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None

    assert rxn.extract_reactants() == [(1, 0, 1)]
    assert rxn.extract_reactants_formula() == "A"

    assert rxn.extract_products() == [(1, 1, 1)]
    assert rxn.extract_products_formula() == "B"
 

    # Another reaction
    rxn = Reaction(chem_data, reactants=[(2, "B", 1)], products=[(5, "C", 1)], forward_rate=9., reverse_rate=7.)

    assert np.allclose(rxn.extract_forward_rate() , 9.)
    assert np.allclose(rxn.extract_reverse_rate() , 7.)
    assert np.allclose(rxn.K , 9./7.)
    assert rxn.extract_reactants() == [(2, 1, 1)]
    assert rxn.extract_products() == [(5, 2, 1)]
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None


    # Add another reaction.  This time, first set the temperature
    chem_data.set_temp(200)
    rxn = Reaction(chem_data, reactants=[(2, "D", 3)], products=[(1, "C", 2)], forward_rate=11., reverse_rate=13.)

    assert np.allclose(rxn.extract_forward_rate() , 11.)
    assert np.allclose(rxn.extract_reverse_rate() , 13.)
    assert np.allclose(rxn.K , 11./13.)
    assert rxn.extract_reactants() == [(2, 3, 3)]
    assert rxn.extract_products() == [(1, 2, 2)]
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert np.allclose(rxn.delta_G, 277.7928942715384)   # - RT log(K)


    # Add a multi-term reaction
    rxn = Reaction(chem_data, reactants=["A", (2, "B", 1)], products=[(3, "C", 2), "D"], forward_rate=5., reverse_rate=1.)

    assert np.allclose(rxn.extract_forward_rate() , 5.)
    assert np.allclose(rxn.extract_reverse_rate() , 1.)
    assert np.allclose(rxn.K , 5./1.)
    assert rxn.extract_reactants() == [(1, 0, 1), (2, 1, 1)]
    assert rxn.extract_products() == [(3, 2, 2), (1, 3, 1)]
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert np.allclose(rxn.delta_G, -2676.321364705849)   # - RT log(K)


    # Add a reaction with thermodynamic data;
    # the reverse reaction rate will get computed from the thermodynamic data
    rxn = Reaction(chem_data, reactants=["A"], products=[(2, "B", 1)], forward_rate=10.,
                   delta_H= 5., delta_S= 0.4)

    assert rxn.extract_reactants() == [(1, 0, 1)]
    assert rxn.extract_products() == [(2, 1, 1)]
    assert rxn.delta_H == 5.
    assert rxn.delta_S == 0.4
    assert np.allclose(rxn.delta_G, -75.0)         # 5 - 200 * 0.4
    assert np.allclose(rxn.K, 1.0461347154679432)  # exp(75/(8.3144598 * 200))
    assert np.allclose(rxn.extract_forward_rate() , 10.)
    assert np.allclose(rxn.extract_reverse_rate() , 9.558998331803693) # 10. / 1.0461347154679432



def test_set_macro_enzyme():
    chem_data = ChemData(names=['A', 'B'])

    rxn = Reaction(chem_data, reactants="A", products="B")
    rxn.set_macro_enzyme("mm", 4)

    assert rxn.macro_enzyme == ("mm", 4)



def test_extract_reactants():
    chem_data = ChemData(names=["CH4", "O2", "CO2", "H2O"])
    rxn = Reaction(chem_data, reactants=["CH4", (2, "O2"), [2, "CH4", 1]],
                              products=["CO2", (2, "H2O")])
    assert rxn.extract_reactants() == [(1, 0, 1), (2, 1, 2), (2, 0, 1)]


def test_extract_reactants_formula():
    chem_data = ChemData(names=["CH4", "O2", "CO2", "H2O"])
    rxn = Reaction(chem_data, reactants=["CH4", (2, "O2"), [2, "CH4", 1]],
                              products=["CO2", (2, "H2O")])
    assert rxn.extract_reactants_formula() == "CH4 + 2 O2 + 2 CH4"



def test_extract_products():
    chem_data = ChemData(names=["CH4", "O2", "CO2", "H2O"])
    rxn = Reaction(chem_data, reactants=["CH4", (2, "O2")],
                              products=["CO2", (2, "H2O"), [2, "CH4", 1]])
    assert rxn.extract_products() == [(1, 2, 1), (2, 3, 2), (2, 0, 1)]


def test_extract_products_formula():
    chem_data = ChemData(names=["CH4", "O2", "CO2", "H2O"])
    rxn = Reaction(chem_data, reactants=["CH4", (2, "O2")],
                              products=["CO2", (2, "H2O"), [2, "CH4", 1]])
    assert rxn.extract_products_formula() == "CO2 + 2 H2O + 2 CH4"



#########  TO DESCRIBE THE DATA  #########

def test_describe():
    chem_data = ChemData(names=["CH4", "O2", "CO2", "H2O"])
    rxn = Reaction(chem_data, reactants=["CH4", (2, "O2", 1)], products=["CO2", (2, "H2O", 1)])

    assert rxn.describe(concise=True) == "CH4 + 2 O2 <-> CO2 + 2 H2O"
    assert rxn.describe(concise=False) == "CH4 + 2 O2 <-> CO2 + 2 H2O  () | 1st order in all reactants & products"   # TODO: eliminate the ()


    # Start over
    chem_data = ChemData(names=["A", "B", "C", "D", "E", "F"])
    chem_data.set_temp(None)

    rxn = Reaction(chem_data, reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    assert rxn.describe(concise=True) == "A <-> B"
    assert rxn.describe(concise=False) == "A <-> B  (kF = 3 / kR = 2 / K = 1.5) | 1st order in all reactants & products"

    rxn = Reaction(chem_data, reactants=[(2, "B", 1)], products=[(5, "C", 1)], forward_rate=9., reverse_rate=7.)
    assert rxn.describe(concise=True) == "2 B <-> 5 C"
    assert rxn.describe(concise=False) == "2 B <-> 5 C  (kF = 9 / kR = 7 / K = 1.2857) | 1st order in all reactants & products"

    chem_data.set_temp(200)
    rxn = Reaction(chem_data, reactants=[(2, "D", 3)], products=[(1, "C", 2)], forward_rate=11., reverse_rate=13.)
    assert rxn.describe(concise=True) == "2 D <-> C"
    assert rxn.describe(concise=False) == "2 D <-> C  (kF = 11 / kR = 13 / delta_G = 277.79 / K = 0.84615) | 3-th order in reactant D | 2-th order in product C"

    rxn = Reaction(chem_data, reactants=["A", (2, "B", 1)], products=[(3, "C", 2), "D"], forward_rate=5., reverse_rate=1.)
    assert rxn.describe(concise=True) == "A + 2 B <-> 3 C + D"
    assert rxn.describe(concise=False) == "A + 2 B <-> 3 C + D  (kF = 5 / kR = 1 / delta_G = -2,676.3 / K = 5) | 2-th order in product C"



def test_extract_rxn_properties():
    chem = ChemData(names=["A", "B"])
    rxn = Reaction(chem, reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn.extract_rxn_properties()
    #print(result)
    assert np.allclose(result["kF"] , 3.)
    assert np.allclose(result["kR"] , 2.)
    assert np.allclose(result["K"] , 1.5)
    assert np.allclose(result["delta_G"] , -1005.1305052750387)



def test_extract_chemicals_in_reaction():
    chem = ChemData(names=["A", "B"])

    rxn = Reaction(chem, reactants="A", products="B")
    assert rxn.extract_chemicals_in_reaction() == {0, 1}

    chem = ChemData(names=["A", "B", "C"])
    rxn = Reaction(chem, reactants=["A", "B"], products=["B", "C"]) # B (1) acts as enzyme
    assert rxn.extract_chemicals_in_reaction(exclude_enzyme=False) == {0, 1, 2}
    assert rxn.extract_chemicals_in_reaction(exclude_enzyme=True) == {0, 2}

    chem = ChemData(names=["A", "B", "C", "D", "E"])
    rxn = Reaction(chem, reactants=[(2, "D"), "C"], products=["C", (3, "E")]) # C (2) acts as enzyme
    assert rxn.extract_chemicals_in_reaction(exclude_enzyme=False) == {2, 3, 4}
    assert rxn.extract_chemicals_in_reaction(exclude_enzyme=True) == {3, 4}




#######  For ANALYSIS  #######

def test_reaction_quotient():
    chem = ChemData(names=["A", "B", "C", "D", "E", "F"])


    # Reaction : A <-> B
    rxn = Reaction(chem, reactants="A", products="B")
    c = {'A': 24., 'B': 36.}
    assert np.allclose(1.5, rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1.5, quotient)
    assert formula == '[B] / [A]'


    # Reaction : A <-> F
    rxn = Reaction(chem, reactants="A", products="F")
    c = {'A': 3., 'F': 33.}
    assert np.allclose(11., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(11., quotient)
    assert formula == '[F] / [A]'


    # Reaction : A <-> 3B
    rxn = Reaction(chem, reactants=["A"], products=[(3,"B",1)])   # 1st order
    c = {'A': 3., 'B': 12.}
    assert np.allclose(4., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(4., quotient)
    assert formula == '[B] / [A]'


    # Reaction 3:  2A <-> 3B
    rxn = Reaction(chem, reactants=[(2,"A",1)], products=[(3,"B",1)])   # 1st order
    c = {'A': 3., 'B': 12.}
    assert np.allclose(4., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(4., quotient)
    assert formula == '[B] / [A]'


    # Reaction :  A + B <-> C , with 1st-order kinetics for each species
    rxn = Reaction(chem, reactants=["A" , "B"], products=["C"])
    c = {'A': 3., 'B': 4., 'C': 12.}
    assert np.allclose(1., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1., quotient)
    assert formula == '[C] / ([A][B])'


    # Reaction :  A <-> 2C + D , with 1st-order kinetics for each species
    rxn = Reaction(chem, reactants=["A"], products=[(2, "C" ,1) , "D"])
    c = {'A': 2., 'C': 4., 'D': 8.}
    assert np.allclose(16., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(16., quotient)
    assert formula == '([C][D]) / [A]'


    # Reaction 6:  2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    rxn = Reaction(chem, reactants=[(2,"A",1) , (5,"B",1)], products=[(4,"C",1) , (3,"D",1)])
    c = {'A': 2., 'B': 1., 'C': 4., 'D': 8.}
    assert np.allclose(16., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(16., quotient)
    assert formula == '([C][D]) / ([A][B])'


    # Reaction :  2A <-> B , with 1st-order kinetics in both directions
    rxn = Reaction(chem, reactants=[(2, "A",1)], products=["B"])
    c = {'A': 4., 'B': 20.}
    assert np.allclose(5., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(5., quotient)
    assert formula == '[B] / [A]'


    # Reaction :  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction
    rxn = Reaction(chem, reactants=[(2, "A", 2)], products="B")
    c = {'A': 4., 'B': 20.}
    assert np.allclose(1.25, rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1.25, quotient)
    assert formula == '[B] / [A]^2 '


    # Reaction :  A <-> 2B , with 2nd-order kinetics in the reverse direction
    rxn = Reaction(chem, reactants="A", products=[(2, "B")])
    c = {'A': 4., 'B': 20.}
    assert np.allclose(100., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(100., quotient)
    assert formula == '[B]^2  / [A]'


    # Reaction :  A + B <-> C + D
    rxn = Reaction(chem, reactants=["A", "B"], products=["C", "D"])

    # with zero concentrations of a reaction product (and non-zero reactants), the reaction quotient will be zero
    c = {'A': 15.2, 'B': 21.3, 'C': 0 , 'D': 4.1}
    assert np.allclose(0., rxn.reaction_quotient(conc=c, explain=False))

    # with zero concentrations of a reactants (and non-zero products), it will be infinite
    c = {'A': 15.2, 'B': 0, 'C': 21.3 , 'D': 4.1}
    assert np.isinf(rxn.reaction_quotient(conc=c, explain=False))

    # with zero concentrations of both a reactant and a reaction product, it will be undefined (nan)
    c = {'A': 15.2, 'B': 0, 'C': 0 , 'D': 4.1}
    assert np.isnan(rxn.reaction_quotient(conc=c, explain=False))




#######  For PRIVATE methods  #######

def test__standard_form_chem_eqn():
    chem_data = ChemData(names=['Fe', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'Cl'])
    rxn = Reaction(chem_data, reactants="A", products="B")  # Won't actually use reactants/products

    assert rxn._standard_form_chem_eqn([(1, 0, 1), (2, 8, 1)]) == "Fe + 2 Cl"

    assert rxn._standard_form_chem_eqn([(3, 0, 1), (5, 7, 2)]) == "3 Fe + 5 G"



def test__parse_reaction_term():
    chem_data = ChemData(names=["A", "B", "C", "D", "E", "F"])
    rxn = Reaction(chem_data, reactants="A", products="B")  # Won't actually use the reactants/products

    with pytest.raises(Exception):
        rxn._parse_reaction_term(5)    # The argument is not a string nor a tuple nor a list

    assert rxn._parse_reaction_term("F") == (1, 5, 1)   # 5 is the index of "F"

    with pytest.raises(Exception):
        rxn._parse_reaction_term("Not the name of any chemical")

    with pytest.raises(Exception):
        rxn._parse_reaction_term( (2, 5) )   # The last item in the pair is not a string

    assert rxn._parse_reaction_term( (2, "F") ) == (2, 5, 2)    # order defaults to stoichiometry
    assert rxn._parse_reaction_term( [2, "F"] ) == (2, 5, 2)

    with pytest.raises(Exception):
        rxn._parse_reaction_term( (2, 5, 1) )   # The mid-item in the triplet is not a string

    assert rxn._parse_reaction_term( (2, "F", 1) ) == (2, 5, 1)
    assert rxn._parse_reaction_term( [2, "F", 1] ) == (2, 5, 1)

    with pytest.raises(Exception):
        rxn._parse_reaction_term( (3, 5, 2, 123) )     # Extra element in tuple



def test__set_kinetic_and_thermodynamic():
    chem_data = ChemData(names=['A', 'B'])

    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=6,  reverse_rate=None,
                                       delta_H=None, delta_S=None, delta_G=None, temp=None)
    assert rxn.kF == 6
    assert rxn.kR is None
    assert rxn.K is None
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=6,  reverse_rate=2,
                                       delta_H=None, delta_S=None, delta_G=None, temp=None)
    assert rxn.kF == 6
    assert rxn.kR == 2
    assert rxn.K == 3
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=6,  reverse_rate=2,
                                       delta_H=None, delta_S=None, delta_G=None, temp=100)
    assert rxn.kF == 6
    assert rxn.kR == 2
    assert rxn.K == 3
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert np.allclose(rxn.delta_G, -913.437080597)

    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    with pytest.raises(Exception):
        # Inconsistent kinetic data and passed thermo data
        rxn._set_kinetic_and_thermodynamic(forward_rate=6,  reverse_rate=2,
                                           delta_H=None, delta_S=None, delta_G=666., temp=100)


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=None)
    assert rxn.kF is None
    assert rxn.kR is None
    assert rxn.K is None
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert rxn.delta_G is None


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=100)
    assert rxn.kF is None
    assert rxn.kR is None
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert rxn.delta_G == 800


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    with pytest.raises(Exception):
        # Inconsistent thermo data
        rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                           delta_H=500, delta_S=-3, delta_G=999, temp=100)


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=10,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert rxn.kR == 10
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert rxn.delta_G == 800


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=None,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert np.allclose(rxn.kR, 10)
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert rxn.delta_G == 800


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    with pytest.raises(Exception):
        # Inconsistent kinetic / thermo
        rxn._set_kinetic_and_thermodynamic(forward_rate=10,  reverse_rate=10,
                                           delta_H=500, delta_S=-3, delta_G=None, temp=100)


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=10,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert rxn.kR == 10
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert np.allclose(rxn.delta_G, 800.)


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                       delta_H=None, delta_S=None, delta_G=800, temp=None)
    assert rxn.kF is None
    assert rxn.kR is None
    assert rxn.K is None
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G == 800


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                       delta_H=None, delta_S=None, delta_G=800, temp=100)
    assert rxn.kF is None
    assert rxn.kR is None
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G == 800


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=10,
                                       delta_H=None, delta_S=None, delta_G=800, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert rxn.kR == 10
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G == 800


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=None,
                                       delta_H=None, delta_S=None, delta_G=800, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert np.allclose(rxn.kR, 10)
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G == 800


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=None,
                                       delta_H=700, delta_S=None, delta_G=800, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert np.allclose(rxn.kR, 10)
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 700
    assert rxn.delta_S == -1
    assert rxn.delta_G == 800


    rxn = Reaction(chem_data, reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=None,
                                       delta_H=None, delta_S=-1, delta_G=800, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert np.allclose(rxn.kR, 10)
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 700
    assert rxn.delta_S == -1
    assert rxn.delta_G == 800
