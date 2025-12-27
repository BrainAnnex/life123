import pytest
import numpy as np
from collections import Counter
from life123.reactions import ReactionCommon, ReactionOneStep, \
        ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition,\
        ReactionEnzyme, ReactionGeneric, ReactionKinetics



#######################   ReactionCommon   ###################################################################

def test_constructor_ReactionCommon():
    rxn = ReactionCommon()
    assert rxn.active == True
    assert rxn.temp is None


    rxn = ReactionCommon(active=False, temp=200)
    assert rxn.active == False
    assert rxn.temp == 200




#######################   ReactionOneStep   ###################################################################

def test_constructor_ReactionOneStep():
    rxn = ReactionOneStep()
    assert rxn.active == True
    assert rxn.reversible == True
    assert rxn.kF is None
    assert rxn.K is None


    rxn = ReactionOneStep(kF=10, kR=2, delta_H=-3000)
    assert rxn.active == True
    assert rxn.reversible == True
    assert rxn.kF == 10
    assert rxn.kR == 2
    assert rxn.delta_H == -3000
    assert rxn.delta_S is None
    assert rxn.K == 5


    with pytest.raises(Exception):
        ReactionOneStep(reversible=False, kF=10, kR=2)   # Irreversible can't have reverse rate


    rxn = ReactionOneStep(active=False, reversible=False, kF=10, kR=0)
    assert rxn.active == False
    assert rxn.reversible == False
    assert rxn.kF == 10
    assert rxn.kR == 0
    assert rxn.K is None




########################   ReactionUnimolecular    #########################################################

def test_constructor_ReactionUnimolecular():
    with pytest.raises(Exception):
        ReactionUnimolecular()      # Missing arguments

    with pytest.raises(Exception):
        ReactionUnimolecular(reactant=123, product="B")     # Bad reactant

    with pytest.raises(Exception):
        ReactionUnimolecular(reactant="A", product=True)    # Bad product

    with pytest.raises(Exception):
        ReactionUnimolecular(reactant="A", product="A")     # Cannot be same


    rxn = ReactionUnimolecular(reactant="A", product="B")
    assert rxn.active == True
    assert rxn.reversible == True
    assert rxn.kF is None
    assert rxn.K is None
    assert rxn.delta_S is None
    assert rxn.reactant == "A"
    assert rxn.product == "B"


    rxn = ReactionUnimolecular(reversible=False, reactant="A", product="B",
                               kF=20)
    assert rxn.active == True
    assert rxn.reversible == False
    assert rxn.kF == 20
    assert rxn.kR is None
    assert rxn.K is None
    assert rxn.delta_S is None
    assert rxn.reactant == "A"
    assert rxn.product == "B"


    rxn = ReactionUnimolecular(active=False, reactant="A", product="B",
                               kF=20, kR=4)
    assert rxn.active == False
    assert rxn.reversible == True
    assert rxn.kF == 20
    assert rxn.kR == 4
    assert rxn.K == 5
    assert rxn.delta_S is None
    assert rxn.reactant == "A"
    assert rxn.product == "B"



def test_extract_reactant_labels_ReactionUnimolecular():
    rxn = ReactionUnimolecular(reactant="A", product="B")
    assert rxn.extract_reactant_labels() == ["A"]


def test_extract_product_labels_ReactionUnimolecular():
    rxn = ReactionUnimolecular(reactant="A", product="B")
    assert rxn.extract_product_labels() == ["B"]


def test_extract_chemicals_in_reaction_ReactionUnimolecular():
    rxn = ReactionUnimolecular(reactant="A", product="B")
    assert rxn.extract_chemicals_in_reaction() == {"A", "B"}




def test_describe_ReactionUnimolecular():
    rxn = ReactionUnimolecular(reactant="A", product="B")

    assert rxn.describe(concise=True) == "A <-> B"
    assert rxn.describe(concise=False) == "A <-> B  (Elementary Unimolecular reaction)"

    # TODO: more tests with temp and kinetic pars



def test_determine_reaction_rate_ReactionUnimolecular():
    rxn = ReactionUnimolecular(reactant="A", product="B",
                               kF=20., kR=2., reversible=True)

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5. - 2. * 8.)  # 84.0

    rxn.reversible = False

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5.)            # 100.0



def test_reaction_quotient_ReactionUnimolecular():
    # Reaction : A <-> B
    rxn = ReactionUnimolecular(reactant="A", product="B")
    c = {'A': 24., 'B': 36.}
    assert np.allclose(1.5, rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1.5, quotient)
    assert formula == '[B] / [A]'


    # Reaction : A <-> F
    rxn = ReactionUnimolecular(reactant="A", product="F")
    c = {'A': 3., 'F': 33.}
    assert np.allclose(11., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(11., quotient)
    assert formula == '[F] / [A]'



def test_step_simulation_ReactionUnimolecular():
    # Reaction : A <-> B
    rxn = ReactionUnimolecular(reactant="A", product="B", kF=3., kR=2.)
    result = rxn.step_simulation(delta_time=0.1, conc_dict={"A": 10, "B": 50})
    assert result[0] == {'A': 7, 'B': -7}
    assert result[1] == -70






########################   ReactionSynthesis    #########################################################

def test_extract_reactant_labels_ReactionSynthesis():
    rxn = ReactionSynthesis(reactants=["A", "B"], product="C")
    assert rxn.extract_reactant_labels() == ["A", "B"]


def test_extract_product_labels_ReactionSynthesis():
    rxn = ReactionSynthesis(reactants=["A", "B"], product="C")
    assert rxn.extract_product_labels() == ["C"]


def test_extract_chemicals_in_reaction_ReactionSynthesis():
    rxn = ReactionSynthesis(reactants=["A", "B"], product="C")
    assert rxn.extract_chemicals_in_reaction() == {"A", "B", "C"}



def test_determine_reaction_rate_ReactionSynthesis():
    # Reaction A + B <-> C
    rxn = ReactionSynthesis(reactants=["A", "B"], product="C",
                            kF=20., reversible=False)

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5. * 8.)

    rxn.reversible = True
    rxn.kR = 2.

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5. * 8. - 2. * 3.)



def test_reaction_quotient_ReactionSynthesis():
    # Reaction :  A + B <-> C
    rxn = ReactionSynthesis(reactants=["A" , "B"], product="C")
    c = {'A': 3., 'B': 4., 'C': 12.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1., quotient)
    assert formula == '[C] / ([A][B])'


    # Reaction :  2A <-> P
    rxn = ReactionSynthesis(reactants=["A" , "A"], product="P")
    c = {'A': 2., 'P': 20.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(5., quotient)
    assert formula == '[P] / ([A][A])'



def test_step_simulation_ReactionSynthesis():
    # Reaction : A + B <-> C
    rxn = ReactionSynthesis(reactants=["A" , "B"], product="C", kF=5., kR=2.)

    result = rxn.step_simulation(delta_time=0.002, conc_dict={"A": 10, "B": 50, "C": 20})
    assert result[0] == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert result[1] == 5*10*50 - 2 * 20        # 2460





########################   ReactionDecomposition    #########################################################

def test_extract_reactant_labels_ReactionDecomposition():
    rxn = ReactionDecomposition(reactant="A", products=["B", "C"])
    assert rxn.extract_reactant_labels() == ["A"]


def test_extract_product_labels_ReactionDecomposition():
    rxn = ReactionDecomposition(reactant="A", products=["B", "C"])
    assert rxn.extract_product_labels() == ["B", "C"]


def test_extract_chemicals_in_reaction_ReactionDecomposition():
    rxn = ReactionDecomposition(reactant="A", products=["B", "C"])
    assert rxn.extract_chemicals_in_reaction() == {"A", "B", "C"}



def test_determine_reaction_rate_ReactionDecomposition():
    # Reaction A <-> B + C
    rxn = ReactionDecomposition(reactant="A", products=["B", "C"],
                                kF=20., reversible=False)

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5.)

    rxn.reversible = True
    rxn.kR = 2.

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5.  - 2. * 8. * 3.)



def test_determine_reaction_quotient_ReactionDecomposition():
    # Reaction :  C <-> A + B
    rxn = ReactionDecomposition(reactant="C", products=["A" , "B"])
    c = {'A': 3., 'B': 4., 'C': 12.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1., quotient)
    assert formula == '([A][B]) / [C]'


    # Reaction :  B <-> 2A
    rxn = ReactionDecomposition(reactant="B", products=["A" , "A"])
    c = {'A': 2., 'B': 20.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1/5., quotient)
    assert formula == '([A][A]) / [B]'



def test_step_simulation_ReactionDecomposition():
    # Reaction : C <-> A + B
    rxn = ReactionDecomposition(reactant="C", products=["A" , "B"], kF=2., kR=5.)

    result = rxn.step_simulation(delta_time=0.002, conc_dict={"A": 10, "B": 50, "C": 20})
    assert result[0] == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert result[1] == -2460




########################   ReactionEnzyme    #########################################################

def test_initialize_ReactionEnzyme():
    rxn = ReactionEnzyme(enzyme="E", substrate="S", product="P",
                         k1_F=10., k1_R=2., k2_F=5.)

    assert rxn.active == True
    assert rxn.temp is None

    assert rxn.enzyme == "E"
    assert rxn.substrate == "S"
    assert rxn.product == "P"
    assert rxn.intermediate == "ES"

    assert np.allclose(rxn.k1_F, 10.)
    assert np.allclose(rxn.k1_R, 2.)
    assert np.allclose(rxn.k2_F, 5.)

    assert np.allclose(rxn.kM, 0.7)     # (rxn.k2_F + rxn.k1_R) / rxn.k1_F)
    assert np.allclose(rxn.kcat, 5.)    # Equal to k2_F

    with pytest.raises(Exception):
        ReactionEnzyme(enzyme="E", substrate="S", product="P",
                       k1_F=10., k1_R=2., k2_F=5., kM=0.71) # Inconsistent

    with pytest.raises(Exception):
        ReactionEnzyme(enzyme="E", substrate="S", product="P",
                       k1_F=10., k1_R=2., k2_F=5., kcat=5.01) # Inconsistent

    # TODO: more testing



def test_extract_reactant_labels_ReactionEnzyme():
    rxn = ReactionEnzyme(enzyme="E", substrate="S", product="P")
    assert rxn.extract_reactant_labels() == ["S", "E"]


def test_extract_product_labels_ReactionEnzyme():
    rxn = ReactionEnzyme(enzyme="E", substrate="S", product="P")
    assert rxn.extract_product_labels() == ["P", "E"]


def test_extract_chemicals_in_reaction_ReactionEnzyme():
    rxn = ReactionEnzyme(enzyme="E", substrate="S", product="P")
    assert rxn.extract_chemicals_in_reaction() == {"S", "P", "E", "ES"}



def test_step_simulation_ReactionEnzyme_1():
    # E + S <-> ES -> E + P
    initial_conc = {"S": 20, "E": 1, "P": 0, "ES": 0}
    dt = 0.002

    rxn = ReactionEnzyme(enzyme="E", substrate="S", product="P",
                         k1_F=18., k1_R=100., k2_F=49.)
    result = rxn.step_simulation(delta_time=dt, conc_dict=initial_conc)
    assert result[0] == {'E': -0.72, 'S': -0.72, 'ES': 0.72, 'P': 0.0}
    assert result[1] == (360, 0)



def test_step_simulation_ReactionEnzyme_2():
    # E + S <-> ES -> E + P
    initial_conc = {"S": 20, "E": 1, "P": 0, "ES": 0}
    dt = 0.002
    print()
    conc = initial_conc

    for _ in range(3):
        # Compare the results given by ReactionEnzyme with those
        # given by simulating the 2 reactions separately
        rxn = ReactionEnzyme(enzyme="E", substrate="S", product="P",
                             k1_F=18., k1_R=100., k2_F=49.)
        result = rxn.step_simulation(delta_time=dt, conc_dict=conc)
        #print(result)

        # Simulate the 2 elementary reactions separately
        rxn1 = ReactionSynthesis(reactants=["E", "S"], product="ES",
                                 kF=18., kR=100.)
        result1 = rxn1.step_simulation(delta_time=dt, conc_dict=conc)

        rxn2 = ReactionDecomposition(reactant="ES", products=["E", "P"],
                                     kF=49., kR=0)
        result2 = rxn2.step_simulation(delta_time=dt, conc_dict=conc)

        # Check the increment dictionary
        c = Counter(result1[0])
        c.update(result2[0])                # Sum up dict values whenever the keys match
        incr_dict = dict(c)
        assert incr_dict == result[0]

         # Check the reaction rates
        assert result[1] == (result1[1], result2[1])

        # Update the system concentrations (advancing the simulation)
        for k, v in incr_dict.items():
            conc[k] += v

        #print("conc:", conc)





########################   ReactionGeneric    #########################################################

def test_initialize():
    rxn = ReactionGeneric(reactants="A", products="B")

    assert rxn.reactants == [(1, "A")]
    assert rxn.products == [(1, "B")]
    assert rxn.kF is None
    assert rxn.kR is None
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None
    assert rxn.K is None


    # More complex scenarios

    # Missing products or reactants
    with pytest.raises(Exception):
        ReactionGeneric(reactants=["A"], products=None)
    with pytest.raises(Exception):
        ReactionGeneric(reactants=None, products="B")

    # Reactants and the products can't be the same
    with pytest.raises(Exception):
        ReactionGeneric(reactants=["A"], products=["A"])
    with pytest.raises(Exception):
        ReactionGeneric(reactants=["A"], products=[("A")])
    with pytest.raises(Exception):
        ReactionGeneric(reactants=["A"], products=[(1, "A")])
    with pytest.raises(Exception):
        ReactionGeneric(reactants=[(2, "B")], products=[(2, "B")])
    with pytest.raises(Exception):
        ReactionGeneric(reactants=["A", "B"], products=["A", "B"])
    with pytest.raises(Exception):
        ReactionGeneric(reactants=["A", (3, "B")], products=["A", (3, "B")])
    with pytest.raises(Exception):
        ReactionGeneric(reactants=["A", "B"], products=["B", "A"])
    with pytest.raises(Exception):
        ReactionGeneric(reactants=[(2, "A"), "B", "C"], products=["B", (1, "C"), (2, "A")])


    rxn = ReactionGeneric(reactants=["A"], products=["B"], kF=3., kR=2.)

    assert np.allclose(rxn.extract_forward_rate() , 3.)
    assert np.allclose(rxn.extract_reverse_rate() , 2.)
    assert np.allclose(rxn.K , 3./2.)
    assert rxn.extract_reactants() == [(1, "A")]
    assert rxn.extract_products() == [(1, "B")]
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None

    assert rxn.extract_reactants() == [(1, "A")]
    assert rxn.extract_reactants_formula() == "A"

    assert rxn.extract_products() == [(1, "B")]
    assert rxn.extract_products_formula() == "B"
 

    # Another reaction
    rxn = ReactionGeneric(reactants=[(2, "B")], products=[(5, "C")], kF=9., kR=7.)

    assert np.allclose(rxn.extract_forward_rate() , 9.)
    assert np.allclose(rxn.extract_reverse_rate() , 7.)
    assert np.allclose(rxn.K , 9./7.)
    assert rxn.extract_reactants() == [(2, "B")]
    assert rxn.extract_products() == [(5, "C")]
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None


    # Add another reaction.  This time, first set the temperature
    rxn = ReactionGeneric(reactants=[(2, "D")], products=[(1, "C")], kF=11., kR=13., temp=200)

    assert np.allclose(rxn.extract_forward_rate() , 11.)
    assert np.allclose(rxn.extract_reverse_rate() , 13.)
    assert np.allclose(rxn.K , 11./13.)
    assert rxn.extract_reactants() == [(2, "D")]
    assert rxn.extract_products() == [(1, "C")]
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert np.allclose(rxn.delta_G, 277.7928942715384)   # - RT log(K)


    # Add a multi-term reaction
    rxn = ReactionGeneric(reactants=["A", (2, "B")], products=[(3, "C"), "D"], kF=5., kR=1., temp=200)

    assert np.allclose(rxn.extract_forward_rate() , 5.)
    assert np.allclose(rxn.extract_reverse_rate() , 1.)
    assert np.allclose(rxn.K , 5./1.)
    assert rxn.extract_reactants() == [(1, "A"), (2, "B")]
    assert rxn.extract_products() == [(3, "C"), (1, "D")]
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert np.allclose(rxn.delta_G, -2676.321364705849)   # - RT log(K)


    # Add a reaction with thermodynamic data;
    # the reverse reaction rate will get computed from the thermodynamic data
    rxn = ReactionGeneric(reactants=["A"], products=[(2, "B")], kF=10.,
                          delta_H= 5., delta_S= 0.4, temp=200)

    assert rxn.extract_reactants() == [(1, "A")]
    assert rxn.extract_products() == [(2, "B")]
    assert rxn.delta_H == 5.
    assert rxn.delta_S == 0.4
    assert np.allclose(rxn.delta_G, -75.0)         # 5 - 200 * 0.4
    assert np.allclose(rxn.K, 1.0461347154679432)  # exp(75/(8.3144598 * 200))
    assert np.allclose(rxn.extract_forward_rate() , 10.)
    assert np.allclose(rxn.extract_reverse_rate() , 9.558998331803693) # 10. / 1.0461347154679432



def test_set_macro_enzyme():
    rxn = ReactionGeneric(reactants="A", products="B")
    rxn.set_macro_enzyme("mm", 4)

    assert rxn.macro_enzyme == ("mm", 4)



def test_extract_reactants():
    rxn = ReactionGeneric(reactants=["CH4", (2, "O2")],
                          products=["CO2", (2, "H2O")])
    assert rxn.extract_reactants() == [(1, "CH4"), (2, "O2")]


def test_extract_reactants_formula():
    rxn = ReactionGeneric(reactants=["CH4", (2, "O2")],
                          products=["CO2", (2, "H2O")])
    assert rxn.extract_reactants_formula() == "CH4 + 2 O2"



def test_extract_products():
    rxn = ReactionGeneric(reactants=["CH4", (2, "O2")],
                          products=["CO2", (2, "H2O")])
    assert rxn.extract_products() == [(1, "CO2"), (2, "H2O")]


def test_extract_products_formula():
    rxn = ReactionGeneric(reactants=["CH4", (2, "O2")],
                          products=["CO2", (2, "H2O")])
    assert rxn.extract_products_formula() == "CO2 + 2 H2O"




#########  TO DESCRIBE THE DATA  #########

def test_describe():
    rxn = ReactionGeneric(reactants=["CH4", (2, "O2")], products=["CO2", (2, "H2O")])

    assert rxn.describe(concise=True) == "CH4 + 2 O2 <-> CO2 + 2 H2O"
    assert rxn.describe(concise=False) == "CH4 + 2 O2 <-> CO2 + 2 H2O | No kinetic rate function provided"

    rxn.set_rate_function(ReactionKinetics.compute_rate_pseudo_elementary)

    assert rxn.describe(concise=False) == "CH4 + 2 O2 <-> CO2 + 2 H2O | Kinetic rate function : `compute_rate_pseudo_elementary()`"


    # Start over

    rxn = ReactionGeneric(reactants=["A"], products=["B"], kF=3., kR=2.)
    assert rxn.describe(concise=True) == "A <-> B"
    assert rxn.describe(concise=False) == "A <-> B  (kF = 3 / kR = 2 / K = 1.5) | No kinetic rate function provided"

    rxn = ReactionGeneric(reactants=[(2, "B")], products=[(5, "C")], kF=9., kR=7.)
    assert rxn.describe(concise=True) == "2 B <-> 5 C"
    assert rxn.describe(concise=False) == "2 B <-> 5 C  (kF = 9 / kR = 7 / K = 1.2857) | No kinetic rate function provided"


    rxn = ReactionGeneric(reactants=[(2, "D")], products=[(1, "C")], kF=11., kR=13., temp=200)
    assert rxn.describe(concise=True) == "2 D <-> C"
    assert rxn.describe(concise=False) == "2 D <-> C  (kF = 11 / kR = 13 / delta_G = 277.79 / K = 0.84615 / Temp = -73.15 C) | No kinetic rate function provided"

    rxn = ReactionGeneric(reactants=["A", (2, "B")], products=[(3, "C"), "D"], kF=5., kR=1., temp=200)
    assert rxn.describe(concise=True) == "A + 2 B <-> 3 C + D"
    assert rxn.describe(concise=False) == "A + 2 B <-> 3 C + D  (kF = 5 / kR = 1 / delta_G = -2,676.3 / K = 5 / Temp = -73.15 C) | No kinetic rate function provided"



def test_extract_rxn_properties():
    rxn = ReactionGeneric(reactants=["A"], products=["B"], kF=3., kR=2., temp = 298.15)

    result = rxn.extract_rxn_properties()
    #print(result)
    assert np.allclose(result["kF"] , 3.)
    assert np.allclose(result["kR"] , 2.)
    assert np.allclose(result["K"] , 1.5)
    assert np.allclose(result["delta_G"] , -1005.1305052750387)



def test_extract_reactant_labels_ReactionGeneric():
    rxn = ReactionGeneric(reactants="A", products="B")
    assert rxn.extract_reactant_labels() == ["A"]

    rxn = ReactionGeneric(reactants=["A", "B"], products="C")
    assert rxn.extract_reactant_labels() == ["A", "B"]

    rxn = ReactionGeneric(reactants=["A", "B"], products=["B", "C"])               # B acts as catalyst
    assert rxn.extract_reactant_labels() == ["A", "B"]

    rxn = ReactionGeneric(reactants=[(2, "D"), "C"], products=["C", (3, "E")])     # C acts as catalyst
    assert rxn.extract_reactant_labels() == ["D", "C"]



def test_extract_product_labels_ReactionGeneric():
    rxn = ReactionGeneric(reactants="A", products="B")
    assert rxn.extract_product_labels() == ["B"]

    rxn = ReactionGeneric(reactants=["A", "B"], products="C")
    assert rxn.extract_product_labels() == ["C"]


    rxn = ReactionGeneric(reactants=["A", "B"], products=["B", "C"])               # B acts as catalyst
    assert rxn.extract_product_labels() == ["B", "C"]

    rxn = ReactionGeneric(reactants=[(2, "D"), "C"], products=["C", (3, "E")])     # C acts as catalyst
    assert rxn.extract_product_labels() == ["C", "E"]



def test_extract_chemicals_in_reaction_ReactionGeneric():
    rxn = ReactionGeneric(reactants="A", products="B")
    assert rxn.extract_chemicals_in_reaction() == {"A", "B"}

    rxn = ReactionGeneric(reactants=["A", "B"], products="C")
    assert rxn.extract_chemicals_in_reaction() == {"A", "B", "C"}

    rxn = ReactionGeneric(reactants=["A", "B"], products=["B", "C"])               # B acts as catalyst
    assert rxn.extract_chemicals_in_reaction() == {"A", "B", "C"}

    rxn = ReactionGeneric(reactants=[(2, "D"), "C"], products=["C", (3, "E")])     # C acts as catalyst
    assert rxn.extract_chemicals_in_reaction() == {"C", "D", "E"}





#############  For ANALYSIS  #############

def test_reaction_quotient():
    # Reaction : A <-> B
    rxn = ReactionGeneric(reactants="A", products="B")
    rxn.set_rate_function(ReactionKinetics.compute_rate_pseudo_elementary)
    c = {'A': 24., 'B': 36.}
    assert np.allclose(1.5, rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1.5, quotient)
    assert formula == '[B] / [A]'


    # Reaction : A <-> F
    rxn = ReactionGeneric(reactants="A", products="F")
    rxn.set_rate_function(ReactionKinetics.compute_rate_pseudo_elementary)
    c = {'A': 3., 'F': 33.}
    assert np.allclose(11., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(11., quotient)
    assert formula == '[F] / [A]'


    # Reaction : A <-> 3B, hypothetically treated as an elementary reaction
    rxn = ReactionGeneric(reactants=["A"], products=[(3, "B")])
    rxn.set_rate_function(ReactionKinetics.compute_rate_pseudo_elementary)
    c = {'A': 5., 'B': 10.}
    assert np.allclose(200., rxn.reaction_quotient(conc=c, explain=False))  # 10^3 / 5
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(200., quotient)    # 10^3 / 5
    assert formula == '[B]^3  / [A]'


    # Reaction :  A + B <-> C , with 1st-order kinetics for each species
    rxn = ReactionGeneric(reactants=["A" , "B"], products=["C"])
    rxn.set_rate_function(ReactionKinetics.compute_rate_pseudo_elementary)
    c = {'A': 3., 'B': 4., 'C': 12.}
    assert np.allclose(1., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1., quotient)
    assert formula == '[C] / ([A][B])'


    # Reaction :  2A <-> B
    rxn = ReactionGeneric(reactants=[(2, "A")], products="B")
    rxn.set_rate_function(ReactionKinetics.compute_rate_pseudo_elementary)
    c = {'A': 4., 'B': 20.}
    assert np.allclose(1.25, rxn.reaction_quotient(conc=c, explain=False))  # 20 / 4^2
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1.25, quotient)  # 20 / 4^2
    assert formula == '[B] / [A]^2 '


    # Reaction :  A <-> 2B , with 2nd-order kinetics in the reverse direction
    rxn = ReactionGeneric(reactants="A", products=[(2, "B")])
    rxn.set_rate_function(ReactionKinetics.compute_rate_pseudo_elementary)
    c = {'A': 4., 'B': 20.}
    assert np.allclose(100., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(100., quotient)
    assert formula == '[B]^2  / [A]'


    # Reaction :  A + B <-> C + D
    rxn = ReactionGeneric(reactants=["A", "B"], products=["C", "D"])
    rxn.set_rate_function(ReactionKinetics.compute_rate_pseudo_elementary)

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
    rxn = ReactionGeneric(reactants="A", products="B")     # Won't actually use reactants/products

    assert rxn._standard_form_chem_eqn([(1, "Fe"), (2, "Cl")]) == "Fe + 2 Cl"

    assert rxn._standard_form_chem_eqn([(3, "Fe"), (5, "G")]) == "3 Fe + 5 G"



def test__parse_reaction_term():
    rxn = ReactionGeneric(reactants="A", products="B")     # Won't actually use the reactants/products

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



def test__set_kinetic_and_thermodynamic():
    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=6,  reverse_rate=None,
                                       delta_H=None, delta_S=None, delta_G=None, temp=None)
    assert rxn.kF == 6
    assert rxn.kR is None
    assert rxn.K is None
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=6,  reverse_rate=2,
                                       delta_H=None, delta_S=None, delta_G=None, temp=None)
    assert rxn.kF == 6
    assert rxn.kR == 2
    assert rxn.K == 3
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G is None


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=6,  reverse_rate=2,
                                       delta_H=None, delta_S=None, delta_G=None, temp=100)
    assert rxn.kF == 6
    assert rxn.kR == 2
    assert rxn.K == 3
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert np.allclose(rxn.delta_G, -913.437080597)

    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    with pytest.raises(Exception):
        # Inconsistent kinetic data and passed thermo data
        rxn._set_kinetic_and_thermodynamic(forward_rate=6,  reverse_rate=2,
                                           delta_H=None, delta_S=None, delta_G=666., temp=100)


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=None)
    assert rxn.kF is None
    assert rxn.kR is None
    assert rxn.K is None
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert rxn.delta_G is None


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=100)
    assert rxn.kF is None
    assert rxn.kR is None
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert rxn.delta_G == 800


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    with pytest.raises(Exception):
        # Inconsistent thermo data
        rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                           delta_H=500, delta_S=-3, delta_G=999, temp=100)


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=10,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert rxn.kR == 10
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert rxn.delta_G == 800


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=None,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert np.allclose(rxn.kR, 10)
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert rxn.delta_G == 800


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    with pytest.raises(Exception):
        # Inconsistent kinetic / thermo
        rxn._set_kinetic_and_thermodynamic(forward_rate=10,  reverse_rate=10,
                                           delta_H=500, delta_S=-3, delta_G=None, temp=100)


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=10,
                                       delta_H=500, delta_S=-3, delta_G=None, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert rxn.kR == 10
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 500
    assert rxn.delta_S == -3
    assert np.allclose(rxn.delta_G, 800.)


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                       delta_H=None, delta_S=None, delta_G=800, temp=None)
    assert rxn.kF is None
    assert rxn.kR is None
    assert rxn.K is None
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G == 800


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=None,
                                       delta_H=None, delta_S=None, delta_G=800, temp=100)
    assert rxn.kF is None
    assert rxn.kR is None
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G == 800


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=None,  reverse_rate=10,
                                       delta_H=None, delta_S=None, delta_G=800, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert rxn.kR == 10
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G == 800


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=None,
                                       delta_H=None, delta_S=None, delta_G=800, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert np.allclose(rxn.kR, 10)
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H is None
    assert rxn.delta_S is None
    assert rxn.delta_G == 800


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=None,
                                       delta_H=700, delta_S=None, delta_G=800, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert np.allclose(rxn.kR, 10)
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 700
    assert rxn.delta_S == -1
    assert rxn.delta_G == 800


    rxn = ReactionGeneric(reactants=["A"], products=["B"])
    rxn._set_kinetic_and_thermodynamic(forward_rate=3.8205953171,  reverse_rate=None,
                                       delta_H=None, delta_S=-1, delta_G=800, temp=100)
    assert np.allclose(rxn.kF, 3.8205953171)
    assert np.allclose(rxn.kR, 10)
    assert np.allclose(rxn.K, 0.38205953171)
    assert rxn.delta_H == 700
    assert rxn.delta_S == -1
    assert rxn.delta_G == 800
