import numpy as np
import pytest
import pandas as pd
from life123.species_registry import Species, SpeciesRegistry, MacroMolecules
from life123.reactions_new import Stoichiometry, Kinetics, ReactionThermodynamics, Reaction
from life123.visualization.colors import Colors
from tests.utilities.comparisons import *



############################  class Stoichiometry  ############################

def test_CONSTRUCTOR_Stoichiometry():

    # Reaction R -> P
    st = Stoichiometry(coefficients={"R": -1, "P": 1})
    assert st.coefficients == {"R": -1, "P": 1}

    # Reaction R -> P + Q
    st = Stoichiometry(coefficients={"R": -1, "P": 1, "Q": 1})
    assert st.coefficients == {"R": -1, "P": 1, "Q": 1}



def test_reaction_pattern():

    # Reaction R -> P
    st = Stoichiometry(coefficients={"R": -1, "P": 1})
    assert st.reaction_pattern() == (1, 1, 0)

    # Reaction R -> P + Q
    st = Stoichiometry(coefficients={"R": -1, "P": 1, "Q": 1})
    assert st.reaction_pattern() == (1, 2, 0)

    # Reaction A + E -> B + E
    st = Stoichiometry(coefficients={"A": -1, "B": 1, "E": 0})
    assert st.reaction_pattern() == (1, 1, 1)



def test_get_reaction_vector():
    # Reaction  2A + 3B + 2 E + F <--> 4C  + 5D  + 2 E + F
    st = Stoichiometry({"A": -2, "B": -3, "E": 0, "F": 0, "C": 4, "D": 5})
    assert st.get_reaction_vector() == {"A": -2, "B": -3, "C": 4, "D": 5}


    sr = SpeciesRegistry(ids=["A", "B", "C", "D", "R", "P", "Q", "S", "E", "F"])

    rxn = Reaction(reactants="R", products="P", species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"R": -1, "P": 1}

    rxn = Reaction(reactants=("R", "S"), products="P", species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"R": -1, "S": -1, "P": 1}

    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"R": -1, "P": 1, "Q": 1}

    rxn = Reaction(reactants=("E", "S"), products=("E", "P"),
                   species_registry=sr, kinetics_type="Michaelis-Menten")
    assert rxn.stoichiometry.get_reaction_vector() == {"S": -1, "P": 1}

    rxn = Reaction(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                   species_registry=sr, kinetics_type="custom")
    assert rxn.stoichiometry.get_reaction_vector() == {"A": -2, "B":- 2, "P": 3, "Q": 1}


    # Reaction  2A + 3B + 2 E + F <--> 4C  + 5D  + 2 E + F
    rxn = Reaction(reactants=["E", "A", "F", (3, "B"), "E", "A"], products=[(4, "C"), "F", (2, "E"), (5, "D")],
                   species_registry=sr, kinetics_type="custom")
    assert rxn.stoichiometry.get_reaction_vector() == {"A": -2, "B": -3, "C": 4, "D": 5}



def test_get_split_reaction_vector():
    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    rxn = Reaction(reactants="R", products="P", species_registry=sr)
    assert rxn.stoichiometry.get_split_reaction_vector() == ({"R": -1}, {"P": 1})

    rxn = Reaction(reactants=("R", "S"), products="P", species_registry=sr)
    assert rxn.stoichiometry.get_split_reaction_vector() == ({"R": -1, "S": -1}, {"P": 1})

    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.stoichiometry.get_split_reaction_vector() == ({"R": -1}, {"P": 1, "Q": 1})

    rxn = Reaction(reactants=("E", "S"), products=("E", "P"),
                   species_registry=sr, kinetics_type="Michaelis-Menten")
    assert rxn.stoichiometry.get_split_reaction_vector() == ({"S": -1}, {"P": 1})

    rxn = Reaction(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                   species_registry=sr, kinetics_type="custom")
    assert rxn.stoichiometry.get_split_reaction_vector() == ({"A": -2, "B":- 2}, {"P": 3, "Q": 1})





############################  class Kinetics  ############################

def test_constructor_Kinetics():
    k = Kinetics(law="mass action")
    assert k.law == "mass action"

    with pytest.raises(Exception):
        Kinetics(law="never heard of this!")

    k = Kinetics(law="mass action", parameters={"kF": 10})
    assert k.law == "mass action"
    assert k.parameters == {"kF": 10, "kR": 0, "reversible": False, "K": None}
    #TODO: more tests



def test_to_dict_Kinetics():
    k = Kinetics(law="mass action")
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 0, "kR": 0, "reversible": False}

    k = Kinetics(law="mass action", parameters={"kF": 10})
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 10, "kR": 0, "reversible": False}

    k = Kinetics(law="mass action", parameters={"kF": 10, "kR": 2})
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 10, "kR": 2, "reversible": True}

    k = Kinetics(law="mass action", parameters={"kF": 10})
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 10, "kR": 0, "reversible": False}

    with pytest.raises(Exception):
        Kinetics(law="random name")



def test_consistency_checker():
    # Reaction  A <-> B
    st = Stoichiometry({"A": -1, "B": 1})

    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 0, "B": 0, "C": 0}, conc_after={"A": 0, "B": 0})

    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 0, "B": 0}, conc_after={"A": 0, "B": 0, "C": 0})

    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 0, "X": 0}, conc_after={"A": 0, "B": 0})

    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 0, "B": 0}, conc_after={"X": 0, "B": 0})

    st.consistency_checker(conc_before={"A": 0, "B": 0}, conc_after={"A": 0, "B": 0})
    st.consistency_checker(conc_before={"A": 0, "B": 50}, conc_after={"A": 10, "B": 40})

    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 0, "B": 50}, conc_after={"A": 10, "B": 39.9})

    st.consistency_checker(conc_before={"A": 100, "B": 0}, conc_after={"A": 90, "B": 10})

    # Reaction 2A <-> B
    st = Stoichiometry({"A": -2, "B": 1})

    st.consistency_checker(conc_before={"A": 0, "B": 0}, conc_after={"A": 0, "B": 0})
    st.consistency_checker(conc_before={"A": 100, "B": 0}, conc_after={"A": 80, "B": 10})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 0}, conc_after={"A": 80, "B": 10.1})

    st.consistency_checker(conc_before={"A": 0, "B": 50}, conc_after={"A": 10, "B": 45})
    
    # Reaction A + B <--> C
    st = Stoichiometry({"A": -1, "B": -1, "C": 1})

    st.consistency_checker(conc_before={"A": 0, "B": 0, "C": 0}, conc_after={"A": 0, "B": 0, "C": 0})
    st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 40, "C": 10})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 40, "C": 9.9})

    # Reaction A + 3B <--> C
    st = Stoichiometry({"A": -1, "B": -3, "C": 1})

    st.consistency_checker(conc_before={"A": 0, "B": 0, "C": 0}, conc_after={"A": 0, "B": 0, "C": 0})
    st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 20, "C": 10})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 20, "C": 10.1})

    # Reaction  A + 3B <--> 4C
    st = Stoichiometry({"A": -1, "B": -3, "C": 4})

    st.consistency_checker(conc_before={"A": 0, "B": 0, "C": 0}, conc_after={"A": 0, "B": 0, "C": 0})
    st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 20, "C": 40})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 20, "C": 39.9})

    # Reaction  2A + 3B <--> 4C + 5D
    st = Stoichiometry({"A": -2, "B": -3, "C": 4, "D": 5})

    st.consistency_checker(conc_before={"A": 0, "B": 0, "C": 0, "D": 0}, conc_after={"A": 0, "B": 0, "C": 0, "D": 0})
    st.consistency_checker(conc_before={"A": 100, "B": 100, "C": 100, "D": 100}, conc_after={"A": 120, "B": 130, "C": 60, "D": 50})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 100, "C": 100, "D": 100}, conc_after={"A": 120.1, "B": 130, "C": 60, "D": 50})

    st.consistency_checker(conc_before={"A": 100, "B": 100, "C": 100, "D": 100}, conc_after={"A": 80, "B": 70, "C": 140, "D": 150})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 100, "C": 100, "D": 100.1}, conc_after={"A": 80, "B": 70, "C": 140, "D": 150})




############################  class ReactionThermodynamics  ############################

def test_CONSTRUCTOR_ReactionThermodynamics():
    rt = ReactionThermodynamics(delta_H=3, delta_S=-200, delta_G=1000, K=0.1)
    assert rt.delta_H == 3
    assert rt.delta_S == -200
    assert rt.delta_G == 1000
    assert np.allclose(rt.K, 0.1)



def test_to_dict():
    rt = ReactionThermodynamics(delta_H=3, delta_S=-200, delta_G=1000, K=0.1)
    assert rt.to_dict() == {'delta_H': 3, 'delta_S': -200, 'delta_G': 1000, 'K': 0.1}

    rt.K = None
    assert rt.to_dict() == {'delta_H': 3, 'delta_S': -200, 'delta_G': 1000}

    rt = ReactionThermodynamics(delta_G=1000)
    assert rt.to_dict() == {'delta_G': 1000}

    rt = ReactionThermodynamics()
    assert rt.to_dict() == {}





########################  class Reaction  ########################

def test_constructor_Reaction_1():

    sr = SpeciesRegistry(ids=["R", "P", "Q"])

    with pytest.raises(Exception):
        Reaction(reactants="X", products="P", species_registry=sr)  # Un-registered reactant

    with pytest.raises(Exception):
        Reaction(reactants="R", products="Y", species_registry=sr)  # Un-registered product

    with pytest.raises(Exception):
        Reaction(reactants="R", products=123, species_registry=sr)  # Bad product

    with pytest.raises(Exception):
        Reaction(reactants="R", products="R", species_registry=sr)     # Cannot be same

    with pytest.raises(Exception):
        Reaction(reactants=("R", (2, "P")), products=[(2, "P"), "R"], species_registry=sr)     # Cannot be same


    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    rxn = Reaction(reactants="R", products="P", species_registry=sr)
    assert rxn.stoichiometry.to_dict() == {"R": -1, "P": 1}

    rxn = Reaction(reactants=("R", "S"), products="P", species_registry=sr)
    assert rxn.stoichiometry.to_dict() == {"R": -1, "S": -1, "P": 1}

    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.stoichiometry.to_dict() == {"R": -1, "P": 1, "Q": 1}

    rxn = Reaction(reactants="R", products=["P", "P"], species_registry=sr)
    assert rxn.stoichiometry.to_dict() == {"R": -1, "P": 2}

    rxn = Reaction(reactants=["R", "R"], products="P", species_registry=sr)
    assert rxn.stoichiometry.to_dict() == {"R": -2, "P": 1}

    rxn = Reaction(reactants=("E", "S"), products=("E", "P"),
                   species_registry=sr, kinetics_type="Michaelis-Menten")
    assert rxn.stoichiometry.to_dict() == {"S": -1, "P": 1, "E": 0}

    rxn = Reaction(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                   species_registry=sr, kinetics_type="custom")
    assert rxn.stoichiometry.to_dict() == {"A": -2, "B":- 2, "P": 3, "Q": 1, "E": 0}



def test_constructor_Reaction_2():

    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    # Reaction R -> P
    rxn = Reaction(reactants="R", products="P", species_registry=sr)

    assert rxn.reactants == [(1, 'R')]
    assert rxn.products == [(1, 'P')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'R': -1, 'P': 1})
    assert rxn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K=None)

    assert rxn.analytic_solution_family == "ONE_TO_ONE"
    assert rxn.reaction_category == "Unimolecular rearrangement/isomerization"
    assert rxn.elementary == True
    assert rxn.active == True


    rxn = Reaction(reactants="R", products="P", species_registry=sr,
                   kinetics_type="mass action", kinetic_parameters={"kF": 20, "kR": 4})

    assert rxn.reactants == [(1, 'R')]
    assert rxn.products == [(1, 'P')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'R': -1, 'P': 1})
    assert rxn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K=5)

    assert rxn.analytic_solution_family == "ONE_TO_ONE"
    assert rxn.reaction_category == "Unimolecular rearrangement/isomerization"
    assert rxn.elementary == True
    assert rxn.active == True

    assert rxn.kinetics.parameters["kF"] == 20
    assert rxn.kinetics.parameters["kR"] == 4
    assert rxn.kinetics.parameters["reversible"] == True



    # Reaction R -> P + Q
    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr)

    assert rxn.reactants == [(1, 'R')]
    assert rxn.products == [(1, 'P'), (1, 'Q')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'R': -1, 'P': 1, 'Q': 1})
    assert rxn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K=None)

    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn.reaction_category == "Unimolecular decomposition"
    assert rxn.elementary == True
    assert rxn.active == True


    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr,
                   delta_H=5, delta_S=-3)

    assert rxn.reactants == [(1, 'R')]
    assert rxn.products == [(1, 'P'), (1, 'Q')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'R': -1, 'P': 1, 'Q': 1})
    assert rxn.thermodynamics == ReactionThermodynamics(delta_H=5, delta_S=-3, delta_G=None, K=None)

    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn.reaction_category == "Unimolecular decomposition"
    assert rxn.elementary == True
    assert rxn.active == True


    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr,
                   delta_H=5, delta_S=-3, temp=100)

    assert rxn.reactants == [(1, 'R')]
    assert rxn.products == [(1, 'P'), (1, 'Q')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'R': -1, 'P': 1, 'Q': 1})

    assert rxn.thermodynamics.delta_H == 5
    assert rxn.thermodynamics.delta_S == -3
    assert np.allclose(rxn.thermodynamics.delta_G, 5.3)
    assert np.allclose(rxn.thermodynamics.K, 0.0017045829244452543)

    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn.reaction_category == "Unimolecular decomposition"
    assert rxn.elementary == True
    assert rxn.active == True


    with pytest.raises(Exception):
        # Inconsistent thermodynamic/kinetic data
        Reaction(reactants="R", products=["P", "Q"], species_registry=sr,
                   delta_H=5, delta_S=-3, temp=100,
                   kinetics_type="mass action", kinetic_parameters={"kF": 10, "kR": 2})


    # The thermodynamic data allows derivation of kinetic parameters not supplied
    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr,
                   delta_H=5, delta_S=-3, temp=100,
                   kinetics_type="mass action", kinetic_parameters={"kF": 10})

    assert rxn.reactants == [(1, 'R')]
    assert rxn.products == [(1, 'P'), (1, 'Q')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'R': -1, 'P': 1, 'Q': 1})

    assert rxn.thermodynamics.delta_H == 5
    assert rxn.thermodynamics.delta_S == -3
    assert np.allclose(rxn.thermodynamics.delta_G, 5.3)
    assert np.allclose(rxn.thermodynamics.K, 0.0017045829244452543)

    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn.reaction_category == "Unimolecular decomposition"
    assert rxn.elementary == True
    assert rxn.active == True

    assert rxn.kinetics.parameters["kF"] == 10
    assert np.allclose(rxn.kinetics.parameters["kR"], 5866.537706433048)
    assert rxn.kinetics.parameters["reversible"] == True


    # Consistent thermodynamic/kinetic data
    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr,
                   delta_H=5, delta_S=-3, temp=100,
                   kinetics_type="mass action", kinetic_parameters={"kF": 10, "kR": 5866.537706433048})

    assert rxn.reactants == [(1, 'R')]
    assert rxn.products == [(1, 'P'), (1, 'Q')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'R': -1, 'P': 1, 'Q': 1})

    assert rxn.thermodynamics.delta_H == 5
    assert rxn.thermodynamics.delta_S == -3
    assert np.allclose(rxn.thermodynamics.delta_G, 5.3)
    assert np.allclose(rxn.thermodynamics.K, 0.0017045829244452543)

    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn.reaction_category == "Unimolecular decomposition"
    assert rxn.elementary == True
    assert rxn.active == True

    assert rxn.kinetics.parameters["kF"] == 10
    assert np.allclose(rxn.kinetics.parameters["kR"], 5866.537706433048)
    assert rxn.kinetics.parameters["reversible"] == True


    # Reaction R -> 2 P
    rxn = Reaction(reactants="R", products=[(2, "P")], species_registry=sr)

    assert rxn.reactants == [(1, 'R')]
    assert rxn.products == [(2, 'P')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'R': -1, 'P': 2})
    assert rxn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K=None)

    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn.reaction_category == "Unimolecular decomposition"
    assert rxn.elementary == True
    assert rxn.active == True


    rxn = Reaction(reactants="R", products=["P", "P"], species_registry=sr)

    assert rxn.reactants == [(1, 'R')]
    assert rxn.products == [(2, 'P')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'R': -1, 'P': 2})
    assert rxn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K=None)
    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn.reaction_category == "Unimolecular decomposition"
    assert rxn.elementary == True
    assert rxn.active == True


def test_constructor_Reaction_3():

    sr = SpeciesRegistry(ids=["S", "P", "E"])

    # Enzymatic reaction
    rxn = Reaction(reactants=["S", "E"], products=["P", "E"], species_registry=sr)

    assert rxn.reactants == [(1, 'S'), (1, 'E')]
    assert rxn.products == [(1, 'P'), (1, 'E')]
    assert rxn.stoichiometry == Stoichiometry(coefficients={'S': -1, 'P': 1, 'E': 0})

    assert rxn.analytic_solution_family is None
    assert rxn.reaction_category == "Enzymatic"
    assert rxn.elementary == False
    assert rxn.active == True



def test_get_signed_stoichiometric_coefficients():
    sr = SpeciesRegistry(["A", "B"])
    rxn = Reaction(reactants="A", products="B", species_registry=sr)

    result = rxn.get_signed_stoichiometric_coefficients(reactants=[(1, "R")], products=[(1, "P")])
    assert result == {"R": -1, "P": 1}

    result = rxn.get_signed_stoichiometric_coefficients(reactants=[(1, "R"), (1, "S")], products=[(1, "P")])
    assert result == {"R": -1, "S": -1, "P": 1}

    result = rxn.get_signed_stoichiometric_coefficients(reactants=[(1, "R"), (1, "R")], products=[(1, "P")])
    assert result == {"R": -2, "P": 1}

    result = rxn.get_signed_stoichiometric_coefficients(reactants=[(1, "R")], products=[(1, "P"), (1, "Q")])
    assert result == {"R": -1, "P": 1, "Q": 1}

    result = rxn.get_signed_stoichiometric_coefficients(reactants=[(1, "E"), (1, "S")], products=[(1, "E"), (1, "P")])
    assert result == {"S": -1, "P": 1, "E": 0}

    result = rxn.get_signed_stoichiometric_coefficients(reactants=[(1, "A"), (2, "B"), (1, "E"), (1, "A")],
                                                        products=[(3, "P"), (1, "Q"), (1, "E")])
    assert result == {"A": -2, "B":- 2, "P": 3, "Q": 1, "E": 0}






############################  class Reaction  ############################

def test_CONSTRUCTOR_Reaction():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   delta_H=-3,
                   kinetic_parameters={"kF":10, "kR":2})
    assert rxn.active == True
    assert rxn.kinetics.parameters == {"kF":10, "kR":2, "K": 5, "reversible": True}
    assert rxn.thermodynamics.delta_H == -3
    assert rxn.thermodynamics.delta_S is None
    assert rxn.thermodynamics.K == 5


    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   active=False,
                   delta_H=-3,
                   kinetic_parameters={"kF":10, "kR":0})
    assert rxn.active == False
    assert rxn.kinetics.parameters == {"kF":10, "kR":0, "K": None, "reversible": False}
    assert rxn.thermodynamics.delta_H == -3
    assert rxn.thermodynamics.delta_S is None
    assert rxn.thermodynamics.K is None



def test_extract_rxn_properties():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   delta_H=-3000,
                   kinetic_parameters={"kF":10, "kR":2})

    assert rxn.extract_rxn_properties() == {'kinetics_type': 'mass action', 'kF': 10, 'kR': 2, 'delta_H': -3000, 'K': 5.0, 'reversible': True}



def test_set_thermodynamic_data():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   kinetic_parameters={"kF":6, "kR":2})
    assert rxn.thermodynamics.delta_H is None
    assert rxn.thermodynamics.delta_S is None
    assert rxn.thermodynamics.delta_G is None
    assert rxn.thermodynamics.K == 3

    rxn.set_thermodynamic_data(temp=100)
    assert np.allclose(rxn.thermodynamics.delta_G, -0.9134370805974775)


    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   delta_H=0.5, delta_S=-3)
    assert np.allclose(rxn.thermodynamics.delta_H, 0.5)
    assert rxn.thermodynamics.delta_S == -3
    assert rxn.thermodynamics.delta_G is None
    assert rxn.kinetics.parameters["kR"] == 0
    assert rxn.kinetics.parameters["kF"] == 0

    rxn.set_thermodynamic_data(temp=100)
    assert np.allclose(rxn.thermodynamics.delta_H, 0.5)
    assert rxn.thermodynamics.delta_S == -3
    assert np.allclose(rxn.thermodynamics.delta_G, 0.8)
    assert np.allclose(rxn.thermodynamics.K, 0.38205953171)
    assert rxn.kinetics.parameters["kR"] == 0
    assert rxn.kinetics.parameters["kF"] == 0



def test_extract_intermediate():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr, delta_H=-3)

    assert rxn.extract_intermediate() is None



def test_describe():
    sr = SpeciesRegistry(ids=["R", "P"])

    rxn = Reaction(reactants="R", products="P", species_registry=sr,
                   delta_H=0.5, delta_S=-3, temp=100,
                   kinetics_type="mass action", kinetic_parameters={"kF": 10})
    assert rxn.describe(concise=True) == "R <-> P"
    print(rxn.describe(concise=False))
    assert rxn.describe(concise=False) == \
            "R <-> P  Elementary Unimolecular rearrangement/isomerization reaction\n" \
            "         (delta_H = 0.5 kJ/mol | delta_S = -3 J/(mol·K) | delta_G = 0.8 kJ/mol | K = 0.38206 | Temp = -173.1 C | kinetics_type = 'mass action' | kF = 10 | kR = 26.174 | reversible = True)"



def test_extract_reactant_ids():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr)
    assert rxn.extract_reactant_ids() == ["A"]


    rxn = Reaction(reactants=["A", "B"], products="C", species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == ["A", "B"]

    rxn = Reaction(reactants=["R", "R"], products="C", species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == ["R"]

    rxn = Reaction(reactants=[(2, "R")], products="C", species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == ["R"]


    rxn = Reaction(reactants="A", products=["B", "C"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == ["A"]



def test_extract_reactants():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr)
    assert rxn.extract_reactants() == [(1, "A")]


def test_extract_reactants_formula():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr)
    assert rxn.extract_reactants_formula() == "A"



def test_extract_product_ids():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr)
    assert rxn.extract_product_ids() == ["B"]

    rxn = Reaction(reactants=["A", "B"], products="C", species_registry=sr, autoregister_species=True)
    assert rxn.extract_product_ids() == ["C"]

    rxn = Reaction(reactants="A", products=["B", "C"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_product_ids() == ["B", "C"]

    rxn = Reaction(reactants="A", products=["F", "F"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_product_ids() == ["F"]


def test_extract_products():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr)
    assert rxn.extract_products() == [(1, "B")]


def test_extract_products_formula():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr)
    assert rxn.extract_products_formula() == "B"



def test_extract_species_in_reaction():
    sr = SpeciesRegistry(ids=["A", "B", "C"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr)
    assert rxn.extract_species_in_reaction() == {"A", "B"}

    rxn = Reaction(reactants=["A", "B"], products="C", species_registry=sr)
    assert rxn.extract_species_in_reaction() == {"A", "B", "C"}

    rxn = Reaction(reactants="A", products=["B", "C"], species_registry=sr)
    assert rxn.extract_species_in_reaction() == {"A", "B", "C"}


def test_reaction_quotient():
    sr = SpeciesRegistry(ids=["A", "B"])

    # Reaction : A <-> B
    rxn = Reaction(reactants="A", products="B", species_registry=sr)
    c = {'A': 24., 'B': 36.}
    assert np.allclose(1.5, rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1.5, quotient)
    assert formula == '[B] / [A]'

    # Reaction : A <-> F
    sr.add_species("F")
    rxn = Reaction(reactants="A", products="F", species_registry=sr)
    c = {'A': 3., 'F': 33.}
    assert np.allclose(11., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(11., quotient)
    assert formula == '[F] / [A]'

    # Reaction :  A + B <-> C
    sr.add_species("C")
    rxn = Reaction(reactants=["A" , "B"], products="C", species_registry=sr)
    c = {'A': 3., 'B': 4., 'C': 12.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1., quotient)
    assert formula == '[C] / ([A][B])'

    # Reaction :  2A <-> P
    rxn = Reaction(reactants=["A" , "A"], products="P", species_registry=sr, autoregister_species=True)
    c = {'A': 2., 'P': 20.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(5., quotient)
    assert formula == '[P] /  [A]^2 '

    # Reaction :  C <-> A + B
    rxn = Reaction(reactants="C", products=["A" , "B"], species_registry=sr, autoregister_species=True)
    c = {'A': 3., 'B': 4., 'C': 12.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1., quotient)
    assert formula == '([A][B]) / [C]'

    # Reaction :  B <-> 2A
    rxn = Reaction(reactants="B", products=["A" , "A"], species_registry=sr, autoregister_species=True)
    c = {'A': 2., 'B': 20.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1/5., quotient)
    assert formula == ' [A]^2  / [B]'


def test_determine_reaction_rate():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   kinetic_parameters={"kF": 20., "kR": 2.})
    assert rxn.kinetics.law == "mass action"

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5. - 2. * 8.)  # 84.0

    # Now just the forward reaction
    rxn.kinetics.set_parameters({"kR": 0})
    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5.)            # 100.0


    # Reaction A + B -> C
    rxn = Reaction(reactants=["A", "B"], products="C",
                   species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 20})
    assert rxn.kinetics.parameters["reversible"] == False

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5. * 8.)

    # Now make reversible
    rxn.kinetics.set_parameters({"kR": 2})
    assert rxn.kinetics.parameters["reversible"] == True

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5. * 8. - 2. * 3.)


    # Reaction A <-> B + C
    rxn = Reaction(reactants="A", products=["B", "C"], species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 20})

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5.)

    # Make reversible
    rxn.kinetics.set_parameters({"kR": 2.})

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5.  - 2. * 8. * 3.)



def test_step_simulation():
    sr = SpeciesRegistry(ids=["A", "B"])

    # Reaction : A <-> B
    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   kinetic_parameters={"kF": 3., "kR": 2.})
    assert rxn.analytic_solution_family == "ONE_TO_ONE"

    # Euler approx
    result = rxn.step_simulation(delta_time=0.1, conc_dict={"A": 10, "B": 50})
    assert result[0] == {'A': 7, 'B': -7}
    assert result[1] == -70     # Rate = 3. * 10. - 2. * 50 .  Reaction is in reverse

    result = rxn.step_simulation(delta_time=0.8, conc_dict={"A": 10, "B": 50})
    assert result[0] == {'A': 56, 'B': -56}              # Note: these increments would make [B] negative!
    assert result[1] == -70

    # Exact solution
    result = rxn.step_simulation(delta_time=0.1, conc_dict={"A": 10, "B": 50}, exact=True)
    assert np.allclose(result[0]['A'],  5.508570764023133)
    assert np.allclose(result[0]['B'], -5.508570764023133)
    assert result[1] == -70

    result = rxn.step_simulation(delta_time=0.8, conc_dict={"A": 10, "B": 50}, exact=True)
    assert np.allclose(result[0]['A'],  13.74358105555772)  # Note: far more sensible than Euler method!
    assert np.allclose(result[0]['B'], -13.74358105555772)
    assert result[1] == -70


    # Reaction : A -> B
    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   kinetic_parameters={"kF": 3.})
    assert rxn.analytic_solution_family == "ONE_TO_ONE"
    assert rxn.kinetics.parameters["reversible"] == False

    # Euler approx
    result = rxn.step_simulation(delta_time=0.1, conc_dict={"A": 10, "B": 50})
    assert result[0] == {'A': -3, 'B': 3}
    assert result[1] == 30    # Rate = 3. * 10.     Reaction is now forward

    result = rxn.step_simulation(delta_time=0.4, conc_dict={"A": 10, "B": 50})
    assert result[0] == {'A': -12, 'B': 12}         # Note: these increments would make [A] negative!
    assert result[1] == 30

    # Exact solution
    result = rxn.step_simulation(delta_time=0.1, conc_dict={"A": 10, "B": 50}, exact=True)
    assert result[0] == {'A': -2.5918177931828215, 'B': 2.5918177931828215}
    assert result[1] == 30    # Rate = 3. * 10.

    result = rxn.step_simulation(delta_time=0.4, conc_dict={"A": 10, "B": 50}, exact=True)
    assert result[0] == {'A': -6.98805788087798, 'B': 6.98805788087798} # Note: far more sensible than Euler method!
    assert result[1] == 30    # Rate = 3. * 10.


    # Reaction : A + B <-> C
    rxn = Reaction(reactants=["A" , "B"], products="C", species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 5., "kR": 2})

    result = rxn.step_simulation(delta_time=0.002, conc_dict={"A": 10, "B": 50, "C": 20})
    assert result[0] == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert result[1] == 5*10*50 - 2 * 20        # 2460


    # Reaction : C <-> A + B
    rxn = Reaction(reactants="C", products=["A" , "B"], species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 2., "kR": 5.})

    result = rxn.step_simulation(delta_time=0.002, conc_dict={"A": 10, "B": 50, "C": 20})
    assert result[0] == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert result[1] == -2460



def test_find_equilibrium_conc():
    sr = SpeciesRegistry()

    # Reaction : A <-> C
    rxn = Reaction(reactants="A", products="C", species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 3, "kR": 2})
    assert rxn.analytic_solution_family == "ONE_TO_ONE"
    assert rxn.kinetics.law == "mass action"
    assert rxn.kinetics.parameters["kF"] == 3
    assert rxn.kinetics.parameters["kR"] == 2

    result = rxn.find_equilibrium_conc(conc_dict={"A":80., "C":10.})
    assert np.allclose(result["A"], 36)
    assert np.allclose(result["C"], 54)
    rxn.stoichiometry.consistency_checker(conc_before={"A":80., "C":10.}, conc_after={"A":36., "C":54.})


    # Only the forward reaction
    rxn = Reaction(reactants="A", products="C", species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 3})
    assert rxn.analytic_solution_family == "ONE_TO_ONE"
    assert rxn.kinetics.parameters["reversible"] == False
    result = rxn.find_equilibrium_conc(conc_dict={"A":80., "C":10.})
    assert np.allclose(result["A"], 0)
    assert np.allclose(result["C"], 90)

    # Only the reverse reaction
    rxn = Reaction(reactants="A", products="C", species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kR": 2})
    assert rxn.analytic_solution_family == "ONE_TO_ONE"
    result = rxn.find_equilibrium_conc(conc_dict={"A":80., "C":10.})
    assert np.allclose(result["A"], 90)
    assert np.allclose(result["C"], 0)

    # Reaction X + Y <-> Z
    rxn = Reaction(reactants=["X", "Y"], products="Z", species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 5, "kR": 2})
    assert rxn.analytic_solution_family == "TWO_TO_ONE"
    result = rxn.find_equilibrium_conc(conc_dict={"X":10., "Y": 50, "Z":20.})
    expected_eq = [0.2948774087575341, 40.294877408757536, 29.705122591242464]
    assert np.allclose(result["X"], expected_eq[0])
    assert np.allclose(result["Y"], expected_eq[1])
    assert np.allclose(result["Z"], expected_eq[2])
    rxn.stoichiometry.consistency_checker(conc_before={"X":10., "Y": 50, "Z": 20.},
                                          conc_after={"X":expected_eq[0], "Y": expected_eq[1], "Z": expected_eq[2]})

    with pytest.raises(Exception):
        rxn.find_equilibrium_conc(conc_dict={"X":10., "Z":20.})     # Missing reactant concentration

    with pytest.raises(Exception):
        rxn.find_equilibrium_conc(conc_dict={"X":10., "Y": 50})     # Missing product concentration


    # 2 A <-> C
    rxn = Reaction(reactants=["A", "A"], products="C", species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 3., "kR": 2.})
    assert rxn.analytic_solution_family == "TWO_TO_ONE"
    result = rxn.find_equilibrium_conc(conc_dict={"A":200., "C": 40.})
    expected_eq = [9.49568869375716, 135.2521556531214]
    assert np.allclose(result["A"], expected_eq[0])
    assert np.allclose(result["C"], expected_eq[1])
    rxn.stoichiometry.consistency_checker(conc_before={"A":200., "C": 40.},
                                          conc_after={"A":expected_eq[0], "C": expected_eq[1]})


    # Reaction Z <-> X + Y
    rxn = Reaction(reactants="Z", products=["X", "Y"], species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 2., "kR": 5.})
    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    result = rxn.find_equilibrium_conc(conc_dict={"X":10., "Y": 50, "Z":20.})
    expected_eq = [0.2948774087575341, 40.294877408757536, 29.705122591242464]
    assert np.allclose(result["X"], expected_eq[0])
    assert np.allclose(result["Y"], expected_eq[1])
    assert np.allclose(result["Z"], expected_eq[2])
    rxn.stoichiometry.consistency_checker(conc_before={"X":10., "Y": 50, "Z": 20.},
                                          conc_after={"X":expected_eq[0], "Y": expected_eq[1], "Z": expected_eq[2]})


    # C <-> 2 A
    rxn = Reaction(reactants="C", products=["A", "A"], species_registry=sr, autoregister_species=True,
                   kinetic_parameters={"kF": 2., "kR": 3.})
    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    result = rxn.find_equilibrium_conc(conc_dict={"C": 40., "A":200.})

    assert np.allclose(result["C"], 135.2521556531214)
    assert np.allclose(result["A"], 9.49568869375716)
    rxn.stoichiometry.consistency_checker(conc_before={"C": 40., "A":200.},
                                          conc_after={"C": 135.2521556531214, "A": 9.49568869375716})
