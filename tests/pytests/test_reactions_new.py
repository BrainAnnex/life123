import numpy as np
import pytest
import pandas as pd
from life123.species_registry import Species, SpeciesRegistry, MacroMolecules
from life123.reactions_new import Stoichiometry, Kinetics, ReactionThermodynamics, Reaction
from life123.visualization.colors import Colors
from tests.utilities.comparisons import *



########################  class Stoichiometry  ########################

def test_constructor_Stoichiometry():

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




########################  class Kinetics  ########################

def test_constructor_Kinetics():
    k = Kinetics(law="mass action")
    assert k.law == "mass action"

    with pytest.raises(Exception):
        Kinetics(law="never heard of this!")

    k = Kinetics(law="mass action", parameters={"kF": 10})
    assert k.law == "mass action"
    assert k.parameters == {"kF": 10, "kR": None, "reversible": True}
    #TODO: more tests



def test_to_dict_Kinetics():
    k = Kinetics(law="mass action")
    assert k.to_dict() == {"kinetics_type": "mass action", 'reversible': True}

    k = Kinetics(law="mass action", parameters={"kF": 10})
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 10, 'reversible': True}

    k = Kinetics(law="mass action", parameters={"kF": 10, "kR": 2})
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 10, "kR": 2, 'reversible': True}

    k = Kinetics(law="mass action", parameters={"kF": 10, "reversible": False})
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 10, "kR": 0, 'reversible': False}

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 10, "kR": 2, "reversible": False})    # Incompatible parameters

    with pytest.raises(Exception):
        Kinetics(law="random name")




########################  class ReactionThermodynamics  ########################

def test_constructor_ReactionThermodynamics():
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
                   kinetics_type="mass action", kinetic_parameters={"kF": 10, "kR": 5866.537706433048, "reversible": True})

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



def test_constructor_Reaction_2():

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
    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    rxn = Reaction(reactants="R", products="P", species_registry=sr)
    assert rxn.get_signed_stoichiometric_coefficients() == {"R": -1, "P": 1}

    rxn = Reaction(reactants=("R", "S"), products="P", species_registry=sr)
    assert rxn.get_signed_stoichiometric_coefficients() == {"R": -1, "S": -1, "P": 1}

    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.get_signed_stoichiometric_coefficients() == {"R": -1, "P": 1, "Q": 1}

    rxn = Reaction(reactants=("E", "S"), products=("E", "P"),
                   species_registry=sr, kinetics_type="Michaelis-Menten")
    assert rxn.get_signed_stoichiometric_coefficients() == {"S": -1, "P": 1, "E": 0}

    rxn = Reaction(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                   species_registry=sr, kinetics_type="custom")
    assert rxn.get_signed_stoichiometric_coefficients() == {"A": -2, "B":- 2, "P": 3, "Q": 1, "E": 0}



def test_get_reaction_vector():
    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    rxn = Reaction(reactants="R", products="P", species_registry=sr)
    assert rxn.get_reaction_vector() == {"R": -1, "P": 1}

    rxn = Reaction(reactants=("R", "S"), products="P", species_registry=sr)
    assert rxn.get_reaction_vector() == {"R": -1, "S": -1, "P": 1}

    rxn = Reaction(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.get_reaction_vector() == {"R": -1, "P": 1, "Q": 1}

    rxn = Reaction(reactants=("E", "S"), products=("E", "P"),
                   species_registry=sr, kinetics_type="Michaelis-Menten")
    assert rxn.get_reaction_vector() == {"S": -1, "P": 1}

    rxn = Reaction(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                   species_registry=sr, kinetics_type="custom")
    assert rxn.get_reaction_vector() == {"A": -2, "B":- 2, "P": 3, "Q": 1}



def test_ReactionElementary():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   delta_H=-3,
                   kinetic_parameters={"kF":10, "kR":2, "reversible": True})
    assert rxn.active == True
    assert rxn.kinetics.parameters == {"kF":10, "kR":2, "K": 5, "reversible": True}
    assert rxn.thermodynamics.delta_H == -3
    assert rxn.thermodynamics.delta_S is None
    assert rxn.thermodynamics.K == 5

    with pytest.raises(Exception):
        # Irreversible mass-action kinetics reactions can't have reverse rate
        Reaction(reactants="A", products="B", species_registry=sr,
                 kinetic_parameters={"kF":10, "kR":2, "reversible": False})


    rxn = Reaction(reactants="A", products="B", species_registry=sr,
                   active=False,
                   delta_H=-3,
                   kinetic_parameters={"kF":10, "kR":0, "reversible": False})
    assert rxn.active == False
    assert rxn.kinetics.parameters == {"kF":10, "kR":0, "reversible": False}
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
    assert rxn.kinetics.parameters["kR"] is None
    assert rxn.kinetics.parameters["kF"] is None

    rxn.set_thermodynamic_data(temp=100)
    assert np.allclose(rxn.thermodynamics.delta_H, 0.5)
    assert rxn.thermodynamics.delta_S == -3
    assert np.allclose(rxn.thermodynamics.delta_G, 0.8)
    assert np.allclose(rxn.thermodynamics.K, 0.38205953171)
    assert rxn.kinetics.parameters["kR"] is None
    assert rxn.kinetics.parameters["kF"] is None



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
            "R <-> P  Elementary Unimolecular rearrangement/isomerization\n" \
            "         (delta_H = 0.5 kJ/mol | delta_S = -3 J/(mol·K) | delta_G = 0.8 kJ/mol | K = 0.38206 | Temp = -173.1 C | kinetics_type = 'mass action' | kF = 10 | kR = 26.174 | reversible = True)"
