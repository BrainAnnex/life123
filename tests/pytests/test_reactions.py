import numpy as np
import pytest
from life123.species_registry import SpeciesRegistry
from life123.reaction_kinetics import ReactionKinetics
from life123.reactions import Stoichiometry, ReactionThermodynamics, \
ReactionDefinition, SimulationReaction, MassAction_Model, MichaelisMenten_Model
from tests.utilities.comparisons import *



def update_concentrations(conc, delta_conc) -> None:
    """
    Update the values of the dict `conc` based on the increments in `delta_conc` for the corresponding keys

    TODO: eventually move to one of the libraries

    :param conc:
    :param delta_conc:
    :return:            None
    """
    for k in conc:
        conc[k] += delta_conc.get(k, 0)     # Missing values default to zero




############################  class Stoichiometry  ############################

def test_CONSTRUCTOR_Stoichiometry():

    # Reaction R -> P
    st = Stoichiometry(vector={"R": -1, "P": 1})
    assert st.vector == {"R": -1, "P": 1}
    assert st.catalysts == []

    # Reaction R -> P + Q
    st = Stoichiometry(vector={"R": -1, "P": 1, "Q": 1})
    assert st.vector == {"R": -1, "P": 1, "Q": 1}
    assert st.catalysts == []

    # Reaction E + S -> P + E
    st = Stoichiometry(vector={"S": -1, "P": 1}, catalysts=["E"])
    assert st.vector == {"S": -1, "P": 1}
    assert st.catalysts == ["E"]


    with pytest.raises(Exception):
        Stoichiometry(vector={"R": -1, "P": 1, "E": 0})    # Zero coefficient



def test_to_dict():
    # Reaction R -> P
    st = Stoichiometry(vector={"R": -1, "P": 1})
    assert st.to_dict() == {"R": -1, "P": 1}

    # Reaction R -> P + Q
    st = Stoichiometry(vector={"R": -1, "P": 1, "Q": 1})
    assert st.to_dict() == {"R": -1, "P": 1, "Q": 1}

    # Reaction A + E -> B + E
    st = Stoichiometry(vector={"A": -1, "B": 1}, catalysts=["E"])
    assert st.to_dict() == {"A": -1, "B": 1, "E": 0}



def test_get_reaction_vector():
    # Reaction  2A + 3B + 2 E + F --> 4C  + 5D  + 2 E + F
    st = Stoichiometry({"A": -2, "B": -3, "C": 4, "D": 5}, catalysts=["E", "F"])
    assert st.get_reaction_vector() == {"A": -2, "B": -3, "C": 4, "D": 5}


    sr = SpeciesRegistry(ids=["A", "B", "C", "D", "R", "P", "Q", "S", "E", "F"])

    rxn = ReactionDefinition(reactants="R", products="P", species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"R": -1, "P": 1}

    rxn = ReactionDefinition(reactants=["R", "S"], products="P", species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"R": -1, "S": -1, "P": 1}

    rxn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"R": -1, "P": 1, "Q": 1}

    rxn = ReactionDefinition(reactants=["E", "S"], products=["E", "P"],
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"S": -1, "P": 1}

    rxn = ReactionDefinition(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"A": -2, "B":- 2, "P": 3, "Q": 1}


    # Reaction  2A + 3B + 2 E + F --> 4C  + 5D  + 2 E + F
    rxn = ReactionDefinition(reactants=["E", "A", "F", (3, "B"), "E", "A"], products=[(4, "C"), "F", (2, "E"), (5, "D")],
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"A": -2, "B": -3, "C": 4, "D": 5}



def test_get_reactant_list():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_reactant_list() == [(2, "A"), (2, "B")]

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_reactant_list() == [(2, "A"), (2, "B"), (1, "E")]

    sr = SpeciesRegistry(ids=["A", "B"])
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn_defn.stoichiometry.get_reactant_list() == [(1, "A")]

    sr = SpeciesRegistry(ids=["CH4", "O2", "CO2", "H2O"])
    rxn_defn = ReactionDefinition(reactants=["CH4", (2, "O2")],
                                  products=["CO2", (2, "H2O")], species_registry=sr)
    assert rxn_defn.stoichiometry.get_reactant_list() == [(1, "CH4"), (2, "O2")]


def test_get_reactant_ids():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_reactant_ids() == {"A", "B"}

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_reactant_ids() == {"A", "B", "E"}

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_reactant_ids(exclude_catalysts=True) == {"A", "B"}



def test_get_product_list():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_product_list() == [(3, "P")]

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_product_list() == [(3, "P"), (1, "E")]

    sr = SpeciesRegistry(ids=["A", "B"])
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn_defn.stoichiometry.get_product_list() == [(1, "B")]


    sr = SpeciesRegistry(ids=["CH4", "O2", "CO2", "H2O"])
    rxn_defn = ReactionDefinition(reactants=["CH4", (2, "O2")],
                                  products=["CO2", (2, "H2O")], species_registry=sr)
    assert rxn_defn.stoichiometry.get_product_list() == [(1, "CO2"), (2, "H2O")]


def test_get_product_ids():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_product_ids() == {"P"}

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_product_ids() == {"P", "E"}

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_product_ids(exclude_catalysts=True) == {"P"}



def test_get_all_species_ids():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_all_species_ids() == {"A", "B", "P"}

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_all_species_ids() == {"A", "B", "P", "E"}



def test_get_reaction_complexes():
    st = Stoichiometry(vector={"A": -1, "P": 2, "Q": 1}, catalysts=["E"])

    assert st.get_reaction_complexes() == ( {"A": 1, "E": 1} ,  {"P": 2, "Q": 1, "E": 1} )


    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    rxn = ReactionDefinition(reactants="R", products="P", species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"R": 1}, {"P": 1})

    rxn = ReactionDefinition(reactants=["R", "S"], products="P", species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"R": 1, "S": 1}, {"P": 1})

    rxn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"R": 1}, {"P": 1, "Q": 1})

    rxn = ReactionDefinition(reactants=["E", "S"], products=["E", "P"],
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"S": 1, "E": 1}, {"P": 1, "E": 1})

    rxn = ReactionDefinition(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"A": 2, "B": 2, "E": 1}, {"P": 3, "Q": 1, "E": 1})



def test__standard_form_complex():
    pass    # TODO


def test_standard_chemical_formula():
    st = Stoichiometry(vector={"A": -1, "P": 2, "Q": 1}, catalysts=["E"])
    assert st.standard_chemical_formula() == "A + E -> 2 P + Q + E"


    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    rxn = ReactionDefinition(reactants="R", products="P", species_registry=sr)
    assert rxn.stoichiometry.standard_chemical_formula() == "R -> P"

    rxn = ReactionDefinition(reactants=["R", "S"], products="P", species_registry=sr)
    assert rxn.stoichiometry.standard_chemical_formula() == "R + S -> P"

    rxn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.stoichiometry.standard_chemical_formula() == "R -> P + Q"

    rxn = ReactionDefinition(reactants=["E", "S"], products=["E", "P"],
                             species_registry=sr)
    assert rxn.stoichiometry.standard_chemical_formula() == "S + E -> P + E"

    rxn = ReactionDefinition(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                             species_registry=sr)
    assert rxn.stoichiometry.standard_chemical_formula(reversible=True) == "2 A + 2 B + E <-> 3 P + Q + E"



def test_reaction_pattern():
    # Reaction R -> P
    st = Stoichiometry(vector={"R": -1, "P": 1})
    assert st.reaction_pattern() == (1, 1, 0)

    # Reaction R -> P + Q
    st = Stoichiometry(vector={"R": -1, "P": 1, "Q": 1})
    assert st.reaction_pattern() == (1, 2, 0)

    # Reaction A + E -> B + E
    st = Stoichiometry(vector={"A": -1, "B": 1}, catalysts=["E"])
    assert st.reaction_pattern() == (1, 1, 1)
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])



def test_consistency_checker():
    # Reaction  A -> B
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

    # Reaction 2A -> B
    st = Stoichiometry({"A": -2, "B": 1})

    st.consistency_checker(conc_before={"A": 0, "B": 0}, conc_after={"A": 0, "B": 0})
    st.consistency_checker(conc_before={"A": 100, "B": 0}, conc_after={"A": 80, "B": 10})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 0}, conc_after={"A": 80, "B": 10.1})

    st.consistency_checker(conc_before={"A": 0, "B": 50}, conc_after={"A": 10, "B": 45})

    # Reaction A + B --> C
    st = Stoichiometry({"A": -1, "B": -1, "C": 1})

    st.consistency_checker(conc_before={"A": 0, "B": 0, "C": 0}, conc_after={"A": 0, "B": 0, "C": 0})
    st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 40, "C": 10})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 40, "C": 9.9})

    # Reaction A + 3B --> C
    st = Stoichiometry({"A": -1, "B": -3, "C": 1})

    st.consistency_checker(conc_before={"A": 0, "B": 0, "C": 0}, conc_after={"A": 0, "B": 0, "C": 0})
    st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 20, "C": 10})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 20, "C": 10.1})

    # Reaction  A + 3B --> 4C
    st = Stoichiometry({"A": -1, "B": -3, "C": 4})

    st.consistency_checker(conc_before={"A": 0, "B": 0, "C": 0}, conc_after={"A": 0, "B": 0, "C": 0})
    st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 20, "C": 40})
    with pytest.raises(Exception):
        st.consistency_checker(conc_before={"A": 100, "B": 50, "C": 0}, conc_after={"A": 90, "B": 20, "C": 39.9})

    # Reaction  2A + 3B --> 4C + 5D
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
    rt = ReactionThermodynamics(delta_H=3, delta_S=-200, delta_G=1000, K_eq=0.1)
    assert rt.delta_H == 3
    assert rt.delta_S == -200
    assert rt.delta_G == 1000
    assert np.allclose(rt.K_eq, 0.1)



def test_to_dict_ReactionThermodynamics():
    rt = ReactionThermodynamics(delta_H=3, delta_S=-200, delta_G=1000, K_eq=0.1)
    assert rt.to_dict() == {'delta_H': 3, 'delta_S': -200, 'delta_G': 1000, 'K_eq': 0.1}

    rt.K_eq = None
    assert rt.to_dict() == {'delta_H': 3, 'delta_S': -200, 'delta_G': 1000}

    rt = ReactionThermodynamics(delta_G=1000)
    assert rt.to_dict() == {'delta_G': 1000}

    rt = ReactionThermodynamics()
    assert rt.to_dict() == {}



def test_get_signed_stoichiometric_coefficients():
    sr = SpeciesRegistry(["A", "B"])
    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr)

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





############################  class SimulationReaction  ############################

def test_CONSTRUCTOR_SimulationReaction():
    pass        # TODO




def test_determine_reaction_rate_1():
    sr = SpeciesRegistry(ids=["A", "B"])

    # Reaction A <-> B , with "mass action"
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF": 20., "kR": 2.})

    sim_rxn = rxn_defn.sim_reactions[0]
    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5. - 2. * 8.)  # 84.0

    # Now just the forward reaction
    sim_rxn.model.set_parameters(parameters={"kR": 0, "K": None})
    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5.)            # 100.0


    # Reaction A + B -> C , with "mass action"
    rxn_defn = ReactionDefinition(reactants=["A", "B"], products="C",
                             species_registry=sr, autoregister_species=True,
                             reaction_model="mass action", kinetic_parameters={"kF": 20})

    sim_rxn = rxn_defn.sim_reactions[0]
    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5. * 8.)

    # Now make reversible
    sim_rxn.model.set_parameters({"kR": 2})
    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5. * 8. - 2. * 3.)


    # Reaction A -> B + C , with "mass action"
    rxn_defn = ReactionDefinition(reactants="A", products=["B", "C"], species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 20})

    sim_rxn = rxn_defn.sim_reactions[0]
    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5.)

    # Make reversible
    sim_rxn.model.set_parameters({"kR": 2.})
    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 3})
    assert np.allclose(result, 20. * 5.  - 2. * 8. * 3.)


    # Reaction  2A <-> B , with "mass action"
    rxn_defn = ReactionDefinition(reactants=[(2, "A")], products="B", species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 5, "kR": 2})

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.stoichiometry == Stoichiometry({"A": -2, "B": 1})

    init_conc = {"A": 4.5, "B": 6.}
    result = sim_rxn.determine_reaction_rate(conc_dict=init_conc)
    assert np.allclose(result, 5. * 4.5 **2 - 2. * 6.)  # 89.25

    # Make irreversible
    sim_rxn.model.set_parameters({"kR": 0, "K": None})     # K was previously automatically set
    result = sim_rxn.determine_reaction_rate(conc_dict=init_conc)
    assert np.allclose(result, 5. * 4.5 **2)            # 101.25

    # Change the kinetic parameters and the concentration values
    sim_rxn.model.set_parameters({"kF": 3, "kR": 2, "K": None})
    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 59.)         # 3. * 5**2 - 2. * 8.


    # Reaction  B <-> 2C , with "mass action"
    rxn_defn = ReactionDefinition(reactants="B", products=[(2, "C")], species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 4, "kR": 2})

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.stoichiometry == Stoichiometry({"B": -1, "C": 2})
    result = sim_rxn.determine_reaction_rate(conc_dict={"B": 5., "C": 4})
    assert np.allclose(result, 4. * 5. - 2. * 4. **2)           # -12.0



def kinetic_constant_rate(stoichiometry, kinetic_parameters, conc_dict):
    return 8

def test_determine_reaction_rate_2_a():
    # Reaction: # A + B -> C + D , with custom reaction model
    sr = SpeciesRegistry(ids=["A", "B", "C", "D"])
    rxn_defn = ReactionDefinition(id=49, reactants=["A", "B"], products=["C", "D"], species_registry=sr,
                                  reaction_model="custom",
                                  kinetic_parameters={"kF": 10, "rate_function": ReactionKinetics.kinetic_rate_first_order})
     #print(rxn_defn.describe(concise=False))

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.model.rate_function.__name__ == "kinetic_rate_first_order"

    initial_conc = {"A": 2, "B": 4, "C": 5, "D": 3}
    result = sim_rxn.determine_reaction_rate(conc_dict=initial_conc)
    assert result == 80       #  10. * 2 * 4   (no reverse reaction)

    # Make reversible
    sim_rxn.model.set_parameters({"kR": 2})
    result = sim_rxn.determine_reaction_rate(conc_dict=initial_conc)
    assert np.allclose(result, 80. - 2. * 5 * 3)


    # Reaction: # A + B -> C + D , with custom reaction model
    sr = SpeciesRegistry(ids=["A", "B", "C", "D"])
    rxn_defn = ReactionDefinition(id=49, reactants=["A", "B"], products=["C", "D"], species_registry=sr,
                             reaction_model="custom",
                             kinetic_parameters={"kF": 10})     # rate_function not provided
    #print(rxn_defn.describe(concise=False))

    sim_rxn = rxn_defn.sim_reactions[0]

    initial_conc = {"A": 2, "B": 4, "C": 5, "D": 3}
    with pytest.raises(Exception):
        sim_rxn.determine_reaction_rate(conc_dict=initial_conc)     # Missing rate_function

    sim_rxn.set_parameters({'rate_function': kinetic_constant_rate})
    result = sim_rxn.determine_reaction_rate(conc_dict=initial_conc)
    assert result == 8


def test_determine_reaction_rate_2_b():
    # Reaction A <-> B , with custom reaction model
    sr = SpeciesRegistry(ids=["A", "B"])
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                                  reaction_model="custom",
                                  kinetic_parameters={"kF": 20, "kR": 2, "rate_function": ReactionKinetics.kinetic_rate_first_order})

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.model.rate_function.__name__ == "kinetic_rate_first_order"

    initial_conc = {"A": 5., "B": 8.}
    result = sim_rxn.determine_reaction_rate(conc_dict=initial_conc)
    assert np.allclose(result, 20. * 5. - 2. * 8.)  # 84.0

    # Make irreversible
    sim_rxn.set_parameters({"kR": 0, "K": None})     # K was previously automatically set to 10.0

    result = sim_rxn.determine_reaction_rate(conc_dict=initial_conc)
    assert np.allclose(result, 20. * 5.)            # 100.0


    # Switch to using the following test function for the reaction rate
    def my_rate_law_1(stoichiometry, kinetic_parameters, conc_dict) -> float:
        return 1234.

    sim_rxn.set_parameters({"rate_function": my_rate_law_1})
    assert sim_rxn.model.rate_function == my_rate_law_1

    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 1234.)


    # Reaction 5A <-> 2B , with custom reaction model (hypothetically with 1st-order kinetics with respect to each species)
    rxn_defn = ReactionDefinition(reactants=[(5, "A")], products=[(2, "B")], species_registry=sr,
                                  reaction_model="custom",
                                  kinetic_parameters={"kF": 20, "kR": 2, "rate_function": ReactionKinetics.kinetic_rate_first_order})

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.model.rate_function.__name__ == "kinetic_rate_first_order"

    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5. - 2. * 8.)      # 84.0


    # Reaction 2B <-> 3C , with custom reaction model (hypothetically with 1st-order kinetics with respect to each species)
    rxn_defn = ReactionDefinition(reactants=[(2, "B")], products=[(3, "C")], species_registry=sr, autoregister_species=True,
                                  reaction_model="custom",
                                  kinetic_parameters={"kF": 10, "kR": 25, "rate_function": ReactionKinetics.kinetic_rate_first_order})

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.model.rate_function.__name__ == "kinetic_rate_first_order"

    result = sim_rxn.determine_reaction_rate(conc_dict={"B": 8., "C": 15.})
    assert np.allclose(result,  10. * 8. - 25. * 15.)   # -295.0


    # Reaction 2A + 5B <-> 4C + 3D , with custom reaction model (hypothetically with 1st-order kinetics with respect to each species)
    rxn_defn = ReactionDefinition(reactants=[(2, "A") , (5, "B")],
                                  products=[(4, "C") , (3, "D")],
                                  species_registry=sr, autoregister_species=True,
                                  reaction_model="custom",
                                  kinetic_parameters={"kF": 5, "kR": 2, "rate_function": ReactionKinetics.kinetic_rate_first_order})

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.model.rate_function.__name__ == "kinetic_rate_first_order"

    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 3.5, "B": 9., "C": 11., "D": 7.})
    assert np.allclose(result,  5. * 3.5 * 9. - 2. * 11. * 7.)  # 3.5

    result = sim_rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8., "C": 15., "D": 7.})
    assert np.allclose(result,  -10.)



def test_determine_reaction_rate_3():
    # E + S -> E + P, with "michaelis menten" model
    sr = SpeciesRegistry(ids=["S", "P", "E"])
    rxn_defn = ReactionDefinition(id=22, reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                                  reaction_model="michaelis menten",
                                  kinetic_parameters={'kM': 2, 'kcat': 5})
    #print(rxn_defn.describe(concise=False))

    sim_rxn = rxn_defn.sim_reactions[0]

    initial_conc = {"E": 0.1, "S": 10, "P": 2}
    result = sim_rxn.determine_reaction_rate(conc_dict=initial_conc)
    assert math.isclose(result, 5/12)       # (5 * 0.1 * 10) / (2 + 10)



def test_step_simulation_1():
    sr = SpeciesRegistry(ids=["A", "B"])

    # Reaction : A <-> B
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                                  reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})
    assert rxn_defn.analytic_solution_family == "ONE_TO_ONE"

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.analytic_solution_family == "ONE_TO_ONE"

    # Euler approx
    result = sim_rxn.step_simulation(delta_time=0.1, conc_dict={"A": 10, "B": 50})
    assert result[0] == {'A': 7, 'B': -7}
    assert result[1] == -70     # Rate = 3. * 10. - 2. * 50 .  Reaction is in reverse

    result = sim_rxn.step_simulation(delta_time=0.8, conc_dict={"A": 10, "B": 50})
    assert result[0] == {'A': 56, 'B': -56}              # Note: these increments would make [B] negative!
    assert result[1] == -70

    # Exact solution
    result = sim_rxn.step_simulation(delta_time=0.1, conc_dict={"A": 10, "B": 50}, exact=True)
    assert np.allclose(result[0]['A'],  5.508570764023133)
    assert np.allclose(result[0]['B'], -5.508570764023133)
    assert result[1] == -70

    result = sim_rxn.step_simulation(delta_time=0.8, conc_dict={"A": 10, "B": 50}, exact=True)
    assert np.allclose(result[0]['A'],  13.74358105555772)  # Note: far more sensible than Euler method!
    assert np.allclose(result[0]['B'], -13.74358105555772)
    assert result[1] == -70


    # Reaction : A -> B  (no reverse reaction)
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                                  reaction_model="mass action", kinetic_parameters={"kF": 3.})
    assert rxn_defn.analytic_solution_family == "ONE_TO_ONE"

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.analytic_solution_family == "ONE_TO_ONE"

    # Euler approx
    result = sim_rxn.step_simulation(delta_time=0.1, conc_dict={"A": 10, "B": 50})
    assert result[0] == {'A': -3, 'B': 3}
    assert result[1] == 30    # Rate = 3. * 10.     Reaction is now forward

    result = sim_rxn.step_simulation(delta_time=0.4, conc_dict={"A": 10, "B": 50})
    assert result[0] == {'A': -12, 'B': 12}         # Note: these increments would make [A] negative!
    assert result[1] == 30

    # Exact solution
    result = sim_rxn.step_simulation(delta_time=0.1, conc_dict={"A": 10, "B": 50}, exact=True)
    assert result[0] == {'A': -2.5918177931828215, 'B': 2.5918177931828215}
    assert result[1] == 30    # Rate = 3. * 10.

    result = sim_rxn.step_simulation(delta_time=0.4, conc_dict={"A": 10, "B": 50}, exact=True)
    assert result[0] == {'A': -6.98805788087798, 'B': 6.98805788087798} # Note: far more sensible than Euler method!
    assert result[1] == 30    # Rate = 3. * 10.


    # Reaction : A + B <-> C
    rxn_defn = ReactionDefinition(reactants=["A" , "B"], products="C", species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 5., "kR": 2})

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.analytic_solution_family == "TWO_TO_ONE"
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"A": -1, "B": -1, "C": 1})

    conc_initial = {"A": 10, "B": 50, "C": 20}
    delta, rate = sim_rxn.step_simulation(delta_time=0.002, conc_dict=conc_initial)
    assert delta == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert rate == 5 * 10 * 50 - 2 * 20        # 2460 , i.e.  kF [A] [B] - kR [C]
    for k in delta:
        assert delta[k] == rate * 0.002 * np.sign(sim_rxn.stoichiometry.to_dict()[k])


    # Reaction : C <-> A + B
    rxn_defn = ReactionDefinition(reactants="C", products=["A" , "B"], species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 2., "kR": 5.})

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.analytic_solution_family == "ONE_TO_TWO"

    result = sim_rxn.step_simulation(delta_time=0.002, conc_dict={"A": 10, "B": 50, "C": 20})
    assert result[0] == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert result[1] == -2460


def test_step_simulation_2():

    # Reaction: # E + S <-> ES* -> E + P, with SingleSubstrateMechanism model
    sr = SpeciesRegistry(ids=["S", "P", "E"])
    rxn_defn = ReactionDefinition(id=85, reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                             reaction_model="single substrate mechanism",
                             kinetic_parameters={"k1_F": 18, "k1_R": 100, "k2_F": 49})

    initial_conc = {"E": 1, "S": 20, "P": 0, "ES*": 0}
    dt = 0.002

    # Look at each derived reaction individually
    rxn_sim_tuple = rxn_defn.sim_reactions
    assert len(rxn_sim_tuple) == 2
    rxn_sim_1, rxn_sim_2 = rxn_sim_tuple

    incr_dict_1, rate_1 = rxn_sim_1.step_simulation(delta_time=dt, conc_dict=initial_conc)
    assert rate_1 == 360       # 18 * 1 * 20 - 100 * 0
    assert incr_dict_1 == {'E': -0.72, 'S': -0.72, 'ES*': 0.72}


    incr_dict_2, rate_2 = rxn_sim_2.step_simulation(delta_time=dt, conc_dict=initial_conc)
    assert rate_2 == 0          # 49 * 0
    assert incr_dict_2 == {'E': 0, 'ES*': 0, 'P': 0}

    
    # Simulating by hand each of the 2 sub-reactions will give the same results
    upstream_rxn_defn = ReactionDefinition(reactants=["E", "S"], products="ES*", species_registry=sr,
                             reaction_model="mass action",
                             kinetic_parameters={"kF": 18, "kR": 100})
    upstream_rxn_sim = upstream_rxn_defn.sim_reactions[0]
    assert upstream_rxn_sim.step_simulation(delta_time=dt, conc_dict=initial_conc) == ( {'E': -0.72, 'S': -0.72, 'ES*': 0.72} , 360 )

    downstream_rxn_defn = ReactionDefinition(reactants="ES*", products=["E", "P"], species_registry=sr,
                             reaction_model="mass action",
                             kinetic_parameters={"kF": 49})
    downstream_rxn_sim = downstream_rxn_defn.sim_reactions[0]
    assert downstream_rxn_sim.step_simulation(delta_time=dt, conc_dict=initial_conc) == ( {'E': 0, 'ES*': 0, 'P': 0} , 0 )


    # Manually advance the reaction simulation by this step, plus one more
    conc = initial_conc

    # Update the system concentrations (thus advancing the simulation)
    update_concentrations(conc, incr_dict_1)
    update_concentrations(conc, incr_dict_2)

    assert conc == {'E': 0.28, 'S': 19.28, 'P': 0.0, 'ES*': 0.72}

    incr_dict_1, rate_1 = upstream_rxn_sim.step_simulation(delta_time=dt, conc_dict=conc)
    assert math.isclose(rate_1, 25.1712)    # 18 * 0.28 * 19.28 - 100 * 0.72
    expected_incr = {'E': -0.0503424, 'S': -0.0503424, 'ES*': 0.0503424}    # Delta_conc = 25.1712 * 0.002 = 0.0503424
    compare_dicts(incr_dict_1, expected_incr)

    incr_dict_2, rate_2 = downstream_rxn_sim.step_simulation(delta_time=dt, conc_dict=initial_conc)
    assert math.isclose(rate_2, 35.28)      # 49 * 0.72
    expected_incr = {'ES*': -0.07056, 'E': 0.07056, 'P': 0.07056}           # Delta_conc = 35.28 * 0.002 = 0.07056
    compare_dicts(incr_dict_2, expected_incr)

    # Update the system concentrations (thus advancing the simulation)
    update_concentrations(conc, incr_dict_1)
    update_concentrations(conc, incr_dict_2)

    expected = {'E': 0.3002176, 'S': 19.2296576, 'P': 0.07056, 'ES*': 0.6997824}
    compare_dicts(conc, expected)


def test_step_simulation_3():

    # Reaction: # A + B -> C + D , with custom reaction model
    sr = SpeciesRegistry(ids=["A", "B", "C", "D"])
    rxn_defn = ReactionDefinition(id=49, reactants=["A", "B"], products=["C", "D"], species_registry=sr,
                             reaction_model="custom",
                             kinetic_parameters={"kF": 10, "rate_function": ReactionKinetics.kinetic_rate_first_order})
    #print(rxn_defn.describe(concise=False))

    sim_rxn = rxn_defn.sim_reactions[0]

    initial_conc = {"A": 2, "B": 4, "C": 5, "D": 3}
    dt = 0.1

    incr_dict, rate = sim_rxn.step_simulation(delta_time=dt, conc_dict=initial_conc)
    assert rate == 80       #  10. * 2 * 4   (no reverse reaction)
    assert incr_dict == {'A': -8, 'B': -8, 'C': 8, 'D': 8}      # 80 * 0.1 = 8

    # Make reversible
    sim_rxn.set_parameters({"kR": 2})
    incr_dict, rate = sim_rxn.step_simulation(delta_time=dt, conc_dict=initial_conc)
    assert rate == 50      # 80 - 2. * 5 * 3
    assert incr_dict == {'A': -5, 'B': -5, 'C': 5, 'D': 5}      # 50 * 0.1 = 5





############################################################################



# TODO: also test the various "Compiler" classes and ReactionModelRegistry




############################  class ReactionDefinition  ############################

def test_CONSTRUCTOR_ReactionDefinition_1():
    # The group of testing below is about parsing the reactant/products,
    # and building the reaction "Stoichiometry" dataclass.
    # No kinetics and no thermodynamics!

    sr = SpeciesRegistry(ids=["R", "P", "Q"])

    with pytest.raises(Exception):
        ReactionDefinition(reactants="X", products="P", species_registry=sr)  # Un-registered reactant

    with pytest.raises(Exception):
        ReactionDefinition(reactants="R", products="Y", species_registry=sr)  # Un-registered product

    # Missing products or reactants
    with pytest.raises(Exception):
        ReactionDefinition(reactants=["R"], products=None, species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants=None, products="P", species_registry=sr)

    # Bad products or reactants
    with pytest.raises(Exception):
        ReactionDefinition(reactants={"k": 666}, products="P", species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants="R", products=123, species_registry=sr)


    # Reactants and the products can't be the same
    with pytest.raises(Exception):
        ReactionDefinition(reactants=["A"], products=["A"], species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants=["A"], products=[("A")], species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants=["A"], products=[(1, "A")], species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants="R", products="R", species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants=[(2, "B")], products=[(2, "B")], species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants=["A", "B"], products=["A", "B"], species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants=["A", (3, "B")], products=["A", (3, "B")], species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants=["A", "B"], products=["B", "A"], species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants=[(2, "A"), "B", "C"], products=["B", (1, "C"), (2, "A")], species_registry=sr)
    with pytest.raises(Exception):
        ReactionDefinition(reactants=["R", (2, "P")], products=[(2, "P"), "R"], species_registry=sr)


    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    # Reaction R -> P
    rxn_defn = ReactionDefinition(reactants="R", products="P", species_registry=sr)
    assert rxn_defn.species_registry == sr
    assert rxn_defn.stoichiometry.to_dict() == {"R": -1, "P": 1}
    assert rxn_defn.reaction_model is None
    assert rxn_defn.sim_reactions == ()
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K_eq=None, derived_pars=set())

    # Reaction R + S -> P
    rxn_defn = ReactionDefinition(reactants=["R", "S"], products="P", species_registry=sr)
    assert rxn_defn.stoichiometry.to_dict() == {"R": -1, "S": -1, "P": 1}

    # Reaction R -> P + Q
    rxn_defn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn_defn.stoichiometry.to_dict() == {"R": -1, "P": 1, "Q": 1}
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K_eq=None, derived_pars=set())
    assert rxn_defn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn_defn.reaction_category == "Unimolecular decomposition"

    # Reaction R -> 2 P
    rxn_defn = ReactionDefinition(reactants="R", products=["P", "P"], species_registry=sr)
    assert rxn_defn.stoichiometry.to_dict() == {"R": -1, "P": 2}
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K_eq=None, derived_pars=set())
    assert rxn_defn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn_defn.reaction_category == "Unimolecular decomposition"

    rxn_defn = ReactionDefinition(reactants="R", products=[(2, "P")], species_registry=sr)
    assert rxn_defn.stoichiometry.to_dict() == {"R": -1, "P": 2}

    # Reaction 2 R -> P
    rxn_defn = ReactionDefinition(reactants=["R", "R"], products="P", species_registry=sr)
    assert rxn_defn.stoichiometry.to_dict() == {"R": -2, "P": 1}

    rxn_defn = ReactionDefinition(reactants=["E", "S"], products=["E", "P"],
                                  species_registry=sr)
    assert rxn_defn.stoichiometry.to_dict() == {"S": -1, "P": 1, "E": 0}

    rxn_defn = ReactionDefinition(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                                  species_registry=sr)
    assert rxn_defn.stoichiometry.to_dict() == {"A": -2, "B":- 2, "P": 3, "Q": 1, "E": 0}


    # Reaction E + S -> P + E
    rxn_defn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr)
    assert rxn_defn.stoichiometry == Stoichiometry(vector={'S': -1, 'P': 1}, catalysts=["E"])
    assert rxn_defn.analytic_solution_family is None
    assert rxn_defn.reaction_category == "Enzymatic"



def test_CONSTRUCTOR_ReactionDefinition_2():
    # No thermodynamic data passed

    # R -> P, with "mass action"
    sr = SpeciesRegistry()
    rxn_defn = ReactionDefinition(id=942, reactants="R", products="P", species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 20, "kR": 4})
    assert rxn_defn.stoichiometry == Stoichiometry(vector={'R': -1, 'P': 1})
    assert rxn_defn.source_kinetic_parameters == {"kF": 20, "kR": 4}
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K_eq=5, derived_pars={'K_eq'})
    assert rxn_defn.analytic_solution_family == "ONE_TO_ONE"
    assert rxn_defn.reaction_category == "Unimolecular rearrangement/isomerization"

    sim_rxn_tuple = rxn_defn.sim_reactions
    assert len(sim_rxn_tuple) == 1
    sim_rxn = sim_rxn_tuple[0]
    assert type(sim_rxn) == SimulationReaction
    assert type(sim_rxn.model) == MassAction_Model
    assert sim_rxn.source_object.id == 942
    assert sim_rxn.derivation == "direct"
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"R": -1, "P": 1})
    assert sim_rxn.model.get_parameters() == {'kF': 20, 'kR': 4, 'K': 5.0, 'reversible': True}
    assert sim_rxn.model.derived_pars == {'K', 'reversible'}

    assert set(sr.get_all_species_ids()) == {"R", "P"}


    # A -> B, with "mass action"
    sr = SpeciesRegistry(ids=["A", "B"])
    rxn_defn = ReactionDefinition(id=8, reactants="A", products="B", species_registry=sr,
                    reaction_model="mass action",
                    kinetic_parameters={"kF": 10, "kR": 2})
    assert rxn_defn.species_registry == sr
    assert rxn_defn.stoichiometry == Stoichiometry(vector={"A": -1, "B": 1})
    assert rxn_defn.reaction_model == "mass action"
    assert rxn_defn.analytic_solution_family == "ONE_TO_ONE"
    assert rxn_defn.reaction_category == "Unimolecular rearrangement/isomerization"
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K_eq=5, derived_pars={'K_eq'})
    assert rxn_defn.source_kinetic_parameters == {"kF": 10, "kR": 2}

    sim_rxn_tuple = rxn_defn.sim_reactions
    assert len(sim_rxn_tuple) == 1
    sim_rxn = sim_rxn_tuple[0]
    assert type(sim_rxn) == SimulationReaction
    assert type(sim_rxn.model) == MassAction_Model
    assert sim_rxn.source_object.id == 8
    assert sim_rxn.derivation == "direct"
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"A": -1, "B": 1})
    assert sim_rxn.model.get_parameters() == {'kF': 10, 'kR': 2, 'K': 5.0, 'reversible': True}


    with pytest.raises(Exception):
        ReactionDefinition(reactants="A", products="B", species_registry=sr,
                           reaction_model="mass action", kinetic_parameters={"intruder":666})


    # E + S -> E + P, with "michaelis menten" model
    sr = SpeciesRegistry(ids=["S", "P", "E"])
    rxn_defn = ReactionDefinition(id=17, reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                                  reaction_model="michaelis menten",
                                  kinetic_parameters={'kM': 2, 'kcat': 5})
    assert rxn_defn.species_registry == sr
    assert rxn_defn.stoichiometry == Stoichiometry(vector={"S": -1, "P": 1}, catalysts=["E"])
    assert rxn_defn.reaction_model == "michaelis menten"
    assert rxn_defn.source_kinetic_parameters == {'kM': 2, 'kcat': 5}

    sim_rxn_tuple = rxn_defn.sim_reactions
    assert len(sim_rxn_tuple) == 1
    sim_rxn = sim_rxn_tuple[0]
    assert type(sim_rxn) == SimulationReaction
    assert sim_rxn.source_object.id == 17
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"S": -1, "P": 1}, catalysts=["E"])
    assert type(sim_rxn.model) == MichaelisMenten_Model
    assert sim_rxn.model.get_parameters() == {'kM': 2, 'kcat': 5, 'Enzyme': 'E', 'Product': 'P', 'Substrate': 'S'}

    with pytest.raises(Exception):
        ReactionDefinition(reactants=["S", "X"], products=["P", "E"],
                           species_registry=sr, autoregister_species=True,
                           reaction_model="michaelis menten")       # Bad stoichiometry

    with pytest.raises(Exception):
        ReactionDefinition(reactants=["S", "E", "E2"], products=["P", "E", "E2"],
                           species_registry=sr, autoregister_species=True,
                           reaction_model="michaelis menten")           # Bad stoichiometry (too many enzymes)

    with pytest.raises(Exception):
        ReactionDefinition(reactants=["A", "B", "E"], products=["P", "E"],
                           species_registry=sr, autoregister_species=True,
                           reaction_model="michaelis menten")       # Bad stoichiometry

    with pytest.raises(Exception):
        ReactionDefinition(reactants=["S", "E"], products=["E"],
                           species_registry=sr, autoregister_species=True,
                           reaction_model="michaelis menten")       # Bad stoichiometry


    # E + S <-> ES -> P + E, with SingleSubstrateMechanism model
    sr = SpeciesRegistry(ids=["S", "P", "E"])
    rxn_defn = ReactionDefinition(id=123, reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                             reaction_model="single substrate mechanism",
                             kinetic_parameters={"k1_F": 10, "k1_R": 2, "k2_F": 3})
    assert rxn_defn.species_registry == sr
    assert rxn_defn.stoichiometry == Stoichiometry(vector={"S": -1, "P": 1}, catalysts=["E"])
    assert rxn_defn.reaction_model == "single substrate mechanism"
    assert rxn_defn.source_kinetic_parameters == {"k1_F": 10, "k1_R": 2, "k2_F": 3}

    sim_rxn_tuple = rxn_defn.sim_reactions
    assert len(sim_rxn_tuple) == 2
    sim_rxn_1, sim_rxn_2 = sim_rxn_tuple

    assert type(sim_rxn_1) == SimulationReaction
    assert type(sim_rxn_1.model) == MassAction_Model
    assert sim_rxn_1.source_object.id == 123
    assert sim_rxn_1.stoichiometry == Stoichiometry(vector={"S": -1, "E": -1, "ES*": 1})
    assert sim_rxn_1.model.get_parameters() == {'kF': 10, 'kR': 2, 'K': 5.0, 'reversible': True}

    assert type(sim_rxn_2) == SimulationReaction
    assert type(sim_rxn_2.model) == MassAction_Model
    assert sim_rxn_2.source_object.id == 123
    assert sim_rxn_2.stoichiometry == Stoichiometry(vector={"ES*": -1, "P": 1, "E": 1})
    assert sim_rxn_2.model.get_parameters() == {'kF': 3, 'kR': None, 'K': None, 'reversible': False}

    assert sr.number_of_species() == 4    # 1 species was automatically added
    set_of_species = set(sr.get_all_species_ids())
    assert set_of_species == {"S", "P", "E", "ES*"}



def test_CONSTRUCTOR_ReactionDefinition_3():
    # Thermodynamic data passed, but no temperature

    # Reaction R -> P
    sr = SpeciesRegistry()
    rxn_defn = ReactionDefinition(reactants="R", products="P", species_registry=sr, autoregister_species=True)

    assert set(rxn_defn.species_registry.get_all_species_ids()) == {"R", "P"}
    assert rxn_defn.stoichiometry == Stoichiometry(vector={'R': -1, 'P': 1})
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K_eq=None)
    assert rxn_defn.analytic_solution_family == "ONE_TO_ONE"
    assert rxn_defn.reaction_category == "Unimolecular rearrangement/isomerization"
    assert rxn_defn.source_kinetic_parameters == {}


    # A -> B, with "mass action" (reversible)
    sr = SpeciesRegistry(ids=["A", "B"])
    rxn_defn = ReactionDefinition(id=41, reactants="A", products="B", species_registry=sr,
                                  thermodynamic_parameters={"delta_H": -3, "K_eq": 5},
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF":10})
    assert rxn_defn.species_registry == sr
    assert rxn_defn.stoichiometry == Stoichiometry(vector={"A": -1, "B": 1})
    assert rxn_defn.reaction_model == "mass action"
    assert rxn_defn.analytic_solution_family == "ONE_TO_ONE"
    assert rxn_defn.reaction_category == "Unimolecular rearrangement/isomerization"
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=-3, delta_S=None, delta_G=None, K_eq=5, derived_pars=set())
    assert rxn_defn.source_kinetic_parameters == {"kF": 10}
    assert rxn_defn.source_thermodynamic_parameters == {"delta_H": -3, "K_eq": 5}

    sim_rxn_tuple = rxn_defn.sim_reactions
    assert len(sim_rxn_tuple) == 1
    sim_rxn = sim_rxn_tuple[0]
    assert type(sim_rxn) == SimulationReaction
    assert type(sim_rxn.model) == MassAction_Model
    assert sim_rxn.source_object.id == 41
    assert sim_rxn.derivation == "direct"
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"A": -1, "B": 1})
    assert sim_rxn.model.get_parameters() == {'kF': 10, 'kR': 2, 'K': 5.0, 'reversible': True}
    assert sim_rxn.model.derived_pars == {'K', 'kR', 'reversible'}


    # A -> B, with "mass action" (irreversible)
    rxn_defn = ReactionDefinition(id=42, reactants="A", products="B", species_registry=sr,
                                  thermodynamic_parameters={"delta_H": -3},
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF":10, "kR":0})

    assert rxn_defn.source_thermodynamic_parameters == {"delta_H": -3}
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=-3, delta_S=None, delta_G=None, K_eq=None, derived_pars=set())
    assert rxn_defn.source_kinetic_parameters == {"kF": 10, "kR":0}

    sim_rxn = rxn_defn.sim_reactions[0]
    assert type(sim_rxn) == SimulationReaction
    assert type(sim_rxn.model) == MassAction_Model
    assert sim_rxn.source_object.id == 42
    assert sim_rxn.derivation == "direct"
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"A": -1, "B": 1})
    assert sim_rxn.model.get_parameters() == {'kF': 10, 'kR': 0, 'K': math.inf, 'reversible': False}
    assert sim_rxn.model.derived_pars == {'K', 'reversible'}


    # Reaction R -> P + Q, with "mass action"
    rxn_defn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr, autoregister_species=True,
                                  thermodynamic_parameters={"delta_H": 5, "delta_S": -3})
    assert rxn_defn.stoichiometry == Stoichiometry(vector={'R': -1, 'P': 1, 'Q': 1})
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=5, delta_S=-3, delta_G=None, K_eq=None, derived_pars=set())
    assert rxn_defn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn_defn.reaction_category == "Unimolecular decomposition"


    # Enz + Sub -> Enz + Prod, with MM model
    sr = SpeciesRegistry(ids=["Sub", "Prod", "Enz"])

    rxn_defn = ReactionDefinition(id=43, reactants=["Sub", "Enz"], products=["Prod", "Enz"], species_registry=sr,
                                  thermodynamic_parameters={"delta_H": -20},
                                  reaction_model="michaelis menten",
                                  kinetic_parameters={'kM': 2, 'kcat': 3})
    assert rxn_defn.source_kinetic_parameters == {'kM': 2, 'kcat': 3}
    assert rxn_defn.source_thermodynamic_parameters == {"delta_H": -20}
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=-20, delta_S=None, delta_G=None, K_eq=None, derived_pars=set())
    assert rxn_defn.reaction_category == "Enzymatic"
    assert rxn_defn.analytic_solution_family is None

    sim_rxn_tuple = rxn_defn.sim_reactions
    assert len(sim_rxn_tuple) == 1
    sim_rxn = sim_rxn_tuple[0]
    assert type(sim_rxn) == SimulationReaction
    assert sim_rxn.source_object.id == 43
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"Sub": -1, "Prod": 1}, catalysts=["Enz"])
    assert type(sim_rxn.model) == MichaelisMenten_Model
    assert sim_rxn.model.get_parameters() == {'Enzyme': 'Enz', 'Product': 'Prod', 'Substrate': 'Sub', 'kM': 2, 'kcat': 3}
    assert sim_rxn.model.derived_pars == set()
    assert sim_rxn.model.E == "Enz"
    assert sim_rxn.model.S == "Sub"
    assert sim_rxn.model.P == "Prod"



    # E + S -> E + P, with MM model (but passing "k1_F", "k1_R" and "k2_F")
    sr = SpeciesRegistry(ids=["S", "P", "E"])

    rxn_defn = ReactionDefinition(id=44, reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                                  thermodynamic_parameters={"delta_S": 100},
                                  reaction_model="michaelis menten",
                                  kinetic_parameters={"k1_F": 10, "k1_R": 2, "k2_F": 5})
    assert rxn_defn.source_kinetic_parameters == {"k1_F": 10, "k1_R": 2, "k2_F": 5}
    assert rxn_defn.source_thermodynamic_parameters == {"delta_S": 100}
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=100, delta_G=None, K_eq=None, derived_pars=set())
    assert rxn_defn.reaction_category == "Enzymatic"
    assert rxn_defn.analytic_solution_family is None

    sim_rxn_tuple = rxn_defn.sim_reactions
    assert len(sim_rxn_tuple) == 1
    sim_rxn = sim_rxn_tuple[0]
    assert type(sim_rxn) == SimulationReaction
    assert sim_rxn.source_object.id == 44
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"S": -1, "P": 1}, catalysts=["E"])
    assert type(sim_rxn.model) == MichaelisMenten_Model
    assert sim_rxn.model.get_parameters() == {'Enzyme': 'E', 'Product': 'P', 'Substrate': 'S', 'kM': 0.7, 'kcat': 5} # (kM = k2_F + k1_R) / k1_F  ; kcat = k2_F)
    assert sim_rxn.model.derived_pars == {'kM', 'kcat'}


    with pytest.raises(Exception):
        ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                           reaction_model="michaelis menten",
                           kinetic_parameters={"k1_F": 10, "k1_R": 2, "k2_F": 5, "kM": 0.71})   # Inconsistent

    with pytest.raises(Exception):
        ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                           reaction_model="michaelis menten",
                           kinetic_parameters={"k1_F": 10, "k1_R": 2, "k2_F": 5, "kcat": 5.01})   # Inconsistent


    # E + S <-> ES -> P + E, with SingleSubstrateMechanism model
    sr = SpeciesRegistry(ids=["S", "P", "E"])

    rxn_defn = ReactionDefinition(id=43, reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                                  thermodynamic_parameters={"delta_S": 100},
                                  reaction_model="single substrate mechanism",
                                  kinetic_parameters={"k1_F": 10, "k1_R": 2, "k2_F": 3})
    assert rxn_defn.source_kinetic_parameters == {"k1_F": 10, "k1_R": 2, "k2_F": 3}
    assert rxn_defn.source_thermodynamic_parameters == {"delta_S": 100}
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=100, delta_G=None, K_eq=None, derived_pars=set())

    sim_rxn_tuple = rxn_defn.sim_reactions
    assert len(sim_rxn_tuple) == 2
    sim_rxn_1, sim_rxn_2 = sim_rxn_tuple

    assert type(sim_rxn_1) == SimulationReaction
    assert type(sim_rxn_1.model) == MassAction_Model
    assert sim_rxn_1.source_object.id == 43
    assert sim_rxn_1.stoichiometry == Stoichiometry(vector={"S": -1, "E": -1, "ES*": 1})
    assert sim_rxn_1.model.get_parameters() == {'kF': 10, 'kR': 2, 'K': 5.0, 'reversible': True}

    assert type(sim_rxn_2) == SimulationReaction
    assert type(sim_rxn_2.model) == MassAction_Model
    assert sim_rxn_2.source_object.id == 43
    assert sim_rxn_2.stoichiometry == Stoichiometry(vector={"ES*": -1, "P": 1, "E": 1})
    assert sim_rxn_2.model.get_parameters() == {'kF': 3, 'kR': None, 'K': None, 'reversible': False}

    assert sr.number_of_species() == 4    # 1 species was automatically added
    set_of_species = set(sr.get_all_species_ids())
    assert set_of_species == {"S", "P", "E", "ES*"}



def test_CONSTRUCTOR_ReactionDefinition_4():
    # Thermodynamic data passed, with temperature

    sr = SpeciesRegistry()

    # Reaction R -> P + Q
    rxn_defn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr, autoregister_species=True,
                                  thermodynamic_parameters={"delta_H": 5, "delta_S": -3, "temp": 100})
    assert rxn_defn.stoichiometry == Stoichiometry(vector={'R': -1, 'P': 1, 'Q': 1})
    assert rxn_defn.source_thermodynamic_parameters == {"delta_H": 5, "delta_S": -3, "temp": 100}
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=5, delta_S=-3, delta_G=5.3,
                                                             K_eq=0.0017045829244452543, temp=100,
                                                             derived_pars={"delta_G", "K_eq"})
    assert rxn_defn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn_defn.reaction_category == "Unimolecular decomposition"


    with pytest.raises(Exception):
        # Inconsistent thermodynamic/kinetic data
        ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr, autoregister_species=True,
                                  thermodynamic_parameters={"delta_H": 5, "delta_S": -3, "temp": 100},
                                  reaction_model="mass action", kinetic_parameters={"kF": 10, "kR": 2})


    # Reaction R -> P + Q, with "mass action"
    # The thermodynamic data allows derivation of kinetic parameters not supplied
    rxn_defn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr, autoregister_species=True,
                                  thermodynamic_parameters={"delta_H": 5, "delta_S": -3, "temp": 100},
                                  reaction_model="mass action", kinetic_parameters={"kF": 10})

    assert rxn_defn.stoichiometry == Stoichiometry(vector={'R': -1, 'P': 1, 'Q': 1})
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=5, delta_S=-3, delta_G=5.3,
                                                             K_eq=0.0017045829244452543, temp=100,
                                                             derived_pars={"delta_G", "K_eq"})
    assert rxn_defn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn_defn.reaction_category == "Unimolecular decomposition"
    assert rxn_defn.source_kinetic_parameters == {"kF": 10}

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.model.get_parameters() == {'kF': 10, 'kR': 5866.537706433048, 'K': 0.0017045829244452543, 'reversible': True}
    assert sim_rxn.model.derived_pars == {'kR', 'K', 'reversible'}


    # Redundant but consistent thermodynamic/kinetic data
    rxn_defn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr, autoregister_species=True,
                                  thermodynamic_parameters={"delta_H": 5, "delta_S": -3, "temp": 100},
                                  reaction_model="mass action", kinetic_parameters={"kF": 10, "kR": 5866.537706433048})

    assert rxn_defn.stoichiometry == Stoichiometry(vector={'R': -1, 'P': 1, 'Q': 1})
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=5, delta_S=-3, delta_G=5.3,
                                                             K_eq=0.0017045829244452543, temp=100,
                                                             derived_pars={"delta_G", "K_eq"})
    assert rxn_defn.analytic_solution_family == "ONE_TO_TWO"
    assert rxn_defn.reaction_category == "Unimolecular decomposition"
    assert rxn_defn.source_kinetic_parameters == {"kF": 10, "kR": 5866.537706433048}

    sim_rxn = rxn_defn.sim_reactions[0]
    assert sim_rxn.model.get_parameters() == {'kF': 10, 'kR': 5866.537706433048, 'K': 0.0017045829244452543, 'reversible': True}
    assert sim_rxn.model.derived_pars == {'K', 'reversible'}




    # Add a reaction with thermodynamic data;
    # the reverse reaction rate will get computed from the thermodynamic data
    rxn_defn = ReactionDefinition(reactants=["A"], products=[(2, "B")], species_registry=sr, autoregister_species=True,
                                  thermodynamic_parameters={"delta_H": 0.005, "delta_S": 0.4, "temp": 200},
                                  reaction_model="mass action", kinetic_parameters={"kF": 10})
    assert rxn_defn.stoichiometry.get_reactant_list() == [(1, "A")]
    assert rxn_defn.stoichiometry.get_product_list()  == [(2, "B")]

    assert rxn_defn.stoichiometry == Stoichiometry(vector={'A': -1, 'B': 2})

    assert math.isclose(rxn_defn.thermodynamics.delta_H, 0.005)
    assert math.isclose(rxn_defn.thermodynamics.delta_S, 0.4)
    assert math.isclose(rxn_defn.thermodynamics.delta_G, -0.075)        # In kJ/mol :  0.005 - 200 * 0.4/1000
    assert math.isclose(rxn_defn.thermodynamics.K_eq, 1.046134699475)   # exp(75/(8.31446261815324 * 200))
    assert rxn_defn.thermodynamics.derived_pars == {"delta_G", "K_eq"}

    sim_rxn = rxn_defn.sim_reactions[0]
    assert math.isclose(sim_rxn.model.K, 1.046134699475)        # From thermodynamics
    assert math.isclose(sim_rxn.model.kF, 10)
    assert math.isclose(sim_rxn.model.kR, 9.5589984779)         # 10. / 1.046134699475
    assert sim_rxn.model.reversible == True
    assert sim_rxn.model.derived_pars == {'K', 'kR', 'reversible'}



def test_extract_rxn_properties():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                             thermodynamic_parameters={"delta_H": -30},
                             reaction_model="mass action", kinetic_parameters={"kF":10, "kR":2})

    assert rxn.extract_rxn_properties() == {'reaction_model': 'mass action', 'delta_H': -30, 'K_eq': 5.0, 'K': 5, 'kF': 10, 'kR': 2, 'reversible': True}



def kinetic_test_func(stoichiometry, kinetic_parameters, conc_dict):
    # Used as place-holder custom kinetic rate function
    pass


def test_describe():
    sr = SpeciesRegistry(ids=["R", "P"])

    rxn_defn = ReactionDefinition(reactants="R", products="P", species_registry=sr,
                             thermodynamic_parameters={"delta_H": 0.5, "delta_S": -3, "temp": 100},
                             reaction_model="mass action", kinetic_parameters={"kF": 10})
    assert rxn_defn.describe(concise=True) == "R <-> P"

    assert rxn_defn.describe(concise=False) ==  \
        'R <-> P\n'  \
        '        Unimolecular rearrangement/isomerization reaction, with Reaction Model: "mass action"\n'  \
        '        Thermodynamics - passed:    (delta_H = 0.5 kJ/mol | delta_S = -3 J/(mol·K) | Temp = -173.1 C)\n'  \
        '        Thermodynamics - derived:   (delta_H = 0.5 kJ/mol | delta_S = -3 J/(mol·K) | delta_G = 0.8 kJ/mol | K_eq = 0.38206 | Temp = -173.1 C)\n'  \
        '        Kinetics - passed:   (kF = 10)\n'  \
        '        Kinetics - derived:  1 derived reaction\n'  \
        '            (1) Type: "mass action"   (kF = 10 | kR = 26.174 | K = 0.38206 | reversible = True)'


    # Reaction: # E + S <-> ES -> E + P, with SingleSubstrateMechanism model
    sr = SpeciesRegistry(ids=["S", "P", "E"])
    rxn_defn = ReactionDefinition(id=85, reactants=["S", "E"], products=["P", "E"],
                                  species_registry=sr, autoregister_species=True,
                                  reaction_model="single substrate mechanism",
                                  kinetic_parameters={"k1_F": 18, "k1_R": 100, "k2_F": 49})
    assert rxn_defn.describe(concise=True) == "S + E -> P + E"
    assert rxn_defn.describe(concise=False) ==  \
        'S + E -> P + E\n'  \
        '        Enzymatic reaction, with Reaction Model: "single substrate mechanism"   (Reaction ID 85)\n'  \
        '        Thermodynamics - passed:  None\n'  \
        '        Thermodynamics - derived: None\n'  \
        '        Kinetics - passed:   (k1_F = 18 | k1_R = 100 | k2_F = 49)\n'  \
        '        Kinetics - derived:  2 derived reactions\n'  \
        '            (1) Type: "mass action"   (kF = 18 | kR = 100 | K = 0.18 | reversible = True)\n'  \
        '            (2) Type: "mass action"   (kF = 49 | reversible = False)'


    # Reaction: # E + S <-> ES -> E + P, with "michaelis menten" model
    sr = SpeciesRegistry(ids=["S", "P", "E"])
    rxn_defn = ReactionDefinition(id=22, reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                                  reaction_model="michaelis menten",
                                  kinetic_parameters={'kM': 2, 'kcat': 5})
    assert rxn_defn.describe(concise=True) == "S + E -> P + E"
    assert rxn_defn.describe(concise=False) ==  \
        'S + E -> P + E\n'  \
        '        Enzymatic reaction, with Reaction Model: "michaelis menten"   (Reaction ID 22)\n'  \
        '        Thermodynamics - passed:  None\n'  \
        '        Thermodynamics - derived: None\n'  \
        '        Kinetics - passed:   (kM = 2 | kcat = 5)\n'  \
        '        Kinetics - derived:  1 derived reaction\n'  \
        '            (1) Type: "Michaelis-Menten"   (kM = 2 | kcat = 5 | Substrate = \'S\' | Enzyme = \'E\' | Product = \'P\')'


    # Reaction: CH4 + 2 O2 <-> CO2 + 2 H2O (no model specified)
    sr = SpeciesRegistry(ids=["CH4", "O2", "CO2", "H2O"])
    rxn_defn = ReactionDefinition(reactants=["CH4", (2, "O2")], products=["CO2", (2, "H2O")], species_registry=sr)
    assert rxn_defn.describe(concise=True) == "CH4 + 2 O2 -> CO2 + 2 H2O"

    assert rxn_defn.describe(concise=False) ==  \
        'CH4 + 2 O2 -> CO2 + 2 H2O\n'  \
        '        General one-step reaction, with no Reaction Model specified\n'  \
        '        Thermodynamics - passed:  None\n'  \
        '        Thermodynamics - derived: None\n'  \
        '        Kinetics - passed: None\n'  \
        '        Kinetics - derived:  None'


    # Same reaction, with more parameters
    rxn_defn = ReactionDefinition(id=777,
                                  reactants=["CH4", (2, "O2")], products=["CO2", (2, "H2O")], species_registry=sr,
                                  reaction_model="custom")
    assert rxn_defn.describe(concise=True) == "CH4 + 2 O2 -> CO2 + 2 H2O"
    assert rxn_defn.describe(concise=False) ==  \
        'CH4 + 2 O2 -> CO2 + 2 H2O\n'  \
        '        General one-step reaction, with Reaction Model: "custom"   (Reaction ID 777)\n'  \
        '        Thermodynamics - passed:  None\n'  \
        '        Thermodynamics - derived: None\n'  \
        '        Kinetics - passed: None\n'  \
        '        Kinetics - derived:  1 derived reaction\n'  \
        '            (1) Type: "custom"   (reversible = False)'


    # Same reaction, with yet more parameters
    rxn_defn = ReactionDefinition(id=999,
                                  reactants=["CH4", (2, "O2")], products=["CO2", (2, "H2O")], species_registry=sr,
                                  reaction_model="custom", kinetic_parameters={"kF": 10, "kR": 2, "rate_function": kinetic_test_func})
    assert rxn_defn.describe(concise=True) == "CH4 + 2 O2 <-> CO2 + 2 H2O"
    assert rxn_defn.describe(concise=False) ==  \
        'CH4 + 2 O2 <-> CO2 + 2 H2O\n'  \
        '        General one-step reaction, with Reaction Model: "custom"   (Reaction ID 999)\n'  \
        '        Thermodynamics - passed:  None\n'  \
        '        Thermodynamics - derived:   (K_eq = 5)\n'  \
        '        Kinetics - passed:   (kF = 10 | kR = 2 | rate_function = kinetic_test_func())\n'  \
        '        Kinetics - derived:  1 derived reaction\n'  \
        '            (1) Type: "custom"   (kF = 10 | kR = 2 | K = 5 | rate_function = kinetic_test_func() | reversible = True)'

    # Same reaction, with yet more parameters
    rxn_defn = ReactionDefinition(id=1001,
                                  reactants=["CH4", (2, "O2")], products=["CO2", (2, "H2O")], species_registry=sr,
                                  reaction_model="custom", kinetic_parameters={"kF": 10, "kR": 2, "rate_function": kinetic_test_func},
                                  thermodynamic_parameters={"temp": 200})
    assert rxn_defn.describe(concise=True) == "CH4 + 2 O2 <-> CO2 + 2 H2O"
    assert rxn_defn.describe(concise=False) ==  \
        'CH4 + 2 O2 <-> CO2 + 2 H2O\n'  \
        '        General one-step reaction, with Reaction Model: "custom"   (Reaction ID 1001)\n'  \
        '        Thermodynamics - passed:    (Temp = -73.15 C)\n'  \
        '        Thermodynamics - derived:   (delta_G = -2.6763 kJ/mol | K_eq = 5 | Temp = -73.15 C)\n'  \
        '        Kinetics - passed:   (kF = 10 | kR = 2 | rate_function = kinetic_test_func())\n'  \
        '        Kinetics - derived:  1 derived reaction\n'  \
        '            (1) Type: "custom"   (kF = 10 | kR = 2 | K = 5 | rate_function = kinetic_test_func() | reversible = True)'



def test_set_thermodynamic_data():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                             reaction_model="mass action", kinetic_parameters={"kF":6, "kR":2})
    assert rxn.thermodynamics.delta_H is None
    assert rxn.thermodynamics.delta_S is None
    assert rxn.thermodynamics.delta_G is None
    assert rxn.thermodynamics.K_eq == 3

    rxn.set_thermodynamic_data(temp=100)
    assert math.isclose(rxn.thermodynamics.delta_G, -0.9134370805974775)


    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                             thermodynamic_parameters={"delta_H": 0.5, "delta_S":-3})
    assert np.allclose(rxn.thermodynamics.delta_H, 0.5)
    assert rxn.thermodynamics.delta_S == -3
    assert rxn.thermodynamics.delta_G is None

    rxn.set_thermodynamic_data(temp=100)
    assert np.allclose(rxn.thermodynamics.delta_H, 0.5)
    assert rxn.thermodynamics.delta_S == -3
    assert np.allclose(rxn.thermodynamics.delta_G, 0.8)
    assert np.allclose(rxn.thermodynamics.K_eq, 0.38205953171)



def test_extract_intermediate():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr)
    assert rxn_defn.extract_intermediate() is None

    rxn_defn = ReactionDefinition(reactants=["R", "S"], products="P", species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action")
    assert rxn_defn.extract_intermediate() is None


    rxn_defn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True,
                                  reaction_model="michaelis menten")
    assert rxn_defn.extract_intermediate() is None


    rxn_defn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True,
                                  reaction_model="single substrate mechanism")
    assert rxn_defn.extract_intermediate() == "ES*"



def test_extract_reactant_ids():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr)
    assert rxn.extract_reactant_ids() == {"A"}

    rxn = ReactionDefinition(reactants=["A", "B"], products="C", species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == {"A", "B"}

    rxn = ReactionDefinition(reactants=["R", "R"], products="C", species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == {"R"}

    rxn = ReactionDefinition(reactants=[(2, "R")], products="C", species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == {"R"}

    rxn = ReactionDefinition(reactants="A", products=["B", "C"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == {"A"}

    rxn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == {"S", "E"}

    rxn = ReactionDefinition(reactants=["A", "B"], products=["B", "C"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == {"A", "B"}

    rxn = ReactionDefinition(reactants=[(2, "D"), "C"], products=["C", (3, "E")], species_registry=sr, autoregister_species=True)
    assert rxn.extract_reactant_ids() == {"D", "C"}


def test_extract_reactants_formula():
    sr = SpeciesRegistry(ids=["A", "B"])
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr)
    assert rxn_defn.extract_reactants_formula() == "A"

    rxn_defn = ReactionDefinition(reactants=["CH4", (2, "O2")],
                                 products=["CO2", (2, "H2O")], species_registry=sr, autoregister_species=True)
    assert rxn_defn.extract_reactants_formula() == "CH4 + 2 O2"


def test_extract_product_ids():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn.extract_product_ids() == {"B"}

    rxn = ReactionDefinition(reactants=["A", "B"], products="C", species_registry=sr, autoregister_species=True)
    assert rxn.extract_product_ids() == {"C"}

    rxn = ReactionDefinition(reactants="A", products=["B", "C"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_product_ids() == {"B", "C"}

    rxn = ReactionDefinition(reactants="A", products=["F", "F"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_product_ids() == {"F"}

    rxn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_product_ids() == {"P", "E"}

    rxn = ReactionDefinition(reactants=["A", "B"], products=["B", "C"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_product_ids() == {"B", "C"}

    rxn = ReactionDefinition(reactants=[(2, "D"), "C"], products=["C", (3, "E")], species_registry=sr, autoregister_species=True)
    assert rxn.extract_product_ids() == {"C", "E"}



def test_extract_products_formula():
    sr = SpeciesRegistry(ids=["A", "B"])
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn_defn.extract_products_formula() == "B"

    rxn_defn = ReactionDefinition(reactants=["CH4", (2, "O2")],
                                  products=["CO2", (2, "H2O")], species_registry=sr, autoregister_species=True)
    assert rxn_defn.extract_products_formula() == "CO2 + 2 H2O"



def test_extract_species_in_reaction():
    sr = SpeciesRegistry(ids=["A", "B", "C"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn.extract_species_in_reaction() == {"A", "B"}

    rxn = ReactionDefinition(reactants=["A", "B"], products="C", species_registry=sr, reaction_model="mass action")
    assert rxn.extract_species_in_reaction() == {"A", "B", "C"}

    rxn = ReactionDefinition(reactants="A", products=["B", "C"], species_registry=sr, reaction_model="mass action")
    assert rxn.extract_species_in_reaction() == {"A", "B", "C"}

    rxn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_species_in_reaction(include_intermediaries=False) == {"S", "P", "E"}
    assert rxn.extract_species_in_reaction(include_intermediaries=True) == {"S", "P", "E"}

    rxn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True,
                             reaction_model="michaelis menten")
    assert rxn.extract_species_in_reaction(include_intermediaries=False) == {"S", "P", "E"}
    assert rxn.extract_species_in_reaction(include_intermediaries=True) == {"S", "P", "E"}

    rxn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True,
                             reaction_model="single substrate mechanism")
    assert rxn.extract_species_in_reaction(include_intermediaries=False) == {"S", "P", "E"}
    assert rxn.extract_species_in_reaction(include_intermediaries=True) == {"S", "P", "E", "ES*"}

    rxn = ReactionDefinition(reactants=[(2, "D"), "C"], products=["C", (3, "E")], species_registry=sr, autoregister_species=True)
    assert rxn.extract_species_in_reaction() == {"C", "D", "E"}



def test_reaction_quotient():
    sr = SpeciesRegistry(ids=["A", "B"])

    # Reaction : A <-> B
    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    c = {'A': 24., 'B': 36.}
    assert np.allclose(1.5, rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1.5, quotient)
    assert formula == '[B] / [A]'

    # Reaction : A <-> F
    sr.add_species("F")
    rxn = ReactionDefinition(reactants="A", products="F", species_registry=sr, reaction_model="mass action")
    c = {'A': 3., 'F': 33.}
    assert np.allclose(11., rxn.reaction_quotient(conc=c, explain=False))
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(11., quotient)
    assert formula == '[F] / [A]'

    # Reaction :  A + B <-> C
    sr.add_species("C")
    rxn = ReactionDefinition(reactants=["A" , "B"], products="C", species_registry=sr, reaction_model="mass action")
    c = {'A': 3., 'B': 4., 'C': 12.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1., quotient)
    assert formula == '[C] / ([A][B])'

    # Reaction :  2A <-> P
    rxn = ReactionDefinition(reactants=["A" , "A"], products="P", species_registry=sr,
                             reaction_model="mass action", autoregister_species=True)
    c = {'A': 2., 'P': 20.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(5., quotient)
    assert formula == '[P] /  [A]^2 '

    # Reaction :  C <-> A + B
    rxn = ReactionDefinition(reactants="C", products=["A" , "B"], species_registry=sr,
                             reaction_model="mass action", autoregister_species=True)
    c = {'A': 3., 'B': 4., 'C': 12.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1., quotient)
    assert formula == '([A][B]) / [C]'

    # Reaction :  B <-> 2A
    rxn = ReactionDefinition(reactants="B", products=["A" , "A"], species_registry=sr,
                             reaction_model="mass action", autoregister_species=True)
    c = {'A': 2., 'B': 20.}
    quotient, formula = rxn.reaction_quotient(conc=c, explain=True)
    assert np.allclose(1/5., quotient)
    assert formula == ' [A]^2  / [B]'

    rxn = ReactionDefinition(reactants=[(2, "B")], products=[(3, "A")], species_registry=sr,
                             reaction_model="custom", autoregister_species=True)
    with pytest.raises(Exception):
        rxn.reaction_quotient(conc=c)       # Not "mass action"



def test_find_equilibrium_conc():
    sr = SpeciesRegistry()

    # Reaction : A <-> C
    rxn = ReactionDefinition(reactants="A", products="C", species_registry=sr, autoregister_species=True,
                             reaction_model="mass action", kinetic_parameters={"kF": 3, "kR": 2})
    assert rxn.analytic_solution_family == "ONE_TO_ONE"

    result = rxn.find_equilibrium_conc(conc_dict={"A":80., "C":10.})
    assert np.allclose(result["A"], 36)
    assert np.allclose(result["C"], 54)
    rxn.stoichiometry.consistency_checker(conc_before={"A":80., "C":10.}, conc_after={"A":36., "C":54.})


    # Only the forward reaction
    rxn = ReactionDefinition(reactants="A", products="C", species_registry=sr, autoregister_species=True,
                             reaction_model="mass action", kinetic_parameters={"kF": 3})
    assert rxn.analytic_solution_family == "ONE_TO_ONE"
    #assert rxn.kinetics.parameters["reversible"] == False
    result = rxn.find_equilibrium_conc(conc_dict={"A":80., "C":10.})
    assert np.allclose(result["A"], 0)
    assert np.allclose(result["C"], 90)

    # Only the reverse reaction
    rxn = ReactionDefinition(reactants="A", products="C", species_registry=sr, autoregister_species=True,
                             reaction_model="mass action", kinetic_parameters={"kR": 2})
    assert rxn.analytic_solution_family == "ONE_TO_ONE"
    result = rxn.find_equilibrium_conc(conc_dict={"A":80., "C":10.})
    assert np.allclose(result["A"], 90)
    assert np.allclose(result["C"], 0)

    # Reaction X + Y <-> Z
    rxn = ReactionDefinition(reactants=["X", "Y"], products="Z", species_registry=sr, autoregister_species=True,
                             reaction_model="mass action", kinetic_parameters={"kF": 5, "kR": 2})
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
    rxn = ReactionDefinition(reactants=["A", "A"], products="C", species_registry=sr, autoregister_species=True,
                             reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})
    assert rxn.analytic_solution_family == "TWO_TO_ONE"
    result = rxn.find_equilibrium_conc(conc_dict={"A":200., "C": 40.})
    expected_eq = [9.49568869375716, 135.2521556531214]
    assert np.allclose(result["A"], expected_eq[0])
    assert np.allclose(result["C"], expected_eq[1])
    rxn.stoichiometry.consistency_checker(conc_before={"A":200., "C": 40.},
                                          conc_after={"A":expected_eq[0], "C": expected_eq[1]})


    # Reaction Z <-> X + Y
    rxn = ReactionDefinition(reactants="Z", products=["X", "Y"], species_registry=sr, autoregister_species=True,
                             reaction_model="mass action", kinetic_parameters={"kF": 2., "kR": 5.})
    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    result = rxn.find_equilibrium_conc(conc_dict={"X":10., "Y": 50, "Z":20.})
    expected_eq = [0.2948774087575341, 40.294877408757536, 29.705122591242464]
    assert np.allclose(result["X"], expected_eq[0])
    assert np.allclose(result["Y"], expected_eq[1])
    assert np.allclose(result["Z"], expected_eq[2])
    rxn.stoichiometry.consistency_checker(conc_before={"X":10., "Y": 50, "Z": 20.},
                                          conc_after={"X":expected_eq[0], "Y": expected_eq[1], "Z": expected_eq[2]})


    # C <-> 2 A
    rxn = ReactionDefinition(reactants="C", products=["A", "A"], species_registry=sr, autoregister_species=True,
                             reaction_model="mass action", kinetic_parameters={"kF": 2., "kR": 3.})
    assert rxn.analytic_solution_family == "ONE_TO_TWO"
    result = rxn.find_equilibrium_conc(conc_dict={"C": 40., "A":200.})

    assert np.allclose(result["C"], 135.2521556531214)
    assert np.allclose(result["A"], 9.49568869375716)
    rxn.stoichiometry.consistency_checker(conc_before={"C": 40., "A":200.},
                                          conc_after={"C": 135.2521556531214, "A": 9.49568869375716})


    rxn = ReactionDefinition(reactants=[(2, "B")], products=[(3, "A")], species_registry=sr,
                             reaction_model="custom", autoregister_species=True)
    with pytest.raises(Exception):
        rxn.find_equilibrium_conc(conc_dict={"A":10., "B": 50})     # Not "mass action"




##########  For PRIVATE methods  ##########

def test__standard_form_chem_eqn():
    sr = SpeciesRegistry()
    rxn = ReactionDefinition(reactants="A", products="B",
                             species_registry=sr, autoregister_species=True)     # Won't actually use reactants/products

    assert rxn._standard_form_chem_eqn([(1, "Fe"), (2, "Cl")]) == "Fe + 2 Cl"

    assert rxn._standard_form_chem_eqn([(3, "Fe"), (5, "G")]) == "3 Fe + 5 G"
