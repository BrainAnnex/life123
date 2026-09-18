import numpy as np
import pytest
import math
from life123.species_registry import Species, SpeciesRegistry, MacroMolecules
from life123.reactions_new import Stoichiometry, ReactionThermodynamics, \
ReactionDefinition, SimulationReaction, MassAction_Model, MichaelisMenten_Model
from tests.utilities.comparisons import *





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

    # Reaction S + E -> P + E
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

    rxn = ReactionDefinition(reactants=("R", "S"), products="P", species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"R": -1, "S": -1, "P": 1}

    rxn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"R": -1, "P": 1, "Q": 1}

    rxn = ReactionDefinition(reactants=("E", "S"), products=("E", "P"),
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"S": -1, "P": 1}

    rxn = ReactionDefinition(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"A": -2, "B":- 2, "P": 3, "Q": 1}


    # Reaction  2A + 3B + 2 E + F --> 4C  + 5D  + 2 E + F
    rxn = ReactionDefinition(reactants=["E", "A", "F", (3, "B"), "E", "A"], products=[(4, "C"), "F", (2, "E"), (5, "D")],
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_vector() == {"A": -2, "B": -3, "C": 4, "D": 5}



def test_get_reaction_complexes():
    st = Stoichiometry(vector={"A": -1, "P": 2, "Q": 1}, catalysts=["E"])

    assert st.get_reaction_complexes() == ( {"A": 1, "E": 1} ,  {"P": 2, "Q": 1, "E": 1} )


    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    rxn = ReactionDefinition(reactants="R", products="P", species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"R": 1}, {"P": 1})

    rxn = ReactionDefinition(reactants=("R", "S"), products="P", species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"R": 1, "S": 1}, {"P": 1})

    rxn = ReactionDefinition(reactants="R", products=["P", "Q"], species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"R": 1}, {"P": 1, "Q": 1})

    rxn = ReactionDefinition(reactants=("E", "S"), products=("E", "P"),
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"S": 1, "E": 1}, {"P": 1, "E": 1})

    rxn = ReactionDefinition(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                             species_registry=sr)
    assert rxn.stoichiometry.get_reaction_complexes() == ({"A": 2, "B": 2, "E": 1}, {"P": 3, "Q": 1, "E": 1})



def test_get_reactant_list():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_reactant_list() == [(2, "A"), (2, "B")]

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_reactant_list() == [(2, "A"), (2, "B"), (1, "E")]

def test_get_reactant_ids():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_reactant_ids() == {"A", "B"}

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_reactant_ids() == {"A", "B", "E"}


def test_get_product_list():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_product_list() == [(3, "P")]

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_product_list() == [(3, "P"), (1, "E")]


def test_get_product_ids():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_product_ids() == {"P"}

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_product_ids() == {"P", "E"}


def test_get_all_species_ids():
    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3})
    assert s.get_all_species_ids() == {"A", "B", "P"}

    s = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert s.get_all_species_ids() == {"A", "B", "P", "E"}



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

    with pytest.raises(Exception):
        ReactionDefinition(reactants="R", products=123, species_registry=sr)  # Bad product

    with pytest.raises(Exception):
        ReactionDefinition(reactants="R", products="R", species_registry=sr, reaction_model="mass action")     # Cannot be same

    with pytest.raises(Exception):
        ReactionDefinition(reactants=("R", (2, "P")), products=[(2, "P"), "R"], species_registry=sr , reaction_model="mass action")     # Cannot be same


    sr = SpeciesRegistry(ids=["A", "B", "R", "P", "Q", "S", "E"])

    # Reaction R -> P
    rxn_defn = ReactionDefinition(reactants="R", products="P", species_registry=sr)
    assert rxn_defn.species_registry == sr
    assert rxn_defn.stoichiometry.to_dict() == {"R": -1, "P": 1}
    assert rxn_defn.reaction_model is None
    assert rxn_defn.sim_reactions == ()
    assert rxn_defn.thermodynamics == ReactionThermodynamics(delta_H=None, delta_S=None, delta_G=None, K_eq=None, derived_pars=set())

    # Reaction R + S -> P
    rxn_defn = ReactionDefinition(reactants=("R", "S"), products="P", species_registry=sr)
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

    rxn_defn = ReactionDefinition(reactants=("E", "S"), products=("E", "P"),
                             species_registry=sr)
    assert rxn_defn.stoichiometry.to_dict() == {"S": -1, "P": 1, "E": 0}

    rxn_defn = ReactionDefinition(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), "Q", "E"],
                             species_registry=sr)
    assert rxn_defn.stoichiometry.to_dict() == {"A": -2, "B":- 2, "P": 3, "Q": 1, "E": 0}


    # Reaction S + E -> P + E
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
    assert sim_rxn.source_definition_id == 942
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
    assert sim_rxn.source_definition_id == 8
    assert sim_rxn.derivation == "direct"
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"A": -1, "B": 1})
    assert sim_rxn.model.get_parameters() == {'kF': 10, 'kR': 2, 'K': 5.0, 'reversible': True}


    with pytest.raises(Exception):
        ReactionDefinition(reactants="A", products="B", species_registry=sr,
                           reaction_model="mass action", kinetic_parameters={"intruder":666})



    # S + E -> P + E , with "michaelis menten" model
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
    assert sim_rxn.source_definition_id == 17
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"S": -1, "P": 1}, catalysts=["E"])
    assert type(sim_rxn.model) == MichaelisMenten_Model
    assert sim_rxn.model.get_parameters() == {'kM': 2, 'kcat': 5}



    # S + E <-> SE -> P + E, with SingleSubstrateMechanism model
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
    assert sim_rxn_1.source_definition_id == 123
    assert sim_rxn_1.stoichiometry == Stoichiometry(vector={"S": -1, "E": -1, "SE*": 1})
    assert sim_rxn_1.model.get_parameters() == {'kF': 10, 'kR': 2, 'K': 5.0, 'reversible': True}

    assert type(sim_rxn_2) == SimulationReaction
    assert type(sim_rxn_2.model) == MassAction_Model
    assert sim_rxn_2.source_definition_id == 123
    assert sim_rxn_2.stoichiometry == Stoichiometry(vector={"SE*": -1, "P": 1, "E": 1})
    assert sim_rxn_2.model.get_parameters() == {'kF': 3, 'kR': None, 'K': None, 'reversible': False}

    assert sr.number_of_species() == 4    # 1 species was automatically added
    set_of_species = set(sr.get_all_species_ids())
    assert set_of_species == {"S", "P", "E", "SE*"}



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
    assert sim_rxn.source_definition_id == 41
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
    assert sim_rxn.source_definition_id == 42
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


    # S + E -> P + E , with MM model
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
    assert sim_rxn.source_definition_id == 44
    assert sim_rxn.stoichiometry == Stoichiometry(vector={"S": -1, "P": 1}, catalysts=["E"])
    assert type(sim_rxn.model) == MichaelisMenten_Model
    assert sim_rxn.model.get_parameters() == {'kM': 0.7, 'kcat': 5} # (kM = k2_F + k1_R) / k1_F  ; kcat = k2_F)
    assert sim_rxn.model.derived_pars == {'kM', 'kcat'}


    with pytest.raises(Exception):
        ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                           reaction_model="michaelis menten",
                           kinetic_parameters={"k1_F": 10, "k1_R": 2, "k2_F": 5, "kM": 0.71})   # Inconsistent

    with pytest.raises(Exception):
        ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                           reaction_model="michaelis menten",
                           kinetic_parameters={"k1_F": 10, "k1_R": 2, "k2_F": 5, "kcat": 5.01})   # Inconsistent


    # S + E <-> SE -> P + E, with SingleSubstrateMechanism model
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
    assert sim_rxn_1.source_definition_id == 43
    assert sim_rxn_1.stoichiometry == Stoichiometry(vector={"S": -1, "E": -1, "SE*": 1})
    assert sim_rxn_1.model.get_parameters() == {'kF': 10, 'kR': 2, 'K': 5.0, 'reversible': True}

    assert type(sim_rxn_2) == SimulationReaction
    assert type(sim_rxn_2.model) == MassAction_Model
    assert sim_rxn_2.source_definition_id == 43
    assert sim_rxn_2.stoichiometry == Stoichiometry(vector={"SE*": -1, "P": 1, "E": 1})
    assert sim_rxn_2.model.get_parameters() == {'kF': 3, 'kR': None, 'K': None, 'reversible': False}

    assert sr.number_of_species() == 4    # 1 species was automatically added
    set_of_species = set(sr.get_all_species_ids())
    assert set_of_species == {"S", "P", "E", "SE*"}



def test_constructor_ReactionDefinition_4():
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



"""
def test_extract_rxn_properties():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                             thermodynamic_parameters={"delta_H": -30},
                             kinetic_parameters={"kF":10, "kR":2})

    assert rxn.extract_rxn_properties() == {'kinetics_type': 'mass action', 'kF': 10, 'kR': 2, 'delta_H': -30, 'K': 5.0, 'reversible': True}


def test_describe():
    sr = SpeciesRegistry(ids=["R", "P"])

    rxn = ReactionDefinition(reactants="R", products="P", species_registry=sr,
                             delta_H=0.5, delta_S=-3, temp=100,
                             reaction_model="mass action", kinetic_parameters={"kF": 10})
    assert rxn.describe(concise=True) == "R <-> P"
    print(rxn.describe(concise=False))
    assert rxn.describe(concise=False) == \
            "R <-> P  Elementary Unimolecular rearrangement/isomerization reaction\n" \
            "         (delta_H = 0.5 kJ/mol | delta_S = -3 J/(mol·K) | delta_G = 0.8 kJ/mol | K = 0.38206 | Temp = -173.1 C | kinetics_type = 'mass action' | kF = 10 | kR = 26.174 | reversible = True)"

"""

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

    rxn_defn = ReactionDefinition(reactants=("R", "S"), products="P", species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action")
    assert rxn_defn.extract_intermediate() is None


    rxn_defn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True,
                                  reaction_model="michaelis menten")
    assert rxn_defn.extract_intermediate() is None


    rxn_defn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True,
                                  reaction_model="single substrate mechanism")
    assert rxn_defn.extract_intermediate() == "SE*"



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



def test_extract_reactants():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn.extract_reactants() == [(1, "A")]


def test_extract_reactants_formula():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn.extract_reactants_formula() == "A"



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



def test_extract_products():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn.extract_products() == [(1, "B")]


def test_extract_products_formula():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn.extract_products_formula() == "B"



def test_extract_species_in_reaction():
    sr = SpeciesRegistry(ids=["A", "B", "C"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr, reaction_model="mass action")
    assert rxn.extract_species_in_reaction() == {"A", "B"}

    rxn = ReactionDefinition(reactants=["A", "B"], products="C", species_registry=sr, reaction_model="mass action")
    assert rxn.extract_species_in_reaction() == {"A", "B", "C"}

    rxn = ReactionDefinition(reactants="A", products=["B", "C"], species_registry=sr, reaction_model="mass action")
    assert rxn.extract_species_in_reaction() == {"A", "B", "C"}

    rxn = ReactionDefinition(reactants=["S", "E"], products=["P", "E"], species_registry=sr, autoregister_species=True)
    assert rxn.extract_species_in_reaction() == {"S", "P", "E"}



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


"""

def test_determine_reaction_rate():
    sr = SpeciesRegistry(ids=["A", "B"])

    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                             kinetic_parameters={"kF": 20., "kR": 2.})
    assert rxn.kinetics.law == "mass action"

    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5. - 2. * 8.)  # 84.0

    # Now just the forward reaction
    rxn.kinetics.set_parameters({"kR": 0})
    result = rxn.determine_reaction_rate(conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5.)            # 100.0


    # Reaction A + B -> C
    rxn = ReactionDefinition(reactants=["A", "B"], products="C",
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
    rxn = ReactionDefinition(reactants="A", products=["B", "C"], species_registry=sr, autoregister_species=True,
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
    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
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
    rxn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
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
    rxn = ReactionDefinition(reactants=["A" , "B"], products="C", species_registry=sr, autoregister_species=True,
                             kinetic_parameters={"kF": 5., "kR": 2})

    result = rxn.step_simulation(delta_time=0.002, conc_dict={"A": 10, "B": 50, "C": 20})
    assert result[0] == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert result[1] == 5*10*50 - 2 * 20        # 2460


    # Reaction : C <-> A + B
    rxn = ReactionDefinition(reactants="C", products=["A" , "B"], species_registry=sr, autoregister_species=True,
                             kinetic_parameters={"kF": 2., "kR": 5.})

    result = rxn.step_simulation(delta_time=0.002, conc_dict={"A": 10, "B": 50, "C": 20})
    assert result[0] == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert result[1] == -2460
"""


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
