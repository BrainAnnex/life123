import numpy as np
import pytest
from life123.species_registry import SpeciesRegistry
from life123.reaction_kinetics import ReactionKinetics
from life123.reactions import Stoichiometry, ReactionThermodynamics, \
ReactionDefinition, SimulationReaction, MassAction_Model, MichaelisMenten_Model
from tests.utilities.comparisons import *



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

    st = Stoichiometry(vector={"A": -2, "B":- 2, "P": 3}, catalysts=["E"])
    assert st.reaction_pattern() == (2, 1, 1)



def test_is_elementary_like():
    assert Stoichiometry(vector={"R": -1, "P": 1}).is_elementary_like()

    assert Stoichiometry(vector={"R": -1, "P": 1, "Q": 1}).is_elementary_like()
    assert Stoichiometry(vector={"R": -1, "P": 2}).is_elementary_like()
    assert not Stoichiometry(vector={"R": -1, "P": 3}).is_elementary_like()

    assert Stoichiometry(vector={"R": -1, "S": -1, "P": 1}).is_elementary_like()
    assert Stoichiometry(vector={"R": -2, "P": 1}).is_elementary_like()
    assert not Stoichiometry(vector={"R": -3, "P": 1}).is_elementary_like()

    assert not Stoichiometry(vector={"R": -2, "P": 2}).is_elementary_like()
    assert not Stoichiometry(vector={"R": -1, "P": 1, "Q": 1, "S": 1}).is_elementary_like()

    assert not Stoichiometry(vector={"S": -1, "P": 1}, catalysts=["E"]).is_elementary_like()



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



def test_from_reactants_products():

    st = Stoichiometry.from_reactants_products(reactants="A", products="B")
    assert st == Stoichiometry(vector={'A': -1, 'B': 1}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants=["A"], products="B")
    assert st == Stoichiometry(vector={'A': -1, 'B': 1}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants="A", products=["B"])
    assert st == Stoichiometry(vector={'A': -1, 'B': 1}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants=["A"], products=["B"])
    assert st == Stoichiometry(vector={'A': -1, 'B': 1}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants=["A", "B"], products=["P"])
    assert st == Stoichiometry(vector={'A': -1, 'B': -1, 'P': 1}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants="R", products=["P", "Q"])
    assert st == Stoichiometry(vector={'R': -1, 'P': 1, 'Q': 1}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants=["R"], products=[(2,"P")])
    assert st == Stoichiometry(vector={'R': -1, 'P': 2}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants=[(2,"A")], products="B")
    assert st == Stoichiometry(vector={'A': -2, 'B': 1}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants=[(2,"A")], products=[(1, "P")])
    assert st == Stoichiometry(vector={'A': -2, 'P': 1}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants=[(2,"A"), "X"], products=[(1, "P")])
    assert st == Stoichiometry(vector={'A': -2, 'X': -1, 'P': 1}, catalysts=[])

    st = Stoichiometry.from_reactants_products(reactants=[(2,"A"), "X"], products=["X", (1, "P")])
    assert st == Stoichiometry(vector={'A': -2, 'P': 1}, catalysts=['X'])

    st = Stoichiometry.from_reactants_products(reactants=["A", (2, "B"), "E", "A"], products=[(3, "P"), (1, "Q"), "E"])  # `A` gets combined
    assert st == Stoichiometry(vector={"A": -2, "B":- 2, "P": 3, "Q": 1}, catalysts=['E'])

    st = Stoichiometry.from_reactants_products(reactants=["S", "E"], products=["P", "E"])
    assert st == Stoichiometry(vector={'S': -1, 'P': 1}, catalysts=["E"])

    # Edge case: 2E + S -> E + P   (E gets consumed, and thus not regarded as a catalyst!)
    st = Stoichiometry.from_reactants_products(reactants=[(2, "E"), (1, "S")], products=[(1, "E"), (1, "P")])
    assert st == Stoichiometry(vector={'E': -1, 'S': -1, 'P': 1}, catalysts=[])

    # Missing products or reactants
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=["R"], products=None)
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=None, products="P")

    # Bad products or reactants type
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants={"k": 666}, products="P")
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants="R", products=123)

    # Bad list elements in products or reactants
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=[1, 2], products="P")
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants="R", products=[("X", "Y", "Z")])
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=[(1, 2)], products="P")
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants="R", products=[("P", "Q")])

    # Reactants and the products can't be the same
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=["A"], products=["A"])
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=["A"], products=[(1, "A")])
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants="R", products="R")
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=[(2, "B")], products=[(2, "B")])
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=["A", "B"], products=["B", "A"])
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=["A", (3, "B")], products=["A", (3, "B")])
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=[(2, "A"), "B", "C"], products=["B", (1, "C"), (2, "A")])
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=["R", (2, "P")], products=[(2, "P"), "R"])
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=["A", "A"], products=[(2, "A")])
    with pytest.raises(Exception):
        Stoichiometry.from_reactants_products(reactants=["A", (2, "B"), "E", "A"], products=[(1, "E"), "B", (2, "A"), "B"])
