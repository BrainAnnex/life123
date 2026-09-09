import numpy as np
import math
import pytest
import pandas as pd
from life123.species_registry import Species, SpeciesRegistry, MacroMolecules
from life123.reactions_new import Stoichiometry, ReactionThermodynamics, Reaction
from life123.kinetics import Kinetics
from tests.utilities.comparisons import *



def test_CONSTRUCTOR_Kinetics():
    with pytest.raises(Exception):
        Kinetics(law="never heard of this!")

    k = Kinetics(law="mass action")
    assert k.law == "mass action"
    assert k.get_parameters() == {'kF': None, 'kR': None, 'reversible': False, 'K': None}

    k = Kinetics(law="mass action", parameters={"kF": 10})
    assert k.law == "mass action"
    assert k.get_parameters() == {"kF": 10, "kR": None, "reversible": False, "K": None}

    k = Kinetics(law="mass action", parameters={"kF": 10, "kR": 2})
    assert k.law == "mass action"
    assert k.get_parameters() == {"kF": 10, "kR": 2, "reversible": True, "K": 5}

    k = Kinetics(law="mass action", parameters={"kF": 10, "K": 5})
    assert k.law == "mass action"
    assert k.get_parameters() == {"kF": 10, "kR": 2, "reversible": True, "K": 5}

    k = Kinetics(law="mass action", parameters={"kR": 2, "K": 5})
    assert k.law == "mass action"
    assert k.get_parameters() == {"kF": 10, "kR": 2, "reversible": True, "K": 5}


    Kinetics(law="mass action", parameters={"kF": 10, "kR": 2, "K": 5})     # Consistent
    assert k.law == "mass action"
    assert k.get_parameters() == {"kF": 10, "kR": 2, "reversible": True, "K": 5}

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 10.01, "kR": 2, "K": 5})       # In-consistent

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 10, "kR": 2.01, "K": 5})       # In-consistent

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 10, "kR": 2, "K": 5.01})  # In-consistent


    k = Kinetics(law="mass action", parameters={"kF": 10, "kR": 0})
    assert k.law == "mass action"
    assert k.get_parameters() == {"kF": 10, "kR": 0, "reversible": False, "K": math.inf}

    k = Kinetics(law="mass action", parameters={"kF": 0, "kR": 0})
    assert k.law == "mass action"
    assert k.get_parameters() == {"kF": 0, "kR": 0, "reversible": False, "K": None}

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 10, "K": 0})

    k = Kinetics(law="mass action", parameters={"kF": 0, "K": 0})
    assert k.law == "mass action"
    assert k.get_parameters() == {"kF": 0, "kR": None, "reversible": False, "K": 0}


    k = Kinetics(law="mass action", parameters={"kF": 10, "K": math.inf})
    assert k.law == "mass action"
    assert k.get_parameters() == {"kF": 10, "kR": 0, "reversible": False, "K": math.inf}

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 0, "K": math.inf})

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": "Yo!", "kR": 2, "K": 5})

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": -0.001, "kR": 2, "K": 5})

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 10, "kR": "Yo!", "K": 5})

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 10, "kR": -0.001, "K": 5})

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 10, "kR": 2, "K": "Yo!"})

    with pytest.raises(Exception):
        Kinetics(law="mass action", parameters={"kF": 10, "kR": 2, "K": -0.001})



def test_to_dict_Kinetics():
    k = Kinetics(law="mass action")
    assert k.to_dict() == {"kinetics_type": "mass action", "reversible": False}

    k = Kinetics(law="mass action", parameters={"kF": 10})
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 10, "reversible": False}

    k = Kinetics(law="mass action", parameters={"kF": 10, "kR": 2})
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 10, "kR": 2, "K": 5, "reversible": True}

    k = Kinetics(law="mass action", parameters={"kF": 10})
    assert k.to_dict() == {"kinetics_type": "mass action", "kF": 10, "reversible": False}

    with pytest.raises(Exception):
        Kinetics(law="random name")
