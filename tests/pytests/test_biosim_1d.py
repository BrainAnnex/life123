# See also "test_system_1d.py" and "test_membranes_1d.py"

import pytest
import numpy as np
from life123 import BioSim1D, ReactionRegistry, UniformCompartment, SpeciesRegistry



#########   TESTS OF INITIALIZATION, SETTING AND VIEWING    #########

def test_constructor():
    with pytest.raises(Exception):
        BioSim1D()              # Missing required arguments

    with pytest.raises(Exception):
        BioSim1D(n_bins="I'm not an integer")

    with pytest.raises(Exception):
        BioSim1D(n_bins=0)      # Must be at least 1

    with pytest.raises(Exception):
        BioSim1D(n_bins=3)      # At least one more arg is needed


    chem_data = SpeciesRegistry(id="A")

    bio = BioSim1D(n_bins=5, species_data=chem_data)
    #bio.describe_state()

    assert bio.n_bins == 5
    assert bio.n_species == 1
    assert bio.species_data == chem_data
    expected = np.zeros((1, 5), dtype=float)
    assert np.allclose(bio.system, expected)


    # New test
    chem_data = SpeciesRegistry(id=["A", "B", "C"])
    bio = BioSim1D(n_bins=15, species_data=chem_data)
    #bio.describe_state()
    assert bio.n_bins == 15
    assert bio.n_species == 3
    assert bio.species_data == chem_data
    expected = np.zeros((3,15), dtype=float)
    assert np.allclose(bio.system, expected)
    assert type(bio.reactions) == ReactionRegistry
    assert type(bio.reaction_dynamics) == UniformCompartment


    # New test
    rxn = UniformCompartment()
    with pytest.raises(Exception):
        BioSim1D(n_bins=5, reaction_handler=rxn)    # No chemicals yet registered


    # New test
    rxn = UniformCompartment(names=["A", "B", "C"])
    bio = BioSim1D(n_bins=5, reaction_handler=rxn)
    #bio.describe_state()
    assert bio.n_bins == 5
    assert bio.n_species == 3
    expected = np.zeros((3, 5), dtype=float)
    assert np.allclose(bio.system, expected)
    assert type(bio.reactions) == ReactionRegistry
    assert type(bio.reaction_dynamics) == UniformCompartment
    assert type(bio.species_data) == SpeciesRegistry



def test_get_reactions():
    pass

def test_get_reaction_handler():
    pass

def test_reaction_in_equilibrium():
    pass

def test_reset_system():
    pass


def test_enable_history():
    pass

def test_capture_snapshot():
    pass

def test_get_bin_history():
    pass





#########################################################################
#                                                                       #
#                               DIFFUSION                               #
#                                                                       #
#########################################################################

# Done in separate file "test_biosim_1d_diffusion.py"



#########################################################################
#                                                                       #
#                               REACTIONS                               #
#                                                                       #
#########################################################################


# Done in separate file "test_biosim_1d_reactions.py"



#########################################################################
#                                                                       #
#                         REACTION-DIFFUSION                            #
#                                                                       #
#########################################################################

def test_react_diffuse():
    pass
