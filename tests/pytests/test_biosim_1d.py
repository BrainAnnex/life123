import pytest
import numpy as np
from life123 import BioSim1D, Reactions, UniformCompartment, ChemData



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


    chem_data = ChemData(names="A")

    bio = BioSim1D(n_bins=5, chem_data=chem_data)
    #bio.describe_state()

    assert bio.n_bins == 5
    assert bio.n_species == 1
    assert bio.chem_data == chem_data
    expected = np.zeros((1, 5), dtype=float)
    assert np.allclose(bio.system, expected)


    # New test
    chem_data = ChemData(names=["A", "B", "C"])
    bio = BioSim1D(n_bins=15, chem_data=chem_data)
    #bio.describe_state()
    assert bio.n_bins == 15
    assert bio.n_species == 3
    assert bio.chem_data == chem_data
    expected = np.zeros((3,15), dtype=float)
    assert np.allclose(bio.system, expected)
    assert type(bio.reactions) == Reactions
    assert type(bio.reaction_dynamics) == UniformCompartment


    # New test
    rxn = UniformCompartment()
    with pytest.raises(Exception):
        BioSim1D(n_bins=5, reaction_handler=rxn)    # No chemicals yet registered


    # New test
    rxn = UniformCompartment(names=["A", "B", "C"])
    bio = BioSim1D(n_bins=5, reaction_handler=rxn)
    bio.describe_state()
    assert bio.n_bins == 5
    assert bio.n_species == 3
    expected = np.zeros((3, 5), dtype=float)
    assert np.allclose(bio.system, expected)
    assert type(bio.reactions) == Reactions
    assert type(bio.reaction_dynamics) == UniformCompartment
    assert type(bio.chem_data) == ChemData



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
