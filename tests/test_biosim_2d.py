import pytest
import numpy as np
from life_2D.bio_sim_2d import BioSim2D
from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions



#########   TESTS OF INITIALIZATION, SETTING AND VIEWING    #########

def test_initialize_system():
    chem_data = chem(names=["A", "B"])
    bio = BioSim2D(n_bins=(3,5), chem_data=chem_data)

    assert bio.n_bins_x == 3
    assert bio.n_bins_y == 5

    bio.describe_state()

    expected = np.zeros((2, 3, 5), dtype=float)
    assert np.allclose(bio.system, expected)



def test_set_bin_conc():
    chem_data = chem(names=["A", "B"])
    bio = BioSim2D(n_bins=(3,4), chem_data=chem_data)

    bio.set_bin_conc(bin_x=0, bin_y=2, species_index=0, conc=0.2)
    bio.describe_state()

    expected_0 = np.array([
                            [0.,  0.,  0.2, 0.],
                            [0. , 0.,  0. , 0.],
                            [0. , 0.,  0. , 0.]
                           ])

    expected_1 = np.zeros((3, 4), dtype=float)

    assert np.allclose(bio.system[0], expected_0)
    assert np.allclose(bio.system[1], expected_1)


    bio.set_bin_conc(bin_x=2, bin_y=3, species_index=1, conc=1.23)
    bio.describe_state()

    expected_1 = np.array([
        [0.,  0.,  0. , 0.],
        [0. , 0.,  0. , 0.],
        [0. , 0.,  0. , 1.23]
    ])
    assert np.allclose(bio.system[0], expected_0)
    assert np.allclose(bio.system[1], expected_1)



def test_set_bin_conc_all_species():
    chem_data = chem(names=["A", "B"])
    bio = BioSim2D(n_bins=(3,4), chem_data=chem_data)

    bio.set_bin_conc_all_species(bin_x=0, bin_y=2, conc_list=[0.02, 1.02])
    bio.describe_state()

    expected_0 = np.array([
        [0.,  0.,  0.02, 0.],
        [0. , 0.,  0.  , 0.],
        [0. , 0.,  0.  , 0.]
    ])

    expected_1 = np.array([
        [0.,  0.,  1.02, 0.],
        [0. , 0.,  0.  , 0.],
        [0. , 0.,  0.  , 0.]
    ])

    assert np.allclose(bio.system[0], expected_0)
    assert np.allclose(bio.system[1], expected_1)




#########################################################################
#                                                                       #
#                               DIFFUSION                               #
#                                                                       #
#########################################################################

# Done in separate file "test_biosim_2d_diffusion.py"




#########################################################################
#                                                                       #
#                               REACTIONS                               #
#                                                                       #
#########################################################################

def test_react():
    chem_data = chem(names=["A", "B"])

    rxn = Reactions(chem_data)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    bio = BioSim2D(n_bins=(3,4), chem_data=chem_data, reactions=rxn)
    bio.set_bin_conc_all_species(bin_x=0, bin_y=0, conc_list=[10.,50.])
    bio.set_bin_conc_all_species(bin_x=0, bin_y=1, conc_list=[20.,35.])
    bio.set_bin_conc_all_species(bin_x=2, bin_y=3, conc_list=[5.,100.])
    bio.describe_state()

    bio.react(time_step=0.1, n_steps=1)
    bio.describe_state()

    assert np.allclose(bio.system[0], [[ 17., 21., 0., 0.] , [ 0., 0., 0., 0.] , [ 0., 0., 0., 23.5]])
    assert np.allclose(bio.system[1], [[43., 34., 0., 0.] , [ 0., 0., 0., 0.] , [ 0., 0., 0., 81.5]])

    # Continue to equilibrium
    bio.react(time_step=0.1, n_steps=20)
    bio.describe_state()

    assert np.allclose(bio.system[1, 0, 0] / bio.system[0, 0, 0], 1.5)  # From the ration forward/back reaction rates
    assert np.allclose(bio.system[1, 0, 1] / bio.system[0, 0, 1], 1.5)
    assert np.allclose(bio.system[1, 2, 3] / bio.system[0, 2, 3], 1.5)



def test_reaction_step():
    chem_data = chem(names=["A", "B"])

    rxn = Reactions(chem_data)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    bio = BioSim2D(n_bins=(3,4), chem_data=chem_data, reactions=rxn)
    bio.set_bin_conc_all_species(bin_x=0, bin_y=0, conc_list=[10.,50.])
    bio.set_bin_conc_all_species(bin_x=0, bin_y=1, conc_list=[20.,35.])
    bio.set_bin_conc_all_species(bin_x=2, bin_y=3, conc_list=[5.,100.])
    #bio.describe_state()

    bio.reaction_step(delta_time=0.1)
    #print("bio.delta_reactions:\n", bio.delta_reactions)

    assert np.allclose(bio.delta_reactions[0], [[ 7., 1., 0., 0.] , [ 0., 0., 0., 0.] , [ 0., 0., 0., 18.5]])
    assert np.allclose(bio.delta_reactions[1], [[-7., -1., 0., 0.] , [ 0., 0., 0., 0.] , [ 0., 0., 0., -18.5]])

    #bio.system += bio.delta_reactions       # Matrix operation to update all the concentrations
    #bio.system_time += 0.1

    #bio.describe_state()

    #conc_dict = {0: 5., 1: 100.}
    #result = rxn.compute_all_rate_deltas(conc_dict=conc_dict, delta_time=0.1)
    #print(result)