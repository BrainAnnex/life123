import pytest
import numpy as np
from life123 import BioSim2D, Reactions, UniformCompartment, ChemData



#########   TESTS OF INITIALIZATION, SETTING AND VIEWING    #########

def test_constructor():
    with pytest.raises(Exception):
        BioSim2D()                  # Missing required arguments

    with pytest.raises(Exception):
        BioSim2D(x_bins=666, y_bins=3.14)        # Wrong data type

    with pytest.raises(Exception):
        BioSim2D(x_bins=0, y_bins=4)      # Must be at least 1 in each dimension

    with pytest.raises(Exception):
        BioSim2D(x_bins=4, y_bins=0)      # Must be at least 1 in each dimension

    with pytest.raises(Exception):
        BioSim2D(x_bins=3, y_bins=5)      # At least one more arg is needed


    chem_data = ChemData(names="A")
    bio = BioSim2D(x_bins=3, y_bins=5, chem_data=chem_data)
    #bio.describe_state()
    assert bio.n_bins_x == 3
    assert bio.n_bins_y == 5
    assert bio.n_species == 1
    assert bio.chem_data == chem_data
    expected = np.zeros((1, 3, 5), dtype=float)
    assert np.allclose(bio.system, expected)


    # New test
    chem_data = ChemData(names=["A", "B", "C", "D"])
    bio = BioSim2D(x_bins=3, y_bins=5, chem_data=chem_data)
    #bio.describe_state()
    assert bio.n_bins_x == 3
    assert bio.n_bins_y == 5
    assert bio.n_species == 4
    assert bio.chem_data == chem_data
    expected = np.zeros((4, 3, 5), dtype=float)
    assert np.allclose(bio.system, expected)
    assert type(bio.reactions) == Reactions
    assert type(bio.reaction_dynamics) == UniformCompartment


    # New test
    rxn = UniformCompartment()
    with pytest.raises(Exception):
        BioSim2D(x_bins=4, y_bins=2, reaction_handler=rxn)    # No chemicals yet registered


    # New test
    rxn = UniformCompartment(names=["A", "B", "C"])
    bio = BioSim2D(x_bins=4, y_bins=2, reaction_handler=rxn)
    #bio.describe_state()
    assert bio.n_bins_x == 4
    assert bio.n_bins_y == 2
    assert bio.n_species == 3
    expected = np.zeros((3, 4, 2), dtype=float)
    assert np.allclose(bio.system, expected)
    assert type(bio.reactions) == Reactions
    assert type(bio.reaction_dynamics) == UniformCompartment
    assert type(bio.chem_data) == ChemData



def test_system_size():
    bio = BioSim2D(x_bins=3, y_bins=5 , chem_data=ChemData(names="A"))
    assert bio.system_size() == (3,5)

    bio = BioSim2D(x_bins=9, y_bins=9 , chem_data=ChemData(names=["A", "B"]))
    assert bio.system_size() == (9,9)



def test_set_bin_conc():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim2D(x_bins=3, y_bins=4, chem_data=chem_data)

    bio.set_bin_conc(bin_address=(0,2), chem_label="A", conc=0.2)
    bio.describe_state()

    expected_0 = np.array([
                            [0.,  0.,  0.2, 0.],
                            [0. , 0.,  0. , 0.],
                            [0. , 0.,  0. , 0.]
                           ])

    expected_1 = np.zeros((3, 4), dtype=float)

    assert np.allclose(bio.system[0], expected_0)
    assert np.allclose(bio.system[1], expected_1)


    bio.set_bin_conc(bin_address=(2,3), chem_label="B", conc=1.23)
    bio.describe_state()

    expected_1 = np.array([
        [0.,  0.,  0. , 0.],
        [0. , 0.,  0. , 0.],
        [0. , 0.,  0. , 1.23]
    ])
    assert np.allclose(bio.system[0], expected_0)
    assert np.allclose(bio.system[1], expected_1)



def test_set_bin_conc_all_species():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim2D(x_bins=3, y_bins=4, chem_data=chem_data)

    bio.set_bin_conc_all_species(bin_address=(0,2), conc_list=[0.02, 1.02])
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



def test_set_species_conc():
    chem_data = ChemData(names=["A"])
    bio = BioSim2D(x_bins=2, y_bins=3, chem_data=chem_data)   # 2 rows and 3 columns
    conc_matrix = [[50, 80, 20], [10, 60, 0]]

    bio.set_species_conc(conc_data=conc_matrix, species_index=0)
    #bio.describe_state()
    assert np.allclose(bio.lookup_species(species_name="A") , conc_matrix)

    bio.set_species_conc(conc_data=conc_matrix, species_name="A")
    assert np.allclose(bio.lookup_species(species_index=0) , conc_matrix)

    bio.set_species_conc(conc_data=((50, 80, 20), (10, 60, 0)), species_name="A")      # Tuple instead of list
    assert np.allclose(bio.lookup_species(species_index=0) , conc_matrix)

    bio.set_species_conc(conc_data=np.array(conc_matrix), species_name="A")
    assert np.allclose(bio.lookup_species(species_index=0) , conc_matrix)

    with pytest.raises(Exception):
        bio.set_species_conc(conc_data=conc_matrix)    # Missing chemical species
        bio.set_species_conc(conc_data="Do I look like a list??", species_name="A")
        bio.set_species_conc(conc_data=[1., 2., 3.], species_name="A")                  # Not a matrix
        bio.set_species_conc(conc_data=[[1., 2., 3.]], species_name="A")                # Wrong number of rows
        bio.set_species_conc(conc_data=[[1., 2.], [3., 4.]], species_name="A")          # Wrong number of columns
        bio.set_species_conc(conc_data=[[50, -80, 20], [10, 60, 0]], species_name="A")  # Negative concentrations
        bio.set_species_conc(conc_data=np.array([[50, 80, 20], [10, 60, -0.01]]), species_name="A")  # Negative concentrations




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
    chem_data = ChemData(names=["A", "B"])

    bio = BioSim2D(x_bins=3, y_bins=4, chem_data=chem_data)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    bio.reactions.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    bio.set_bin_conc_all_species(bin_address=(0,0), conc_list=[10.,50.])
    bio.set_bin_conc_all_species(bin_address=(0,1), conc_list=[20.,35.])
    bio.set_bin_conc_all_species(bin_address=(2,3), conc_list=[5.,100.])
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
    chem_data = ChemData(names=["A", "B"])

    bio = BioSim2D(x_bins=3, y_bins=4, chem_data=chem_data)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    bio.reactions.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    bio.set_bin_conc_all_species(bin_address=(0,0), conc_list=[10.,50.])
    bio.set_bin_conc_all_species(bin_address=(0,1), conc_list=[20.,35.])
    bio.set_bin_conc_all_species(bin_address=(2,3), conc_list=[5.,100.])
    #bio.describe_state()

    bio.reaction_step(delta_time=0.1)
    #print("bio.delta_reactions:\n", bio.delta_reactions)

    assert np.allclose(bio.delta_reactions[0], [[ 7., 1., 0., 0.] , [ 0., 0., 0., 0.] , [ 0., 0., 0., 18.5]])
    assert np.allclose(bio.delta_reactions[1], [[-7., -1., 0., 0.] , [ 0., 0., 0., 0.] , [ 0., 0., 0., -18.5]])

    #bio.system += bio.delta_reactions       # Matrix operation to update all the concentrations
    #bio.system_time += 0.1

    #bio.describe_state()

    #conc_dict = {0: 5., 1: 100.}
    #result = chem_data.compute_all_rate_deltas(conc_dict=conc_dict, delta_time=0.1)
    #print(result)
