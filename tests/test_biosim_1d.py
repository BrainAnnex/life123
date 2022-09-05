import pytest
import numpy as np
from life_1D.bio_sim_1d import BioSim1D as bio
from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions


# Do an initialization operation for every test that uses it as argument
@pytest.fixture()       # By default, the scope is the function
def biomsim1D():
    bio.reset_system()




#########   TESTS OF INITIALIZATION, SETTING AND VIEWING    #########

def test_initialize_system(biomsim1D):
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=10, chem_data=chem_data)

    assert bio.n_bins == 10
    assert bio.n_species == 1

    bio.describe_state()

    expected = np.zeros((1, 10), dtype=float)
    assert np.allclose(bio.system, expected)

    # New test
    chem_data = chem(names=["A", "B", "C"])
    bio.initialize_system(n_bins=15, chem_data=chem_data)
    bio.describe_state()
    expected = np.zeros((3,15), dtype=float)
    assert np.allclose(bio.system, expected)


def test_reset_system(biomsim1D):
    pass

def test_replace_system(biomsim1D):
    pass



def test_set_uniform_concentration(biomsim1D):
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=0.3)
    expected = np.full(5, 0.3, dtype=float)
    assert np.allclose(bio.lookup_species(0), expected)

    # New test
    bio.reset_system()
    chem_data = chem(names=["A", "B", "C"])
    bio.initialize_system(n_bins=15, chem_data=chem_data)
    bio.set_uniform_concentration(species_name="A", conc=10.)
    bio.set_uniform_concentration(species_name="B", conc=11.)
    bio.set_uniform_concentration(species_name="C", conc=12.)
    #bio.describe_state()
    assert np.allclose(bio.lookup_species(0), np.full(15, 10., dtype=float))
    assert np.allclose(bio.lookup_species(1), np.full(15, 11., dtype=float))
    assert np.allclose(bio.lookup_species(2), np.full(15, 12., dtype=float))

    # New test, with membranes
    bio.reset_system()
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)
    bio.set_membranes(membrane_pos=[1])   # A single membrane, passing thru bin 1

    bio.set_uniform_concentration(species_name="A", conc=8.8)
    expected = np.full(5, 8.8, dtype=float)
    assert np.allclose(bio.lookup_species(0), expected)
    #print(bio.system_B)
    expected = np.array([0, 8.8, 0, 0, 0])
    assert np.allclose(bio.lookup_species(species_name="A", trans_membrane=True), expected)




def test_set_all_uniform_concentrations(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.set_all_uniform_concentrations(conc_list=[3., 7.])
    assert np.allclose(bio.lookup_species(species_name="A"), np.full(5, 3., dtype=float))
    assert np.allclose(bio.lookup_species(species_name="B"), np.full(5, 7., dtype=float))

    #print(bio.system)



def test_set_bin_conc(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.set_bin_conc(bin_address=2, conc=3.14, species_name="A")
    assert np.allclose(bio.lookup_species(species_name="A"), [0, 0, 3.14, 0, 0])

    bio.set_bin_conc(bin_address=2, conc=2.78, species_index=0)
    assert np.allclose(bio.lookup_species(species_name="A"), [0, 0, 2.78, 0, 0])
    #print(bio.system)




def test_set_species_conc(biomsim1D):
    pass



def test_inject_conc_to_bin(biomsim1D):
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    with pytest.raises(Exception):
        # cell_index out of bounds
        bio.inject_conc_to_bin(bin_address=5, species_index=0, delta_conc=10.)

    with pytest.raises(Exception):
        # species_index out of bounds
        bio.inject_conc_to_bin(bin_address=5, species_index=1, delta_conc=10.)

    bio.inject_conc_to_bin(bin_address=1, species_index=0, delta_conc=10.)

    bio.describe_state()




########  DIMENSION-RELATED  ################

def test_set_dimensions(biomsim1D):
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=4, chem_data=chem_data)
    bio.set_dimensions(21.)
    assert np.allclose(bio.system_length, 21.)
    assert np.allclose(bio.global_Dx, 7.)

    with pytest.raises(Exception):
        bio.set_dimensions("I'm no number")
        bio.set_dimensions(0)


def test_x_coord(biomsim1D):
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=4, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.x_coord(0)      # must first call set_dimensions()

    bio.set_dimensions(21.)
    assert np.allclose(bio.x_coord(0), 0.)
    assert np.allclose(bio.x_coord(1), 7.)
    assert np.allclose(bio.x_coord(3), 21.)




########  MEMBRANE-RELATED  ################

def test_uses_membranes(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    assert not bio.uses_membranes()
    bio.set_membranes(membrane_pos=[2, 4])
    assert bio.uses_membranes()



def test_bins_with_membranes(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.set_membranes(membrane_pos=[2, 4])
    assert bio.bins_with_membranes() == [2, 4]



def test_set_membranes(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.set_membranes(membrane_pos=[2, 4])
    expected = np.array([False, False, True, False, True])
    assert (bio.membranes == expected).all()
    assert np.allclose(bio.A_fraction, [0, 0, 0.5, 0, 0.5])

    bio.set_membranes(membrane_pos=(0,2))
    expected = np.array([True, False, True, False, False])
    assert (bio.membranes == expected).all()
    assert np.allclose(bio.A_fraction, [0.5, 0, 0.5, 0, 0])

    with pytest.raises(Exception):  # Invalid bin numbers
        bio.set_membranes(membrane_pos=[5])
        bio.set_membranes(membrane_pos=[-1])
        bio.set_membranes(membrane_pos=[2, 4, 6])


    # Now, specify the "A fraction" of (some of) the membrane-containing bins
    bio.set_membranes(membrane_pos=[1, [2, .8], (4, .3)])
    expected = np.array([False, True, True, False, True])
    assert (bio.membranes == expected).all()
    assert np.allclose(bio.A_fraction, [0, 0.5, 0.8, 0, 0.3])

    # Re-do the last setting of membrane, this time to a system that has concentration values set
    bio.reset_system()
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.set_uniform_concentration(conc=8., species_name="A")
    bio.set_uniform_concentration(conc=5., species_name="B")
    assert np.allclose(bio.lookup_species(species_name="A"), [8, 8, 8, 8, 8])
    assert np.allclose(bio.lookup_species(species_name="B"), [5, 5, 5, 5, 5])
    assert bio.membranes is None

    bio.set_membranes(membrane_pos=[1, [2, .8], (4, .3)])
    expected = np.array([False, True, True, False, True])
    assert (bio.membranes == expected).all()
    assert np.allclose(bio.A_fraction, [0, 0.5, 0.8, 0, 0.3])
    assert np.allclose(bio.lookup_species(species_name="A"), [8, 8, 8, 8, 8])
    assert np.allclose(bio.lookup_species(species_name="B"), [5, 5, 5, 5, 5])
    assert np.allclose(bio.lookup_species(species_name="A", trans_membrane=True), [0, 8, 8, 0, 8])
    assert np.allclose(bio.lookup_species(species_name="B", trans_membrane=True), [0, 5, 5, 0, 5])

    #bio.show_membranes()




    #########################################################################
    #                                                                       #
    #                              TO VIEW                                  #
    #                                                                       #
    #########################################################################


def test_assert_valid_bin(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=3, chem_data=chem_data)
    with pytest.raises(Exception):
        bio.assert_valid_bin(-1)
        bio.assert_valid_bin(3)
        bio.assert_valid_bin("2")
        bio.assert_valid_bin((0, 2))

    bio.assert_valid_bin(0)
    bio.assert_valid_bin(2)



def test_lookup_species(biomsim1D):
    pass    # TODO


def test_bin_concentration(biomsim1D):
    pass    # TODO


def test_bin_snapshot(biomsim1D):
    pass    # TODO


def test_system_snapshot(biomsim1D):
    pass    # TODO



def test_show_membranes(biomsim1D):
    print("cls.membranes: ", bio.membranes)
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    result = bio.show_membranes()
    assert result == ""

    bio.set_membranes(membrane_pos=[2, 4])
    result = bio.show_membranes()
    assert result == "\n_____________________\n|   |   |0.5|   |0.5|\n---------------------"

    bio.set_membranes(membrane_pos=[1, [2, .8], (4, .3)])
    result = bio.show_membranes()
    assert result == "\n_____________________\n|   |0.5|0.8|   |0.3|\n---------------------"



#######  CHANGE RESOLUTIONS #########

def test_increase_spacial_resolution(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_list=[11., 12., 13.])
    bio.set_species_conc(species_index=1, conc_list=[5., 15., 25.])

    bio.increase_spacial_resolution(2)

    assert np.allclose(bio.lookup_species(species_index=0), [11., 11., 12., 12., 13., 13.])
    assert np.allclose(bio.lookup_species(species_name="B"), [ 5., 5.,  15., 15., 25., 25.])



def test_decrease_spacial_resolution(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=6, chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_list=[10., 20., 30., 40., 50., 60.])
    bio.set_species_conc(species_index=1, conc_list=[ 2., 8.,   5., 15., 4.,   2.])
    original_state = bio.system

    # Group by pairs
    bio.decrease_spacial_resolution(2)
    assert np.allclose(bio.system, [[15., 35., 55.],
                                    [ 5., 10.,  3.]])

    # Group by triplets
    bio.replace_system(original_state)
    bio.decrease_spacial_resolution(3)
    assert np.allclose(bio.system, [[20., 50.],
                                    [ 5., 7.]])

    # Group into just 1 single bin
    bio.replace_system(original_state)
    bio.decrease_spacial_resolution(6)
    assert np.allclose(bio.system, [[35.],
                                    [ 6.]])



def test_varying_spacial_resolution(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(species_name="A", conc_list=[11., 12., 13.])
    bio.set_species_conc(species_name="B", conc_list=[5., 15., 25.])

    # First, increase the spacial resolution, after saving the original matrix
    original_state = bio.system
    bio.increase_spacial_resolution(3)

    # Then, decrease the spacial resolution, by the same factor
    bio.decrease_spacial_resolution(3)
    assert np.allclose(bio.system, original_state)



def test_smooth_spacial_resolution(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(species_name="A", conc_list=[10., 20., 30.])
    bio.set_species_conc(species_name="B", conc_list=[2.,   8., 4.])

    bio.smooth_spacial_resolution()
    expected = np.array([[10., 15., 20., 25., 30.],
                         [ 2.,  5.,  8.,  6.,  4.]])
    assert np.allclose(bio.system, expected)
    assert bio.n_bins == 5
    assert bio.n_species == 2




#########   TESTS OF DIFFUSION : single species, one step    #########

def test_diffuse_step_single_species_1(biomsim1D):
    chem_data = chem(diffusion_rates=[10.])

    bio.initialize_system(n_bins=2, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, species_index=0, delta_conc=100.)
    bio.describe_state()

    with pytest.raises(Exception):
        #Excessive time step
        bio.diffuse_step_single_species(time_step=0.034)

    # Diffuse by a single step
    increment_vector = bio.diffuse_step_single_species(time_step=0.02)
    print("Increment vector is: ", increment_vector)
    assert np.allclose(increment_vector, [-20., 20.])



def test_diffuse_step_single_species_1b(biomsim1D):
    chem_data = chem(diffusion_rates=[10.])

    bio.initialize_system(n_bins=2, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, species_index=0, delta_conc=100.)
    bio.describe_state()

    # Diffuse by a single step
    bio.diffuse(time_step=0.02, n_steps=1)
    print(bio.system)
    assert np.allclose(bio.lookup_species(0), [80, 20])

    # Another single step
    bio.diffuse(time_step=0.01, n_steps=1)
    print(bio.system)
    assert np.allclose(bio.lookup_species(0), [74, 26])



def test_diffuse_step_single_species_2(biomsim1D):
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=1, chem_data=chem_data)
    bio.set_uniform_concentration(species_index=0, conc=8.0)
    with pytest.raises(Exception):
        bio.diffuse_step_single_species(time_step=.001)    # Must set the diffusion rates first

    chem_data.set_diffusion_rates([20.])
    bio.describe_state()    # 1 bins and 1 species:  [[8.]]

    increment_vector = bio.diffuse_step_single_species(time_step=3)    # With just 1 bin, nothing happens
    print("Increment vector is: ", increment_vector)
    assert np.allclose(increment_vector, [0.])



def test_diffuse_step_single_species_2b(biomsim1D):
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=1, chem_data=chem_data)
    bio.set_uniform_concentration(species_index=0, conc=8.0)

    chem_data.set_diffusion_rates([20.])
    bio.describe_state()    # 1 bins and 1 species:  [[8.]]

    bio.diffuse(time_step=3, n_steps=1)    # With just 1 bin, nothing happens
    bio.describe_state()
    assert np.allclose(bio.lookup_species(0), [8])



def test_diffuse_step_single_species_3(biomsim1D):
    bio.system = None
    with pytest.raises(Exception):
        bio.diffuse_step_single_species(time_step=.001)    # Must first initialize the system

    chem_data = chem(diffusion_rates=[2.78, 3.14])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=22.2)
    bio.set_uniform_concentration(species_index=1, conc=66.6)

    bio.describe_state()
    """
    5 bins and 2 species:
    Species 0. Diff rate: 2.78. Conc:  [22.2 22.2 22.2 22.2 22.2]
    Species 1. Diff rate: 3.14. Conc:  [66.6 66.6 66.6 66.6 66.6]
    """

    # Diffusing a uniform distribution won't change it
    bio.diffuse(time_step=0.08, n_steps=1)

    assert np.allclose(bio.lookup_species(0), np.full(5, 22.2, dtype=float))
    assert np.allclose(bio.lookup_species(1), np.full(5, 66.6, dtype=float))



def test_diffuse_step_4(biomsim1D):
    # Multiple diffusion steps, with 2 bins
    chem_data = chem(diffusion_rates=[1.])
    bio.initialize_system(n_bins=2, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, species_index=0, delta_conc=10.)
    bio.describe_state()
    """
    2 bins and 1 species:
    Species 0. Diff rate: 1.0. Conc:  [10.  0.]
    """

    for i in range(4):
        bio.diffuse_step(time_step=.3)
        bio.system += bio.delta_diffusion
        print(f"At end of step {i+1}:")
        bio.describe_state()

    assert np.allclose(bio.lookup_species(0), [5.128, 4.872])



def test_diffuse_step_5(biomsim1D):
    # Multiple diffusion steps, with 3 bins, and a large time step
    chem_data = chem(diffusion_rates=[.5])
    bio.initialize_system(n_bins=3, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=1, species_index=0, delta_conc=10.)
    bio.describe_state()

    # The time step is so large that the system immediately equilibrates
    print("The default max allowed time step is: ", bio.max_time_step(.5))
    for i in range(3):
        bio.diffuse_step(time_step=0.6666)
        bio.system += bio.delta_diffusion
        print(f"At end of step {i+1}:")
        print(bio.system)

    assert np.allclose(bio.lookup_species(0), np.full(3, 3.3333333, dtype=float))



def test_diffuse_step_6(biomsim1D):
    # Multiple diffusion steps, with 5 bins, and a large time step
    chem_data = chem(diffusion_rates=[.5])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, species_index=0, delta_conc=10.)
    bio.describe_state()

    print("The default max allowed time step is: ", bio.max_time_step(.5))
    for i in range(20):
        bio.diffuse_step(time_step=0.6666)
        bio.system += bio.delta_diffusion
        print(f"At end of step {i+1}:")
        print(bio.lookup_species(0))

    assert np.allclose(bio.lookup_species(0), [2.23752027, 2.14678423, 1.99998594, 1.85320708, 1.76250248])



def test_diffuse_step_7(biomsim1D):
    # Multiple diffusion steps, with 5 bins,
    # 1/2 the time step of the previous test, and double the duration
    chem_data = chem(diffusion_rates=[.5])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, species_index=0, delta_conc=10.)
    bio.describe_state()

    for i in range(20*2):
        bio.diffuse_step(time_step=0.6666/2)
        bio.system += bio.delta_diffusion
        if i<10 or i >35:
            print(f"At end of step {i+1}:")
            print(bio.lookup_species(0))

    assert np.allclose(bio.lookup_species(0), [2.26063875, 2.16100935, 1.99990821, 1.83893393, 1.73950977])



def test_diffuse_step_8(biomsim1D):
    # Many diffusion steps that the system equilibrates, no matter the starting point
    chem_data = chem(diffusion_rates=[.3])
    bio.initialize_system(n_bins=15, chem_data=chem_data)

    np.random.seed(18)
    bio.set_species_conc(species_index=0, conc_list= 100*np.random.rand(15))
    print()
    bio.describe_state()

    avg_conc = sum(bio.lookup_species(0)) / 15.
    print("Avg of concentrations: ", avg_conc)

    for i in range(2000):
        bio.diffuse_step(time_step=1)
        bio.system += bio.delta_diffusion
        if i<4:
            print(f"At end of step {i+1}:")
            print(bio.lookup_species(0))

    print(f"At end of ALL steps:")
    print(bio.lookup_species(0))
    # With such a large number of steps, all concentrations will
    # equilibrate to their average
    assert np.allclose(bio.lookup_species(0), np.full(15, avg_conc, dtype=float))



#########   TESTS OF DIFFUSION : all species, one step    #########

def test_diffuse_step_1(biomsim1D):
    chem_data = chem(diffusion_rates=[5., 20.])
    bio.initialize_system(n_bins=3, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.set_bin_conc(bin_address=3, species_index=0, conc=100.) # Bin number out of range

    bio.set_bin_conc(species_index=0, bin_address=1, conc=100.)

    with pytest.raises(Exception):
        bio.set_species_conc(species_index=2, conc_list=[10, 20, 50])   # species_index out of range

    bio.set_species_conc(species_index=1, conc_list=[10, 20, 50])

    bio.describe_state()
    """
    3 bins and 2 species:
    Species 0. Diff rate: 5.0. Conc:  [  0. 100.   0.]
    Species 1. Diff rate: 20.0. Conc:  [10. 20. 50.]
    """

    bio.diffuse_step(0.01)
    bio.system += bio.delta_diffusion
    bio.describe_state()
    """
    3 bins and 2 species:
     [[ 5. 90.  5.]
     [12. 24. 44.]]
    """
    assert np.allclose(bio.system, [[5., 90., 5.] , [12., 24., 44.]])



#########   TESTS OF DIFFUSION : all species, multiple steps    #########

def test_diffuse_1(biomsim1D):
    chem_data = chem(diffusion_rates=[5., 20.])
    bio.initialize_system(n_bins=3, chem_data=chem_data)
    bio.set_bin_conc(species_index=0, bin_address=1, conc=100.)
    bio.set_species_conc(species_index=1, conc_list=[10, 20, 50])
    bio.describe_state()
    """
    Species 0. Diff rate: 5.0. Conc:  [  0. 100.   0.]
    Species 1. Diff rate: 20.0. Conc:  [10. 20. 50.]   
    """

    with pytest.raises(Exception):
        bio.diffuse()               # Is not passing any arguments

    with pytest.raises(Exception):
        bio.diffuse(total_duration= 5)      # Is not passing enough arguments

    with pytest.raises(Exception):
        bio.diffuse(time_step = 0.2)        # Is not passing enough arguments

    with pytest.raises(Exception):
        bio.diffuse(n_steps=3)              # Is not passing enough arguments

    with pytest.raises(Exception):
        bio.diffuse(total_duration= 5, time_step = 0.2, n_steps=3)  # Too many args


    # Do 1 step
    bio.diffuse(time_step = 0.01, n_steps = 1)
    bio.describe_state(concise=True)
    """
    DIFFUSION STEP 0:
    [[95.  5.  0.]
     [12. 24. 44.]]
    """
    assert np.allclose(bio.system, [[5., 90., 5.] , [12., 24., 44.]])

    # Repeat 1 step
    bio.diffuse(time_step = 0.01, n_steps = 1)
    bio.describe_state(concise=True)
    """
    DIFFUSION STEP 0:
    [[ 9.25 81.5   9.25]
     [14.4  25.6  40.  ]]
    """
    assert np.allclose(bio.system, [[9.25, 81.5, 9.25] , [14.4, 25.6, 40.]])


    # Reset the system
    bio.set_species_conc(species_index=0, conc_list=[0, 100, 0])
    bio.set_species_conc(species_index=1, conc_list=[10, 20, 50])
    bio.describe_state()
    # Re-take the 2 steps
    bio.diffuse(time_step = 0.01, n_steps = 2)
    assert np.allclose(bio.system, [[9.25, 81.5, 9.25] , [14.4, 25.6, 40.]])



def test_diffuse_2(biomsim1D):
    chem_data = chem(diffusion_rates=[0.1])
    # Diffuse 1 species almost to equilibrium, starting from a single concentration pulse
    bio.initialize_system(n_bins=10, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=0.)
    bio.inject_conc_to_bin(species_index=0, bin_address=2, delta_conc=10.)

    status = bio.diffuse(total_duration=800, time_step=0.1)
    assert status['steps'] == 8000
    assert np.allclose(bio.lookup_species(0),
                                [1.00055275, 1.00049864, 1.00039572, 1.00025407, 1.00008755, 0.99991245,
                                 0.99974593, 0.99960428, 0.99950136, 0.99944725])




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

def test_react_diffuse(biomsim1D):
    pass




def manual_test_compare_states(biomsim1D):     # MANUAL TEST
    state1 = np.array([10, 20, 30])
    state2 = np.array([10.3, 19.9, 30.2])

    print()
    bio.compare_states(state1, state2, verbose=True)
    """
    Max of unsigned absolute differences:  0.3000000000000007
    Relative differences:  [-0.03        0.005      -0.00666667]
    Max of unsigned relative differences:  0.030000000000000072
    Mean of relative differences:  -0.010555555555555547
    Median of relative differences:  -0.006666666666666643
    Standard deviation of relative differences:  0.014550889837454275
    np.allclose with lax tolerance?  (rtol=1e-01, atol=1e-01) :  True
    np.allclose with mid tolerance?  (rtol=1e-02, atol=1e-03) :  False
    np.allclose with tight tolerance?  (rtol=1e-03, atol=1e-05) :  False
    np.allclose with extra-tight tolerance?  (rtol=1e-05, atol=1e-08) :  False
    """


    state1 = np.array([[10, 20, 30],
                       [100, 200, 300]])
    state2 = np.array([[10.3, 19.9, 30.2],
                       [103, 199, 302]])

    print()
    bio.compare_states(state1, state2, verbose=True)
    """
    Max of unsigned absolute differences:  3.0
    Relative differences:  [[-0.03        0.005      -0.00666667]
     [-0.03        0.005      -0.00666667]]
    Max of unsigned relative differences:  0.030000000000000072
    Mean of relative differences:  -0.010555555555555552
    Median of relative differences:  -0.006666666666666655
    Standard deviation of relative differences:  0.014550889837454246
    np.allclose with lax tolerance?  (rtol=1e-01, atol=1e-01) :  True
    np.allclose with mid tolerance?  (rtol=1e-02, atol=1e-03) :  False
    np.allclose with tight tolerance?  (rtol=1e-03, atol=1e-05) :  False
    np.allclose with extra-tight tolerance?  (rtol=1e-05, atol=1e-08) :  False
    """


    state1 = np.array([[10, 20, 30],
                       [100, 200, 300]])
    state2 = np.array([[10.0001, 19.9999, 30.0001],
                       [100.0004, 199.9985, 300.0003]])

    print()
    bio.compare_states(state1, state2, verbose=True)
    """
    Max of unsigned absolute differences:  0.0014999999999929514
    Relative differences:  [[-1.00000000e-05  5.00000000e-06 -3.33333333e-06]
     [-4.00000000e-06  7.50000000e-06 -1.00000000e-06]]
    Max of unsigned relative differences:  9.999999999976694e-06
    Mean of relative differences:  -9.722222222130482e-07
    Median of relative differences:  -2.166666666632011e-06
    Standard deviation of relative differences:  5.82651718172416e-06
    np.allclose with lax tolerance?  (rtol=1e-01, atol=1e-01) :  True
    np.allclose with mid tolerance?  (rtol=1e-02, atol=1e-03) :  True
    np.allclose with tight tolerance?  (rtol=1e-03, atol=1e-05) :  True
    np.allclose with extra-tight tolerance?  (rtol=1e-05, atol=1e-08) :  True
    """
