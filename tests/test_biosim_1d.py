import pytest
import numpy as np
from life_1D.bio_sim_1d import BioSim1D as bio
from modules.chemicals.chemicals import Chemicals as chem



#########   TESTS OF INITIALIZATION, SETTING AND VIEWING    #########

def test_initialize_universe():
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



def test_set_diffusion_rates():
    chem_data = chem()
    chem_data.set_diffusion_rates([.2])
    print("diffusion_rates: ", chem_data.diffusion_rates)
    assert np.allclose(chem_data.diffusion_rates, [.2])

    # New test
    chem_data = chem()
    chem_data.set_diffusion_rates([.1, .7, .2, .4])
    print("diffusion_rates: ", chem_data.diffusion_rates)
    assert np.allclose(chem_data.diffusion_rates, [.1, .7, .2, .4])



def test_set_uniform_concentration():
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=8, chem_data=chem_data)
    bio.set_uniform_concentration(species_index=0, conc=0.3)
    bio.describe_state()
    expected = np.full(8, 0.3, dtype=float)
    assert np.allclose(bio.lookup_species(0), expected)

    # New test
    chem_data = chem(names=["A", "B", "C"])
    bio.initialize_system(n_bins=15, chem_data=chem_data)
    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=11.)
    bio.set_uniform_concentration(species_index=2, conc=12.)
    bio.describe_state()
    assert np.allclose(bio.lookup_species(0), np.full(15, 10., dtype=float))
    assert np.allclose(bio.lookup_species(1), np.full(15, 11., dtype=float))
    assert np.allclose(bio.lookup_species(2), np.full(15, 12., dtype=float))



def test_inject_conc_to_cell():
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    with pytest.raises(Exception):
        # cell_index out of bounds
        bio.inject_conc_to_bin(bin=5, species_index=0, delta_conc=10.)

    with pytest.raises(Exception):
        # species_index out of bounds
        bio.inject_conc_to_bin(bin=5, species_index=1, delta_conc=10.)

    bio.inject_conc_to_bin(bin=1, species_index=0, delta_conc=10.)

    bio.describe_state()



def test_increase_spacial_resolution():
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(0, [11., 12., 13.])
    bio.set_species_conc(1, [5., 15., 25.])
    #bio.describe_state()

    result = bio.increase_spacial_resolution(2)
    #print(result)
    assert np.allclose(result[0], [11., 11., 12., 12., 13., 13.])
    assert np.allclose(result[1], [ 5., 5.,  15., 15., 25., 25.])


def test_decrease_spacial_resolution():
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=6, chem_data=chem_data)
    bio.set_species_conc(0, [10., 20., 30., 40., 50., 60.])
    bio.set_species_conc(1, [ 2., 8.,   5., 15., 4.,   2.])
    bio.describe_state()

    # Group by pairs
    result = bio.decrease_spacial_resolution(2)
    #print(result)
    assert np.allclose(result, [[15., 35., 55.],
                                [ 5., 10.,  3.]])

    # Group by triplets
    result = bio.decrease_spacial_resolution(3)
    #print(result)
    assert np.allclose(result, [[20., 50.],
                                [ 5., 7.]])

    # Group into just 1 single bin
    result = bio.decrease_spacial_resolution(6)
    #print(result)
    assert np.allclose(result, [[35.],
                                [ 6.]])


def test_varying_spacial_resolution():
    chem_data = chem(names=["A", "B"])
    bio.initialize_system(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(0, [11., 12., 13.])
    bio.set_species_conc(1, [5., 15., 25.])

    # First, increase the spacial resolution, after saving the original matrix
    original_state = bio.system
    result = bio.increase_spacial_resolution(3)
    #print(result)
    bio.replace_system(result)

    # Then, decrease the spacial resolution, by the same factor
    result = bio.decrease_spacial_resolution(3)
    #print(result)
    assert np.allclose(result, original_state)




#########   TESTS OF DIFFUSION : single species, one step    #########

def test_diffuse_step_single_species_1():
    chem_data = chem(diffusion_rates=[10.])

    bio.initialize_system(n_bins=2, chem_data=chem_data)

    bio.inject_conc_to_bin(bin=0, species_index=0, delta_conc=100.)
    bio.describe_state()

    with pytest.raises(Exception):
        #Excessive time step
        bio.diffuse_step_single_species(time_step=0.034)

    # Diffuse by a single step
    increment_vector = bio.diffuse_step_single_species(time_step=0.02)
    print("Increment vector is: ", increment_vector)
    assert np.allclose(increment_vector, [-20., 20.])



def test_diffuse_step_single_species_1b():
    chem_data = chem(diffusion_rates=[10.])

    bio.initialize_system(n_bins=2, chem_data=chem_data)

    bio.inject_conc_to_bin(bin=0, species_index=0, delta_conc=100.)
    bio.describe_state()

    # Diffuse by a single step
    bio.diffuse(time_step=0.02, n_steps=1)
    print(bio.system)
    assert np.allclose(bio.lookup_species(0), [80, 20])

    # Another single step
    bio.diffuse(time_step=0.01, n_steps=1)
    print(bio.system)
    assert np.allclose(bio.lookup_species(0), [74, 26])



def test_diffuse_step_single_species_2():
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



def test_diffuse_step_single_species_2b():
    chem_data = chem(names=["A"])
    bio.initialize_system(n_bins=1, chem_data=chem_data)
    bio.set_uniform_concentration(species_index=0, conc=8.0)

    chem_data.set_diffusion_rates([20.])
    bio.describe_state()    # 1 bins and 1 species:  [[8.]]

    bio.diffuse(time_step=3, n_steps=1)    # With just 1 bin, nothing happens
    bio.describe_state()
    assert np.allclose(bio.lookup_species(0), [8])



def test_diffuse_step_single_species_3():
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



def test_diffuse_step_4():
    # Multiple diffusion steps, with 2 bins
    chem_data = chem(diffusion_rates=[1.])
    bio.initialize_system(n_bins=2, chem_data=chem_data)

    bio.inject_conc_to_bin(bin=0, species_index=0, delta_conc=10.)
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



def test_diffuse_step_5():
    # Multiple diffusion steps, with 3 bins, and a large time step
    chem_data = chem(diffusion_rates=[.5])
    bio.initialize_system(n_bins=3, chem_data=chem_data)

    bio.inject_conc_to_bin(bin=1, species_index=0, delta_conc=10.)
    bio.describe_state()

    # The time step is so large that the system immediately equilibrates
    print("The default max allowed time step is: ", bio.max_time_step(.5))
    for i in range(3):
        bio.diffuse_step(time_step=0.6666)
        bio.system += bio.delta_diffusion
        print(f"At end of step {i+1}:")
        print(bio.system)

    assert np.allclose(bio.lookup_species(0), np.full(3, 3.3333333, dtype=float))



def test_diffuse_step_6():
    # Multiple diffusion steps, with 5 bins, and a large time step
    chem_data = chem(diffusion_rates=[.5])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.inject_conc_to_bin(bin=0, species_index=0, delta_conc=10.)
    bio.describe_state()

    print("The default max allowed time step is: ", bio.max_time_step(.5))
    for i in range(20):
        bio.diffuse_step(time_step=0.6666)
        bio.system += bio.delta_diffusion
        print(f"At end of step {i+1}:")
        print(bio.lookup_species(0))

    assert np.allclose(bio.lookup_species(0), [2.23752027, 2.14678423, 1.99998594, 1.85320708, 1.76250248])



def test_diffuse_step_7():
    # Multiple diffusion steps, with 5 bins,
    # 1/2 the time step of the previous test, and double the duration
    chem_data = chem(diffusion_rates=[.5])
    bio.initialize_system(n_bins=5, chem_data=chem_data)

    bio.inject_conc_to_bin(bin=0, species_index=0, delta_conc=10.)
    bio.describe_state()

    for i in range(20*2):
        bio.diffuse_step(time_step=0.6666/2)
        bio.system += bio.delta_diffusion
        if i<10 or i >35:
            print(f"At end of step {i+1}:")
            print(bio.lookup_species(0))

    assert np.allclose(bio.lookup_species(0), [2.26063875, 2.16100935, 1.99990821, 1.83893393, 1.73950977])



def test_diffuse_step_8():
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

def test_diffuse_step_1():
    chem_data = chem(diffusion_rates=[5., 20.])
    bio.initialize_system(n_bins=3, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.set_bin_conc(bin=3, species_index=0, conc=100.) # Bin number out of range

    bio.set_bin_conc(species_index=0, bin=1, conc=100.)

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

def test_diffuse_1():
    chem_data = chem(diffusion_rates=[5., 20.])
    bio.initialize_system(n_bins=3, chem_data=chem_data)
    bio.set_bin_conc(species_index=0, bin=1, conc=100.)
    bio.set_species_conc(species_index=1, conc_list=[10, 20, 50])
    bio.describe_state()
    """
    Species 0. Diff rate: 5.0. Conc:  [  0. 100.   0.]
    Species 1. Diff rate: 20.0. Conc:  [10. 20. 50.]   
    """

    with pytest.raises(Exception):
        bio.diffuse()               # Is not passing any arguments

    with pytest.raises(Exception):
        bio.diffuse(time_duration = 5)      # Is not passing enough arguments

    with pytest.raises(Exception):
        bio.diffuse(time_step = 0.2)        # Is not passing enough arguments

    with pytest.raises(Exception):
        bio.diffuse(n_steps=3)              # Is not passing enough arguments

    with pytest.raises(Exception):
        bio.diffuse(time_duration = 5, time_step = 0.2, n_steps=3)  # Too many args


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



def test_diffuse_2():
    chem_data = chem(diffusion_rates=[0.1])
    # Diffuse 1 species almost to equilibrium, starting from a single concentration pulse
    bio.initialize_system(n_bins=10, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=0.)
    bio.inject_conc_to_bin(species_index=0, bin=2, delta_conc=10.)

    status = bio.diffuse(time_duration=800, time_step=0.1)
    assert status['steps'] == 8000
    assert np.allclose(bio.lookup_species(0),
                                [1.00055275, 1.00049864, 1.00039572, 1.00025407, 1.00008755, 0.99991245,
                                 0.99974593, 0.99960428, 0.99950136, 0.99944725])
