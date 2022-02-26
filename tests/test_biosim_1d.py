import numpy as np
import pytest
from life_1D.bio_sim_1d import BioSim1D as bio


def test_initialize_universe():
    bio.initialize_universe(n_bins=10, n_species=1)

    assert bio.n_bins == 10
    assert bio.n_species == 1

    bio.describe_state()

    expected = np.zeros((1, 10), dtype=float)
    assert np.allclose(bio.univ, expected)

    # New test
    bio.initialize_universe(n_bins=15, n_species=3)
    bio.describe_state()
    expected = np.zeros((3,15), dtype=float)
    assert np.allclose(bio.univ, expected)



def test_set_diffusion_rates():
    bio.initialize_universe(n_bins=5, n_species=1)
    bio.set_diffusion_rates([.2])
    print("diffusion_rates: ", bio.diffusion_rates)
    assert np.allclose(bio.diffusion_rates, [.2])

    # New test
    bio.initialize_universe(n_bins=25, n_species=4)
    bio.set_diffusion_rates([.1, .7, .2, .4])
    print("diffusion_rates: ", bio.diffusion_rates)
    assert np.allclose(bio.diffusion_rates, [.1, .7, .2, .4])



def test_set_uniform_concentration():
    bio.initialize_universe(n_bins=8, n_species=1)
    bio.set_uniform_concentration(species_index=0, conc=0.3)
    bio.describe_state()
    expected = np.full(8, 0.3, dtype=float)
    assert np.allclose(bio.lookup_species(0), expected)

    # New test
    bio.initialize_universe(n_bins=15, n_species=3)
    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=11.)
    bio.set_uniform_concentration(species_index=2, conc=12.)
    bio.describe_state()
    assert np.allclose(bio.lookup_species(0), np.full(15, 10., dtype=float))
    assert np.allclose(bio.lookup_species(1), np.full(15, 11., dtype=float))
    assert np.allclose(bio.lookup_species(2), np.full(15, 12., dtype=float))



def test_inject_conc_to_cell():
    bio.initialize_universe(n_bins=5, n_species=1)

    with pytest.raises(Exception):
        # cell_index out of bounds
        bio.inject_conc_to_cell(bin=5, species_index=0, delta_conc=10.)

    with pytest.raises(Exception):
        # species_index out of bounds
        bio.inject_conc_to_cell(bin=5, species_index=1, delta_conc=10.)

    bio.inject_conc_to_cell(bin=1, species_index=0, delta_conc=10.)

    bio.describe_state()



#########   TESTS OF DIFFUSION : single species, one step    #########

def test_diffuse_step_single_species_1():
    bio.initialize_universe(n_bins=2, n_species=1)

    bio.inject_conc_to_cell(bin=0, species_index=0, delta_conc=100.)
    bio.describe_state()

    bio.set_diffusion_rates([10.])

    with pytest.raises(Exception):
        #Excessive time step
        bio.diffuse_step_single_species(time_step=0.034)

    # Diffuse by a single step
    bio.diffuse_step_single_species(time_step=0.02)
    print(bio.univ)
    assert np.allclose(bio.lookup_species(0), [80, 20])

    # Another single step
    bio.diffuse_step_single_species(time_step=0.01)
    print(bio.univ)
    assert np.allclose(bio.lookup_species(0), [74, 26])



def test_diffuse_step_single_species_2():
    bio.initialize_universe(n_bins=1, n_species=1)
    bio.set_uniform_concentration(species_index=0, conc=8.0)
    with pytest.raises(Exception):
        bio.diffuse_step_single_species(time_step=.001)    # Must set the diffusion rates first

    bio.set_diffusion_rates([20.])
    bio.describe_state()    # 1 bins and 1 species:  [[8.]]

    bio.diffuse_step_single_species(time_step=3)    # With just 1 bin, nothing happens
    bio.describe_state()
    assert np.allclose(bio.lookup_species(0), [8])



def test_diffuse_step_single_species_3():
    bio.univ = None
    with pytest.raises(Exception):
        bio.diffuse_step_single_species(time_step=.001)    # Must first initialize the system

    bio.initialize_universe(n_bins=5, n_species=2)

    bio.set_uniform_concentration(species_index=0, conc=22.2)
    bio.set_uniform_concentration(species_index=1, conc=66.6)
    bio.set_diffusion_rates([2.78, 3.14])
    bio.describe_state(show_diffusion_rates=True)
    """
    5 bins and 2 species:
    Species 0. Diff rate: 2.78. Conc:  [22.2 22.2 22.2 22.2 22.2]
    Species 1. Diff rate: 3.14. Conc:  [66.6 66.6 66.6 66.6 66.6]
    """

    # Diffusing a uniform distribution won't change it
    bio.diffuse_step_single_species(time_step=0.08, species_index=0)
    bio.diffuse_step_single_species(time_step=0.08, species_index=1)

    assert np.allclose(bio.lookup_species(0), np.full(5, 22.2, dtype=float))
    assert np.allclose(bio.lookup_species(1), np.full(5, 66.6, dtype=float))



def test_diffuse_step_single_species_4():
    # Multiple diffusion steps, with 2 bins
    bio.initialize_universe(n_bins=2, n_species=1)

    bio.inject_conc_to_cell(bin=0, species_index=0, delta_conc=10.)
    bio.set_diffusion_rates([1.])
    bio.describe_state(show_diffusion_rates=True)
    """
    2 bins and 1 species:
    Species 0. Diff rate: 1.0. Conc:  [10.  0.]
    """

    for i in range(4):
        bio.diffuse_step_single_species(time_step=.3)
        print(f"At end of step {i+1}:")
        bio.describe_state()

    assert np.allclose(bio.lookup_species(0), [5.128, 4.872])



def test_diffuse_step_single_species_5():
    # Multiple diffusion steps, with 3 bins, and a large time step
    bio.initialize_universe(n_bins=3, n_species=1)

    bio.inject_conc_to_cell(bin=1, species_index=0, delta_conc=10.)
    bio.set_diffusion_rates([.5])
    bio.describe_state(show_diffusion_rates=True)

    # The time step is so large that the system immediately equilibrates
    print("The default max allowed time step is: ", bio.max_time_step(.5))
    for i in range(3):
        bio.diffuse_step_single_species(time_step=0.6666)
        print(f"At end of step {i+1}:")
        print(bio.univ)

    assert np.allclose(bio.lookup_species(0), np.full(3, 3.3333333, dtype=float))



def test_diffuse_step_single_species_6():
    # Multiple diffusion steps, with 5 bins, and a large time step
    bio.initialize_universe(n_bins=5, n_species=1)

    bio.inject_conc_to_cell(bin=0, species_index=0, delta_conc=10.)
    bio.set_diffusion_rates([.5])
    bio.describe_state(show_diffusion_rates=True)

    print("The default max allowed time step is: ", bio.max_time_step(.5))
    for i in range(20):
        bio.diffuse_step_single_species(time_step=0.6666)
        print(f"At end of step {i+1}:")
        print(bio.lookup_species(0))

    assert np.allclose(bio.lookup_species(0), [2.23752027, 2.14678423, 1.99998594, 1.85320708, 1.76250248])


def test_diffuse_step_single_species_7():
    # Multiple diffusion steps, with 5 bins,
    # 1/2 the time step of the previous test, and double the duration
    bio.initialize_universe(n_bins=5, n_species=1)

    bio.inject_conc_to_cell(bin=0, species_index=0, delta_conc=10.)
    bio.set_diffusion_rates([.5])
    bio.describe_state(show_diffusion_rates=True)

    for i in range(20*2):
        bio.diffuse_step_single_species(time_step=0.6666/2)
        if i<10 or i >35:
            print(f"At end of step {i+1}:")
            print(bio.lookup_species(0))

    assert np.allclose(bio.lookup_species(0), [2.26063875, 2.16100935, 1.99990821, 1.83893393, 1.73950977])



def test_diffuse_step_single_species_8():
    # Many diffusion steps that the system equilibrates, no matter the starting point
    bio.initialize_universe(n_bins=15, n_species=1)

    np.random.seed(18)
    bio.set_species_conc(species_index=0, conc_list= 100*np.random.rand(15))
    bio.set_diffusion_rates([.3])
    print()
    bio.describe_state(show_diffusion_rates=True)

    avg_conc = sum(bio.lookup_species(0)) / 15.
    print("Avg of concentrations: ", avg_conc)

    for i in range(2000):
        bio.diffuse_step_single_species(time_step=1)
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
    bio.initialize_universe(n_bins=3, n_species=2)

    with pytest.raises(Exception):
        bio.set_bin_conc(bin=3, species_index=0, conc=100.) # Bin number out of range

    bio.set_bin_conc(species_index=0, bin=1, conc=100.)

    with pytest.raises(Exception):
        bio.set_species_conc(species_index=2, conc_list=[10, 20, 50])   # species_index out of range

    bio.set_species_conc(species_index=1, conc_list=[10, 20, 50])

    bio.set_diffusion_rates([5., 20.])
    bio.describe_state(show_diffusion_rates=True)
    """
    3 bins and 2 species:
    Species 0. Diff rate: 5.0. Conc:  [  0. 100.   0.]
    Species 1. Diff rate: 20.0. Conc:  [10. 20. 50.]
    """

    bio.diffuse_step(0.01)
    bio.describe_state()
    """
    3 bins and 2 species:
     [[ 5. 90.  5.]
     [12. 24. 44.]]
    """
    assert np.allclose(bio.univ, [[ 5., 90.,  5.] , [12., 24., 44.]])



#########   TESTS OF DIFFUSION : all species, multiple steps    #########

def test_diffuse_1():
    bio.initialize_universe(n_bins=3, n_species=2)
    bio.set_bin_conc(species_index=0, bin=1, conc=100.)
    bio.set_species_conc(species_index=1, conc_list=[10, 20, 50])
    bio.set_diffusion_rates([5., 20.])
    bio.describe_state(show_diffusion_rates=True)
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
    assert np.allclose(bio.univ, [[ 5., 90.,  5.] , [12., 24., 44.]])

    # Repeat 1 step
    bio.diffuse(time_step = 0.01, n_steps = 1)
    bio.describe_state(concise=True)
    """
    DIFFUSION STEP 0:
    [[ 9.25 81.5   9.25]
     [14.4  25.6  40.  ]]
    """
    assert np.allclose(bio.univ, [[ 9.25, 81.5,  9.25] , [14.4,  25.6, 40.]])


    # Reset the system
    bio.set_species_conc(species_index=0, conc_list=[0, 100, 0])
    bio.set_species_conc(species_index=1, conc_list=[10, 20, 50])
    bio.describe_state()
    # Re-take the 2 steps
    bio.diffuse(time_step = 0.01, n_steps = 2)
    assert np.allclose(bio.univ, [[ 9.25, 81.5,  9.25] , [14.4,  25.6, 40.]])



def test_diffuse_2():
    # Diffuse 1 species almost to equilibrium, starting from a single concentration pulse
    bio.initialize_universe(n_bins=10, n_species=1)

    bio.set_uniform_concentration(species_index=0, conc=0.)
    bio.inject_conc_to_cell(species_index=0, bin=2, delta_conc=10.)
    bio.set_diffusion_rates([0.1])

    status = bio.diffuse(time_duration=800, time_step=0.1)
    assert status['steps'] == 8000
    assert np.allclose(bio.lookup_species(0),
                                [1.00055275, 1.00049864, 1.00039572, 1.00025407, 1.00008755, 0.99991245,
                                 0.99974593, 0.99960428, 0.99950136, 0.99944725])