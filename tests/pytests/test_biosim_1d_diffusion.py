# These are tests specifically for diffusion in 1D;
# for general tests of 1D system, see test_biosim_1d.py


import pytest
import numpy as np
from scipy.ndimage import shift
#from life123 import BioSim1D
from life123.bio_sim_1d_experimental import BioSim1D_NEW, Diffusion1D
from life123 import ChemData as chem
from life123 import Numerical as num
from life123 import CollectionArray



##############################################################################################

def test_diffuse_step_1():
    # Test with just 1 bin
    chem_data = chem(names="A")
    bio = BioSim1D_NEW(n_bins=1, chem_data=chem_data)

    bio.set_uniform_concentration(chem_index=0, conc=8.0)

    chem_data.set_diffusion_rate(label="A", diff_rate = 20.)
    #bio.describe_state()    # 1 bins and 1 species:  [[8.]]

    bio.diffuse_step(time_step=3)
    assert np.allclose(bio.delta_diffusion, [[0.]])     # With just 1 bin, nothing happens!

    bio.diffuse_step(time_step=77., algorithm="5_1_explicit")
    assert np.allclose(bio.delta_diffusion, [[0.]])     # With just 1 bin, nothing happens!



def test_diffuse_step_4():
    # Multiple diffusion steps, with 2 bins
    chem_data = chem(diffusion_rates=[1.])
    bio = BioSim1D_NEW(n_bins=2, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=10.)
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
    bio = BioSim1D_NEW(n_bins=3, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=1, chem_index=0, delta_conc=10.)
    #bio.describe_state()

    # The time step is so large that the system immediately equilibrates
    #diff_obj = Diffusion1D(n_bins=3)
    #print("The default max allowed time step is: ", diff_obj.max_time_step(.5, delta_x=1))

    for i in range(3):
        bio.diffuse_step(time_step=0.6666)
        bio.system += bio.delta_diffusion
        #print(f"At end of step {i+1}:")
        #print(bio.system)

    assert np.allclose(bio.lookup_species(0), np.full(3, 3.3333333, dtype=float))



def test_diffuse_step_6():
    # Multiple diffusion steps, with 5 bins, and a large time step
    chem_data = chem(diffusion_rates=0.5)
    bio = BioSim1D_NEW(n_bins=5, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=10.)
    #bio.describe_state()

    # The time step is its max allowed value
    #diff_obj = Diffusion1D(n_bins=3)
    #print("The default max allowed time step is: ", diff_obj.max_time_step(.5, delta_x=1))
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
    bio = BioSim1D_NEW(n_bins=5, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=10.)
    #bio.describe_state()

    for i in range(20*2):
        bio.diffuse_step(time_step=0.6666/2)
        bio.system += bio.delta_diffusion
        """
        if i<10 or i >35:
            print(f"At end of step {i+1}:")
            print(bio.lookup_species(0))
        """

    assert np.allclose(bio.lookup_species(0), [2.26063875, 2.16100935, 1.99990821, 1.83893393, 1.73950977])



def test_diffuse_step_8():
    # Many diffusion steps that the system equilibrates, no matter the starting point
    chem_data = chem(diffusion_rates=[.3])
    bio = BioSim1D_NEW(n_bins=15, chem_data=chem_data)

    np.random.seed(18)
    bio.set_species_conc(chem_index=0, conc_list=100 * np.random.rand(15))
    #bio.describe_state()

    avg_conc = sum(bio.lookup_species(0)) / 15.
    print("Avg of concentrations: ", avg_conc)

    for i in range(2000):
        bio.diffuse_step(time_step=1)
        bio.system += bio.delta_diffusion
        '''
        if i<4:
            print(f"At end of step {i+1}:")
            print(bio.lookup_species(0))
        '''

    #print(f"At end of ALL steps:")
    #print(bio.lookup_species(0))
    # With such a large number of steps, all concentrations will
    # equilibrate to their average
    assert np.allclose(bio.lookup_species(0), np.full(15, avg_conc, dtype=float))




#########   TESTS OF DIFFUSION : all species, one step    #########

def test_diffuse_step_1_b():
    chem_data = chem(diffusion_rates=[5., 20.])
    bio = BioSim1D_NEW(n_bins=3, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.set_bin_conc(bin_address=3, chem_index=0, conc=100.) # Bin number out of range

    bio.set_bin_conc(chem_index=0, bin_address=1, conc=100.)

    with pytest.raises(Exception):
        bio.set_species_conc(chem_index=2, conc_list=[10, 20, 50])   # species_index out of range

    bio.set_species_conc(chem_index=1, conc_list=[10, 20, 50])

    #bio.describe_state()
    """
    3 bins and 2 species:
    Species 0. Diff rate: 5.0. Conc:  [  0. 100.   0.]
    Species 1. Diff rate: 20.0. Conc:  [10. 20. 50.]
    """

    bio.diffuse_step(0.01)
    bio.system += bio.delta_diffusion
    #bio.describe_state()
    """
    3 bins and 2 species:
     [[ 5. 90.  5.]
     [12. 24. 44.]]
    """
    assert np.allclose(bio.system, [[5., 90., 5.] , [12., 24., 44.]])



#########   TESTS OF DIFFUSION : all species, multiple steps    #########

def test_diffuse_1():
    chem_data = chem(diffusion_rates=[10.])

    bio = BioSim1D_NEW(n_bins=2, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=100.)
    #bio.describe_state()

    # Diffuse by a single step
    bio.diffuse(time_step=0.02, n_steps=1)
    #print(bio.system)
    assert np.allclose(bio.lookup_species(0), [80, 20])

    # Another single step
    bio.diffuse(time_step=0.01, n_steps=1)
    #print(bio.system)
    assert np.allclose(bio.lookup_species(0), [74, 26])



def test_diffuse_2():
    chem_data = chem(names=["A"])
    bio = BioSim1D_NEW(n_bins=1, chem_data=chem_data)
    bio.set_uniform_concentration(chem_index=0, conc=8.0)

    chem_data.set_diffusion_rate(label="A", diff_rate = 20.)
    #bio.describe_state()    # 1 bins and 1 species:  [[8.]]

    bio.diffuse(time_step=3, n_steps=1)    # With just 1 bin, nothing happens
    #bio.describe_state()
    assert np.allclose(bio.lookup_species(0), [8])



def test_diffuse_3():
    n_bins=5

    chem_data = chem(names="A")
    bio = BioSim1D_NEW(n_bins=n_bins, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.diffuse(time_step=0.08, n_steps=1)    # Must first initialize the system


    chem_data = chem(diffusion_rates=[2.78, 3.14])
    bio = BioSim1D_NEW(n_bins=n_bins, chem_data=chem_data)

    bio.set_uniform_concentration(chem_index=0, conc=22.2)
    bio.set_uniform_concentration(chem_index=1, conc=66.6)

    #bio.describe_state()
    """
    5 bins and 2 species:
    Species 0. Diff rate: 2.78. Conc:  [22.2 22.2 22.2 22.2 22.2]
    Species 1. Diff rate: 3.14. Conc:  [66.6 66.6 66.6 66.6 66.6]
    """

    # Diffusing a uniform distribution won't change it
    bio.diffuse(time_step=0.08, n_steps=1)

    assert np.allclose(bio.lookup_species(0), np.full(5, 22.2, dtype=float))
    assert np.allclose(bio.lookup_species(1), np.full(5, 66.6, dtype=float))



def test_diffuse_4():
    chem_data = chem(diffusion_rates=[5., 20.])
    bio = BioSim1D_NEW(n_bins=3, chem_data=chem_data)
    bio.set_bin_conc(chem_index=0, bin_address=1, conc=100.)
    bio.set_species_conc(chem_index=1, conc_list=[10, 20, 50])
    #bio.describe_state()
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
    #bio.describe_state(concise=True)
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
    bio.set_species_conc(chem_index=0, conc_list=[0, 100, 0])
    bio.set_species_conc(chem_index=1, conc_list=[10, 20, 50])
    bio.describe_state()
    # Re-take the 2 steps
    bio.diffuse(time_step = 0.01, n_steps = 2)
    assert np.allclose(bio.system, [[9.25, 81.5, 9.25] , [14.4, 25.6, 40.]])



def test_diffuse_5():
    chem_data = chem(diffusion_rates=[0.1])
    # Diffuse 1 species almost to equilibrium, starting from a single concentration pulse
    bio = BioSim1D_NEW(n_bins=10, chem_data=chem_data)

    bio.set_uniform_concentration(chem_index=0, conc=0.)
    bio.inject_conc_to_bin(chem_index=0, bin_address=2, delta_conc=10.)

    status = bio.diffuse(total_duration=800, time_step=0.1)
    assert status['steps'] == 8000
    assert np.allclose(bio.lookup_species(0),
                                [1.00055275, 1.00049864, 1.00039572, 1.00025407, 1.00008755, 0.99991245,
                                 0.99974593, 0.99960428, 0.99950136, 0.99944725])



def test_diffuse_6():
    # Compare a low-level and a higher-level function for one-step diffusion
    n_bins=5
    delta_t = 0.01
    delta_x = 2
    diff = 10.

    chem_data = chem(diffusion_rates=diff)
    initial_concs = np.array([50, 80, 40, 100, 120])

    bio = BioSim1D_NEW(n_bins=n_bins, chem_data=chem_data)
    bio.set_species_conc(chem_index=0, conc_list=initial_concs)
    #bio.describe_state()

    # Compute the increments to the concentrations, from a single diffusion step,
    # using a low-level function
    diffusion_obj = Diffusion1D(n_bins=n_bins)
    increment_vector = diffusion_obj.diffuse_step_5_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs, delta_x=delta_x)


    # Redo computations on an identical system
    bio2 = BioSim1D_NEW(n_bins=n_bins, chem_data=chem_data)
    bio2.set_species_conc(chem_index=0, conc_list=initial_concs)

    status = bio2.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x, algorithm="5_1_explicit")
    assert status["steps"] == 1
    bio2.describe_state()

    assert np.allclose(initial_concs + increment_vector, bio2.lookup_species(chem_index=0))



def test_diffuse_7():
    # Based on experiment 1D/diffusion/validate_diffusion_3

    # Parameters of the simulation run.  We'll be considering just 1 chemical species, "A"
    diffusion_rate = 10.
    delta_t = 0.01
    n_bins = 5000
    delta_x = 2
    algorithm = None    # This corresponds to a 3+1 stencil, explicit method

    chem_data = chem(diffusion_rates=[diffusion_rate], names=["A"])

    bio = BioSim1D_NEW(n_bins=n_bins, chem_data=chem_data)

    bio.inject_sine_conc(chem_label="A", number_cycles=1, amplitude=12, bias=40)
    bio.inject_sine_conc(chem_label="A", number_cycles=2, amplitude=10)
    bio.inject_sine_conc(chem_label="A", number_cycles=16, amplitude=5)

    history = CollectionArray()   # All the system state will get collected in this object

    # Store the initial state
    arr = bio.lookup_species(chem_index=0, copy=True)
    history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

    # Do the 4 rounds of single-step diffusion; accumulate all data in the history object
    for _ in range(4):
        status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)
        assert status["steps"] == 1

        arr = bio.lookup_species(chem_index=0, copy=True)
        history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")


    # Now, let's examine the data collected at the 5 time points
    all_history = history.get_array()
    assert all_history.shape == (5, n_bins)

    # Compute time derivatives (for each bin), using 5-point stencils
    df_dt_all_bins = np.apply_along_axis(num.gradient_order4_1d, 0, all_history, delta_t)

    # Let's consider the state at the midpoint in time (t2)
    f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
    assert f_at_t2.shape == (n_bins, )

    # Computer the second spatial derivative, using 5-point stencils
    gradient_x_at_t2 = num.gradient_order4_1d(arr=f_at_t2, dx=delta_x)
    second_gradient_x_at_t2 = num.gradient_order4_1d(arr=gradient_x_at_t2, dx=delta_x)
    assert second_gradient_x_at_t2.shape == (n_bins, )

    # Compare the left and right hand sides of the diffusion equation
    lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
    rhs = diffusion_rate*second_gradient_x_at_t2

    dist = num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end
    assert np.allclose(dist, 0.0017647994920801059)



def test_diffuse_8():
    # Based on experiment 1D/diffusion/validate_diffusion_3
    # Identical to test_diffuse_4(), but for a different diffusion-computing algorithm

    # Parameters of the simulation run.  We'll be considering just 1 chemical species, "A"
    diffusion_rate = 10.
    delta_t = 0.01
    n_bins = 5000
    delta_x = 2
    algorithm = "5_1_explicit"  # A 5+1 stencil, explicit method

    chem_data = chem(diffusion_rates=[diffusion_rate], names=["A"])

    bio = BioSim1D_NEW(n_bins=n_bins, chem_data=chem_data)

    bio.inject_sine_conc(chem_label="A", number_cycles=1, amplitude=12, bias=40)
    bio.inject_sine_conc(chem_label="A", number_cycles=2, amplitude=10)
    bio.inject_sine_conc(chem_label="A", number_cycles=16, amplitude=5)

    history = CollectionArray()   # All the system state will get collected in this object

    # Store the initial state
    arr = bio.lookup_species(chem_index=0, copy=True)
    history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

    # Do the 4 rounds of single-step diffusion; accumulate all data in the history object
    for _ in range(4):
        status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)
        assert status["steps"] == 1

        arr = bio.lookup_species(chem_index=0, copy=True)
        history.store(par=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")


    # Now, let's examine the data collected at the 5 time points
    all_history = history.get_array()
    assert all_history.shape == (5, n_bins)

    # Compute time derivatives (for each bin), using 5-point stencils
    df_dt_all_bins = np.apply_along_axis(num.gradient_order4_1d, 0, all_history, delta_t)

    # Let's consider the state at the midpoint in time (t2)
    f_at_t2 = all_history[2]     # The middle of the 5 time snapshots
    assert f_at_t2.shape == (n_bins, )

    # Computer the second spatial derivative, using 5-point stencils
    gradient_x_at_t2 = num.gradient_order4_1d(arr=f_at_t2, dx=delta_x)
    second_gradient_x_at_t2 = num.gradient_order4_1d(arr=gradient_x_at_t2, dx=delta_x)
    assert second_gradient_x_at_t2.shape == (n_bins, )

    # Compare the left and right hand sides of the diffusion equation
    lhs = df_dt_all_bins[2]   # t2 is the middle point of the 5
    rhs = diffusion_rate*second_gradient_x_at_t2

    dist = num.compare_vectors(lhs, rhs, trim_edges=2)  # Euclidean distance, ignoring 2 edge points at each end
    assert np.allclose(dist, 0.003517310789846865)
