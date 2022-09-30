import pytest
import numpy as np
import pandas as pd
from scipy.ndimage import shift
from pandas.testing import assert_frame_equal
from life_1D.bio_sim_1d import BioSim1D
from modules.chemicals.chemicals import Chemicals as chem
from modules.numerical.numerical import Numerical as num


# Do an initialization operation for every test that uses it as argument
@pytest.fixture()       # By default, the scope is the function
def biomsim1D():
    pass
    #bio.reset_system()




#########   TESTS OF INITIALIZATION, SETTING AND VIEWING    #########

def test_initialize_system(biomsim1D):
    chem_data = chem(names=["A"])

    bio = BioSim1D(n_bins=10, chem_data=chem_data)

    assert bio.n_bins == 10
    assert bio.n_species == 1

    bio.describe_state()

    expected = np.zeros((1, 10), dtype=float)
    assert np.allclose(bio.system, expected)

    # New test
    chem_data = chem(names=["A", "B", "C"])
    bio = BioSim1D(n_bins=15, chem_data=chem_data)
    bio.describe_state()
    expected = np.zeros((3,15), dtype=float)
    assert np.allclose(bio.system, expected)



def test_reset_system(biomsim1D):
    pass

def test_replace_system(biomsim1D):
    pass



def test_set_uniform_concentration(biomsim1D):
    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=0.3)
    expected = np.full(5, 0.3, dtype=float)
    assert np.allclose(bio.lookup_species(0), expected)

    # New test
    bio.reset_system()
    chem_data = chem(names=["A", "B", "C"])
    bio = BioSim1D(n_bins=15, chem_data=chem_data)
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
    bio = BioSim1D(n_bins=5, chem_data=chem_data)
    bio.set_membranes(membrane_pos=[1])   # A single membrane, passing thru bin 1

    bio.set_uniform_concentration(species_name="A", conc=8.8)
    expected = np.full(5, 8.8, dtype=float)
    assert np.allclose(bio.lookup_species(0), expected)
    #print(bio.system_B)
    expected = np.array([0, 8.8, 0, 0, 0])
    assert np.allclose(bio.lookup_species(species_name="A", trans_membrane=True), expected)




def test_set_all_uniform_concentrations(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.set_all_uniform_concentrations(conc_list=[3., 7.])
    assert np.allclose(bio.lookup_species(species_name="A"), np.full(5, 3., dtype=float))
    assert np.allclose(bio.lookup_species(species_name="B"), np.full(5, 7., dtype=float))

    #print(bio.system)



def test_set_bin_conc(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.set_bin_conc(bin_address=2, conc=3.14, species_name="A")
    assert np.allclose(bio.lookup_species(species_name="A"), [0, 0, 3.14, 0, 0])

    bio.set_bin_conc(bin_address=2, conc=2.78, species_index=0)
    assert np.allclose(bio.lookup_species(species_name="A"), [0, 0, 2.78, 0, 0])
    #print(bio.system)




def test_set_species_conc(biomsim1D):
    pass



def test_inject_conc_to_bin(biomsim1D):
    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    with pytest.raises(Exception):
        # cell_index out of bounds
        bio.inject_conc_to_bin(bin_address=5, species_index=0, delta_conc=10.)

    with pytest.raises(Exception):
        # species_index out of bounds
        bio.inject_conc_to_bin(bin_address=5, species_index=1, delta_conc=10.)

    bio.inject_conc_to_bin(bin_address=1, species_index=0, delta_conc=10.)

    bio.describe_state()



def test_inject_sine_conc(biomsim1D):
    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=8, chem_data=chem_data)

    bio.inject_sine_conc(species_name="A", amplitude=5, bias=20, frequency=1, phase=0)
    assert np.allclose(bio.lookup_species(species_name="A"),
                       [20., 23.53553391, 25., 23.53553391, 20., 16.46446609, 15., 16.46446609])

    # Adding a 2nd wave in opposite phase will erase everything except for the biases
    bio.inject_sine_conc(species_name="A", amplitude=5, bias=10, frequency=1, phase=180)
    assert np.allclose(bio.lookup_species(species_name="A"),
                       [30., 30., 30., 30., 30., 30., 30., 30.])


    bio.set_uniform_concentration(conc=0, species_name="A")     # RESET concentrations

    bio.inject_sine_conc(species_name="A", amplitude=5, bias=20, frequency=1, phase=90)
    assert np.allclose(bio.lookup_species(species_name="A"),
                       [15., 16.46446609, 20., 23.53553391, 25.,  23.53553391, 20., 16.46446609])

    # Adding a 2nd wave in opposite phase will erase everything except for the biases
    bio.inject_sine_conc(species_name="A", amplitude=5, bias=10, frequency=1, phase=-90)
    assert np.allclose(bio.lookup_species(species_name="A"),
                       [30., 30., 30., 30., 30., 30., 30., 30.])


    bio.set_uniform_concentration(conc=0, species_name="A")     # RESET concentrations

    with pytest.raises(Exception):
        bio.inject_sine_conc(species_name="A", amplitude=10, bias=0, frequency=1/2, phase=180)  # Would lead to negative values


    bio.set_uniform_concentration(conc=0, species_name="A")     # RESET concentrations

    bio.inject_sine_conc(species_name="A", amplitude=10, bias=0, frequency=1/2, phase=0)

    assert np.allclose(bio.lookup_species(species_name="A"),
                       [0., 3.82683432, 7.07106781, 9.23879533, 10., 9.23879533, 7.07106781, 3.82683432])

    # Adding a 2nd wave in opposite phase will erase everything except for the biases
    bio.inject_sine_conc(species_name="A", amplitude=10, bias=30, frequency=1/2, phase=180)
    assert np.allclose(bio.lookup_species(species_name="A"),
                       [30., 30., 30., 30., 30., 30., 30., 30.])


    bio.set_uniform_concentration(conc=0, species_name="A")     # RESET concentrations
    bio.inject_sine_conc(species_name="A", amplitude=10, bias=30, frequency=2, phase=0)
    assert np.allclose(bio.lookup_species(species_name="A"),
                       [30., 40., 30., 20., 30., 40., 30., 20.])

    bio.set_uniform_concentration(conc=0, species_name="A")     # RESET concentrations
    bio.inject_sine_conc(species_name="A", amplitude=10, bias=30, frequency=2, phase=-90)
    assert np.allclose(bio.lookup_species(species_name="A"),
                       [40., 30., 20., 30., 40., 30., 20., 30.])

    bio.set_uniform_concentration(conc=0, species_name="A")     # RESET concentrations
    bio.inject_sine_conc(species_name="A", amplitude=10, bias=30, frequency=2, phase=90)
    assert np.allclose(bio.lookup_species(species_name="A"),
                       [20., 30., 40., 30., 20., 30., 40., 30.])

    bio.set_uniform_concentration(conc=0, species_name="A")     # RESET concentrations
    bio.inject_sine_conc(species_name="A", amplitude=10, bias=0, frequency=2, phase=0, zero_clip=True)
    assert np.allclose(bio.lookup_species(species_name="A"),
                       [0., 10., 0., 0., 0., 10., 0., 0.])

    #print(bio.lookup_species(species_name="A"))

    # Curves could be visualized as follows:
    #import plotly.express as px
    #fig = px.line(y=bio.lookup_species(species_name="A"))
    #fig.show()



def test_frequency_analysis(biomsim1D):
    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=100, chem_data=chem_data)

    bio.inject_sine_conc(species_name="A", frequency=2, amplitude=1, bias=3)

    result = bio.frequency_analysis(species_name="A")

    expected = pd.DataFrame({"Frequency": [0,2], "Relative Amplitude": [3.0,1.0]})
    assert_frame_equal(result, expected, check_dtype=False)

    bio.inject_sine_conc(species_name="A", frequency=4, amplitude=0.5)
    bio.inject_sine_conc(species_name="A", frequency=8, amplitude=0.2)

    #import plotly.express as px
    #fig = px.line(y=bio.lookup_species(species_name="A"))
    #fig.show()

    result = bio.frequency_analysis(species_name="A")
    #print(result)
    expected = pd.DataFrame({"Frequency": [0,2,4,8], "Relative Amplitude": [3.0,1.0,0.5,0.2]})
    assert_frame_equal(result, expected, check_dtype=False)

    result = bio.frequency_analysis(species_name="A", n_largest=3)  # ditch the smallest amplitude
    #print(result)
    expected = pd.DataFrame({"Frequency": [0,2,4], "Relative Amplitude": [3.0,1.0,0.5]})
    assert_frame_equal(result, expected, check_dtype=False)




########  DIMENSION-RELATED  ################

def test_set_dimensions(biomsim1D):
    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=4, chem_data=chem_data)
    bio.set_dimensions(21.)
    assert np.allclose(bio.system_length, 21.)
    assert np.allclose(bio.global_Dx, 7.)

    with pytest.raises(Exception):
        bio.set_dimensions("I'm no number")
        bio.set_dimensions(0)


def test_x_coord(biomsim1D):
    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=4, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.x_coord(0)      # must first call set_dimensions()

    bio.set_dimensions(21.)
    assert np.allclose(bio.x_coord(0), 0.)
    assert np.allclose(bio.x_coord(1), 7.)
    assert np.allclose(bio.x_coord(3), 21.)




########  MEMBRANE-RELATED  ################

def test_uses_membranes(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    assert not bio.uses_membranes()
    bio.set_membranes(membrane_pos=[2, 4])
    assert bio.uses_membranes()



def test_bins_with_membranes(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.set_membranes(membrane_pos=[2, 4])
    assert bio.bins_with_membranes() == [2, 4]



def test_set_membranes(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
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
    bio = BioSim1D()
    print("cls.membranes: ", bio.membranes)
    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_list=[11., 12., 13.])
    bio.set_species_conc(species_index=1, conc_list=[5., 15., 25.])

    bio.increase_spacial_resolution(2)

    assert np.allclose(bio.lookup_species(species_index=0), [11., 11., 12., 12., 13., 13.])
    assert np.allclose(bio.lookup_species(species_name="B"), [ 5., 5.,  15., 15., 25., 25.])



def test_double_spacial_resolution_linear_inter(biomsim1D):
    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=2, chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_list=[1., 2.])
    bio.double_spacial_resolution_linear()
    assert np.allclose(bio.system, [[1., 1.5, 2.]])

    chem_data = chem(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_list=[11., 12., 13.])
    bio.set_species_conc(species_index=1, conc_list=[5., 15., 25.])
    bio.double_spacial_resolution_linear()
    assert np.allclose(bio.system, [[11.,11.5,12.,12.5,13.],
                                    [ 5.,10. ,15.,20. ,25.]])

    bio = BioSim1D(n_bins=1, chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_list=[3.14])
    with pytest.raises(Exception):
        bio.double_spacial_resolution_linear()



def test_decrease_spacial_resolution(biomsim1D):
    chem_data = chem(names=["A", "B"])
    bio = BioSim1D(n_bins=6, chem_data=chem_data)
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
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
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
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
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
    delta_t = 0.01
    delta_x = 2
    diff = 10.

    # The coefficients for the 2nd derivative of order 2
    C1 = 1
    C0 = -2
    coeffs = np.array([C1, C0, C1])
    assert np.allclose(np.sum(coeffs), 0)   # The coefficients should add up to zero

    chem_data = chem(diffusion_rates=[diff])

    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    initial_concs = np.array([50, 80, 20])
    bio.set_species_conc(species_index=0, conc_list=initial_concs)
    #bio.describe_state()

    result = bio.diffuse_step_single_species(time_step=delta_t, species_index=0, delta_x=delta_x)
    #print(result)

    # Manually computes the expected increments at each bin
    concs = shift(initial_concs, 1, cval=initial_concs[0])      # [50, 50, 80]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, result[0])

    concs = initial_concs                                       # [50, 80, 20]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, result[1])

    concs = shift(initial_concs, -1, cval=initial_concs[-1])    # [80, 20, 20]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, result[2])



def test_diffuse_step_single_species_1a(biomsim1D):
    chem_data = chem(diffusion_rates=[10.])

    bio = BioSim1D(n_bins=2, chem_data=chem_data)

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

    bio = BioSim1D(n_bins=2, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=1, chem_data=chem_data)
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
    bio = BioSim1D(n_bins=1, chem_data=chem_data)
    bio.set_uniform_concentration(species_index=0, conc=8.0)

    chem_data.set_diffusion_rates([20.])
    bio.describe_state()    # 1 bins and 1 species:  [[8.]]

    bio.diffuse(time_step=3, n_steps=1)    # With just 1 bin, nothing happens
    bio.describe_state()
    assert np.allclose(bio.lookup_species(0), [8])



def test_diffuse_step_single_species_3(biomsim1D):
    bio = BioSim1D()
    bio.system = None
    with pytest.raises(Exception):
        bio.diffuse_step_single_species(time_step=.001)    # Must first initialize the system

    chem_data = chem(diffusion_rates=[2.78, 3.14])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=2, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=1, species_index=0, delta_conc=10.)
    bio.describe_state()

    # The time step is so large that the system immediately equilibrates
    print("The default max allowed time step is: ", bio.max_time_step(.5, delta_x=1))
    for i in range(3):
        bio.diffuse_step(time_step=0.6666)
        bio.system += bio.delta_diffusion
        print(f"At end of step {i+1}:")
        print(bio.system)

    assert np.allclose(bio.lookup_species(0), np.full(3, 3.3333333, dtype=float))



def test_diffuse_step_6(biomsim1D):
    # Multiple diffusion steps, with 5 bins, and a large time step
    chem_data = chem(diffusion_rates=[.5])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, species_index=0, delta_conc=10.)
    bio.describe_state()

    print("The default max allowed time step is: ", bio.max_time_step(.5, delta_x=1))
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
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=15, chem_data=chem_data)

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



###     Alternate methods to compute single-step diffusion    ###

def test_diffuse_step_single_species_5_1_stencils(biomsim1D):
    delta_t = 0.01
    delta_x = 2
    diff = 10.

    # The coefficients for the 2nd derivative of order 4
    C2 = -1/12
    C1 = 4/3
    C0 = -5/2
    coeffs = np.array([C2, C1, C0, C1, C2])
    assert np.allclose(np.sum(coeffs), 0)   # The coefficients should add up to zero

    chem_data = chem(diffusion_rates=[diff])

    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    initial_concs = np.array([50, 80, 40, 100, 120])
    bio.set_species_conc(species_index=0, conc_list=initial_concs)

    #bio.describe_state()

    result = bio.diffuse_step_single_species_5_1_stencils(time_step=delta_t, species_index=0, delta_x=delta_x)
    #print(result)

    # Manually computes the expected increments at each bin
    concs = shift(initial_concs, 2, cval=initial_concs[0])      # [50, 50, 50, 80, 40]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, result[0])

    concs = shift(initial_concs, 1, cval=initial_concs[0])      # [50, 50, 80, 40, 100]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, result[1])

    concs = initial_concs                                       # [50, 80, 40, 100, 120]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, result[2])

    concs = shift(initial_concs, -1, cval=initial_concs[-1])    # [80, 40, 100, 120, 120]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, result[3])

    concs = shift(initial_concs, -2, cval=initial_concs[-1])    # [40, 100, 120, 120, 120])
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, result[4])




#########   TESTS OF DIFFUSION : all species, one step    #########

def test_diffuse_step_1(biomsim1D):
    chem_data = chem(diffusion_rates=[5., 20.])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
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
    bio = BioSim1D(n_bins=10, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=0.)
    bio.inject_conc_to_bin(species_index=0, bin_address=2, delta_conc=10.)

    status = bio.diffuse(total_duration=800, time_step=0.1)
    assert status['steps'] == 8000
    assert np.allclose(bio.lookup_species(0),
                                [1.00055275, 1.00049864, 1.00039572, 1.00025407, 1.00008755, 0.99991245,
                                 0.99974593, 0.99960428, 0.99950136, 0.99944725])



def test_diffuse_3(biomsim1D):
    # Compare a low-level and a higher-level function for one-step diffusion
    delta_t = 0.01
    delta_x = 2
    diff = 10.

    chem_data = chem(diffusion_rates=[diff])

    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    initial_concs = np.array([50, 80, 40, 100, 120])
    bio.set_species_conc(species_index=0, conc_list=initial_concs)

    bio.describe_state()
    # Compute the increments to the concentrations, from a single diffusion step
    incs = bio.diffuse_step_single_species_5_1_stencils(time_step=delta_t, species_index=0, delta_x=delta_x)

    # Redo computations on an identical system
    bio2 = BioSim1D(n_bins=5, chem_data=chem_data)
    bio2.set_species_conc(species_index=0, conc_list=initial_concs)

    status = bio2.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x, algorithm="5_1_explicit")
    assert status["steps"] == 1
    bio2.describe_state()

    assert np.allclose(initial_concs + incs , bio2.lookup_species(species_index=0))



def test_diffuse_4(biomsim1D):
    # TODO: doesn't work!
    delta_t = 0.001
    delta_x = 2
    diff = 10.
    i = 3
    f_t = []

    chem_data = chem(diffusion_rates=[diff])

    bio = BioSim1D(n_bins=17, chem_data=chem_data)

    initial_concs = np.array([50, 58, 65, 72, 80, 70, 60, 50, 40, 55, 70, 85, 100, 105, 110, 115, 120])
    bio.set_species_conc(species_index=0, conc_list=initial_concs)
    bio.describe_state()
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))

    status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x) # , algorithm="5_1_explicit"
    assert status["steps"] == 1
    bio.describe_state()    # [ 51.02083333 , 77.6875 , 43.14583333 , 98.66666667 , 119.5]
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))

    f = bio.lookup_species(species_index=0)
    gradient_x = np.gradient(f, delta_x)
    #print("gradient_x: ", gradient_x)   # [13.33333333, -1.96875   ,  5.24479167, 19.08854167, 10.41666667]
    second_gradient_x = np.gradient(gradient_x, delta_x)
    print("second_gradient_x: ", second_gradient_x) # [-7.65104167, -2.02213542,  5.26432292,  1.29296875, -4.3359375 ]

    # A second diffusion step
    status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x) # , algorithm="5_1_explicit"
    assert status["steps"] == 1
    bio.describe_state()    # [ 51.92612847  75.65907118  45.97243924  97.51072049 118.96462674]
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))

    print(f"\nf_t at x{i}:", f_t)
    print(f"df_t at x{i}:", np.gradient(f_t, delta_t))
    print(f"df_t at x{i}, t1:", np.gradient(f_t, delta_t)[1])
    print(f"D d2/dx2 at x{i}:", diff*second_gradient_x[i])



def test_diffuse_5(biomsim1D):
    delta_t = 0.015
    delta_x = 2
    diff = 10.
    i = 34
    f_t = []
    algorithm = None # "5_1_explicit"    # "5_1_explicit"

    chem_data = chem(diffusion_rates=[diff], names=["A"])

    bio = BioSim1D(n_bins=33, chem_data=chem_data)

    initial_concs = np.array([50, 54, 58, 60, 65, 69, 72, 76, 80, 75, 70, 65, 60, 55, 50, 45, 40, 48, 55, 62, 70, 77, 85,
                              93, 100, 102.5, 105, 108, 110, 112, 115, 118, 120])
    '''
    initial_concs = np.array([50.  ,       53.61241666, 56.73695644, 58.95163291, 59.95734176, 59.61825643,
     57.98017227, 55.26432163, 51.83749518, 48.16250482, 44.73567837, 42.01982773,
     40.38174357, 40.04265824, 41.04836709, 43.26304356, 46.38758334])
     '''
    bio.set_species_conc(species_name="A", conc_list=initial_concs)

    bio.double_spacial_resolution_linear()
    bio.double_spacial_resolution_linear()


    #bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
    #bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)
    #bio.inject_sine_conc(species_name="A", frequency=3, amplitude=6)


    status = bio.diffuse(time_step=delta_t, n_steps=60, delta_x=delta_x , algorithm=algorithm)


    bio.describe_state()
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))

    #status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x) # , algorithm="5_1_explicit"
    status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)
    assert status["steps"] == 1
    bio.describe_state()    # [ 51.02083333 , 77.6875 , 43.14583333 , 98.66666667 , 119.5]
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))

    f = bio.lookup_species(species_index=0)
    gradient_x = np.gradient(f, delta_x)
    #print("gradient_x: ", gradient_x)   # [13.33333333, -1.96875   ,  5.24479167, 19.08854167, 10.41666667]
    second_gradient_x = np.gradient(gradient_x, delta_x)
    print("second_gradient_x: ", second_gradient_x) # [-7.65104167, -2.02213542,  5.26432292,  1.29296875, -4.3359375 ]

    # A second diffusion step
    #status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x) # , algorithm="5_1_explicit"
    status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)
    assert status["steps"] == 1
    bio.describe_state()    # [ 51.92612847  75.65907118  45.97243924  97.51072049 118.96462674]
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))

    print(f"\nf_t at x{i}:", f_t)
    print(f"df_t at x{i}:", np.gradient(f_t, delta_t))
    print(f"df_t at x{i}, t1:", np.gradient(f_t, delta_t)[1])
    print(f"D d2/dx2 at x{i}:", diff*second_gradient_x[i])




def test_diffuse_6(biomsim1D):      #TODO: ditch
    from modules.movies.movies import Movie
    delta_t = 0.015
    delta_x = 2
    diff = 10.
    i = 10
    algorithm = None # "5_1_explicit"    # "5_1_explicit"
    history = Movie(tabular=False)
    f_t = []

    chem_data = chem(diffusion_rates=[diff], names=["A"])

    bio = BioSim1D(n_bins=100, chem_data=chem_data)

    bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
    bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)
    #bio.inject_sine_conc(species_name="A", frequency=3, amplitude=6)

    bio.describe_state()
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))
    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

    status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)
    assert status["steps"] == 1
    bio.describe_state()
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))
    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

    f = bio.lookup_species(species_index=0)
    gradient_x = np.gradient(f, delta_x)
    #print("gradient_x: ", gradient_x)   # [13.33333333, -1.96875   ,  5.24479167, 19.08854167, 10.41666667]
    second_gradient_x = np.gradient(gradient_x, delta_x)
    #print("second_gradient_x: ", second_gradient_x) # [-7.65104167, -2.02213542,  5.26432292,  1.29296875, -4.3359375 ]

    # A second diffusion step
    #status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x) # , algorithm="5_1_explicit"
    status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)
    assert status["steps"] == 1
    bio.describe_state()    # [ 51.92612847  75.65907118  45.97243924  97.51072049 118.96462674]
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))
    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

    print(f"\nf_t at x{i}:", f_t)
    print(f"df_t at x{i}:", np.gradient(f_t, delta_t))
    print(f"df_t at (x{i}, t1):", np.gradient(f_t, delta_t)[1])
    print(f"D d2/dx2 at x{i}:", diff*second_gradient_x[i], "\n")

    all_history = history.get()
    #print(all_history)

    H2 = all_history[1]  # The 2nd point in the 3-point series
    f = H2[2]
    #print(f)
    gradient_x = np.gradient(f, delta_x)
    #print("gradient_x: ", gradient_x)   # [13.33333333, -1.96875   ,  5.24479167, 19.08854167, 10.41666667]
    second_gradient_x = np.gradient(gradient_x, delta_x)
    #print("second_gradient_x: ", second_gradient_x) # [-7.65104167, -2.02213542,  5.26432292,  1.29296875, -4.3359375 ]

    for j in range(100):
        time_history_at_bin = np.array([all_history[0][2][j] , all_history[1][2][j], all_history[2][2][j]])
        #print(time_history_at_bin)
        print(f"\ntime_history_at_bin x{j}:", time_history_at_bin)
        print(f"df_t at x{j}:", np.gradient(time_history_at_bin, delta_t))
        lhs = np.gradient(time_history_at_bin, delta_t)[1]
        rhs = diff*second_gradient_x[j]
        print(f"df_t at (x{j}, t1):", lhs)
        print(f"D d2/dx2 at (x{j}, t1):", rhs)
        print("ABS(diff): ", abs(lhs-rhs))
        print(f"ABS(diff) as %: {100 * (abs(lhs - rhs)) / lhs}%")



def test_diffuse_7(biomsim1D):
    from modules.movies.movies import MovieArray
    delta_t = 0.015
    delta_x = 2
    diffusion_rate = 10.
    i = 4
    algorithm = None # "5_1_explicit"    # "5_1_explicit"
    history = MovieArray()
    f_t = []

    chem_data = chem(diffusion_rates=[diffusion_rate], names=["A"])

    bio = BioSim1D(n_bins=150, chem_data=chem_data)

    bio.inject_sine_conc(species_name="A", frequency=1, amplitude=10, bias=50)
    bio.inject_sine_conc(species_name="A", frequency=2, amplitude=8)
    #bio.inject_sine_conc(species_name="A", frequency=3, amplitude=6)

    bio.describe_state()
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))
    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

    status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)
    assert status["steps"] == 1
    bio.describe_state()
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))
    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

    f = bio.lookup_species(species_index=0)
    gradient_x = np.gradient(f, delta_x)
    #print("gradient_x: ", gradient_x)   # [13.33333333, -1.96875   ,  5.24479167, 19.08854167, 10.41666667]
    second_gradient_x = np.gradient(gradient_x, delta_x)
    #print("second_gradient_x: ", second_gradient_x) # [-7.65104167, -2.02213542,  5.26432292,  1.29296875, -4.3359375 ]

    # A second diffusion step
    #status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x) # , algorithm="5_1_explicit"
    status = bio.diffuse(time_step=delta_t, n_steps=1, delta_x=delta_x , algorithm=algorithm)
    assert status["steps"] == 1
    bio.describe_state()    # [ 51.92612847  75.65907118  45.97243924  97.51072049 118.96462674]
    f_t.append(bio.bin_concentration(bin_address=i, species_index=0))
    arr = bio.lookup_species(species_index=0, copy=True)
    history.store(pars=bio.system_time, data_snapshot=arr, caption=f"State at time {bio.system_time}")

    #print(f"\nf_t at x{i}:", f_t)
    #print(f"df_t at x{i}:", np.gradient(f_t, delta_t))
    #print(f"df_t at (x{i}, t1):", np.gradient(f_t, delta_t)[1])
    #print(f"D * d2/dx2 at x{i}:", diffusion_rate*second_gradient_x[i], "\n")

    all_history = history.get_array()
    #print(all_history)

    df_dt_all_bins = np.apply_along_axis(np.gradient, 0, all_history, delta_t)
    #print("df_dt_all_bins: \n", df_dt_all_bins)

    print(f"df_t at t1, across all bins:", df_dt_all_bins[1])    # t1 is the middle point of 3

    f = all_history[1]  # The middle of the 3 time snapshots
    #print("State at Middle time: ", f)
    gradient_x = np.gradient(f, delta_x)
    second_gradient_x = np.gradient(gradient_x, delta_x)
    #print("second_gradient_x at time t1: ", second_gradient_x)
    print("D * second_gradient_x at time t1: ", diffusion_rate*second_gradient_x)

    lhs = df_dt_all_bins[1]
    #print("lhs: ", lhs)

    rhs = diffusion_rate*second_gradient_x
    #print("rhs: ", rhs)

    print("ABS(diff): ", abs(lhs-rhs))
    print(f"ABS(diff) as %: {100 * (abs(lhs - rhs)) / lhs}")

    num.compare_states(lhs, rhs, verbose=True)




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
