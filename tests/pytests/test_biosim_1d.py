import pytest
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
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



def test_system_size():
    bio = BioSim1D(n_bins=8 , chem_data=ChemData(names="A"))
    assert bio.system_size() == 8

    bio = BioSim1D(n_bins=1 , chem_data=ChemData(names=["A", "B"]))
    assert bio.system_size() == 1



def test_reset_system():
    pass

def test_replace_system():
    pass



def test_set_uniform_concentration():
    chem_data = ChemData(names=["A"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.set_uniform_concentration(chem_index=0, conc=0.3)
    expected = np.full(5, 0.3, dtype=float)
    assert np.allclose(bio.lookup_species(0), expected)

    # New test
    bio.reset_system()
    chem_data = ChemData(names=["A", "B", "C"])
    bio = BioSim1D(n_bins=15, chem_data=chem_data)
    bio.set_uniform_concentration(chem_label="A", conc=10.)
    bio.set_uniform_concentration(chem_label="B", conc=11.)
    bio.set_uniform_concentration(chem_label="C", conc=12.)
    #bio.describe_state()
    assert np.allclose(bio.lookup_species(0), np.full(15, 10., dtype=float))
    assert np.allclose(bio.lookup_species(1), np.full(15, 11., dtype=float))
    assert np.allclose(bio.lookup_species(2), np.full(15, 12., dtype=float))



def test_set_all_uniform_concentrations():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.set_all_uniform_concentrations(conc_list=[3., 7.])
    assert np.allclose(bio.lookup_species(chem_label="A"), np.full(5, 3., dtype=float))
    assert np.allclose(bio.lookup_species(chem_label="B"), np.full(5, 7., dtype=float))

    #print(bio.system)



def test_set_bin_conc():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.set_bin_conc(bin_address=2, conc=3.14, chem_label="A")
    assert np.allclose(bio.lookup_species(chem_label="A"), [0, 0, 3.14, 0, 0])

    bio.set_bin_conc(bin_address=2, conc=2.78, chem_index=0)
    assert np.allclose(bio.lookup_species(chem_label="A"), [0, 0, 2.78, 0, 0])
    #print(bio.system)



def test_set_species_conc():
    chem_data = ChemData(names=["A"])
    bio = BioSim1D(n_bins=4, chem_data=chem_data)
    bio.set_species_conc(conc_list=[1., 2., 3., 4.], chem_index=0)
    assert np.allclose(bio.lookup_species(chem_label="A"), [1., 2., 3., 4.])

    bio.set_species_conc(conc_list=[1., 2., 3., 4.], chem_label="A")
    assert np.allclose(bio.lookup_species(chem_index=0), [1., 2., 3., 4.])

    bio.set_species_conc(conc_list=(1., 2., 3., 4.), chem_label="A")      # Tuple instead of list
    assert np.allclose(bio.lookup_species(chem_index=0), [1., 2., 3., 4.])

    bio.set_species_conc(conc_list=np.array([1., 2., 3., 4.]), chem_label="A")
    assert np.allclose(bio.lookup_species(chem_index=0), [1., 2., 3., 4.])

    with pytest.raises(Exception):
        bio.set_species_conc(conc_list=[1., 2., 3., 4.])     # Missing chemical species

    with pytest.raises(Exception):
        bio.set_species_conc(conc_list="Do I look like a list??", chem_label="A")

    with pytest.raises(Exception):
        bio.set_species_conc(conc_list=[1., 2., 3.], chem_label="A")                  # Wrong size

    with pytest.raises(Exception):
        bio.set_species_conc(conc_list=np.array([1., 2.]), chem_label="A")            # Wrong size

    with pytest.raises(Exception):
        bio.set_species_conc(conc_list=[1., -0.09, 3., 4.], chem_label="A")           # Negative concentrations

    with pytest.raises(Exception):
        bio.set_species_conc(conc_list=np.array([1., -0.08, 3., 4.]), chem_label="A") # Negative concentrations



def test_inject_conc_to_bin():
    chem_data = ChemData(names=["A"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    with pytest.raises(Exception):
        # cell_index out of bounds
        bio.inject_conc_to_bin(bin_address=5, chem_index=0, delta_conc=10.)

    with pytest.raises(Exception):
        # species_index out of bounds
        bio.inject_conc_to_bin(bin_address=5, chem_index=1, delta_conc=10.)

    bio.inject_conc_to_bin(bin_address=1, chem_index=0, delta_conc=10.)

    bio.describe_state()



def test_inject_gradient():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=8, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.inject_gradient("B", conc_left = 10., conc_right = 20.)     # Non-existent chemical species

        bio.inject_gradient("A", conc_left = -10., conc_right = 20.)    # Negative concentration
        bio.inject_gradient("A", conc_left = 10., conc_right = -20.)    # Negative concentration

    bio.inject_gradient("A", conc_left = 11., conc_right = 18.)
    assert np.allclose(bio.lookup_species(chem_label="A"), [11., 12., 13., 14., 15., 16., 17., 18.])

    bio.inject_gradient("A", conc_left = 18., conc_right = 11.)
    assert np.allclose(bio.lookup_species(chem_label="A"), [29., 29., 29., 29., 29., 29., 29., 29.])

    bio = BioSim1D(n_bins=1, chem_data=chem_data)
    with pytest.raises(Exception):
        bio.inject_gradient("A", conc_left = 11., conc_right = 18.)     # Too few bins in system

    #print(bio.lookup_species(species_name="A"))
    # Curves could be visualized as follows:
    #import plotly.express as px
    #fig = px.line(y=bio.lookup_species(species_name="A"))
    #fig.show()



def test_inject_sine_conc():
    chem_data = ChemData(names=["A"])
    bio = BioSim1D(n_bins=8, chem_data=chem_data)

    bio.inject_sine_conc(chem_label="A", amplitude=5, bias=20, number_cycles=1, phase=0)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [20., 23.53553391, 25., 23.53553391, 20., 16.46446609, 15., 16.46446609])

    # Adding a 2nd wave in opposite phase will erase everything except for the biases
    bio.inject_sine_conc(chem_label="A", amplitude=5, bias=10, number_cycles=1, phase=180)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [30., 30., 30., 30., 30., 30., 30., 30.])


    bio.set_uniform_concentration(conc=0, chem_label="A")     # RESET concentrations

    bio.inject_sine_conc(chem_label="A", amplitude=5, bias=20, number_cycles=1, phase=90)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [15., 16.46446609, 20., 23.53553391, 25.,  23.53553391, 20., 16.46446609])

    # Adding a 2nd wave in opposite phase will erase everything except for the biases
    bio.inject_sine_conc(chem_label="A", amplitude=5, bias=10, number_cycles=1, phase=-90)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [30., 30., 30., 30., 30., 30., 30., 30.])


    bio.set_uniform_concentration(conc=0, chem_label="A")     # RESET concentrations

    with pytest.raises(Exception):
        bio.inject_sine_conc(chem_label="A", amplitude=10, bias=0, number_cycles=1 / 2, phase=180)  # Would lead to negative values


    bio.set_uniform_concentration(conc=0, chem_label="A")     # RESET concentrations

    bio.inject_sine_conc(chem_label="A", amplitude=10, bias=0, number_cycles=1 / 2, phase=0)

    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [0., 3.82683432, 7.07106781, 9.23879533, 10., 9.23879533, 7.07106781, 3.82683432])

    # Adding a 2nd wave in opposite phase will erase everything except for the biases
    bio.inject_sine_conc(chem_label="A", amplitude=10, bias=30, number_cycles=1 / 2, phase=180)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [30., 30., 30., 30., 30., 30., 30., 30.])


    bio.set_uniform_concentration(conc=0, chem_label="A")     # RESET concentrations
    bio.inject_sine_conc(chem_label="A", amplitude=10, bias=30, number_cycles=2, phase=0)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [30., 40., 30., 20., 30., 40., 30., 20.])

    bio.set_uniform_concentration(conc=0, chem_label="A")     # RESET concentrations
    bio.inject_sine_conc(chem_label="A", amplitude=10, bias=30, number_cycles=2, phase=-90)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [40., 30., 20., 30., 40., 30., 20., 30.])

    bio.set_uniform_concentration(conc=0, chem_label="A")     # RESET concentrations
    bio.inject_sine_conc(chem_label="A", amplitude=10, bias=30, number_cycles=2, phase=90)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [20., 30., 40., 30., 20., 30., 40., 30.])

    bio.set_uniform_concentration(conc=0, chem_label="A")     # RESET concentrations
    bio.inject_sine_conc(chem_label="A", amplitude=10, bias=0, number_cycles=2, phase=0, zero_clip=True)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [0., 10., 0., 0., 0., 10., 0., 0.])

    #print(bio.lookup_species(species_name="A"))

    # Curves could be visualized as follows:
    #import plotly.express as px
    #fig = px.line(y=bio.lookup_species(species_name="A"))
    #fig.show()



def test_inject_bell_curve():
    chem_data = ChemData(names=["A"])

    bio = BioSim1D(n_bins=10, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.inject_bell_curve(chem_label="A", bias=-100)      # Negative bias
        bio.inject_bell_curve(chem_label="A", amplitude=-1)   # Negative amplitude

    bio.inject_bell_curve(chem_label="A", mean=0.5, sd=0.15, amplitude=1., bias=0)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [0.01028186, 0.09231148, 0.47878856, 1.43461641, 2.48331496, 2.48331496,
                        1.43461641, 0.47878856, 0.09231148, 0.01028186])

    bio = BioSim1D(n_bins=10, chem_data=chem_data)
    bio.inject_bell_curve(chem_label="A", mean=0.5, sd=0.05, amplitude=2., bias=10)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [10., 10., 10.00000317, 10.06169116, 18.60769844, 18.60769844,
                         10.06169116, 10.00000317, 10., 10.])

    bio = BioSim1D(n_bins=20, chem_data=chem_data)
    bio.inject_bell_curve(chem_label="A", mean=0.2, sd=0.05, amplitude=2., bias=10)
    assert np.allclose(bio.lookup_species(chem_label="A"),
                       [10.00535321, 10.20730806, 12.65097387, 21.19391452, 25.60794776, 17.18616024,
                        11.09253477, 10.05484802, 10.00090923, 10.00000498, 10.00000001, 10.,
                        10., 10.,  10.,  10.,  10.,  10., 10.,  10.])

    # print(bio.lookup_species(species_name="A"))
    # Curves could be visualized as follows:
    #import plotly.express as px
    #fig = px.line(y=bio.lookup_species(species_name="A"))
    #fig.show()



def test_frequency_analysis():
    chem_data = ChemData(names=["A"])
    bio = BioSim1D(n_bins=100, chem_data=chem_data)

    bio.inject_sine_conc(chem_label="A", number_cycles=2, amplitude=1, bias=3)

    result = bio.frequency_analysis(chem_label="A")

    expected = pd.DataFrame({"Frequency": [0,2], "Relative Amplitude": [3.0,1.0]})
    assert_frame_equal(result, expected, check_dtype=False)

    bio.inject_sine_conc(chem_label="A", number_cycles=4, amplitude=0.5)
    bio.inject_sine_conc(chem_label="A", number_cycles=8, amplitude=0.2)

    #import plotly.express as px
    #fig = px.line(y=bio.lookup_species(species_name="A"))
    #fig.show()

    result = bio.frequency_analysis(chem_label="A")
    #print(result)
    expected = pd.DataFrame({"Frequency": [0,2,4,8], "Relative Amplitude": [3.0,1.0,0.5,0.2]})
    assert_frame_equal(result, expected, check_dtype=False)

    result = bio.frequency_analysis(chem_label="A", n_largest=3)  # ditch the smallest amplitude
    #print(result)
    expected = pd.DataFrame({"Frequency": [0,2,4], "Relative Amplitude": [3.0,1.0,0.5]})
    assert_frame_equal(result, expected, check_dtype=False)




########  DIMENSION-RELATED  ################

def test_set_dimensions():
    chem_data = ChemData(names=["A"])
    bio = BioSim1D(n_bins=4, chem_data=chem_data)
    bio.set_dimensions(21.)
    assert np.allclose(bio.system_length, 21.)
    assert np.allclose(bio.global_Dx, 7.)

    with pytest.raises(Exception):
        bio.set_dimensions("I'm no number")
        bio.set_dimensions(0)


def test_x_coord():
    chem_data = ChemData(names=["A"])
    bio = BioSim1D(n_bins=4, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.x_coord(0)      # must first call set_dimensions()

    bio.set_dimensions(21.)
    assert np.allclose(bio.x_coord(0), 0.)
    assert np.allclose(bio.x_coord(1), 7.)
    assert np.allclose(bio.x_coord(3), 21.)




########  MEMBRANE-RELATED  ################

def test_uses_membranes():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    assert not bio.uses_membranes()
    bio.set_membranes(membranes=[(2, 4)])
    assert bio.uses_membranes()



def test_assert_valid_bin():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
    with pytest.raises(Exception):
        bio.assert_valid_bin(-1)

    with pytest.raises(Exception):
        bio.assert_valid_bin(3)

    with pytest.raises(Exception):
        bio.assert_valid_bin("2")

    with pytest.raises(Exception):
        bio.assert_valid_bin((0, 2))

    bio.assert_valid_bin(0)
    bio.assert_valid_bin(2)



def test_lookup_species():
    pass    # TODO


def test_bin_concentration():
    pass    # TODO


def test_bin_snapshot():
    pass    # TODO


def test_system_snapshot():
    pass    # TODO





##################  MEMBRANES  ##################


def test_set_membranes():
    bio = BioSim1D(n_bins=10, chem_data=ChemData("A")) # 10 bins, numbered 0 thru 9

    with pytest.raises(Exception):
        bio.set_membranes()     # Missing required arg

    with pytest.raises(Exception):
        bio.set_membranes(membranes=123)     # Wrong arg type

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[])      # Empty arg

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[8])    # Element isn't a pair

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(1,2,3)])  # Element isn't a pair

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(1,"x")])  # Element isn't a pair of integers

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(8,3)])    # Value not in sorted order

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(-1,3)])   # Bad left value

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(2,11)])   # Bad right value

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(2,4) , (4,9)])    # Touching

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(2,7) , (3,9)])    # Overlapping

    bio.set_membranes(membranes=[(2,4) , (5,9)])
    assert bio.membranes == [(2,4) , (5,9)]
    assert bio.permeability is None

    with pytest.raises(Exception):
        bio.set_membranes(membranes=[(5,9), (2,4)])     # Not in sorted order

    bio.set_membranes(membranes=[(0,3) , (4,6), (8,10)] , permeability=123)
    assert bio.membranes == [(0,3) , (4,6), (8,10)]
    assert bio.permeability == 123






##################  CHANGE RESOLUTIONS  ##################

def test_increase_spacial_resolution():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(chem_index=0, conc_list=[11., 12., 13.])
    bio.set_species_conc(chem_index=1, conc_list=[5., 15., 25.])

    bio.increase_spatial_resolution(2)

    assert np.allclose(bio.lookup_species(chem_index=0), [11., 11., 12., 12., 13., 13.])
    assert np.allclose(bio.lookup_species(chem_label="B"), [5., 5., 15., 15., 25., 25.])



def test_double_spacial_resolution_linear_inter():
    chem_data = ChemData(names=["A"])
    bio = BioSim1D(n_bins=2, chem_data=chem_data)
    bio.set_species_conc(chem_index=0, conc_list=[1., 2.])
    bio.double_spatial_resolution_linear()
    assert np.allclose(bio.system, [[1., 1.5, 2.]])

    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(chem_index=0, conc_list=[11., 12., 13.])
    bio.set_species_conc(chem_index=1, conc_list=[5., 15., 25.])
    bio.double_spatial_resolution_linear()
    assert np.allclose(bio.system, [[11.,11.5,12.,12.5,13.],
                                    [ 5.,10. ,15.,20. ,25.]])

    bio = BioSim1D(n_bins=1, chem_data=chem_data)
    bio.set_species_conc(chem_index=0, conc_list=[3.14])
    with pytest.raises(Exception):
        bio.double_spatial_resolution_linear()



def test_decrease_spacial_resolution():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=6, chem_data=chem_data)
    bio.set_species_conc(chem_index=0, conc_list=[10., 20., 30., 40., 50., 60.])
    bio.set_species_conc(chem_index=1, conc_list=[2., 8., 5., 15., 4., 2.])
    original_state = bio.system

    # Group by pairs
    bio.decrease_spatial_resolution(2)
    assert np.allclose(bio.system, [[15., 35., 55.],
                                    [ 5., 10.,  3.]])

    # Group by triplets
    bio.replace_system(original_state)
    bio.decrease_spatial_resolution(3)
    assert np.allclose(bio.system, [[20., 50.],
                                    [ 5., 7.]])

    # Group into just 1 single bin
    bio.replace_system(original_state)
    bio.decrease_spatial_resolution(6)
    assert np.allclose(bio.system, [[35.],
                                    [ 6.]])



def test_varying_spacial_resolution():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(chem_label="A", conc_list=[11., 12., 13.])
    bio.set_species_conc(chem_label="B", conc_list=[5., 15., 25.])

    # First, increase the spacial resolution, after saving the original matrix
    original_state = bio.system
    bio.increase_spatial_resolution(3)

    # Then, decrease the spacial resolution, by the same factor
    bio.decrease_spatial_resolution(3)
    assert np.allclose(bio.system, original_state)



def test_smooth_spatial_resolution():
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
    bio.set_species_conc(chem_label="A", conc_list=[10., 20., 30.])
    bio.set_species_conc(chem_label="B", conc_list=[2., 8., 4.])

    bio.smooth_spatial_resolution()
    expected = np.array([[10., 15., 20., 25., 30.],
                         [ 2.,  5.,  8.,  6.,  4.]])
    assert np.allclose(bio.system, expected)
    assert bio.n_bins == 5
    assert bio.n_species == 2




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
