# These are tests specifically for diffusion in 1D;
# for general tests of 1D system, see test_biosim_1d.py


import pytest
import numpy as np
from scipy.ndimage import shift
from life123 import BioSim1D
from life123 import ChemData as chem
from life123 import Numerical as num
from life123 import CollectionArray



#########   TESTS OF DIFFUSION : single species, one step    #########


def test__diffuse_step_single_chem_3_1_stencil_A():
    n_bins = 2
    diff = 10.

    chem_data = chem(diffusion_rates=[diff])    # Just 1 chemical species

    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=100.)
    #bio.describe_state()

    increment_vector = np.zeros(n_bins, dtype=float)

    with pytest.raises(Exception):
        #Excessive time step
        bio._diffuse_step_single_chem_3_1_stencil(time_step=0.034, diff=diff, increment_vector=increment_vector)

    # Diffuse by a single step
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector, [-20., 20.])

    # No difference when running it a second time, because the system concentrations weren't updated
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [-20., 20.])



def test__diffuse_step_single_chem_3_1_stencil_B():
    n_bins = 3
    delta_t = 0.01
    delta_x = 2
    diff = 10.

    # Initialize the system
    chem_data = chem(diffusion_rates=[diff])    # Just 1 chemical species
    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)
    initial_concs = np.array([50, 80, 20])
    bio.set_species_conc(chem_index=0, conc_list=initial_concs)
    #bio.describe_state()

    # The coefficients for the 2nd derivative of order 2
    C1 = 1
    C0 = -2
    coeffs = np.array([C1, C0, C1])
    assert np.allclose(np.sum(coeffs), 0)   # The coefficients should add up to zero

    increment_vector = np.zeros(n_bins, dtype=float)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector, chem_index=0, delta_x=delta_x)
    #print("increment_vector: ", increment_vector)

    # MANUALLY computes the expected increments at each bin (i.e. manually do a convolution operation)
    concs = shift(initial_concs, 1, cval=initial_concs[0])      # [50, 50, 80]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, increment_vector[0])

    concs = initial_concs                                       # [50, 80, 20]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, increment_vector[1])

    concs = shift(initial_concs, -1, cval=initial_concs[-1])    # [80, 20, 20]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, increment_vector[2])



def test__diffuse_step_single_chem_3_1_stencil_C():
    n_bins = 2
    diff = 10.

    chem_data = chem(labels="A", diffusion_rates=diff)    # Just 1 chemical species

    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=100.)

    increment_vector = np.zeros(n_bins, dtype=float)

    # Diffuse by a single step (no membranes)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [-20., 20.])

    # Now introduce membranes
    bio.set_membranes([ (0,2) ])    # Membranes at the outer edges of the system
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [-20., 20.])      # Nothing has changes from earlier (w/o membranes)

    # Membranes around the leftmost bin
    bio.set_membranes([ (0,1) ])    # By default, impermeable
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [0., 0.])     # Diffusion is blocked by membrane

    bio.change_permeability("A", 0.5)   # Make the membranes permeable to "A"
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [-1., 1.])     # Passive transport thru membrane, out of leftmost bin

    bio.change_permeability("A", 1.)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [-2., 2.])     # Passive transport thru membrane, out of leftmost bin

    bio.change_permeability("A", diff)     # Permeability identical to diffusion value
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [-20., 20.])  # Passive transport thru membrane, out of leftmost bin

    # Membranes around the rightmost bin
    bio.set_membranes([ (1,2) ])
    bio.change_permeability("A", 0)    # Make impermeable
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [0., 0.])     # Diffusion is blocked by membrane

    bio.change_permeability("A", 0.5)   # Make the membranes permeable to "A"
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [-1., 1.])    # Passive transport thru membrane, out of leftmost bin

    bio.change_permeability("A", 1.)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [-2., 2.])    # Passive transport thru membrane, out of leftmost bin

    bio.change_permeability("A", diff)     # Permeability identical to diffusion value
    bio._diffuse_step_single_chem_3_1_stencil(time_step=0.02, diff=diff, increment_vector=increment_vector)
    assert np.allclose(increment_vector,  [-20., 20.])   # Passive transport thru membrane, out of leftmost bin



def test__diffuse_step_single_chem_3_1_stencil_D():
    n_bins = 3
    delta_t = 0.01
    diff = 5.

    # Initialize the system
    chem_data = chem(labels="A", diffusion_rates=diff)    # Just 1 chemical species
    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

    bio.set_species_conc(chem_index=0, conc_list=[50, 80, 20])

    increment_vector = np.zeros(n_bins, dtype=float)

    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)

    assert np.allclose(increment_vector, [1.5, -4.5, 3.])
    assert np.allclose(np.sum(increment_vector), 0.)


    # Now introduce impermeable membranes at the outer edges of the system
    bio.set_membranes([ (0,3) ])    # By default impermeable

    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.5, -4.5, 3.])      # Identical to membrane-free diffusion


    # Place membranes around the leftmost bin (still impermeable)
    bio.set_membranes([ (0,1) ])

    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0, -3., 3.])  # Diffusion localized to just bins 1 and 2

    bio.change_permeability("A", 1.)   # Make the membranes permeable to "A"
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0.3, -3.3, 3.])  # Diffusion localized to just bins 1 and 2


    # Impermeable membranes around the rightmost bin
    bio.set_membranes([ (2,3) ])
    bio.change_permeability("A", 0)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.5, -1.5, 0.])  # Diffusion localized to just bins 0 and 1

    # Somewhat permeable membranes around the rightmost bin
    bio.change_permeability("A", 2.5)   # Make the membranes permeable to "A"
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.5, -3.0, 1.5])  # Diffusion localized to just bins 0 and 1


    # Impermeable membranes around the middle bin
    bio.set_membranes([ (1,2) ])
    bio.change_permeability("A", 0)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0., 0., 0.])      # Diffusion is blocked by the membranes

    # Permeable membranes around the middle bin, with permeability identical to diffusion rate
    bio.change_permeability("A", diff)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.5, -4.5, 3.])      # Identical to membrane-free diffusion

    # Slightly permeable membranes around the middle bin
    bio.change_permeability("A", diff/3)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0.5, -1.5, 1.])


    # Impermeable membranes around the leftmost and the rightmost bins
    bio.set_membranes([ (0,1), (2,3) ])
    bio.change_permeability("A", 0)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0., 0., 0.])      # Diffusion is blocked by the membranes

    # Permeable membranes around the leftmost and the rightmost bins, with permeability identical to diffusion rate
    bio.change_permeability("A", diff)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.5, -4.5, 3.])      # Identical to membrane-free diffusion

    # Slightly permeable membranes around the leftmost and the rightmost bins
    bio.change_permeability("A", diff/3)
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [.5, -1.5, 1.])



def test__diffuse_step_single_chem_3_1_stencil_E():
    n_bins = 4
    delta_t = 0.01
    diff = 5.

    # Initialize the system
    chem_data = chem(labels="A", diffusion_rates=diff)    # Just 1 chemical species
    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)
    bio.set_species_conc(chem_index=0, conc_list=[50, 100, 20, 60])

    increment_vector = np.zeros(n_bins, dtype=float)

    # Diffusion without membranes
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)

    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])
    assert np.allclose(np.sum(increment_vector), 0.)


    # Now introduce membranes at the outer edges of the system
    bio.set_membranes([ (0,4) ])    # By default impermeable

    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.change_permeability("A", diff/2)    # Make the outer membrane permeable
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion


    # Try out a series of internal membranes with permeability identical to the diffusion rate
    bio.change_permeability("A", diff)     # Permeability identical to diffusion value

    bio.set_membranes([ (0,1) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (1,2) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (2,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (3,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (0,2) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (1,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (2,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (0,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (1, 4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (0,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (0,1), (2,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (0,1), (2,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    bio.set_membranes([ (0,2), (3,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion



    # Try out a series of IMPERMEABLE membranes
    bio.change_permeability("A", 0)

    bio.set_membranes([ (0,1) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0, -4., 6., -2.]) # Diffusion limited to 3 rightmost bins

    bio.set_membranes([ (1,2) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0, 0, 2., -2.])   # Diffusion limited to 2 rightmost bins

    bio.set_membranes([ (2,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -2.5, 0, 0])  # Diffusion limited to 2 leftmost bins

    bio.set_membranes([ (3,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 4., 0])  # Diffusion limited to 3 leftmost bins

    bio.set_membranes([ (0,2) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -2.5, 2., -2.])   # Two sets of 2-bin diffusion

    bio.set_membranes([ (1,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0, -4., 4., 0])      # 2-bin diffusion in middle bins

    bio.set_membranes([ (2,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -2.5, 2., -2.])  # Two sets of 2-bin diffusion

    bio.set_membranes([ (0,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 4., 0])  # Diffusion limited to 3 leftmost bins

    bio.set_membranes([ (1, 4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0, -4., 6., -2.]) # Diffusion limited to 3 rightmost bins

    bio.set_membranes([ (0,1), (2,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [0, 0, 2., -2.])   # Diffusion limited to 2 rightmost bins

    bio.set_membranes([ (0,2), (3,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -2.5, 0, 0])  # Diffusion limited to 2 leftmost bins



    # Try out a series of membranes with a medium amount of permeability
    bio.change_permeability("A", diff/2)

    bio.set_membranes([ (0,1) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.25, -5.25, 6., -2.])    # 2 rightmost bins same as membrane-free scenario

    bio.set_membranes([ (1,2) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.25, -3.25, 4., -2.])

    bio.set_membranes([ (2,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -4.5, 3., -1.])

    bio.set_membranes([ (3,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 5., -1.])

    bio.set_membranes([ (0,2) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -4.5, 4., -2.])

    bio.set_membranes([ (1,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.25, -5.25, 5., -1.])

    bio.set_membranes([ (2,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -4.5, 4., -2.])

    bio.set_membranes([ (0,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 5., -1.])

    bio.set_membranes([ (1, 4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.25, -5.25, 6., -2.])

    bio.set_membranes([ (0,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])     # Identical to membrane-free diffusion

    bio.set_membranes([ (0,1), (2,3) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.25, -3.25, 3., -1.])

    bio.set_membranes([ (0,1), (2,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [1.25, -3.25, 4., -2.])

    bio.set_membranes([ (0,2), (3,4) ])
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(increment_vector, [2.5, -4.5, 3., -1.])



def test__diffuse_step_single_chem_3_1_stencil_F():
    n_bins = 40
    delta_t = 0.015
    diff = 4.2

    # Initialize the system, with a complex pattern of initial concentrations
    chem_data = chem(labels="A", diffusion_rates=diff)    # Just 1 chemical species
    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)
    bio.inject_gradient(chem_label="A", conc_left = 10., conc_right = 1000.)
    bio.inject_sine_conc(chem_label="A", number_cycles=4, amplitude=350, bias=20, phase=60, zero_clip = True)
    bio.inject_bell_curve(chem_label="A", mean=0.8, sd=0.15, amplitude=2.3, bias=18)
    #print(bio.system)

    increment_vector = np.zeros(n_bins, dtype=float)

    # Diffusion without membranes
    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    #print(increment_vector)
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation

    # Now add various complex sets of membranes
    bio.set_membranes([ (0,1), (2,4), (18, 25), (28,29), (38,40) ])

    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(np.sum(increment_vector), 0.)


    bio.set_membranes([ (1,2), (3,8), (38,40) ])

    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(np.sum(increment_vector), 0.)


    bio.set_membranes([ (11,19), (22,28) , (31,34), (39,40) ])

    bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                              increment_vector=increment_vector)
    assert np.allclose(np.sum(increment_vector), 0.)

    for _ in range(5):
        bio._diffuse_step_single_chem_3_1_stencil(time_step=delta_t, diff=diff,
                                                  increment_vector=increment_vector)
        assert np.allclose(np.sum(increment_vector), 0.)



def test__diffuse_step_single_chem_3_1_stencil_G():
    n_bins = 60
    delta_t = 0.025
    diff = [3.2, 2.7]

    # Initialize the system, with a complex pattern of initial concentrations
    chem_data = chem(labels=["A", "B"], diffusion_rates=diff)
    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

    bio.inject_gradient(chem_label="A", conc_left = 20., conc_right = 800.)
    bio.inject_sine_conc(chem_label="A", number_cycles=5, amplitude=450, bias=10, phase=90, zero_clip = True)
    bio.inject_bell_curve(chem_label="A", mean=0.7, sd=0.25, amplitude=13.3, bias=18)

    bio.inject_gradient(chem_label="B", conc_left = 500., conc_right = 18.)
    bio.inject_bell_curve(chem_label="B", mean=0.2, sd=0.45, amplitude=53.3, bias=68)
    bio.inject_sine_conc(chem_label="B", number_cycles=7, amplitude=350, bias=2, phase=123, zero_clip = True)
    #print(bio.system)

    increment_vector = np.zeros(n_bins, dtype=float)

    # Diffusion without membranes
    bio._diffuse_step_single_chem_3_1_stencil(chem_index=0, time_step=delta_t, diff=diff[0],
                                              increment_vector=increment_vector)
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "A"

    bio._diffuse_step_single_chem_3_1_stencil(chem_index=1, time_step=delta_t, diff=diff[1],
                                              increment_vector=increment_vector)
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "B"


    # Now add a complex sets of membranes, initially impermeable
    bio.set_membranes([ (0,1), (2,4) , (11, 19), (20,37), (40,41), (58,60) ])

    bio._diffuse_step_single_chem_3_1_stencil(chem_index=0, time_step=delta_t, diff=diff[0],
                                              increment_vector=increment_vector)    # Diffuse "A"
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "A"

    bio._diffuse_step_single_chem_3_1_stencil(chem_index=1, time_step=delta_t, diff=diff[1],
                                              increment_vector=increment_vector)    # Diffuse "B"
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "B"


    # Make the membranes permeable
    bio.set_permeability({"A": 2.1, "B": 1.4})

    for _ in range(5):
        bio._diffuse_step_single_chem_3_1_stencil(chem_index=0, time_step=delta_t, diff=diff[0],
                                                  increment_vector=increment_vector)    # Diffuse "A"
        assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "A"

        bio._diffuse_step_single_chem_3_1_stencil(chem_index=1, time_step=delta_t, diff=diff[1],
                                                  increment_vector=increment_vector)    # Diffuse "B"
        assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "B"



def test__diffuse_step_single_chem_3_1_stencil_H():
    # Verify that membranes whose permeability is identical to the diffusion rate
    # produce a diffusion identical to what would happens in the absence of membranes
    n_bins = 50
    delta_t = 0.02
    diff = [3.1, 1.7]

    # Initialize the system, with a complex pattern of initial concentrations
    chem_data = chem(labels=["A", "B"], diffusion_rates=diff)
    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

    bio.inject_gradient(chem_label="A", conc_left = 120., conc_right = 800.)
    bio.inject_sine_conc(chem_label="A", number_cycles=3, amplitude=450, bias=20, phase=210, zero_clip = True)
    bio.inject_bell_curve(chem_label="A", mean=0.2, sd=0.5, amplitude=53.3, bias=18)

    bio.inject_gradient(chem_label="B", conc_left = 500., conc_right = 28.)
    bio.inject_bell_curve(chem_label="B", mean=0.4, sd=0.35, amplitude=23.3, bias=68)
    bio.inject_sine_conc(chem_label="B", number_cycles=8, amplitude=350, bias=2, phase=123, zero_clip = True)
    #print(bio.system)
    initial_system = bio.system.copy()

    increment_vector = np.zeros(n_bins, dtype=float)

    # Diffusion without membranes
    bio._diffuse_step_single_chem_3_1_stencil(chem_index=0, time_step=delta_t, diff=diff[0],
                                              increment_vector=increment_vector)    # Diffuse A
    incr_vector_A_no_membranes = increment_vector.copy()

    bio._diffuse_step_single_chem_3_1_stencil(chem_index=1, time_step=delta_t, diff=diff[1],
                                              increment_vector=increment_vector)    # Diffuse B
    incr_vector_B_no_membranes = increment_vector.copy()


    # Now restore the earlier initial state, and add a complex set of membranes,
    # with identical permeabilities to their respective diffusion rates
    bio.set_membranes([ (0,1), (2,4) , (12, 19), (20,37), (40,41), (47,49) ])
    bio.set_permeability({"A":diff[0], "B": diff[1]})

    bio.system = initial_system

    bio._diffuse_step_single_chem_3_1_stencil(chem_index=0, time_step=delta_t, diff=diff[0],
                                              increment_vector=increment_vector)    # Diffuse A
    assert np.allclose(increment_vector, incr_vector_A_no_membranes)

    bio._diffuse_step_single_chem_3_1_stencil(chem_index=1, time_step=delta_t, diff=diff[1],
                                              increment_vector=increment_vector)    # Diffuse B
    assert np.allclose(increment_vector, incr_vector_B_no_membranes)




def test_diffuse_step_single_species_5_1_stencils():
    n_bins=5
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

    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

    initial_concs = np.array([50, 80, 40, 100, 120])
    bio.set_species_conc(chem_index=0, conc_list=initial_concs)

    #bio.describe_state()

    increment_vector = np.zeros(n_bins, dtype=float)
    bio._diffuse_step_single_chem_5_1_stencil(time_step=delta_t, diff=diff,
                                             increment_vector=increment_vector, chem_index=0, delta_x=delta_x)
    #print(increment_vector)

    # Manually computes the expected increments at each bin
    concs = shift(initial_concs, 2, cval=initial_concs[0])      # [50, 50, 50, 80, 40]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, increment_vector[0])

    concs = shift(initial_concs, 1, cval=initial_concs[0])      # [50, 50, 80, 40, 100]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, increment_vector[1])

    concs = initial_concs                                       # [50, 80, 40, 100, 120]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, increment_vector[2])

    concs = shift(initial_concs, -1, cval=initial_concs[-1])    # [80, 40, 100, 120, 120]
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, increment_vector[3])

    concs = shift(initial_concs, -2, cval=initial_concs[-1])    # [40, 100, 120, 120, 120])
    incr = np.dot(concs, coeffs) * diff * delta_t / (delta_x)**2
    assert np.allclose(incr, increment_vector[4])



##############################################################################################

def test_diffuse_step_1():
    # Test with just 1 bin
    chem_data = chem(names="A")
    bio = BioSim1D(n_bins=1, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=2, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=1, chem_index=0, delta_conc=10.)
    bio.describe_state()

    # The time step is so large that the system immediately equilibrates
    print("The default max allowed time step is: ", bio.max_time_step(.5, delta_x=1))
    for i in range(3):
        bio.diffuse_step(time_step=0.6666)
        bio.system += bio.delta_diffusion
        print(f"At end of step {i+1}:")
        print(bio.system)

    assert np.allclose(bio.lookup_species(0), np.full(3, 3.3333333, dtype=float))



def test_diffuse_step_6():
    # Multiple diffusion steps, with 5 bins, and a large time step
    chem_data = chem(diffusion_rates=[.5])
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=10.)
    bio.describe_state()

    print("The default max allowed time step is: ", bio.max_time_step(.5, delta_x=1))
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
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=10.)
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
    bio = BioSim1D(n_bins=15, chem_data=chem_data)

    np.random.seed(18)
    bio.set_species_conc(chem_index=0, conc_list=100 * np.random.rand(15))
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

    #print(f"At end of ALL steps:")
    #print(bio.lookup_species(0))
    # With such a large number of steps, all concentrations will
    # equilibrate to their average
    assert np.allclose(bio.lookup_species(0), np.full(15, avg_conc, dtype=float))




#########   TESTS OF DIFFUSION : all species, one step    #########

def test_diffuse_step_1():
    chem_data = chem(diffusion_rates=[5., 20.])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.set_bin_conc(bin_address=3, chem_index=0, conc=100.) # Bin number out of range

    bio.set_bin_conc(chem_index=0, bin_address=1, conc=100.)

    with pytest.raises(Exception):
        bio.set_species_conc(chem_index=2, conc_list=[10, 20, 50])   # species_index out of range

    bio.set_species_conc(chem_index=1, conc_list=[10, 20, 50])

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
    chem_data = chem(diffusion_rates=[10.])

    bio = BioSim1D(n_bins=2, chem_data=chem_data)

    bio.inject_conc_to_bin(bin_address=0, chem_index=0, delta_conc=100.)
    bio.describe_state()

    # Diffuse by a single step
    bio.diffuse(time_step=0.02, n_steps=1)
    print(bio.system)
    assert np.allclose(bio.lookup_species(0), [80, 20])

    # Another single step
    bio.diffuse(time_step=0.01, n_steps=1)
    print(bio.system)
    assert np.allclose(bio.lookup_species(0), [74, 26])



def test_diffuse_2():
    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=1, chem_data=chem_data)
    bio.set_uniform_concentration(chem_index=0, conc=8.0)

    chem_data.set_diffusion_rate(label="A", diff_rate = 20.)
    #bio.describe_state()    # 1 bins and 1 species:  [[8.]]

    bio.diffuse(time_step=3, n_steps=1)    # With just 1 bin, nothing happens
    #bio.describe_state()
    assert np.allclose(bio.lookup_species(0), [8])



def test_diffuse_3():
    n_bins=5

    chem_data = chem(names=["A"])
    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

    with pytest.raises(Exception):
        bio.diffuse(time_step=0.08, n_steps=1)    # Must first initialize the system


    chem_data = chem(diffusion_rates=[2.78, 3.14])
    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

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
    bio = BioSim1D(n_bins=3, chem_data=chem_data)
    bio.set_bin_conc(chem_index=0, bin_address=1, conc=100.)
    bio.set_species_conc(chem_index=1, conc_list=[10, 20, 50])
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
    bio.set_species_conc(chem_index=0, conc_list=[0, 100, 0])
    bio.set_species_conc(chem_index=1, conc_list=[10, 20, 50])
    bio.describe_state()
    # Re-take the 2 steps
    bio.diffuse(time_step = 0.01, n_steps = 2)
    assert np.allclose(bio.system, [[9.25, 81.5, 9.25] , [14.4, 25.6, 40.]])



def test_diffuse_5():
    chem_data = chem(diffusion_rates=[0.1])
    # Diffuse 1 species almost to equilibrium, starting from a single concentration pulse
    bio = BioSim1D(n_bins=10, chem_data=chem_data)

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

    chem_data = chem(diffusion_rates=[diff])

    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

    initial_concs = np.array([50, 80, 40, 100, 120])
    bio.set_species_conc(chem_index=0, conc_list=initial_concs)

    #bio.describe_state()

    increment_vector = np.zeros(n_bins, dtype=float)

    # Compute the increments to the concentrations, from a single diffusion step
    bio._diffuse_step_single_chem_5_1_stencil(time_step=delta_t, diff=diff,
                                             increment_vector=increment_vector, chem_index=0, delta_x=delta_x)


    # Redo computations on an identical system
    bio2 = BioSim1D(n_bins=n_bins, chem_data=chem_data)
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

    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

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

    bio = BioSim1D(n_bins=n_bins, chem_data=chem_data)

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




#########   TESTS OF MEMBRANES    #########

def test_uses_membranes():
    bio = BioSim1D(n_bins=123, chem_data=chem(labels="A"))
    assert not bio.uses_membranes()

    bio.set_membranes(membranes=[(0, 4)])
    assert bio.uses_membranes()

    bio.set_membranes(membranes=[])
    assert not bio.uses_membranes()



def test_set_membranes():
    bio = BioSim1D(n_bins=40, chem_data=chem(labels="A"))

    bio.set_membranes([ (0, 8) , (17, 31) ])
    assert bio.membranes == [ (0, 8) , (17, 31) ]
    assert bio.permeability == {}

    with pytest.raises(Exception):
        bio.set_membranes(123)

    with pytest.raises(Exception):
        bio.set_membranes([1, 2, 3])

    with pytest.raises(Exception):
        bio.set_membranes([(1, 2, 3)])

    with pytest.raises(Exception):
        bio.set_membranes([(33, -22)])  # Not in sorted order

    with pytest.raises(Exception):
        bio.set_membranes([(-22, 33)])

    with pytest.raises(Exception):
        bio.set_membranes([(0, 66)])

    with pytest.raises(Exception):
        bio.set_membranes([ (3,4), (1,2) ])

    with pytest.raises(Exception):
        bio.set_membranes([ (1,2), (2,4) ])



def test_membranes_list():
    bio = BioSim1D(n_bins=40, chem_data=chem(labels="A"))

    assert bio.membranes_list() == []

    bio.set_membranes([ (8,10) ])
    assert bio.membranes_list() == [8,10]

    bio.set_membranes([ (1,2), (3,4) ])
    assert bio.membranes_list() == [1,2,3,4]



def test_membrane_on_left():
    bio = BioSim1D(n_bins=40, chem_data=chem(labels="A"))

    assert not bio.membrane_on_left(10)

    bio.set_membranes([ (6,10) ])

    assert not bio.membrane_on_left(5)
    assert bio.membrane_on_left(6)
    assert not bio.membrane_on_left(7)

    assert not bio.membrane_on_left(9)
    assert bio.membrane_on_left(10)
    assert not bio.membrane_on_left(11)

    bio.set_membranes([ (6,10), (11,18) ])

    assert not bio.membrane_on_left(9)
    assert bio.membrane_on_left(10)
    assert bio.membrane_on_left(11)
    assert not bio.membrane_on_left(12)
    assert not bio.membrane_on_left(17)
    assert bio.membrane_on_left(18)
    assert not bio.membrane_on_left(19)



def test_membrane_on_right():
    bio = BioSim1D(n_bins=40, chem_data=chem(labels="A"))

    assert not bio.membrane_on_right(10)

    bio.set_membranes([ (6,10) ])

    assert not bio.membrane_on_right(4)
    assert bio.membrane_on_right(5)
    assert not bio.membrane_on_right(6)

    assert not bio.membrane_on_right(8)
    assert bio.membrane_on_right(9)
    assert not bio.membrane_on_right(10)

    bio.set_membranes([ (6,10), (11,18) ])

    assert not bio.membrane_on_right(8)
    assert bio.membrane_on_right(9)
    assert bio.membrane_on_right(10)
    assert not bio.membrane_on_right(11)
    assert not bio.membrane_on_right(16)
    assert bio.membrane_on_right(17)
    assert not bio.membrane_on_right(18)
