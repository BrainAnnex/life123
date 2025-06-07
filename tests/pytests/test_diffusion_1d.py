# These are tests specifically for the class "Diffusion1D";
# for general tests of diffusion in 1D system, see "test_biosim_1d.py"


import pytest
import numpy as np
from scipy.ndimage import shift
from life123 import ChemData, Diffusion1D, Membranes1D, System1D




def test_max_time_step():
    bio = Diffusion1D(n_bins=3)
    result = bio.max_time_step(diff_rate=0.5, delta_x=1)
    assert np.allclose(result, 2/3)



def test_diffuse_step_3_1_stencil_A():
    n_bins = 2
    diff = 10.

    bio = Diffusion1D(n_bins=n_bins)

    initial_concs = np.array([100, 0])      # 2 bins

    with pytest.raises(Exception):
        #Excessive time step
        bio.diffuse_step_3_1_stencil(time_step=0.034, diff=diff, conc_array=initial_concs)
       
    # Diffuse by a single step
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs)
    assert np.allclose(increment_vector, [-20., 20.])



def test_diffuse_step_3_1_stencil_B():
    n_bins = 3
    delta_t = 0.01
    delta_x = 2
    diff = 10.

    # Initialize the system
    bio = Diffusion1D(n_bins=n_bins)
    initial_concs = np.array([50, 80, 20])

    # The coefficients for the 2nd derivative of order 2
    C1 = 1
    C0 = -2
    coeffs = np.array([C1, C0, C1])
    assert np.allclose(np.sum(coeffs), 0)   # The coefficients should add up to zero


    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                   conc_array=initial_concs, delta_x=delta_x)
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



def test_diffuse_step_3_1_stencil_C():
    n_bins = 2
    diff = 10.

    membranes = Membranes1D(n_bins=n_bins)
    bio = Diffusion1D(n_bins=n_bins, membranes=membranes)
    
    initial_concs = np.array([100, 0])      # 2 bins

    # Diffuse by a single step (no membranes)
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs)
    assert np.allclose(increment_vector,  [-20., 20.])

    # Now introduce membranes at the outer edges of the system
    membranes.set_membranes([ (0,2) ])    # By default, impermeable

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs)
    assert np.allclose(increment_vector,  [-20., 20.])      # Nothing has changes from earlier (w/o membranes)

    # Membranes around the leftmost bin
    membranes.set_membranes([ (0,1) ])    # By default, impermeable

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs)
    assert np.allclose(increment_vector,  [0., 0.])     # Diffusion is blocked by membrane

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs,
                                                    permeability=0.5)  # Make the membranes permeable to "A"
    assert np.allclose(increment_vector,  [-1., 1.])     # Passive transport thru membrane, out of leftmost bin

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs,
                                                    permeability=1.)
    assert np.allclose(increment_vector,  [-2., 2.])     # Passive transport thru membrane, out of leftmost bin

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs,
                                                    permeability=diff)     # Permeability identical to diffusion value
    assert np.allclose(increment_vector,  [-20., 20.])  # Passive transport thru membrane, out of leftmost bin


    # Membranes around the rightmost bin
    membranes.set_membranes([ (1,2) ])

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs,
                                                    permeability=0)  # Make impermeable
    assert np.allclose(increment_vector,  [0., 0.])     # Diffusion is blocked by membrane

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs,
                                                    permeability=0.5)   # Make the membranes permeable to "A"
    assert np.allclose(increment_vector,  [-1., 1.])    # Passive transport thru membrane, out of leftmost bin

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs,
                                                    permeability=1.0)
    assert np.allclose(increment_vector,  [-2., 2.])    # Passive transport thru membrane, out of leftmost bin

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=0.02, diff=diff, conc_array=initial_concs,
                                                    permeability=diff)     # Permeability identical to diffusion value
    assert np.allclose(increment_vector,  [-20., 20.])   # Passive transport thru membrane, out of leftmost bin



def test_diffuse_step_3_1_stencil_D():
    n_bins = 3
    delta_t = 0.01
    diff = 5.

    membranes = Membranes1D(n_bins=n_bins)
    bio = Diffusion1D(n_bins=n_bins, membranes=membranes)

    initial_concs = np.array([50, 80, 20])      # 3 bins

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)

    assert np.allclose(increment_vector, [1.5, -4.5, 3.])
    assert np.allclose(np.sum(increment_vector), 0.)


    # Now introduce impermeable membranes at the outer edges of the system
    membranes.set_membranes([ (0,3) ])    # By default impermeable

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [1.5, -4.5, 3.])      # Identical to membrane-free diffusion


    # Place membranes around the leftmost bin (still impermeable)
    membranes.set_membranes([ (0,1) ])

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                   conc_array=initial_concs)
    assert np.allclose(increment_vector, [0, -3., 3.])  # Diffusion localized to just bins 1 and 2

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=1.)   # Make the membranes permeable to "A"
    assert np.allclose(increment_vector, [0.3, -3.3, 3.])  # Diffusion localized to just bins 1 and 2


    # Impermeable membranes around the rightmost bin
    membranes.set_membranes([ (2,3) ])

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [1.5, -1.5, 0.])  # Diffusion localized to just bins 0 and 1

    # Somewhat permeable membranes around the rightmost bin
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs, permeability=2.5)  # Make the membranes permeable
    assert np.allclose(increment_vector, [1.5, -3.0, 1.5])  # Diffusion localized to just bins 0 and 1


    # Impermeable membranes around the middle bin
    membranes.set_membranes([ (1,2) ])

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [0., 0., 0.])      # Diffusion is blocked by the membranes

    # Permeable membranes around the middle bin, with permeability identical to diffusion rate
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)
    assert np.allclose(increment_vector, [1.5, -4.5, 3.])      # Identical to membrane-free diffusion

    # Slightly permeable membranes around the middle bin
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/3)
    assert np.allclose(increment_vector, [0.5, -1.5, 1.])


    # Impermeable membranes around the leftmost and the rightmost bins
    membranes.set_membranes([ (0,1), (2,3) ])

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs)
    assert np.allclose(increment_vector, [0., 0., 0.])      # Diffusion is blocked by the membranes

    # Permeable membranes around the leftmost and the rightmost bins, with permeability identical to diffusion rate
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)
    assert np.allclose(increment_vector, [1.5, -4.5, 3.])      # Identical to membrane-free diffusion

    # Slightly permeable membranes around the leftmost and the rightmost bins
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/3)
    assert np.allclose(increment_vector, [.5, -1.5, 1.])



def test_diffuse_step_3_1_stencil_E():
    n_bins = 4
    delta_t = 0.01
    diff = 5.

    membranes = Membranes1D(n_bins=n_bins)
    bio = Diffusion1D(n_bins=n_bins, membranes=membranes)
    initial_concs = np.array([50, 100, 20, 60])


    # Diffusion without membranes
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)

    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])
    assert np.allclose(np.sum(increment_vector), 0.)


    # Now introduce membranes at the outer edges of the system
    membranes.set_membranes([ (0,4) ])    # By default impermeable

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)    # Make the outer membrane permeable
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion


    # Try out a series of internal membranes with permeability identical to the diffusion rate
    membranes.set_membranes([ (0,1) ])

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (1,2) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (2,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (3,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (0,2) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (1,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (2,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (0,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (1, 4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (0,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (0,1), (2,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (0,1), (2,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion

    membranes.set_membranes([ (0,2), (3,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff)    # Permeability identical to diffusion value)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])      # Identical to membrane-free diffusion


    # Try out a series of IMPERMEABLE membranes
    membranes.set_membranes([ (0,1) ])

    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [0, -4., 6., -2.]) # Diffusion limited to 3 rightmost bins

    membranes.set_membranes([ (1,2) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [0, 0, 2., -2.])   # Diffusion limited to 2 rightmost bins

    membranes.set_membranes([ (2,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [2.5, -2.5, 0, 0])  # Diffusion limited to 2 leftmost bins

    membranes.set_membranes([ (3,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [2.5, -6.5, 4., 0])  # Diffusion limited to 3 leftmost bins

    membranes.set_membranes([ (0,2) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [2.5, -2.5, 2., -2.])   # Two sets of 2-bin diffusion

    membranes.set_membranes([ (1,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [0, -4., 4., 0])      # 2-bin diffusion in middle bins

    membranes.set_membranes([ (2,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [2.5, -2.5, 2., -2.])  # Two sets of 2-bin diffusion

    membranes.set_membranes([ (0,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [2.5, -6.5, 4., 0])  # Diffusion limited to 3 leftmost bins

    membranes.set_membranes([ (1, 4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [0, -4., 6., -2.]) # Diffusion limited to 3 rightmost bins

    membranes.set_membranes([ (0,1), (2,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [0, 0, 2., -2.])   # Diffusion limited to 2 rightmost bins

    membranes.set_membranes([ (0,2), (3,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                    conc_array=initial_concs)
    assert np.allclose(increment_vector, [2.5, -2.5, 0, 0])  # Diffusion limited to 2 leftmost bins



    # Try out a series of membranes with a medium amount of permeability
    membranes.set_membranes([ (0,1) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [1.25, -5.25, 6., -2.])    # 2 rightmost bins same as membrane-free scenario

    membranes.set_membranes([ (1,2) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [1.25, -3.25, 4., -2.])

    membranes.set_membranes([ (2,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [2.5, -4.5, 3., -1.])

    membranes.set_membranes([ (3,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [2.5, -6.5, 5., -1.])

    membranes.set_membranes([ (0,2) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [2.5, -4.5, 4., -2.])

    membranes.set_membranes([ (1,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [1.25, -5.25, 5., -1.])

    membranes.set_membranes([ (2,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [2.5, -4.5, 4., -2.])

    membranes.set_membranes([ (0,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [2.5, -6.5, 5., -1.])

    membranes.set_membranes([ (1, 4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [1.25, -5.25, 6., -2.])

    membranes.set_membranes([ (0,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [2.5, -6.5, 6., -2.])     # Identical to membrane-free diffusion

    membranes.set_membranes([ (0,1), (2,3) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [1.25, -3.25, 3., -1.])

    membranes.set_membranes([ (0,1), (2,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [1.25, -3.25, 4., -2.])

    membranes.set_membranes([ (0,2), (3,4) ])
    increment_vector = bio.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs, permeability=diff/2)
    assert np.allclose(increment_vector, [2.5, -4.5, 3., -1.])



def test_diffuse_step_3_1_stencil_F():
    n_bins = 40
    delta_t = 0.015
    diff = 4.2

    # Initialize the system, with a complex pattern of initial concentrations
    membranes = Membranes1D(n_bins=n_bins)
    diff_obj = Diffusion1D(n_bins=n_bins, membranes=membranes)
    bio = System1D(n_bins=n_bins, chem_data=ChemData("A"))
    
    bio.inject_gradient(chem_label="A", conc_left = 10., conc_right = 1000.)
    bio.inject_sine_conc(chem_label="A", number_cycles=4, amplitude=350, bias=20, phase=60, zero_clip = True)
    bio.inject_bell_curve(chem_label="A", mean=0.8, sd=0.15, amplitude=2.3, bias=18)
    #print(bio.system)
    initial_concs = bio.lookup_species(chem_index=0)

    # Diffusion without membranes
    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                         conc_array=initial_concs)
    #print(increment_vector)
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation


    # Now add various complex sets of membranes (all impermeable)
    membranes.set_membranes([ (0,1), (2,4), (18, 25), (28,29), (38,40) ])

    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs)
    assert np.allclose(np.sum(increment_vector), 0.)


    membranes.set_membranes([ (1,2), (3,8), (38,40) ])

    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs)
    assert np.allclose(np.sum(increment_vector), 0.)


    membranes.set_membranes([ (11,19), (22,28) , (31,34), (39,40) ])

    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                              conc_array=initial_concs)
    assert np.allclose(np.sum(increment_vector), 0.)

    for _ in range(5):
        increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff,
                                                  conc_array=initial_concs)
        assert np.allclose(np.sum(increment_vector), 0.)



def test_diffuse_step_3_1_stencil_G():
    n_bins = 60
    delta_t = 0.025
    diff = [3.2, 2.7]

    # Initialize the system, with a complex pattern of initial concentrations
    membranes = Membranes1D(n_bins=n_bins)
    diff_obj = Diffusion1D(n_bins=n_bins, membranes=membranes)
    bio = System1D(n_bins=n_bins, chem_data=ChemData(["A", "B"]))

    bio.inject_gradient(chem_label="A", conc_left = 20., conc_right = 800.)
    bio.inject_sine_conc(chem_label="A", number_cycles=5, amplitude=450, bias=10, phase=90, zero_clip = True)
    bio.inject_bell_curve(chem_label="A", mean=0.7, sd=0.25, amplitude=13.3, bias=18)
    initial_concs_A = bio.lookup_species(chem_index=0)

    bio.inject_gradient(chem_label="B", conc_left = 500., conc_right = 18.)
    bio.inject_bell_curve(chem_label="B", mean=0.2, sd=0.45, amplitude=53.3, bias=68)
    bio.inject_sine_conc(chem_label="B", number_cycles=7, amplitude=350, bias=2, phase=123, zero_clip = True)
    #print(bio.system)
    initial_concs_B = bio.lookup_species(chem_index=1)

    # Diffusion without membranes
    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[0],
                                                         conc_array=initial_concs_A)
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "A"

    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[1],
                                                    conc_array=initial_concs_B)
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "B"


    # Now add a complex sets of membranes, initially impermeable
    membranes.set_membranes([ (0,1), (2,4) , (11, 19), (20,37), (40,41), (58,60) ])

    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[0],
                                              conc_array=initial_concs_A)    # Diffuse "A"
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "A"

    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[1],
                                              conc_array=initial_concs_B)    # Diffuse "B"
    assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "B"


    # Make the membranes permeable

    for _ in range(5):
        increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[0], permeability=2.1,
                                                  conc_array=initial_concs_A)    # Diffuse "A"
        assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "A"

        increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[1], permeability=1.4,
                                                  conc_array=initial_concs_B)    # Diffuse "B"
        assert np.allclose(np.sum(increment_vector), 0.)    # Check mass conservation for "B"



def test_diffuse_step_3_1_stencil_H():
    # Verify that membranes whose permeability is identical to the diffusion rate
    # produce a diffusion identical to what would happens in the absence of membranes
    n_bins = 50
    delta_t = 0.02
    diff = [3.1, 1.7]

    # Initialize the system, with a complex pattern of initial concentrations
    membranes = Membranes1D(n_bins=n_bins)
    diff_obj = Diffusion1D(n_bins=n_bins, membranes=membranes)
    bio = System1D(n_bins=n_bins, chem_data=ChemData(["A", "B"]))

    bio.inject_gradient(chem_label="A", conc_left = 120., conc_right = 800.)
    bio.inject_sine_conc(chem_label="A", number_cycles=3, amplitude=450, bias=20, phase=210, zero_clip = True)
    bio.inject_bell_curve(chem_label="A", mean=0.2, sd=0.5, amplitude=53.3, bias=18)

    bio.inject_gradient(chem_label="B", conc_left = 500., conc_right = 28.)
    bio.inject_bell_curve(chem_label="B", mean=0.4, sd=0.35, amplitude=23.3, bias=68)
    bio.inject_sine_conc(chem_label="B", number_cycles=8, amplitude=350, bias=2, phase=123, zero_clip = True)
    #print(bio.system)


    # Diffusion without membranes
    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[0],
                                              conc_array=bio.lookup_species(chem_index=0))    # Diffuse A
    incr_vector_A_no_membranes = increment_vector.copy()

    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[1],
                                              conc_array=bio.lookup_species(chem_index=1))    # Diffuse B
    incr_vector_B_no_membranes = increment_vector.copy()


    # Now restore the earlier initial state, and add a complex set of membranes,
    # with identical permeabilities to their respective diffusion rates
    membranes.set_membranes([ (0,1), (2,4) , (12, 19), (20,37), (40,41), (47,49) ])

    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[0], permeability=diff[0],
                                              conc_array=bio.lookup_species(chem_index=0))    # Diffuse A
    assert np.allclose(increment_vector, incr_vector_A_no_membranes)

    increment_vector = diff_obj.diffuse_step_3_1_stencil(time_step=delta_t, diff=diff[1], permeability=diff[1],
                                              conc_array=bio.lookup_species(chem_index=1))    # Diffuse B
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

    bio = Diffusion1D(n_bins=n_bins)

    initial_concs = np.array([50, 80, 40, 100, 120])

    increment_vector = bio.diffuse_step_5_1_stencil(time_step=delta_t, diff=diff,
                                                   conc_array=initial_concs, delta_x=delta_x)
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
