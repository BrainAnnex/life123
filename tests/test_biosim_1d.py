import numpy as np
import pytest
from life_1D.bio_sim_1d import BioSim1D as bio


def test_initialize_universe():
    bio.initialize_universe(n_cells=10, n_species=1)

    assert bio.n_cells == 10
    assert bio.n_species == 1

    bio.describe_state()

    expected = np.zeros((1, 10), dtype=float)
    assert np.allclose(bio.univ, expected)

    # New test
    bio.initialize_universe(n_cells=15, n_species=3)
    bio.describe_state()
    expected = np.zeros((3,15), dtype=float)
    assert np.allclose(bio.univ, expected)



def test_set_diffusion_rates():
    bio.initialize_universe(n_cells=5, n_species=1)
    bio.set_diffusion_rates([.2])
    print("diffusion_rates: ", bio.diffusion_rates)
    assert np.allclose(bio.diffusion_rates, [.2])

    # New test
    bio.initialize_universe(n_cells=25, n_species=4)
    bio.set_diffusion_rates([.1, .7, .2, .4])
    print("diffusion_rates: ", bio.diffusion_rates)
    assert np.allclose(bio.diffusion_rates, [.1, .7, .2, .4])



def test_set_uniform_concentration():
    bio.initialize_universe(n_cells=8, n_species=1)
    bio.set_uniform_concentration(species_index=0, conc=0.3)
    bio.describe_state()
    expected = np.full(8, 0.3, dtype=float)
    assert np.allclose(bio.univ[0], expected)

    # New test
    bio.initialize_universe(n_cells=15, n_species=3)
    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=11.)
    bio.set_uniform_concentration(species_index=2, conc=12.)
    bio.describe_state()
    assert np.allclose(bio.univ[0], np.full(15, 10., dtype=float))
    assert np.allclose(bio.univ[1], np.full(15, 11., dtype=float))
    assert np.allclose(bio.univ[2], np.full(15, 12., dtype=float))



def test_inject_conc_to_cell():
    bio.initialize_universe(n_cells=5, n_species=1)

    with pytest.raises(Exception):
        # cell_index out of bounds
        bio.inject_conc_to_cell(bin=5, species_index=0, delta_conc=10.)

    with pytest.raises(Exception):
        # species_index out of bounds
        bio.inject_conc_to_cell(bin=5, species_index=1, delta_conc=10.)

    bio.inject_conc_to_cell(bin=1, species_index=0, delta_conc=10.)

    bio.describe_state()



#########   TESTS OF DIFFUSION    #########

def test_diffuse_step_single_species_1():
    bio.initialize_universe(n_cells=2, n_species=1)

    bio.inject_conc_to_cell(bin=0, species_index=0, delta_conc=100.)
    bio.describe_state()

    bio.set_diffusion_rates([10.])

    with pytest.raises(Exception):
        #Excessive time step
        bio.diffuse_step_single_species(time_step=0.034)

    # Diffuse by a single step
    bio.diffuse_step_single_species(time_step=0.02)
    print(bio.univ)
    assert np.allclose(bio.univ[0], [80, 20])

    # Another single step
    bio.diffuse_step_single_species(time_step=0.01)
    print(bio.univ)
    assert np.allclose(bio.univ[0], [74, 26])



def test_diffuse_step_single_species2():
    bio.initialize_universe(n_cells=1, n_species=1)
    bio.set_uniform_concentration(species_index=0, conc=8.0)
    with pytest.raises(Exception):
        bio.diffuse_step_single_species(time_step=.001)    # Must set the diffusion rates first

    bio.set_diffusion_rates([20.])
    bio.describe_state()    # 1 bins and 1 species:  [[8.]]

    bio.diffuse_step_single_species(time_step=3)    # With just 1 bin, nothing happens
    bio.describe_state()
    assert np.allclose(bio.univ[0], [8])



def test_diffuse_step_single_species3():
    with pytest.raises(Exception):
        bio.diffuse_step_single_species(time_step=.001)    # Must first initialize the system

    bio.initialize_universe(n_cells=5, n_species=2)

    bio.set_uniform_concentration(species_index=0, conc=22.2)
    bio.set_uniform_concentration(species_index=1, conc=66.6)
    bio.set_diffusion_rates([2.78, 3.14])
    bio.describe_state(show_diff=True)
    """
    5 bins and 2 species:
    Species 0. Diff rate: 2.78. Conc:  [22.2 22.2 22.2 22.2 22.2]
    Species 1. Diff rate: 3.14. Conc:  [66.6 66.6 66.6 66.6 66.6]
    """

    # Diffusing a uniform distribution won't change it
    bio.diffuse_step_single_species(time_step=0.08, species_index=0)
    bio.diffuse_step_single_species(time_step=0.08, species_index=1)

    assert np.allclose(bio.univ[0], np.full(5, 22.2, dtype=float))
    assert np.allclose(bio.univ[1], np.full(5, 66.6, dtype=float))



def test_diffuse_step_1():
    bio.initialize_universe(n_cells=3, n_species=2)

    with pytest.raises(Exception):
        bio.set_bin_conc(bin=3, species_index=0, conc=100.) # Bin number out of range

    bio.set_bin_conc(bin=1, species_index=0, conc=100.)

    with pytest.raises(Exception):
        bio.set_species_conc(species_index=2, conc_list=[10, 20, 50])   # species_index out of range

    bio.set_species_conc(species_index=1, conc_list=[10, 20, 50])

    bio.set_diffusion_rates([5., 20.])
    bio.describe_state(show_diff=True)
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




######### TODO: fix all the tests below

def test_diffuse_step_single_species_3x():   # TODO: fix
    bio.initialize_universe(n_cells=2, n_species=1)

    bio.inject_conc_to_cell(bin=0, species_index=0, delta_conc=10.)
    print(bio.univ)

    bio.set_diffusion_rates([1.])

    for i in range(4):
        bio.diffuse_step_single_species(time_step=.4)
        print(f"At end of step {i+1}:")
        print(bio.univ)

def test_diffuse_step_single_species_4x():   # TODO: fix
    bio.initialize_universe(n_cells=3, n_species=1)

    bio.inject_conc_to_cell(bin=1, species_index=0, delta_conc=10.)
    print(bio.univ)

    bio.set_diffusion_rates([.5])

    for i in range(4):
        bio.diffuse_step_single_species(time_step=0.6666)
        print(f"At end of step {i+1}:")
        print(bio.univ)

def test_diffuse_step_single_species_5():   # TODO: fix
    bio.initialize_universe(n_cells=5, n_species=1)

    bio.inject_conc_to_cell(bin=0, species_index=0, delta_conc=10.)
    print(bio.univ)

    bio.set_diffusion_rates([.5])

    for i in range(20):
        bio.diffuse_step_single_species(time_step=0.6666)
        print(f"At end of step {i+1}:")
        print(bio.univ)



def test_diffuse_step_single_species_6():   # TODO: fix
    # Exploring effect of time resolution on accuracy
    bio.initialize_universe(n_cells=10, n_species=1)

    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
    print(bio.univ)

    bio.set_diffusion_rates([0.1])

    # The time advance (tf * n_steps) always remains constant at 33.3
    # TODO: provide a function that accepts the time advance and n_steps, and determines tf
    # TODO: provide a function that accepts the time advance and tf, and determines n_steps

    n_steps = 10
    tf = 3.33
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_step=tf)
        #print(f"At end of step {i+1}:")

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 20
    tf = 3.33/2
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_step=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 30
    tf = 3.33/3
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_step=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 50
    tf = 3.33/5
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_step=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 100
    tf = 3.33/10
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_step=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 1000
    tf = 3.33/100
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_step=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 10000
    tf = 3.33/1000
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_step=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 100000
    tf = 3.33/10000
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_step=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)