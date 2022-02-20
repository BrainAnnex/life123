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
        bio.inject_conc_to_cell(cell_index=5, species_index=0, delta_conc=10.)

    with pytest.raises(Exception):
        # species_index out of bounds
        bio.inject_conc_to_cell(cell_index=5, species_index=1, delta_conc=10.)

    bio.inject_conc_to_cell(cell_index=1, species_index=0, delta_conc=10.)

    bio.describe_state()



def test_diffuse_step_single_species():
    bio.initialize_universe(n_cells=5, n_species=1)
    bio.set_uniform_concentration(species_index=0, conc=1.0)
    with pytest.raises(Exception):
        bio.diffuse_step_single_species()    # Must set the diffusion rates first

    bio.set_diffusion_rates([10.])

    old_univ = bio.univ
    result = bio.diffuse_step_single_species(time_fraction=0.01)
    bio.describe_state()
    assert np.allclose(old_univ, result)    # Diffusing a uniform distribution won't change it



def test_diffuse_step_single_species_2():
    bio.initialize_universe(n_cells=2, n_species=1)

    bio.inject_conc_to_cell(cell_index=0, species_index=0, delta_conc=10.)
    print(bio.univ)

    bio.set_diffusion_rates([1.])

    bio.diffuse_step_single_species(time_fraction=.5)
    print(bio.univ)

    bio.diffuse_step_single_species(time_fraction=.5)
    print(bio.univ)


def test_diffuse_step_single_species_3():
    bio.initialize_universe(n_cells=2, n_species=1)

    bio.inject_conc_to_cell(cell_index=0, delta_conc=10.)
    print(bio.univ)

    bio.set_diffusion_rates([1.])

    for i in range(4):
        bio.diffuse_step_single_species(time_fraction=.4)
        print(f"At end of step {i+1}:")
        print(bio.univ)

def test_diffuse_step_single_species_4():
    bio.initialize_universe(n_cells=3, n_species=1)

    bio.inject_conc_to_cell(cell_index=1, delta_conc=10.)
    print(bio.univ)

    bio.set_diffusion_rates([.5])

    for i in range(4):
        bio.diffuse_step_single_species(time_fraction=0.6666)
        print(f"At end of step {i+1}:")
        print(bio.univ)

def test_diffuse_step_single_species_5():
    bio.initialize_universe(n_cells=5, n_species=1)

    bio.inject_conc_to_cell(cell_index=0, delta_conc=10.)
    print(bio.univ)

    bio.set_diffusion_rates([.5])

    for i in range(20):
        bio.diffuse_step_single_species(time_fraction=0.6666)
        print(f"At end of step {i+1}:")
        print(bio.univ)



def test_diffuse_step_single_species_6():
    # Exploring effect of time resolution on accuracy
    bio.initialize_universe(n_cells=10, n_species=1)

    bio.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    print(bio.univ)

    bio.set_diffusion_rates([0.1])

    # The time advance (tf * n_steps) always remains constant at 33.3
    # TODO: provide a function that accepts the time advance and n_steps, and determines tf
    # TODO: provide a function that accepts the time advance and tf, and determines n_steps

    n_steps = 10
    tf = 3.33
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_fraction=tf)
        #print(f"At end of step {i+1}:")

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 20
    tf = 3.33/2
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 30
    tf = 3.33/3
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 50
    tf = 3.33/5
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 100
    tf = 3.33/10
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 1000
    tf = 3.33/100
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 10000
    tf = 3.33/1000
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)


    bio.initialize_universe(n_cells=10, n_species=1)
    bio.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 100000
    tf = 3.33/10000
    for i in range(n_steps):
        bio.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(bio.univ)