import numpy as np
import pytest
from bio_sim import BioSim


def test_initialize_universe():
    BioSim.initialize_universe(n_cells=10, n_species=1)

    assert BioSim.n_cells == 10
    assert BioSim.n_species == 1

    print(BioSim.univ)



def test_set_diffusion_rates():
    BioSim.initialize_universe(n_cells=5, n_species=1)
    BioSim.set_diffusion_rates([.2])



def test_set_uniform_concentration():
    BioSim.initialize_universe(n_cells=5, n_species=1)
    BioSim.set_uniform_concentration(c=0.1)

    print(BioSim.univ)


def test_inject_conc_to_cell():
    BioSim.initialize_universe(n_cells=5, n_species=1)

    with pytest.raises(Exception):
        # cell_index out-of-bounds
        BioSim.inject_conc_to_cell(cell_index=5, delta_conc=10.)

    BioSim.inject_conc_to_cell(cell_index=1, delta_conc=10.)

    print(BioSim.univ)



def test_diffuse_step_single_species():
    BioSim.initialize_universe(n_cells=5, n_species=1)
    BioSim.set_uniform_concentration(c=1.0)
    with pytest.raises(Exception):
        BioSim.diffuse_step_single_species()    # Must set the diffusion rates first

    BioSim.set_diffusion_rates([10.])

    old_univ = BioSim.univ
    result = BioSim.diffuse_step_single_species()
    print(BioSim.univ)
    assert np.allclose(old_univ, result)    # Diffusing a uniform distribution won't change it


def test_diffuse_step_single_species_2():
    BioSim.initialize_universe(n_cells=2, n_species=1)

    BioSim.inject_conc_to_cell(cell_index=0, delta_conc=10.)
    print(BioSim.univ)

    BioSim.set_diffusion_rates([1.])

    BioSim.diffuse_step_single_species(time_fraction=.5)
    print(BioSim.univ)

    BioSim.diffuse_step_single_species(time_fraction=.5)
    print(BioSim.univ)


def test_diffuse_step_single_species_3():
    BioSim.initialize_universe(n_cells=2, n_species=1)

    BioSim.inject_conc_to_cell(cell_index=0, delta_conc=10.)
    print(BioSim.univ)

    BioSim.set_diffusion_rates([1.])

    for i in range(4):
        BioSim.diffuse_step_single_species(time_fraction=.4)
        print(f"At end of step {i+1}:")
        print(BioSim.univ)

def test_diffuse_step_single_species_4():
    BioSim.initialize_universe(n_cells=3, n_species=1)

    BioSim.inject_conc_to_cell(cell_index=1, delta_conc=10.)
    print(BioSim.univ)

    BioSim.set_diffusion_rates([.5])

    for i in range(4):
        BioSim.diffuse_step_single_species(time_fraction=0.6666)
        print(f"At end of step {i+1}:")
        print(BioSim.univ)

def test_diffuse_step_single_species_5():
    BioSim.initialize_universe(n_cells=5, n_species=1)

    BioSim.inject_conc_to_cell(cell_index=0, delta_conc=10.)
    print(BioSim.univ)

    BioSim.set_diffusion_rates([.5])

    for i in range(20):
        BioSim.diffuse_step_single_species(time_fraction=0.6666)
        print(f"At end of step {i+1}:")
        print(BioSim.univ)



def test_diffuse_step_single_species_6():
    # Exploring effect of time resolution on accuracy
    BioSim.initialize_universe(n_cells=10, n_species=1)

    BioSim.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    print(BioSim.univ)

    BioSim.set_diffusion_rates([0.1])

    # The time advance (tf * n_steps) always remains constant at 33.3
    # TODO: provide a function that accepts the time advance and n_steps, and determines tf
    # TODO: provide a function that accepts the time advance and tf, and determines n_steps

    n_steps = 10
    tf = 3.33
    for i in range(n_steps):
        BioSim.diffuse_step_single_species(time_fraction=tf)
        #print(f"At end of step {i+1}:")

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim.univ)


    BioSim.initialize_universe(n_cells=10, n_species=1)
    BioSim.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 20
    tf = 3.33/2
    for i in range(n_steps):
        BioSim.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim.univ)


    BioSim.initialize_universe(n_cells=10, n_species=1)
    BioSim.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 30
    tf = 3.33/3
    for i in range(n_steps):
        BioSim.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim.univ)


    BioSim.initialize_universe(n_cells=10, n_species=1)
    BioSim.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 50
    tf = 3.33/5
    for i in range(n_steps):
        BioSim.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim.univ)


    BioSim.initialize_universe(n_cells=10, n_species=1)
    BioSim.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 100
    tf = 3.33/10
    for i in range(n_steps):
        BioSim.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim.univ)


    BioSim.initialize_universe(n_cells=10, n_species=1)
    BioSim.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 1000
    tf = 3.33/100
    for i in range(n_steps):
        BioSim.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim.univ)


    BioSim.initialize_universe(n_cells=10, n_species=1)
    BioSim.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 10000
    tf = 3.33/1000
    for i in range(n_steps):
        BioSim.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim.univ)


    BioSim.initialize_universe(n_cells=10, n_species=1)
    BioSim.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 100000
    tf = 3.33/10000
    for i in range(n_steps):
        BioSim.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim.univ)