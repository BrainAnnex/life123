import numpy as np
import pytest
from Life_1D.bio_sim import BioSim1D


def test_initialize_universe():
    BioSim1D.initialize_universe(n_cells=10, n_species=1)

    assert BioSim1D.n_cells == 10
    assert BioSim1D.n_species == 1

    print(BioSim1D.univ)



def test_set_diffusion_rates():
    BioSim1D.initialize_universe(n_cells=5, n_species=1)
    BioSim1D.set_diffusion_rates([.2])



def test_set_uniform_concentration():
    BioSim1D.initialize_universe(n_cells=5, n_species=1)
    BioSim1D.set_uniform_concentration(c=0.1)

    print(BioSim1D.univ)


def test_inject_conc_to_cell():
    BioSim1D.initialize_universe(n_cells=5, n_species=1)

    with pytest.raises(Exception):
        # cell_index out-of-bounds
        BioSim1D.inject_conc_to_cell(cell_index=5, delta_conc=10.)

    BioSim1D.inject_conc_to_cell(cell_index=1, delta_conc=10.)

    print(BioSim1D.univ)



def test_diffuse_step_single_species():
    BioSim1D.initialize_universe(n_cells=5, n_species=1)
    BioSim1D.set_uniform_concentration(c=1.0)
    with pytest.raises(Exception):
        BioSim1D.diffuse_step_single_species()    # Must set the diffusion rates first

    BioSim1D.set_diffusion_rates([10.])

    old_univ = BioSim1D.univ
    result = BioSim1D.diffuse_step_single_species()
    print(BioSim1D.univ)
    assert np.allclose(old_univ, result)    # Diffusing a uniform distribution won't change it


def test_diffuse_step_single_species_2():
    BioSim1D.initialize_universe(n_cells=2, n_species=1)

    BioSim1D.inject_conc_to_cell(cell_index=0, delta_conc=10.)
    print(BioSim1D.univ)

    BioSim1D.set_diffusion_rates([1.])

    BioSim1D.diffuse_step_single_species(time_fraction=.5)
    print(BioSim1D.univ)

    BioSim1D.diffuse_step_single_species(time_fraction=.5)
    print(BioSim1D.univ)


def test_diffuse_step_single_species_3():
    BioSim1D.initialize_universe(n_cells=2, n_species=1)

    BioSim1D.inject_conc_to_cell(cell_index=0, delta_conc=10.)
    print(BioSim1D.univ)

    BioSim1D.set_diffusion_rates([1.])

    for i in range(4):
        BioSim1D.diffuse_step_single_species(time_fraction=.4)
        print(f"At end of step {i+1}:")
        print(BioSim1D.univ)

def test_diffuse_step_single_species_4():
    BioSim1D.initialize_universe(n_cells=3, n_species=1)

    BioSim1D.inject_conc_to_cell(cell_index=1, delta_conc=10.)
    print(BioSim1D.univ)

    BioSim1D.set_diffusion_rates([.5])

    for i in range(4):
        BioSim1D.diffuse_step_single_species(time_fraction=0.6666)
        print(f"At end of step {i+1}:")
        print(BioSim1D.univ)

def test_diffuse_step_single_species_5():
    BioSim1D.initialize_universe(n_cells=5, n_species=1)

    BioSim1D.inject_conc_to_cell(cell_index=0, delta_conc=10.)
    print(BioSim1D.univ)

    BioSim1D.set_diffusion_rates([.5])

    for i in range(20):
        BioSim1D.diffuse_step_single_species(time_fraction=0.6666)
        print(f"At end of step {i+1}:")
        print(BioSim1D.univ)



def test_diffuse_step_single_species_6():
    # Exploring effect of time resolution on accuracy
    BioSim1D.initialize_universe(n_cells=10, n_species=1)

    BioSim1D.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    print(BioSim1D.univ)

    BioSim1D.set_diffusion_rates([0.1])

    # The time advance (tf * n_steps) always remains constant at 33.3
    # TODO: provide a function that accepts the time advance and n_steps, and determines tf
    # TODO: provide a function that accepts the time advance and tf, and determines n_steps

    n_steps = 10
    tf = 3.33
    for i in range(n_steps):
        BioSim1D.diffuse_step_single_species(time_fraction=tf)
        #print(f"At end of step {i+1}:")

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim1D.univ)


    BioSim1D.initialize_universe(n_cells=10, n_species=1)
    BioSim1D.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 20
    tf = 3.33/2
    for i in range(n_steps):
        BioSim1D.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim1D.univ)


    BioSim1D.initialize_universe(n_cells=10, n_species=1)
    BioSim1D.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 30
    tf = 3.33/3
    for i in range(n_steps):
        BioSim1D.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim1D.univ)


    BioSim1D.initialize_universe(n_cells=10, n_species=1)
    BioSim1D.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 50
    tf = 3.33/5
    for i in range(n_steps):
        BioSim1D.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim1D.univ)


    BioSim1D.initialize_universe(n_cells=10, n_species=1)
    BioSim1D.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 100
    tf = 3.33/10
    for i in range(n_steps):
        BioSim1D.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim1D.univ)


    BioSim1D.initialize_universe(n_cells=10, n_species=1)
    BioSim1D.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 1000
    tf = 3.33/100
    for i in range(n_steps):
        BioSim1D.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim1D.univ)


    BioSim1D.initialize_universe(n_cells=10, n_species=1)
    BioSim1D.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 10000
    tf = 3.33/1000
    for i in range(n_steps):
        BioSim1D.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim1D.univ)


    BioSim1D.initialize_universe(n_cells=10, n_species=1)
    BioSim1D.inject_conc_to_cell(cell_index=2, delta_conc=10.)
    #print(BioSim.univ)
    n_steps = 100000
    tf = 3.33/10000
    for i in range(n_steps):
        BioSim1D.diffuse_step_single_species(time_fraction=tf)

    print(f"At end of {n_steps} steps of size {tf}:")
    print(BioSim1D.univ)