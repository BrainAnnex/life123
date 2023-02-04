# These are tests specifically for diffusion in 2D;
# for general tests of 2D system, see test_biosim_2d.py


import numpy as np
from src.life_2D.bio_sim_2d import BioSim2D
from src.modules.reactions.reaction_data import ReactionData as chem
from src.modules.numerical.numerical import Numerical as num



def evaluate_diffusion(bio, initial_concs, diff, delta_t, delta_s):
    # Use a straightforward computation with a convolution matrix,
    # to determine what the diffusion ought to be

    conv_matrix = np.array([[0, 1, 0],
                            [1,-4, 1],
                            [0, 1, 0]])     # TODO: just using 5-point stencils for now

    m_expanded = num.expand_matrix_boundary(initial_concs)
    #print(m_expanded)

    dims_x, dims_y = bio.system_size()

    expected_m = np.zeros((dims_x, dims_y), dtype=float)

    for i in range(dims_x):
        for j in range(dims_y):
            #print("i,j :", i, j)
            tile = m_expanded[i:i+3, j:j+3]    # 3x3 is the size of the convolution matrix
            #print(tile)
            #print("---")
            prod = tile * conv_matrix
            #print(prod)
            #print(prod.sum())
            #print(prod.sum() * diff * delta_t / (delta_s ** 2))
            #print("~~~~~~~~~~~")
            expected_m[i, j] = prod.sum() * diff * delta_t / (delta_s ** 2)

    #print(expected_m)
    return expected_m



#########   TESTS OF DIFFUSION : all species, multiple steps    #########

def test_diffuse_1():
    delta_t = 0.01
    delta_s = 2
    diff = 10.
    initial_concs = np.array([[50, 80, 20], [10, 60, 0], [30, 15, 100]])

    # Initialize the system
    chem_data = chem(diffusion_rates=[diff])    # Just 1 chemical species
    bio = BioSim2D(n_bins=(3,3), chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_data=initial_concs)
    #bio.describe_state()

    # 1st diffusion step
    bio.diffuse_step(time_step=delta_t, h=delta_s)
    #print(bio.delta_diffusion)
    bio.system += bio.delta_diffusion
    #bio.describe_state()

    # 2nd diffusion step
    bio.diffuse_step(time_step=delta_t, h=delta_s)
    #print(bio.delta_diffusion)
    bio.system += bio.delta_diffusion
    #bio.describe_state()
    snapshot_of_system_state = bio.system.copy()

    # Now, re-initialize the system, and do both steps with a single function call
    chem_data = chem(diffusion_rates=[diff])    # Just 1 chemical species
    bio = BioSim2D(n_bins=(3,3), chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_data=initial_concs)
    #bio.describe_state()

    bio.diffuse(time_step=delta_t, n_steps=2, h=delta_s)    # 2 diffusion steps
    #bio.describe_state()
    assert np.allclose(bio.system, snapshot_of_system_state)



def test_diffuse_2():
    delta_t = 0.01
    delta_s = 2
    diff = (10., 4.)    # 2 chemical species

    initial_concs_0 = np.array([[50, 80, 34.4, 20], [5.34, 10, 60, 0], [30, 15, 100, 9.24], [40, 45.3, 64.4, 80], [30, 9.4, 39.4, 10.5]])
    initial_concs_1 = np.array([[35, 110, 1.5, 13], [53, 60, 6.7, 40], [54, 3.5, 0, 19.8], [120, 25.3, 87, 11.5], [10, 19, 79.14, 50.8]])

    # Initialize the system
    chem_data = chem(diffusion_rates=diff)
    bio = BioSim2D(n_bins=(5,4), chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_data=initial_concs_0)
    bio.set_species_conc(species_index=1, conc_data=initial_concs_1)
    #bio.describe_state()

    # 3 diffusion steps
    for _ in range(3):
        bio.diffuse_step(time_step=delta_t, h=delta_s)
        #print(bio.delta_diffusion)
        bio.system += bio.delta_diffusion
        #bio.describe_state()

    snapshot_of_system_state = bio.system.copy()


    # Now, re-initialize the system, and do all diffusion steps with a single function call
    chem_data = chem(diffusion_rates=diff)
    bio = BioSim2D(n_bins=(5,4), chem_data=chem_data)
    bio.set_species_conc(species_index=0, conc_data=initial_concs_0)
    bio.set_species_conc(species_index=1, conc_data=initial_concs_1)
    #bio.describe_state()

    bio.diffuse(time_step=delta_t, n_steps=3, h=delta_s)    # 3 diffusion steps
    #bio.describe_state()
    assert np.allclose(bio.system, snapshot_of_system_state)



#########   TESTS OF DIFFUSION : all species, one step    #########

def test_diffuse_step_1():
    delta_t = 0.01
    delta_s = 3
    diff_rates = [8., 6.]    # 2 chemical species

    # Initialize the system
    chem_data = chem(diffusion_rates=diff_rates)
    bio = BioSim2D(n_bins=(5,4), chem_data=chem_data)

    initial_concs_0 = np.array([[50, 80, 34.4, 20], [5.34, 10, 60, 0], [30, 15, 100, 9.24], [40, 45.3, 64.4, 80], [30, 9.4, 39.4, 10.5]])
    initial_concs_1 = np.array([[30, 110, 1.5, 20], [53, 60, 6.7, 40], [54, 3.5, 0, 19.8], [120, 25.3, 87, 11.5], [10, 19, 79.14, 50.8]])

    bio.set_species_conc(species_index=0, conc_data=initial_concs_0)
    bio.set_species_conc(species_index=1, conc_data=initial_concs_1)
    #bio.describe_state()

    # For now, just diffuse one chemical at a time
    delta_conc_0 = bio.diffuse_step_single_species(time_step=delta_t,
                                                   h=delta_s, species_index=0)
    #print(delta_conc_0)

    expected_m = evaluate_diffusion(bio, initial_concs_0, diff_rates[0], delta_t, delta_s)
    assert np.allclose(delta_conc_0, expected_m)
    assert np.allclose(delta_conc_0.sum(), 0)     # Verify mass conservation

    delta_conc_1 = bio.diffuse_step_single_species(time_step=delta_t,
                                                   h=delta_s, species_index=1)
    #print(delta_conc_1)

    expected_m = evaluate_diffusion(bio, initial_concs_1, diff_rates[1], delta_t, delta_s)
    assert np.allclose(delta_conc_1, expected_m)
    assert np.allclose(delta_conc_1.sum(), 0)     # Verify mass conservation


    # Note: the system is still in its initial value.  Now diffuse both chemicals at once
    #bio.describe_state()
    bio.diffuse_step(time_step=delta_t, h=delta_s)

    #print(bio.delta_diffusion)
    assert np.allclose(bio.delta_diffusion[0], delta_conc_0)
    assert np.allclose(bio.delta_diffusion[1], delta_conc_1)




#########   TESTS OF DIFFUSION : single species, one step    #########

def test_diffuse_step_single_species_1():
    delta_t = 0.01
    delta_s = 2
    diff = 10.

    # Initialize the system
    chem_data = chem(diffusion_rates=[diff])    # Just 1 chemical species
    bio = BioSim2D(n_bins=(3,3), chem_data=chem_data)
    initial_concs = np.array([[50, 80, 20], [10, 60, 0], [30, 15, 100]])
    bio.set_species_conc(species_index=0, conc_data=initial_concs)
    #bio.describe_state()

    delta_conc = bio.diffuse_step_single_species(time_step=delta_t,
                                                 h=delta_s, species_index=0)
    #print(delta_conc)

    expected_m = evaluate_diffusion(bio, initial_concs, diff, delta_t, delta_s)
    assert np.allclose(delta_conc, expected_m)
    assert np.allclose(delta_conc.sum(), 0)     # Verify mass conservation



def test_diffuse_step_single_species_2():
    delta_t = 0.01
    delta_s = 3
    diff = 8.

    # Initialize the system
    chem_data = chem(diffusion_rates=[diff])    # Just 1 chemical species
    bio = BioSim2D(n_bins=(5,4), chem_data=chem_data)
    initial_concs = np.array([[50, 80, 34.4, 20], [5.34, 10, 60, 0], [30, 15, 100, 9.24], [40, 45.3, 64.4, 80], [30, 9.4, 39.4, 10.5]])
    bio.set_species_conc(species_index=0, conc_data=initial_concs)
    #bio.describe_state()

    delta_conc = bio.diffuse_step_single_species(time_step=delta_t,
                                                 h=delta_s, species_index=0)
    #print(delta_conc)

    expected_m = evaluate_diffusion(bio, initial_concs, diff, delta_t, delta_s)
    assert np.allclose(delta_conc, expected_m)
    assert np.allclose(delta_conc.sum(), 0)     # Verify mass conservation
