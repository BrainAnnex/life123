# These are tests specifically for diffusion in 2D;
# for general tests of 2D system, see test_biosim_2d.py


import pytest
import numpy as np
from life_2D.bio_sim_2d import BioSim2D
from modules.chemicals.chemicals import Chemicals as chem
from modules.numerical.numerical import Numerical as num


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
    bio.describe_state()

    delta_conc = bio.diffuse_step_single_species(time_step=delta_t,
                                                 h=delta_s, species_index=0)
    print(delta_conc)

    # Convolution matrix
    conv_matrix = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]])

    m_expanded = num.expand_matrix_boundary(initial_concs)
    print(m_expanded)

    print()
    for i in range(3):
        for j in range(3):
            print("i,j :", i, j)
            tile = m_expanded[i:i+3, j:j+3]
            print(tile)
            print("---")
            prod = tile * conv_matrix
            print(prod)
            print(prod.sum())
            print(prod.sum() * diff * delta_t / (delta_s ** 2))
            print("~~~~~~~~~~~")




