# These are tests specifically for diffusion in 2D;
# for general tests of 2D system, see test_biosim_2d.py


import pytest
import numpy as np
from life_2D.bio_sim_2d import BioSim2D
from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions


#########   TESTS OF DIFFUSION : single species, one step    #########

def test_diffuse_step_single_species_1():
    delta_t = 0.01
    delta_x = 2
    diff = 10.

    # Initialize the system
    chem_data = chem(diffusion_rates=[diff])    # Just 1 chemical species
    bio = BioSim2D(n_bins=(3,3), chem_data=chem_data)
    initial_concs = np.array([[50, 80, 20], [10, 60, 0], [30, 15, 100]])
    bio.set_species_conc(species_index=0, conc_list=initial_concs)
    #bio.describe_state()
