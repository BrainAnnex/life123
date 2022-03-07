"""
Exploring effect of time resolution on accuracy.

In the examples below, the time advance always remains constant,
but the number of steps used to arrive there vary
"""

from modules.chemicals.chemicals import Chemicals as chem
from life_1D.bio_sim_1d import BioSim1D as bio


def set_initial_condition():
    # Set or reset the initial concentrations
    bio.set_uniform_concentration(species_index=0, conc=0.)
    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)



chem_data = chem(diffusion_rates=[0.1])
bio.initialize_universe(n_bins=10, chem_data=chem_data)

set_initial_condition()

bio.describe_state(show_diffusion_rates=True)


t_final = 33.3

bio.diffuse(time_duration=t_final, n_steps=10, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(time_duration=t_final, n_steps=20, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(time_duration=t_final, n_steps=30, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(time_duration=t_final, n_steps=50, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(time_duration=t_final, n_steps=100, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(time_duration=t_final, n_steps=1000, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(time_duration=t_final, n_steps=10000, verbose=True)

set_initial_condition()     # Reset the concentrations
bio.diffuse(time_duration=t_final, n_steps=100000, verbose=True)
