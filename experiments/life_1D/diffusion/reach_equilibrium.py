"""
Exploring reaching equilibrium, first on a shorter timescale and then a longer one (but
both with identical time steps.)

The system starts out with a pulse in bin 2 (the 3rd bin from the left)

Notice the diffusing pulse "bouncing" off the left wall after total time 30
"""

from modules.chemicals.chemicals import Chemicals as chem
from life_1D.bio_sim_1d import BioSim1D as bio


chem_data = chem(diffusion_rates=[0.1])
bio.initialize_universe(n_bins=10, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=0.)
bio.inject_conc_to_cell(species_index=0, bin=2, delta_conc=10.)

bio.describe_state(show_diffusion_rates=True)


print("\n\nSTARTING on SHORTER time scales.  Dtime=10, with time steps of 0.1 ...")

total_time = 0.
for i in range(10):
    delta_time = 10.
    status = bio.diffuse(time_duration=delta_time, time_step=0.1)
    total_time += delta_time
    print(f"\nAfter Delta time {delta_time}.  TOTAL TIME {total_time}  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)



print("\n\nREPEATING to LONGER time scales.  Dtime=100, again with time steps of 0.1 ...")

# Reset the concentrations
bio.set_uniform_concentration(species_index=0, conc=0.)
bio.inject_conc_to_cell(species_index=0, bin=2, delta_conc=10.)

total_time = 0.
for i in range(20):
    delta_time = 100.
    status = bio.diffuse(time_duration=delta_time, time_step=0.1)
    total_time += delta_time
    print(f"\nAfter Delta time {delta_time}.  TOTAL TIME {total_time}  ({status['steps']} steps taken):")
    bio.describe_state(concise=True)