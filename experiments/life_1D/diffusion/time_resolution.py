"""
Exploring effect of time resolution on accuracy
"""

from life_1D.bio_sim_1d import BioSim1D as bio


def set_initial_condition():
    # Set or reset the initial concentrations
    bio.set_uniform_concentration(species_index=0, conc=0.)
    bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)


def diffuse_single_species_from_give_initial(t_final, n_steps):
    set_initial_condition()     # Reset the concentrations
    #bio.describe_state(show_diff=True)
    time_step = t_final / n_steps

    for i in range(n_steps):
        bio.diffuse_step_single_species(time_step=time_step)

    print(f"System at time {t_final}, at end of {n_steps} steps of size {time_step}:")
    print(bio.univ[0])



bio.initialize_universe(n_bins=10, n_species=1)

bio.set_diffusion_rates([0.1])

set_initial_condition()

bio.describe_state(show_diff=True)


# The time advance (tf * n_steps) always remains constant at 33.3
# TODO: provide a function that accepts the time advance and n_steps, and determines tf
# TODO: provide a function that accepts the time advance and tf, and determines n_steps


diffuse_single_species_from_give_initial(33.3, 10)

diffuse_single_species_from_give_initial(33.3, 20)

diffuse_single_species_from_give_initial(33.3, 30)

diffuse_single_species_from_give_initial(33.3, 50)

diffuse_single_species_from_give_initial(33.3, 100)

diffuse_single_species_from_give_initial(33.3, 1000)

diffuse_single_species_from_give_initial(33.3, 10000)

diffuse_single_species_from_give_initial(33.3, 100000)
