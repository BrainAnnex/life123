"""
THE DANGER OF EXCESSIVELY LARGE SINGLE TIME STEPS IN DIFFUSION

(We're assuming bins of size 1 ; TODO: adjust for any size)

->  When the time step is (0.5 / diffusion rate),
    a 2-bin system equilibrates in a single step,
    and some 3-bin systems can over-shoot equilibrium!

->  When the time step is (0.33333 / diffusion rate),
    some 3-bin systems equilibrate in a single step

So, (0.33333 / diffusion rate) is a - rather lax - upper bound for
sensible single time steps!

That value of 0.33333 is saved in the Class variable "time_step_threshold"
"""

from life_1D.bio_sim_1d import BioSim1D as bio


###########################
# Simulate a 2-bin system
###########################
bio.initialize_universe(n_cells=2, n_species=1)
bio.inject_conc_to_cell(bin=0, delta_conc=100., species_index=0)
bio.describe_state()
# 2 cells and 1 species:  [[100.   0.]]

bio.set_diffusion_rates([10.])

bio.time_step_threshold = 0.51  # To bypass limits that are typically in place

# When the time step is (0.5 / diffusion rate),
# a 2-bin system equilibrates in a single step!
bio.diffuse_step_single_species(time_step=0.05)
print(bio.univ)     # [[50. 50.]]  : the two bins have equilibrated!



######################################
# Simulate a 3-bin system
# with an excessive single time step
######################################
bio.initialize_universe(n_cells=3, n_species=1)
bio.inject_conc_to_cell(bin=1, delta_conc=100., species_index=0)
bio.describe_state()
#3 bins and 1 species:   [[  0. 100.   0.]]

bio.set_diffusion_rates([10.])

bio.time_step_threshold = 0.51      # To bypass a limit that is typically in place

# When the time step is (0.5 / diffusion rate),
# a 3-bin system can overshoot equilibrium!
bio.diffuse_step_single_species(time_step=0.05)
print(bio.univ)     # [[50.  0. 50.]] : the diffusion has over-shot equilibrium!!!



#############################################
# Re-Simulate the 3-bin system,
# with a somewhat smaller single time step
#############################################
bio.initialize_universe(n_cells=3, n_species=1)
bio.inject_conc_to_cell(bin=1, delta_conc=100., species_index=0)
bio.describe_state()
#3 bins and 1 species:   [[  0. 100.   0.]]

bio.set_diffusion_rates([10.])

bio.time_step_threshold = 0.34      # To bypass a limit that is typically in place

# When the time step is (0.33333 / diffusion rate),
# a 3-bin system, configured as above, equilibrates in a single step!
bio.diffuse_step_single_species(time_step=0.033333)
print(bio.univ)     # [[33.333 33.334 33.333]] : the three bins have equilibrated!
