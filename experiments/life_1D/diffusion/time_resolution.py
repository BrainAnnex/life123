"""
Exploring effect of time resolution on accuracy
"""

from life_1D.bio_sim_1d import BioSim1D as bio


bio.initialize_universe(n_bins=10, n_species=1)

bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
print(bio.univ)

bio.set_diffusion_rates([0.1])

# The time advance (tf * n_steps) always remains constant at 33.3
# TODO: provide a function that accepts the time advance and n_steps, and determines tf
# TODO: provide a function that accepts the time advance and tf, and determines n_steps

n_steps = 10
tf = 3.33
for i in range(n_steps):
    bio.diffuse_step_single_species(time_step=tf)
    #print(f"At end of step {i+1}:")

print(f"At end of {n_steps} steps of size {tf}:")
print(bio.univ)


bio.initialize_universe(n_bins=10, n_species=1)
bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
#print(BioSim.univ)
n_steps = 20
tf = 3.33/2
for i in range(n_steps):
    bio.diffuse_step_single_species(time_step=tf)

print(f"At end of {n_steps} steps of size {tf}:")
print(bio.univ)


bio.initialize_universe(n_bins=10, n_species=1)
bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
#print(BioSim.univ)
n_steps = 30
tf = 3.33/3
for i in range(n_steps):
    bio.diffuse_step_single_species(time_step=tf)

print(f"At end of {n_steps} steps of size {tf}:")
print(bio.univ)


bio.initialize_universe(n_bins=10, n_species=1)
bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
#print(BioSim.univ)
n_steps = 50
tf = 3.33/5
for i in range(n_steps):
    bio.diffuse_step_single_species(time_step=tf)

print(f"At end of {n_steps} steps of size {tf}:")
print(bio.univ)


bio.initialize_universe(n_bins=10, n_species=1)
bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
#print(BioSim.univ)
n_steps = 100
tf = 3.33/10
for i in range(n_steps):
    bio.diffuse_step_single_species(time_step=tf)

print(f"At end of {n_steps} steps of size {tf}:")
print(bio.univ)


bio.initialize_universe(n_bins=10, n_species=1)
bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
#print(BioSim.univ)
n_steps = 1000
tf = 3.33/100
for i in range(n_steps):
    bio.diffuse_step_single_species(time_step=tf)

print(f"At end of {n_steps} steps of size {tf}:")
print(bio.univ)


bio.initialize_universe(n_bins=10, n_species=1)
bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
#print(BioSim.univ)
n_steps = 10000
tf = 3.33/1000
for i in range(n_steps):
    bio.diffuse_step_single_species(time_step=tf)

print(f"At end of {n_steps} steps of size {tf}:")
print(bio.univ)


bio.initialize_universe(n_bins=10, n_species=1)
bio.inject_conc_to_cell(bin=2, species_index=0, delta_conc=10.)
#print(BioSim.univ)
n_steps = 100000
tf = 3.33/10000
for i in range(n_steps):
    bio.diffuse_step_single_species(time_step=tf)

print(f"At end of {n_steps} steps of size {tf}:")
print(bio.univ)