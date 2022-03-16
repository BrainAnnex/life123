"""
Reactions  A + B <-> C, with 1st-order kinetics for each species,
taken to equilibrium
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(names=["A", "B", "C"], diffusion_rates=[50., 50., 1.])

rxn = Reactions(chem_data)

# Reactions A + B <-> C , with 1st-order kinetics for each species
rxn.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=20., reverse_rate=2.)

bio.initialize_universe(n_bins=7, chem_data=chem_data, reactions=rxn)

bio.set_bin_conc(bin=0, species_index=0, conc=20.)
bio.set_bin_conc(bin=6, species_index=1, conc=20.)

bio.describe_state(show_diffusion_rates=True)

rxn.describe_reactions()

# Low-level view of the reactions data
for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")


delta_t = 0.002


# First step
bio.react_diffuse(time_step=delta_t, n_steps=1)
bio.describe_state()
"""

"""

for _ in range(7):
    bio.react_diffuse(time_step=delta_t, n_steps=1)
    bio.describe_state()

for _ in range(4):
    print("\n- 10 steps later...")
    bio.react_diffuse(time_step=delta_t, n_steps=10)
    bio.describe_state()

for _ in range(4):
    print("\n--- 30 steps later...")
    bio.react_diffuse(time_step=delta_t, n_steps=30)
    bio.describe_state()

for _ in range(4):
    print("\n------ 50 steps later...")
    bio.react_diffuse(time_step=delta_t, n_steps=50)
    bio.describe_state()

for _ in range(4):
    print("\n------------------ 150 steps later...")
    bio.react_diffuse(time_step=delta_t, n_steps=150)
    bio.describe_state()

for _ in range(2):
    print("\n---------------------------------------- 1,000 steps later...")
    bio.react_diffuse(time_step=delta_t, n_steps=1000)
    bio.describe_state()


# Verify equilibrium concentrations (sampled in the 1st bin; at this point, all bins have equilibrated)
A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
print(f"\nRatio of equilibrium concentrations ((C_eq) / (A_eq * B_eq)) : {(C_eq) / (A_eq * B_eq)}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")
# Both are essentially equal, as expected