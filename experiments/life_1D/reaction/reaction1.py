"""
A simple A <-> B reaction between 2 species with initial uniform concentrations,
with 1st-order kinetics in both directions, taken to equilibrium.
Diffusion (non-applicable) not done
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(diffusion_rates=[0.1, 0.1], names=["A", "B"])  # NOTE: diffusion_rates not used for now
bio.initialize_universe(n_bins=3, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)

bio.describe_state(show_diffusion_rates=True)


rxn = Reactions(chem_data)

# Reaction A -> B , with 1st-order kinetics in both directions
rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
#rxn.add_reaction(reactants=[(1,0,1)], products=[(1,1,1)], forward_rate=3., reverse_rate=2.)

bio.all_reactions = rxn


print("Number of reactions: ", rxn.number_of_reactions())

for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")

rxn.describe_reactions()


# First step
bio.reaction_step(0.1)
bio.describe_state()
"""
3 bins and 2 species:
 [[17. 17. 17.]
 [43. 43. 43.]]
"""


# Numerous more steps
for i in range(10):
    bio.reaction_step(0.1)

bio.describe_state()

# Consistent with the 3/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:
"""
3 bins and 2 species:
 [[23.99316406 23.99316406 23.99316406]
 [36.00683594 36.00683594 36.00683594]]
"""

A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
print(f"Ratio of equilibrium concentrations: {B_eq / A_eq}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")
