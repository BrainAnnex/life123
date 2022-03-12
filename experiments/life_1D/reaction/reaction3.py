"""
A one-bin 2A <-> 3B reaction, with 1st-order kinetics in both directions,
taken to equilibrium
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(diffusion_rates=[0.1, 0.1], names=["A", "B"])  # NOTE: diffusion_rates not used for now

rxn = Reactions(chem_data)

# Reaction A -> 3B , with 1st-order kinetics in both directions
rxn.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)

bio.describe_state(show_diffusion_rates=True)

rxn.describe_reactions()

# Low-level view of the reactions data
for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")


# First step
bio.reaction_step(0.1)
bio.describe_state()
"""
1 bins and 2 species:
 [[20.]
 [35.]]
"""


# Numerous more steps
bio.react(time_step=0.1, n_steps=100)

bio.describe_state()

# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:
"""
1 bins and 2 species:
 [[16.25]
 [40.625]]
"""

A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
print(f"Ratio of equilibrium concentrations: {B_eq / A_eq}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")
