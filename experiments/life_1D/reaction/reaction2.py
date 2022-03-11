"""
A one-bin A <-> 3B reaction, taken to equilibrium
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(diffusion_rates=[0.1, 0.2], names=["A", "B"])
bio.initialize_universe(n_bins=1, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)

bio.describe_state(show_diffusion_rates=True)


rxn = Reactions(chem_data)

# Reaction A -> 3B , with 1st-order kinetics in both directions
rxn.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

bio.all_reactions = rxn


print("Number of reactions: ", rxn.number_of_reactions())

for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")

rxn.describe_reactions()


# First step
bio.reaction_step(0.1)
bio.describe_state()
"""
1 bins and 2 species:
 [[15.]
 [35.]]
"""


# Numerous more steps
for i in range(20):
    bio.reaction_step(0.1)

bio.describe_state()

# As expected by the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:
"""
1 bins and 2 species:
 [[14.54545455]
 [36.36363636]]
"""