"""
A reaction between 2 species with initial uniform concentrations,
taken to equilibrium.  Diffusion (non-applicable) disregarded
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(diffusion_rates=[0.1, 0.2], names=["A", "B"])
bio.initialize_universe(n_bins=3, chem_data=chem_data)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)

bio.describe_state(show_diffusion_rates=True)


rxn = Reactions(chem_data)

# Reaction A -> B , with 1st-order kinetics in both directions
rxn.add_reaction(reactants=[(1,0,1)], products=[(1,1,1)], forward_rate=3., reverse_rate=2.)

bio.all_reactions = rxn


print("Number of reactions: ", rxn.number_of_reactions())

for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")

rxn.describe_reactions()


bio.reaction_step(0.1)
bio.describe_state()

for i in range(10):
    bio.reaction_step(0.1)

bio.describe_state()