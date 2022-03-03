"""
A reaction between 2 species with initial uniform concentrations,
taken to equilibrium.  Diffusion (non-applicable) disregarded
"""

from life_1D.bio_sim_1d import BioSim1D as bio
from modules.reactions.reactions import Reactions

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
bio.initialize_universe(n_species=2, n_bins=3)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)

bio.set_diffusion_rates([0.1, 0.2])
bio.set_names(["A", "B"])

bio.describe_state(show_diffusion_rates=True)


rxn = Reactions()

rxn.add_reaction(reactants=[(0, 1)], products=[(1, 1)], forward_rate=3., back_rate=2.)  # A -> B

bio.all_reactions = rxn


print("Number of reactions: ", rxn.number_of_reactions())

for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} -> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_back_rate(i)}")
    print(f"{i}: {rxn.get_reactant_names(i)} -> {rxn.get_product_names(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_back_rate(i)}")



bio.reaction_step(0.1)
bio.describe_state()

for i in range(10):
    bio.reaction_step(0.1)

bio.describe_state()