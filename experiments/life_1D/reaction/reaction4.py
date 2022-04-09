"""
A one-bin A + B <-> C, with 1st-order kinetics for each species,
taken to equilibrium
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(diffusion_rates=[0.1, 0.1, 0.1], names=["A", "B", "C"])    # NOTE: diffusion_rates not used for now

rxn = Reactions(chem_data)

# Reaction A + B <-> C , with 1st-order kinetics for each species
rxn.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                 forward_rate=5., reverse_rate=2.)

bio.initialize_system(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_uniform_concentration(species_index=0, conc=10.)
bio.set_uniform_concentration(species_index=1, conc=50.)
bio.set_uniform_concentration(species_index=2, conc=20.)

bio.describe_state(show_diffusion_rates=True)

rxn.describe_reactions()

# Low-level view of the reactions data
for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")


#bio.verbose = True

# First step
bio.reaction_step(0.002)
bio.describe_state()
"""
1 bins and 3 species:
 [[ 5.08]
 [45.08]
 [24.92]]
"""


# Numerous more steps
bio.react(time_step=0.002, n_steps=40)

bio.describe_state()

# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:
"""
1 bins and 3 species:
 [[ 0.29487741]
 [40.29487741]
 [29.70512259]]
"""
# Note: "A" is largely the limiting reagent

A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
print(f"Ratio of equilibrium concentrations (C_eq / (A_eq * B_eq)) : {C_eq / (A_eq * B_eq)}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")
