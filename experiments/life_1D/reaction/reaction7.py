"""     TODO: verify
A one-bin 2A <-> B reaction, comparing 1st-order and 2nd-order kinetics in forward direction;
reverse direction 1-st order
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(diffusion_rates=[0.1, 0.1], names=["A", "B"])  # NOTE: diffusion_rates not used for now

rxn = Reactions(chem_data)

# Reaction  2A <-> B , for now with 1st-order kinetics in both directions
rxn.add_reaction(reactants=[(2, "A")], products=["B"], forward_rate=5., reverse_rate=2.)

bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_all_uniform_concentrations( [3., 5.] )

bio.describe_state(show_diffusion_rates=True)

rxn.describe_reactions()

# Low-level view of the reactions data
for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")


# First step
bio.reaction_step(0.02)
bio.describe_state()
"""
1 bins and 2 species:
 [[2.8]
 [5.1]]
"""


# Numerous more steps
bio.react(time_step=0.02, n_steps=20)

bio.describe_state()

# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:
"""
1 bins and 2 species:
 [[2.16928427]
 [5.41535786]]
"""

A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
print(f"Ratio of equilibrium concentrations: {B_eq / A_eq}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")


# Start over again
print("\nSTARTING OVER, this time with 2nd-order kinetics in the forward reaction")
rxn.clear_reactions()

# Reaction  2A <-> B , now with 2nd-order kinetics in the forward directions
rxn.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)

bio.set_all_uniform_concentrations( [3., 5.] )

bio.describe_state(show_diffusion_rates=True)

rxn.describe_reactions()
# Low-level view of the reactions data
for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")


# First step
bio.reaction_step(0.02)
bio.describe_state()
"""
1 bins and 2 species:
 [[1.6]
 [5.7]]
"""

# Numerous more steps
bio.react(time_step=0.02, n_steps=20)

bio.describe_state()
"""
1 bins and 2 species:
 [[1.51554944]
 [5.74222528]]
"""