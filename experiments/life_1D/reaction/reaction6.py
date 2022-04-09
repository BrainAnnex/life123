"""
A one-bin 2A + 5B <-> 4C + 3D reaction, with 1st-order kinetics for each species,
taken to equilibrium
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(diffusion_rates=[0.1, 0.1, 0.1, 0.1], names=["A", "B", "C", "D"])    # NOTE: diffusion_rates not used for now

rxn = Reactions(chem_data)

# Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
rxn.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                 forward_rate=5., reverse_rate=2.)

bio.initialize_system(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_all_uniform_concentrations( [4., 7., 5., 2.] )

bio.describe_state(show_diffusion_rates=True)

rxn.describe_reactions()

# Low-level view of the reactions data
for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")


# First step
bio.reaction_step(0.001)
bio.describe_state()
"""
1 bins and 4 species:
 [[3.76]
 [6.4 ]
 [5.48]
 [2.36]]
"""


# Numerous more steps
bio.react(time_step=0.001, n_steps=40)

bio.describe_state()

# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:
"""
1 bins and 4 species:
 [[2.80284552]
 [4.00711381]
 [7.39430896]
 [3.79573172]]
"""

A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
D_eq = bio.bin_concentration(0, 3)
print(f"Ratio of equilibrium concentrations ((C_eq * D_eq) / (A_eq * B_eq)) : {(C_eq * D_eq) / (A_eq * B_eq)}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")
