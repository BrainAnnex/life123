"""
Two coupled reactions:  A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species,
taken to equilibrium
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(diffusion_rates=[.1, .1, .1, .1, .1], names=["A", "B", "C", "D", "E"])  # NOTE: diffusion_rates not used for now

rxn = Reactions(chem_data)

# Reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species
rxn.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
rxn.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)

bio.initialize_system(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_all_uniform_concentrations( [3., 5., 1., 0.4, 0.1] )

bio.describe_state(show_diffusion_rates=True)

rxn.describe_reactions()

# Low-level view of the reactions data
for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")


# First step
bio.reaction_step(0.01)
bio.describe_state()
"""
1 bins and 5 species:
 [[2.27 ]
 [4.27 ]
 [1.702]
 [0.372]
 [0.128]]
"""

# 2nd step
bio.reaction_step(0.01)
bio.describe_state()
"""
1 bins and 5 species:
 [[1.819395  ]
 [3.819395  ]
 [2.10707348]
 [0.32646848]
 [0.17353152]]
"""

# Numerous more steps
bio.react(time_step=0.01, n_steps=200)
bio.describe_state()
"""
1 bins and 5 species:
 [[0.50508029]
 [2.50508029]
 [3.16316668]
 [0.06824696]
 [0.43175304]]
"""


# Do a consistent check with the equilibrium concentrations:

A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
D_eq = bio.bin_concentration(0, 3)
E_eq = bio.bin_concentration(0, 4)

Rf0 = rxn.get_forward_rate(0)
Rb0 = rxn.get_reverse_rate(0)

Rf1 = rxn.get_forward_rate(1)
Rb1 = rxn.get_reverse_rate(1)

equil = -(Rf0 * A_eq * B_eq - Rf1 * C_eq * D_eq) + (Rb0 * C_eq - Rb1 * E_eq)

print("\nAt equilibrium: ", equil, " (this should be close to 0 at equilibrium)")
