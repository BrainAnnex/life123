"""     TODO: verify
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

# Reaction  2A <-> B , for now with 1st-order kinetics in both directions
rxn.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
rxn.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)

bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_all_uniform_concentrations( [3., 5., 1., 0.4, 0.1] )

bio.describe_state(show_diffusion_rates=True)

rxn.describe_reactions()

# Low-level view of the reactions data
for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")


# First step
bio.react(time_step=0.00002, n_steps=1000)
bio.describe_state()
"""
1 bins and 5 species:
 [[1.54 ]
 [3.54 ]
 [2.404]
 [0.344]
 [0.156]]
 
vs.
 [[1.819395  ]
 [3.819395  ]
 [2.10707348]
 [0.32646848]
 [0.17353152]]
 
vs.
 [[1.90738687]
 [3.90738687]
 [2.01639306]
 [0.32377993]
 [0.17622007]]
 
vs.
 [[1.96490304]
 [3.96490304]
 [1.95778584]
 [0.32268888]
 [0.17731112]]
 
vs.
 [[1.97663547]
 [3.97663547]
 [1.94589685]
 [0.32253232]
 [0.17746768]]
 
vs.
 [[1.97765931]
 [3.97765931]
 [1.94486042]
 [0.32251973]
 [0.17748027]]
"""


# Numerous more steps
bio.react(time_step=0.02, n_steps=50)

bio.describe_state()

# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:
"""
1 bins and 5 species:
 [[0.50508029]
 [2.50508029]
 [3.16316668]
 [0.06824696]
 [0.43175304]]
"""

A_eq = bio.bin_concentration(0, 0)
B_eq = bio.bin_concentration(0, 1)
C_eq = bio.bin_concentration(0, 2)
D_eq = bio.bin_concentration(0, 3)
E_eq = bio.bin_concentration(0, 4)
print(f"Ratio of equilibrium concentrations: {B_eq / A_eq}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")