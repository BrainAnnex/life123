"""
A one-bin A <-> 2C + D, with 1st-order kinetics for each species,
taken to equilibrium
"""

from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio

from modules.html_log.html_log import HtmlLog as log



# Initialize the system
chem_data = chem(diffusion_rates=[0.1, 0.1, 0.1], names=["A", "C", "D"])    # NOTE: diffusion_rates not used for now

rxn = Reactions(chem_data)

# Reaction A <-> 2C + D , with 1st-order kinetics for each species
rxn.add_reaction(reactants=[("A")], products=[(2, "C") , ("D")],
                 forward_rate=5., reverse_rate=2.)

bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

bio.set_all_uniform_concentrations( [4., 7., 2.] )

bio.describe_state(show_diffusion_rates=True)

rxn.describe_reactions()

# Low-level view of the reactions data
for i in range(rxn.number_of_reactions()):
    print(f"{i}: {rxn.get_reactants(i)} <-> {rxn.get_products(i)}   ; Fwd: {rxn.get_forward_rate(i)} / Back: {rxn.get_reverse_rate(i)}")


#bio.verbose = True

# First step
bio.reaction_step(.2)
bio.describe_state()
"""
1 bins and 3 species:
 [[5.6]
 [3.8]
 [0.4]]
"""
# Note: the above values are quite inaccurate because of the large time step 0.2
#       For example, the last concentration is a wild overshot from the initial 2.0 to the equilibrium value of 1.68941267
#       A more precise calculation with bio.react(time_step=.1, n_steps=2) gives conc_D(0.2) = 2.304
#       An even more precise calculation with bio.react(time_step=.05, n_steps=4) gives conc_D(0.2) = 1.69037202
#       I.e. the system is almost at equilibrium already at t=0.2 !
#       TODO: explore the early dynamics of the system in a separate experiment


# Numerous more steps
bio.react(time_step=.05, n_steps=30)

bio.describe_state()

# Consistent with the 5/2 ratio of forward/reverse rates (and the 1st order reactions),
# the systems settles in the following equilibrium:
"""
1 bins and 3 species:
 [[4.31058733]
 [6.37882534]
 [1.68941267]]
"""

A_eq = bio.bin_concentration(0, 0)
C_eq = bio.bin_concentration(0, 1)
D_eq = bio.bin_concentration(0, 2)
print(f"Ratio of equilibrium concentrations (C_eq * D_eq / A_eq) : {C_eq * D_eq / A_eq}")
print(f"Ratio of forward/reverse rates: {rxn.get_forward_rate(0) / rxn.get_reverse_rate(0)}")
