import pytest
import numpy as np
from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio



# Provide an instantiated object that can be used by the various tests that need it
@pytest.fixture(scope="module")
def rxn():
    pass
    #chem_data = chem(names=["A", "B", "C", "D", "E", "F"])
    #bio.initialize_system(n_bins=10, chem_data=chem_data)
    #rnx_obj = Reactions(chem_data)
    #yield rnx_obj



def test_single_bin_reaction_step_1(rxn):
    chem_data = chem(names=["A", "B"])

    rxn = Reactions(chem_data)

    conc_dict = {0: 10., 1: 50.}

    # Reaction A <-> B , with 1st-order kinetics in both directions.  Based on experiment "reaction1"
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=0.1)
    assert np.allclose(result, [ 7. , -7.])
    assert result[0] == - result[1]   # From the stoichiometry


    rxn.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , with 1st-order kinetics in both directions.  Based on experiment "reaction2"
    rxn.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=0.1)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]   # From the stoichiometry


    rxn.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions.  Based on experiment "reaction3"
    rxn.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=0.1)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test_single_bin_reaction_step_2(rxn):
    chem_data = chem(names=["A", "B", "C"])

    rxn = Reactions(chem_data)

    conc_dict = {0: 10., 1: 50., 2: 20.}

    # Reaction A <-> B , with 1st-order kinetics in both directions.  Based on experiment "reaction1"
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=0.1)
    assert np.allclose(result, [ 7. , -7. , 0.])
    assert result[0] == - result[1]   # From the stoichiometry


    rxn.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction A + B <-> C , with 1st-order kinetics for each species.  Based on experiment "reaction4"
    rxn.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                     forward_rate=5., reverse_rate=2.)

    result = rxn.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=0.002)
    assert np.allclose(result, [-4.92, -4.92, 4.92])
    assert result[0] == result[1]    # From the stoichiometry
    assert result[1] == - result[2]   # From the stoichiometry



def test_single_bin_reaction_step_3(rxn):
    chem_data = chem(names=["A", "C", "D"])

    rxn = Reactions(chem_data)

    conc_dict = {0: 4., 1: 7., 2: 2.}

    # Reaction A <-> 2C + D , with 1st-order kinetics for each species.  Based on experiment "reaction5"
    rxn.add_reaction(reactants=[("A")], products=[(2, "C") , ("D")],
                     forward_rate=5., reverse_rate=2.)

    result = rxn.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=0.05)
    assert np.allclose(result, [0.4 , -0.8 , -0.4])
    assert result[0] == - result[1] /2    # From the stoichiometry
    assert result[0] == - result[2]       # From the stoichiometry



def test_single_bin_reaction_step_4(rxn):
    chem_data = chem(names=["A", "B", "C", "D"])

    rxn = Reactions(chem_data)

    conc_dict = {0: 4., 1: 7., 2: 5., 3: 2.}

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species.  Based on experiment "reaction6"
    rxn.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                     forward_rate=5., reverse_rate=2.)

    result = rxn.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=0.001)
    assert np.allclose(result, [-0.24 , -0.6 , 0.48, 0.36])
    assert  np.allclose(result[0] /2 , result[1] /5)    # From the stoichiometry
    assert  np.allclose(result[1] /5 , -result[2] /4)   # From the stoichiometry
    assert  np.allclose(result[2] /4 , result[3] /3)    # From the stoichiometry



def test_single_bin_reaction_step_5(rxn):
    chem_data = chem(names=["A", "B"])

    rxn = Reactions(chem_data)

    conc_dict = {0: 3., 1: 5.}

    # Reaction  2A <-> B , with 2nd-order kinetics in forward reaction, and 1st-order in reverse.  Based on experiment "reaction7"
    rxn.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)

    result = rxn.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=0.02)
    assert np.allclose(result, [-1.4 , 0.7])
    assert  np.allclose(result[0] /2 , -result[1])    # From the stoichiometry



def test_single_bin_reaction_step_6(rxn):
    chem_data = chem(names=["A", "B", "C", "D", "E"])

    rxn = Reactions(chem_data)

    conc_dict = {0: 3., 1: 5., 2: 1., 3: 0.4, 4: 0.1}

    # Coupled reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species.  Based on experiment "reaction8"
    rxn.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
    rxn.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)
    assert rxn.number_of_reactions() == 2

    result = rxn.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=0.02)
    assert np.allclose(result, [-1.46 , -1.46  , 1.404 , -0.056 ,  0.056])
    assert  np.allclose(result[0] , result[1])              # From the stoichiometry
    assert  np.allclose(result[3] , -result[4])             # From the stoichiometry
    assert  np.allclose(result[0] + result[4], -result[2])  # From the stoichiometry
                                                            # The increase in [A] and [E] combined
                                                            # must match the decrease in [C]


def test_compute_all_rate_deltas():
    chem_data = chem(names=["A", "B", "C", "D"])
    rxn = Reactions(chem_data)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=20., reverse_rate=2.)
    conc_dict = {0: 5., 1: 8.}
    result = rxn.compute_all_rate_deltas(conc_dict=conc_dict, delta_time=0.5)
    assert np.allclose(result, [42.0])

    # Reaction 2B <-> 3C , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=[(2, "B")], products=[(3, "C")], forward_rate=10., reverse_rate=25.)
    conc_dict = {0: 5., 1: 8., 2: 15.}
    result = rxn.compute_all_rate_deltas(conc_dict=conc_dict, delta_time=0.5)
    assert np.allclose(result, [42.0, -147.5])

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    rxn.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                     forward_rate=5., reverse_rate=2.)
    conc_dict = {0: 5., 1: 8., 2: 15., 3: 7.}
    result = rxn.compute_all_rate_deltas(conc_dict=conc_dict, delta_time=0.5)
    assert np.allclose(result, [42.0, -147.5, -5.0])

    # Reaction  2A <-> B , with 2nd-order kinetics in the forward direction
    rxn.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=3., reverse_rate=2.)
    conc_dict = {0: 5., 1: 8., 2: 15., 3: 7.}
    result = rxn.compute_all_rate_deltas(conc_dict=conc_dict, delta_time=0.5)
    assert np.allclose(result, [42.0, -147.5, -5.0, 29.5])

    # FLUSH OUT ALL REACTIONS
    rxn.clear_reactions()
    # Reaction A <-> B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=20., reverse_rate=2.)
    conc_dict = {0: 5., 1: 8.}
    result = rxn.compute_all_rate_deltas(conc_dict=conc_dict, delta_time=0.25)
    assert np.allclose(result, [21.0])

    #print(result)
    """
    bio.initialize_system(n_bins=1, chem_data=chem_data, reactions=rxn)
    bio.set_all_uniform_concentrations([5., 8., 0.])
    old_result = bio.compute_rates(bin_n=0, delta_time=0.5, number_reactions=1)
    print(old_result)
    """



def test_compute_rate_diff():
    chem_data = chem(names=["A", "B", "C", "D"])
    rxn = Reactions(chem_data)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=20., reverse_rate=2.)
    conc_dict = {0: 5., 1: 8.}
    result = rxn.compute_rate_delta(rxn_index=0, conc_dict=conc_dict, delta_time=0.5)
    assert np.allclose(result, 42.0)

    # Reaction 2B <-> 3C , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=[(2, "B")], products=[(3, "C")], forward_rate=10., reverse_rate=25.)
    conc_dict = {1: 8., 2: 15.}
    result = rxn.compute_rate_delta(rxn_index=1, conc_dict=conc_dict, delta_time=1.5)
    assert np.allclose(result, -442.5)

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    rxn.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                    forward_rate=5., reverse_rate=2.)
    conc_dict = {0: 3.5, 1: 9., 2: 11., 3: 7.}
    result = rxn.compute_rate_delta(rxn_index=2, conc_dict=conc_dict, delta_time=0.5)
    assert np.allclose(result, 1.75)

    # Reaction  2A <-> B , with 2nd-order kinetics in the forward direction
    rxn.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    conc_dict = {0: 4.5, 1: 6.}
    result = rxn.compute_rate_delta(rxn_index=3, conc_dict=conc_dict, delta_time=2.0)
    assert np.allclose(result, 178.5)

    # Reaction  B <-> 2C , with 2nd-order kinetics in the reverse direction
    rxn.add_reaction(reactants=[("B")], products=[(2, "C", 2)], forward_rate=4., reverse_rate=2.)
    conc_dict = {1: 5., 2: 4.}
    result = rxn.compute_rate_delta(rxn_index=4, conc_dict=conc_dict, delta_time=0.5)
    assert np.allclose(result, -6.0)



def test_is_in_equilibrium():
    chem_data = chem(names=["A", "B", "C", "D", "E", "F"])

    rxn = Reactions(chem_data)

    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    c = {'A': 23.9931640625, 'B': 36.0068359375}
    assert rxn.is_in_equilibrium(rxn_index = 0, conc = c)

    rxn.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    c = {'A': 14.54545455, 'B': 36.36363636}
    assert rxn.is_in_equilibrium(rxn_index = 1, conc = c)

    rxn.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    c = {'A': 16.25, 'B': 40.625}
    assert rxn.is_in_equilibrium(rxn_index = 2, conc = c)

    # Reaction A + B <-> C , with 1st-order kinetics for each species
    rxn.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 0.29487741, 'B': 40.29487741, 'C': 29.70512259}
    assert rxn.is_in_equilibrium(rxn_index = 3, conc = c)

    # Reaction A <-> 2C + D , with 1st-order kinetics for each species
    rxn.add_reaction(reactants=[("A")], products=[(2, "C") , ("D")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 4.31058733, 'C': 6.37882534, 'D': 1.68941267}
    assert rxn.is_in_equilibrium(rxn_index = 4, conc = c)

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    rxn.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 2.80284552, 'B': 4.00711381, 'C': 7.39430896, 'D': 3.79573172}
    assert rxn.is_in_equilibrium(rxn_index = 5, conc = c)


    rxn.clear_reactions()   # This will reset the reaction count to 0

    # Reaction  2A <-> B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=[(2, "A")], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 2.16928427, 'B': 5.41535786}
    assert rxn.is_in_equilibrium(rxn_index = 0, conc = c)

    # Reaction  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction
    rxn.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 1.51554944, 'B': 5.74222528}
    assert rxn.is_in_equilibrium(rxn_index = 1, conc = c)
