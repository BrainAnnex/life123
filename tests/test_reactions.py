import pytest
import numpy as np
from modules.chemicals.chemicals import Chemicals as chem
from modules.reactions.reactions import Reactions
from life_1D.bio_sim_1d import BioSim1D as bio



# Provide an instantiated object that can be used by the various tests that need it
@pytest.fixture(scope="module")
def rxn():
    chem_data = chem()
    chem_data.set_names(["A", "B", "C", "D", "E", "F"])

    bio.initialize_universe(n_bins=10, chem_data=chem_data)

    rnx_obj = Reactions(chem_data)
    yield rnx_obj



def test_chem_data(rxn):
    assert rxn.chem_data.names == ["A", "B", "C", "D", "E", "F"]
    assert rxn.chem_data.name_dict == {"A": 0, "B": 1, "C": 2, "D": 3, "E": 4, "F": 5}



def test_parse_reaction_term(rxn):
    assert rxn._parse_reaction_term(5) == (1, 5, 1)
    assert rxn._parse_reaction_term("F") == (1, 5, 1)

    assert rxn._parse_reaction_term( (3, 5) ) == (3, 5, 1)
    assert rxn._parse_reaction_term( (3, "F") ) == (3, 5, 1)

    assert rxn._parse_reaction_term( (3, 5, 2) ) == (3, 5, 2)
    assert rxn._parse_reaction_term( (3, "F", 2) ) == (3, 5, 2)



def test_add_reaction(rxn):
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert rxn.number_of_reactions() == 1
    assert rxn.get_reaction(0) == {'reactants': [(1, 0, 1)], 'products': [(1, 1, 1)], 'Rf': 3.0, 'Rb': 2.0}

    assert rxn.get_reactants(0) == [(1, 0, 1)]
    assert rxn.get_reactants_formula(0) == "A"

    assert rxn.get_products(0) == [(1, 1, 1)]
    assert rxn.get_products_formula(0) == "B"

    assert rxn.get_forward_rate(0) == 3.
    assert rxn.get_reverse_rate(0) == 2.

    # A 2nd reaction
    rxn.add_reaction(reactants=[(2, "B")], products=[(5, "C")], forward_rate=9., reverse_rate=7.)

    assert rxn.number_of_reactions() == 2
    assert rxn.get_reaction(0) == {'reactants': [(1, 0, 1)], 'products': [(1, 1, 1)], 'Rf': 3.0, 'Rb': 2.0}
    assert rxn.get_reaction(1) == {'reactants': [(2, 1, 1)], 'products': [(5, 2, 1)], 'Rf': 9.0, 'Rb': 7.0}

    # A 3rd reaction
    rxn.add_reaction(reactants=[(2, "D", 3)], products=[(1, "C", 2)], forward_rate=11., reverse_rate=13.)
    assert rxn.number_of_reactions() == 3
    assert rxn.get_reaction(2) == {'reactants': [(2, 3, 3)], 'products': [(1, 2, 2)], 'Rf': 11.0, 'Rb': 13.0}

    # A multi-term reaction
    rxn.add_reaction(reactants=["A", (2, "B")], products=[(3, "C", 2), "D"], forward_rate=5., reverse_rate=1.)
    assert rxn.number_of_reactions() == 4
    assert rxn.get_reaction(3) == {'reactants': [(1, 0, 1), (2, 1, 1)], 'products': [(3, 2, 2), (1, 3, 1)], 'Rf': 5.0, 'Rb': 1.0}

    rxn_list = rxn.describe_reactions()
    assert rxn_list[0] == '0: A <-> B  (Rf = 3.0 / Rb = 2.0)'
    assert rxn_list[1] == '1: 2 B <-> 5 C  (Rf = 9.0 / Rb = 7.0)'
    assert rxn_list[2] == '2: 2 D <-> C  (Rf = 11.0 / Rb = 13.0)'
    assert rxn_list[3] == '3: A + 2 B <-> 3 C + D  (Rf = 5.0 / Rb = 1.0)'



def test_reaction_step_1(rxn):
    # Based on experiment "reaction1"
    chem_data = chem(diffusion_rates=[0.1, 0.1], names=["A", "B"])   # NOTE: diffusion_rates not used
    bio.initialize_universe(n_bins=3, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)

    rxn = Reactions(chem_data)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    bio.set_reactions(rxn)

    assert rxn.number_of_reactions() == 1
    assert np.allclose(bio.univ, [[ 10., 10., 10.] , [50., 50., 50.]])

    # Just the first step, with the lower-level call to reaction_step()
    bio.reaction_step(0.1)
    assert np.allclose(bio.delta_reactions, [[ 7., 7., 7.] , [-7., -7., -7.]])



def test_reaction_step_1b(rxn):
    # Based on experiment "reaction1"
    chem_data = chem(diffusion_rates=[0.1, 0.1], names=["A", "B"])   # NOTE: diffusion_rates not used
    bio.initialize_universe(n_bins=3, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)

    rxn = Reactions(chem_data)

    # Reaction A -> B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    bio.set_reactions(rxn)

    assert rxn.number_of_reactions() == 1
    assert rxn.describe_reactions() == ["0: A <-> B  (Rf = 3.0 / Rb = 2.0)"]
    assert np.allclose(bio.univ, [[ 10., 10., 10.] , [50., 50., 50.]])

    # First step
    bio.react(time_step=0.1, n_steps=1)
    assert np.allclose(bio.delta_reactions, [[ 7., 7., 7.] , [-7., -7., -7.]])
    assert np.allclose(bio.univ, [[ 17., 17., 17.] , [43., 43., 43.]])

    # More steps
    for i in range(10):
        bio.reaction_step(0.1)
        bio.univ += bio.delta_reactions     # Matrix operation

    assert np.allclose(bio.univ,
                       [[ 23.99316406, 23.99316406, 23.99316406] ,
                        [ 36.00683594, 36.00683594, 36.00683594]])

    # Taken to equilibrium
    bio.react(time_step=0.1, n_steps=20)
    assert np.allclose(bio.univ, [[ 24., 24., 24.] , [36., 36., 36.]])



def test_react_2(rxn):
    # Based on experiment "reaction2"
    chem_data = chem(diffusion_rates=[0.1, 0.1], names=["A", "B"])   # NOTE: diffusion_rates not used

    rxn = Reactions(chem_data)

    # Reaction A <-> 3B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    assert rxn.number_of_reactions() == 1

    bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)

    # Large number of steps
    bio.react(time_step=0.1, n_steps=15)
    assert np.allclose(bio.univ, [[14.54545455] , [36.36363636]])
    assert bio.n_bins == 1



def test_react_3(rxn):
    # Based on experiment "reaction3"
    chem_data = chem(diffusion_rates=[0.1, 0.1], names=["A", "B"])   # NOTE: diffusion_rates not used

    rxn = Reactions(chem_data)

    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions
    rxn.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    assert rxn.number_of_reactions() == 1

    bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)

    # First step
    bio.react(time_step=0.1, n_steps=1)
    assert np.allclose(bio.univ, [[20.] , [35.]])

    # Large number of steps
    bio.react(time_step=0.1, n_steps=100)
    assert np.allclose(bio.univ, [[16.25] , [40.625]])
    assert bio.n_bins == 1



def test_react_4(rxn):
    # Based on experiment "reaction4"
    chem_data = chem(diffusion_rates=[0.1, 0.1, 0.1], names=["A", "B", "C"])   # NOTE: diffusion_rates not used

    rxn = Reactions(chem_data)

    # Reaction A + B <-> C , with 1st-order kinetics for each species
    rxn.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                     forward_rate=5., reverse_rate=2.)
    assert rxn.number_of_reactions() == 1

    bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)
    bio.set_uniform_concentration(species_index=2, conc=20.)

    # First step
    bio.react(time_step=0.002, n_steps=1)
    assert np.allclose(bio.univ, [[ 5.08], [45.08], [24.92]])

    # Large number of steps
    bio.react(time_step=0.002, n_steps=40)
    assert np.allclose(bio.univ, [[ 0.29487741], [40.29487741], [29.70512259]])
    assert bio.n_bins == 1


def test_react_5(rxn):
    # Based on experiment "reaction5"
    chem_data = chem(diffusion_rates=[0.1, 0.1, 0.1], names=["A", "C", "D"])   # NOTE: diffusion_rates not used

    rxn = Reactions(chem_data)

    # Reaction A <-> 2C + D , with 1st-order kinetics for each species
    rxn.add_reaction(reactants=[("A")], products=[(2, "C") , ("D")],
                     forward_rate=5., reverse_rate=2.)
    assert rxn.number_of_reactions() == 1

    bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

    bio.set_all_uniform_concentrations( [4., 7., 2.] )

    # First step
    bio.react(time_step=.2, n_steps=1)
    assert np.allclose(bio.univ, [[5.6], [3.8], [0.4]])

    # Numerous more steps
    bio.react(time_step=.05, n_steps=30)
    assert np.allclose(bio.univ, [[4.31058733], [6.37882534], [1.68941267]])
    assert bio.n_bins == 1



def test_react_6(rxn):
    # Based on experiment "reaction6"
    chem_data = chem(diffusion_rates=[0.1, 0.1, 0.1, 0.1], names=["A", "B", "C", "D"])   # NOTE: diffusion_rates not used

    rxn = Reactions(chem_data)

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    rxn.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                     forward_rate=5., reverse_rate=2.)
    assert rxn.number_of_reactions() == 1

    bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

    bio.set_all_uniform_concentrations( [4., 7., 5., 2.] )

    # First step
    bio.react(time_step=0.001, n_steps=1)
    assert np.allclose(bio.univ, [[3.76],
                                  [6.4 ],
                                  [5.48],
                                  [2.36]])

    # Numerous more steps
    bio.react(time_step=0.001, n_steps=40)
    assert np.allclose(bio.univ, [[2.80284552],
                                  [4.00711381],
                                  [7.39430896],
                                  [3.79573172]])
    assert bio.n_bins == 1



def test_react_7(rxn):
    # Based on experiment "reaction7"
    chem_data = chem(diffusion_rates=[0.1, 0.1], names=["A", "B"])   # NOTE: diffusion_rates not used

    rxn = Reactions(chem_data)

    # Reaction  2A <-> B , with 2nd-order kinetics in forward reaction, and 1st-order in reverse
    rxn.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    assert rxn.number_of_reactions() == 1

    bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

    bio.set_all_uniform_concentrations( [3., 5.] )

    # First step
    bio.react(time_step=0.02, n_steps=1)
    assert np.allclose(bio.univ, [[1.6],
                                  [5.7]])

    # Numerous more steps
    bio.react(time_step=0.02, n_steps=20)
    assert np.allclose(bio.univ, [[1.51554944],
                                  [5.74222528]])
    assert bio.n_bins == 1



def test_react_8(rxn):
    # Based on experiment "reaction8"
    chem_data = chem(diffusion_rates=[.1, .1, .1, .1, .1], names=["A", "B", "C", "D", "E"])   # NOTE: diffusion_rates not used

    rxn = Reactions(chem_data)

    # Reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species
    rxn.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
    rxn.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)
    assert rxn.number_of_reactions() == 2

    bio.initialize_universe(n_bins=1, chem_data=chem_data, reactions=rxn)

    bio.set_all_uniform_concentrations( [3., 5., 1., 0.4, 0.1] )

    # First step
    bio.react(time_step=0.01, n_steps=1)
    assert np.allclose(bio.univ, [[2.27 ],
                                  [4.27 ],
                                  [1.702],
                                  [0.372],
                                  [0.128]])

    # 2nd step
    bio.react(time_step=0.01, n_steps=1)
    assert np.allclose(bio.univ, [[1.819395  ],
                                  [3.819395  ],
                                  [2.10707348],
                                  [0.32646848],
                                  [0.17353152]])

    # Numerous more steps
    bio.react(time_step=0.01, n_steps=200)
    assert np.allclose(bio.univ, [[0.50508029],
                                  [2.50508029],
                                  [3.16316668],
                                  [0.06824696],
                                  [0.43175304]])
    assert bio.n_bins == 1
