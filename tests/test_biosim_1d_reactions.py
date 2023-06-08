# These are tests specifically for reactions in 1D;
# for general tests of 1D system, see test_biosim_1d.py

import numpy as np
from src.modules.reactions.reaction_data import ChemData
from src.life_1D.bio_sim_1d import BioSim1D



def test_initialization():
    chem_data = ChemData(names=["A", "B", "C"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    chem_data.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=8., reverse_rate=2.)

    assert bio.reaction_dynamics.reaction_data == chem_data



def test_reaction_step_1():
    # Based on experiment "reaction1"
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)

    #rxn = ReactionDynamics(chem_data)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert chem_data.number_of_reactions() == 1
    assert np.allclose(bio.system, [[10., 10., 10.] , [50., 50., 50.]])

    # Just the first reaction step, with the lower-level call to reaction_step()
    bio.reaction_step(0.1)
    assert np.allclose(bio.delta_reactions, [[ 7., 7., 7.] , [-7., -7., -7.]])



def test_reaction_step_1b():
    # Based on experiment "reaction1"
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)


    # Reaction A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert chem_data.number_of_reactions() == 1
    assert chem_data.multiple_reactions_describe() == ["0: A <-> B  (kF = 3 / kR = 2 / Delta_G = -1,005.13 / K = 1.5) | 1st order in all reactants & products"]
    assert np.allclose(bio.system, [[10., 10., 10.] , [50., 50., 50.]])

    # First step
    bio.react(time_step=0.1, n_steps=1)
    assert np.allclose(bio.delta_reactions, [[ 7., 7., 7.] , [-7., -7., -7.]])
    assert np.allclose(bio.system, [[17., 17., 17.] , [43., 43., 43.]])

    # More steps
    for i in range(10):
        bio.reaction_step(0.1)
        bio.system += bio.delta_reactions     # Matrix operation

    assert np.allclose(bio.system,
                       [[ 23.99316406, 23.99316406, 23.99316406] ,
                        [ 36.00683594, 36.00683594, 36.00683594]])

    # Taken to equilibrium
    bio.react(time_step=0.1, n_steps=20)
    assert np.allclose(bio.system, [[24., 24., 24.] , [36., 36., 36.]])



def test_react_1():
    # Based on experiment "reaction2"
    chem_data = ChemData(names=["A", "B"])

    # Reaction A <-> 3B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    assert chem_data.number_of_reactions() == 1

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)

    # Large number of steps
    bio.react(time_step=0.1, n_steps=15)
    assert np.allclose(bio.system, [[14.54545455] , [36.36363636]])
    assert bio.n_bins == 1



def test_react_2():
    # Based on experiment "reaction3"
    chem_data = ChemData(names=["A", "B"])

    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    assert chem_data.number_of_reactions() == 1

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)

    # First step
    bio.react(time_step=0.1, n_steps=1)
    assert np.allclose(bio.system, [[20.] , [35.]])

    # Large number of steps
    bio.react(time_step=0.1, n_steps=100)
    assert np.allclose(bio.system, [[16.25] , [40.625]])
    assert bio.n_bins == 1



def test_react_3():
    # Based on experiment "reaction4"
    chem_data = ChemData(names=["A", "B", "C"])

    # Reaction A + B <-> C , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                     forward_rate=5., reverse_rate=2.)
    assert chem_data.number_of_reactions() == 1

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    bio.set_uniform_concentration(species_index=0, conc=10.)
    bio.set_uniform_concentration(species_index=1, conc=50.)
    bio.set_uniform_concentration(species_index=2, conc=20.)

    # First step
    bio.react(time_step=0.002, n_steps=1)
    assert np.allclose(bio.system, [[5.08], [45.08], [24.92]])

    # Large number of steps
    bio.react(time_step=0.002, n_steps=40)
    assert np.allclose(bio.system, [[0.29487741], [40.29487741], [29.70512259]])
    assert bio.n_bins == 1



def test_react_4():
    # Based on experiment "reaction5"
    chem_data = ChemData(names=["A", "C", "D"])

    # Reaction A <-> 2C + D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[("A")], products=[(2, "C") , ("D")],
                           forward_rate=5., reverse_rate=2.)
    assert chem_data.number_of_reactions() == 1

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    bio.set_all_uniform_concentrations( [4., 7., 2.] )

    # First step
    bio.react(time_step=.2, n_steps=1)
    assert np.allclose(bio.system, [[5.6], [3.8], [0.4]])

    # Numerous more steps
    bio.react(time_step=.05, n_steps=30)
    assert np.allclose(bio.system, [[4.31058733], [6.37882534], [1.68941267]])
    assert bio.n_bins == 1



def test_react_5():
    # Based on experiment "reaction6"
    chem_data = ChemData(names=["A", "B", "C", "D"])

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                           forward_rate=5., reverse_rate=2.)
    assert chem_data.number_of_reactions() == 1

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    bio.set_all_uniform_concentrations( [4., 7., 5., 2.] )

    # First step
    bio.react(time_step=0.001, n_steps=1)
    assert np.allclose(bio.system, [[3.76],
                                    [6.4 ],
                                    [5.48],
                                    [2.36]])

    # Numerous more steps
    bio.react(time_step=0.001, n_steps=40)
    assert np.allclose(bio.system, [[2.80284552],
                                    [4.00711381],
                                    [7.39430896],
                                    [3.79573172]])
    assert bio.n_bins == 1



def test_react_6():
    # Based on experiment "reaction7"
    chem_data = ChemData(names=["A", "B"])

    # Reaction  2A <-> B , with 2nd-order kinetics in forward reaction, and 1st-order in reverse
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    assert chem_data.number_of_reactions() == 1

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    bio.set_all_uniform_concentrations( [3., 5.] )

    # First step
    bio.react(time_step=0.02, n_steps=1)
    assert np.allclose(bio.system, [[1.6],
                                    [5.7]])

    # Numerous more steps
    bio.react(time_step=0.02, n_steps=20)
    assert np.allclose(bio.system, [[1.51554944],
                                    [5.74222528]])
    assert bio.n_bins == 1



def test_react_7():
    # Based on experiment "reaction8"
    chem_data = ChemData(names=["A", "B", "C", "D", "E"])

    # Coupled reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)
    assert chem_data.number_of_reactions() == 2

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    bio.set_all_uniform_concentrations( [3., 5., 1., 0.4, 0.1] )

    # First step
    bio.react(time_step=0.01, n_steps=1)
    assert np.allclose(bio.system, [[2.27],
                                    [4.27 ],
                                    [1.702],
                                    [0.372],
                                    [0.128]])

    # 2nd step
    bio.react(time_step=0.01, n_steps=1)
    assert np.allclose(bio.system, [[1.819395],
                                    [3.819395  ],
                                    [2.10707348],
                                    [0.32646848],
                                    [0.17353152]])

    # Numerous more steps
    bio.react(time_step=0.01, n_steps=200)
    assert np.allclose(bio.system, [[0.50508029],
                                    [2.50508029],
                                    [3.16316668],
                                    [0.06824696],
                                    [0.43175304]])
    assert bio.n_bins == 1



def test_react_with_membrane():
    # Based on experiment "reaction/membranes_1"
    chem_data = ChemData(names=["A", "B", "C"])     # NOTE: Diffusion not done
    bio = BioSim1D(n_bins=5, chem_data=chem_data)

    bio.set_membranes(membrane_pos=[1])   # A single membrane, passing thru bin 1

    bio.set_all_uniform_concentrations(conc_list=[4., 8., 12.])
    bio.set_bin_conc(bin_address=1, species_name="A", conc=10.)
    bio.set_bin_conc(bin_address=1, species_name="A", conc=55., across_membrane=True)
    bio.set_bin_conc(bin_address=1, species_name="B", conc=20.)
    bio.set_bin_conc(bin_address=1, species_name="C", conc=30., both_sides=True)

    # Make the last bin match all the concentrations of the "post-membrane" section of bin 1
    bio.set_bin_conc(bin_address=4, species_name="A", conc=55.)
    bio.set_bin_conc(bin_address=4, species_name="C", conc=30.)

    assert np.allclose(bio.lookup_species(species_name="A"), [4, 10, 4, 4, 55])
    assert np.allclose(bio.lookup_species(species_name="A", trans_membrane=True), [0, 55, 0, 0, 0])

    assert np.allclose(bio.lookup_species(species_name="B"), [8, 20, 8, 8, 8])
    assert np.allclose(bio.lookup_species(species_name="B", trans_membrane=True), [0, 8, 0, 0, 0])

    assert np.allclose(bio.lookup_species(species_name="C"), [12, 30, 12, 12, 30])
    assert np.allclose(bio.lookup_species(species_name="C", trans_membrane=True), [0, 30, 0, 0, 0])

    # Reaction A + B <-> C , with 1st-order kinetics in both directions, mostly forward
    chem_data.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=8., reverse_rate=2.)

    bio.react(time_step=0.002, n_steps=1)

    bin1_trans = bio.bin_concentration(bin_address=1, species_name="A", trans_membrane=True)
    bin4 = bio.bin_concentration(bin_address=4, species_name="A")
    assert np.allclose(bin1_trans, bin4)
    assert np.allclose(bin4, 48.08)

    bin1_trans = bio.bin_concentration(bin_address=1, species_name="B", trans_membrane=True)
    bin4 = bio.bin_concentration(bin_address=4, species_name="B")
    assert np.allclose(bin1_trans, bin4)
    assert np.allclose(bin4, 1.08)

    bin1_trans = bio.bin_concentration(bin_address=1, species_name="C", trans_membrane=True)
    bin4 = bio.bin_concentration(bin_address=4, species_name="C")
    assert np.allclose(bin1_trans, bin4)
    assert np.allclose(bin4, 36.92)
