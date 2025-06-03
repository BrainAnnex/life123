# These are tests specifically for reactions in 1D;
# for general tests of 1D system, see "test_biosim_1d.py"

import numpy as np
from life123 import ChemData
from life123 import BioSim1D



def test_initialization():
    chem_data = ChemData(names=["A", "B", "C"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    bio.reactions.add_reaction(reactants=["A", "B"], products="C", forward_rate=8., reverse_rate=2.)

    assert bio.reaction_dynamics.chem_data == chem_data



def test_reaction_step_1():
    # Based on experiment "reaction1"
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    bio.set_uniform_concentration(chem_index=0, conc=10.)
    bio.set_uniform_concentration(chem_index=1, conc=50.)


    # Reaction A <-> B , with 1st-order kinetics in both directions
    bio.reactions.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert bio.reactions.number_of_reactions() == 1
    assert np.allclose(bio.system, [[10., 10., 10.] , [50., 50., 50.]])

    # Just the first reaction step, with the lower-level call to reaction_step()
    bio.reaction_step(0.1)
    assert np.allclose(bio.delta_reactions, [[ 7., 7., 7.] , [-7., -7., -7.]])



def test_reaction_step_1b():
    # Based on experiment "reaction1"
    chem_data = ChemData(names=["A", "B"])
    bio = BioSim1D(n_bins=3, chem_data=chem_data)

    bio.set_uniform_concentration(chem_index=0, conc=10.)
    bio.set_uniform_concentration(chem_index=1, conc=50.)


    # Reaction A <-> B , with 1st-order kinetics in both directions
    bio.reactions.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    assert bio.reactions.number_of_reactions() == 1
    assert bio.reactions.multiple_reactions_describe() == ["0: A <-> B  (kF = 3 / kR = 2 / delta_G = -1,005.1 / K = 1.5) | 1st order in all reactants & products"]
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

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    # Reaction A <-> 3B , with 1st-order kinetics in both directions
    bio.reactions.add_reaction(reactants=["A"], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)
    assert bio.reactions.number_of_reactions() == 1

    bio.set_uniform_concentration(chem_index=0, conc=10.)
    bio.set_uniform_concentration(chem_index=1, conc=50.)

    # Large number of steps
    bio.react(time_step=0.1, n_steps=15)
    assert np.allclose(bio.system, [[14.54545455] , [36.36363636]])
    assert bio.n_bins == 1



def test_react_2():
    # Based on experiment "reaction3"
    chem_data = ChemData(names=["A", "B"])

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions
    bio.reactions.add_reaction(reactants=[(2,"A",1)], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)
    assert bio.reactions.number_of_reactions() == 1

    bio.set_uniform_concentration(chem_index=0, conc=10.)
    bio.set_uniform_concentration(chem_index=1, conc=50.)

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

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    # Reaction A + B <-> C , with 1st-order kinetics for each species
    bio.reactions.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                               forward_rate=5., reverse_rate=2.)
    assert bio.reactions.number_of_reactions() == 1

    bio.set_uniform_concentration(chem_index=0, conc=10.)
    bio.set_uniform_concentration(chem_index=1, conc=50.)
    bio.set_uniform_concentration(chem_index=2, conc=20.)

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

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    # Reaction A <-> 2C + D , with 1st-order kinetics for each species
    bio.reactions.add_reaction(reactants=[("A")], products=[(2, "C",1) , ("D")],
                           forward_rate=5., reverse_rate=2.)
    assert bio.reactions.number_of_reactions() == 1

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

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    bio.reactions.add_reaction(reactants=[(2,"A",1) , (5,"B",1)], products=[(4,"C",1) , (3,"D",1)],
                           forward_rate=5., reverse_rate=2.)
    assert bio.reactions.number_of_reactions() == 1

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

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    # Reaction  2A <-> B , with 2nd-order kinetics in forward reaction, and 1st-order in reverse
    bio.reactions.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    assert bio.reactions.number_of_reactions() == 1

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

    bio = BioSim1D(n_bins=1, chem_data=chem_data)

    # Coupled reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species
    bio.reactions.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
    bio.reactions.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)
    assert bio.reactions.number_of_reactions() == 2

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
