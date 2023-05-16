import pytest
import numpy as np
from src.modules.reactions.reaction_data import ReactionData
from src.modules.reactions.reaction_dynamics import ReactionDynamics



def test_set_conc():
    #TODO: test the snapshot argument
    chem_data = ReactionData(names=["A", "B", "C"])
    rxn = ReactionDynamics(chem_data)

    with pytest.raises(Exception):
        rxn.set_conc(conc=[1, 2, 3, 4])         # Wrong number of entries
    with pytest.raises(Exception):
        rxn.set_conc(conc=[1., -2., -3.])       # Negative values
    with pytest.raises(Exception):
        rxn.set_conc(conc=(10., 20., -0.01))    # Negative values

    rxn.set_conc(conc=[1., 2., 3.])
    assert np.allclose(rxn.system, [1., 2., 3.])

    rxn.set_conc(conc=(10., 20., 30.))
    assert np.allclose(rxn.system, [10., 20., 30.])


def test_get_conc():
    chem_data = ReactionData(names=["A", "B", "C"])
    rxn = ReactionDynamics(chem_data)
    rxn.set_conc(conc=(10., 20., 30.))

    result = rxn.get_system_conc()
    assert np.allclose(result, [10., 20., 30.])


def test_get_conc_dict():
    chem_data = ReactionData(names=["A", "B", "C", "D"])
    rxn = ReactionDynamics(chem_data)
    rxn.set_conc(conc=(100, 200, 300, 400))

    result = rxn.get_conc_dict()
    assert result == {"A": 100, "B": 200, "C": 300, "D": 400}

    result = rxn.get_conc_dict(species=["D", "A"])
    assert result == {"A": 100, "D": 400}

    result = rxn.get_conc_dict(species=("C",))      # Tuple with 1 element
    assert result == {"C": 300}

    with pytest.raises(Exception):
        rxn.get_conc_dict(species="C")                      # Wrong data type
        rxn.get_conc_dict(system_data=np.array([1, 2]))     # Wrong number of entries

    result = rxn.get_conc_dict(system_data=np.array([1, 2, 3, 4]))
    assert result == {"A": 1, "B": 2, "C": 3, "D": 4}

    result = rxn.get_conc_dict(species=["B"], system_data=np.array([1, 2, 3, 4]))
    assert result == {"B": 2}



def test_specify_steps():
    rxn = ReactionDynamics(None)

    with pytest.raises(Exception):
        # Too few arguments
        rxn.specify_steps()
        rxn.specify_steps(total_duration=15.1)
        rxn.specify_steps(time_step=0.2)
        rxn.specify_steps(n_steps=30)
        # Too many arguments
        rxn.specify_steps(total_duration=15.1, time_step=0.2, n_steps=30)

    assert rxn.specify_steps(time_step=0.5, n_steps=24) == (0.5, 24)
    assert rxn.specify_steps(total_duration=12.0, time_step=0.5) == (0.5, 24)
    assert rxn.specify_steps(total_duration=12.0, n_steps=24) == (0.5, 24)



def test_single_compartment_react():

    # Test based on experiment "cycles_1"
    chem_data = ReactionData(names=["A", "B", "C", "E_high", "E_low"])

    # Reaction A <-> B, mostly in forward direction (favored energetically)
    chem_data.add_reaction(reactants="A", products="B",
                           forward_rate=9., reverse_rate=3.)

    # Reaction B <-> C, also favored energetically
    chem_data.add_reaction(reactants="B", products="C",
                           forward_rate=8., reverse_rate=4.)

    # Reaction C + E_High <-> A + E_Low, also favored energetically, but kinetically slow
    chem_data.add_reaction(reactants=["C" , "E_high"], products=["A", "E_low"],
                           forward_rate=1., reverse_rate=0.2)

    initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}

    dynamics = ReactionDynamics(reaction_data=chem_data)
    dynamics.set_conc(conc=initial_conc, snapshot=False)

    dynamics.set_diagnostics()
    dynamics.single_compartment_react(time_step=0.0010, stop_time=0.0035)

    run1 = dynamics.get_system_conc()
    assert np.allclose(run1, [9.69124339e+01, 3.06982519e+00, 1.77408783e-02, 9.99985757e+02, 1.42431633e-02])
    assert np.allclose(dynamics.system_time, 0.0035)
    assert dynamics.explain_time_advance(return_times=True, silent=True) == \
               ([0.0, 0.002, 0.0025, 0.0035],
                [0.001, 0.0005, 0.001])

    # The above computation automatically took 2 normal steps, 1 half-size step and 1 normal step;
    # now repeat the process manually

    dynamics2 = ReactionDynamics(reaction_data=chem_data)
    dynamics2.set_conc(conc=initial_conc, snapshot=False)

    dynamics2.set_diagnostics()
    dynamics2.single_compartment_react(time_step=0.0010, n_steps=2)
    dynamics2.single_compartment_react(time_step=0.0005, n_steps=1)
    dynamics2.single_compartment_react(time_step=0.0010, n_steps=1)
    run2 = dynamics.get_system_conc()
    assert np.allclose(run2, run1)      # Same result as before
    assert np.allclose(dynamics2.system_time, 0.0035)
    assert dynamics2.explain_time_advance(return_times=True) == \
           ([0.0, 0.002, 0.0025, 0.0035],
            [0.001, 0.0005, 0.001])
    # The time advance is now different, because multiple calls to single_compartment_react() break up the counting
    # (just a convention)



def test_single_reaction_fixed_step():
    chem_data = ReactionData(names=["A", "B"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([10., 50.])

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [ 7. , -7.])
    assert result[0] == - result[1]         # From the stoichiometry


    rxn.clear_reactions()       # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction2"
    chem_data.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]      # From the stoichiometry


    rxn.clear_reactions()       # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction3"
    chem_data.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test_single_reaction_variable_step_1():
    chem_data = ReactionData(names=["A", "B"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([10., 50.])

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [ 7. , -7.])
    assert result[0] == - result[1]         # From the stoichiometry


    rxn.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction2"
    chem_data.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]      # From the stoichiometry


    rxn.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction3"
    chem_data.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test_single_reaction_step_2():
    chem_data = ReactionData(names=["A", "B", "C"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([10., 50., 20.])

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # # Based on experiment "reactions_single_compartment/react_1"
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [ 7. , -7. , 0.])    # Chemical "C" not participating in this reaction; its delta conc. is 0
    assert result[0] == - result[1]         # From the stoichiometry


    chem_data.clear_reactions_data()   # Re-start with a blank slate of reactions
    # Reaction A + B <-> C , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction4"
    chem_data.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                     forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.002, conc_array=conc_array)
    assert np.allclose(result, [-4.92, -4.92, 4.92])
    assert result[0] == result[1]           # From the stoichiometry
    assert result[1] == - result[2]         # From the stoichiometry



def test_single_reaction_step_3():
    chem_data = ReactionData(names=["A", "C", "D"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([4., 7., 2.])

    # Reaction A <-> 2C + D , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction5"
    chem_data.add_reaction(reactants=[("A")], products=[(2, "C") , ("D")],
                     forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.05, conc_array=conc_array)
    assert np.allclose(result, [0.4 , -0.8 , -0.4])
    assert result[0] == - result[1] /2    # From the stoichiometry
    assert result[0] == - result[2]       # From the stoichiometry



def test_single_reaction_step_4():
    chem_data = ReactionData(names=["A", "B", "C", "D"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([4., 7., 5., 2.])

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction6"
    chem_data.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                     forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.001, conc_array=conc_array)
    assert np.allclose(result, [-0.24 , -0.6 , 0.48, 0.36])
    assert  np.allclose(result[0] /2 , result[1] /5)    # From the stoichiometry
    assert  np.allclose(result[1] /5 , -result[2] /4)   # From the stoichiometry
    assert  np.allclose(result[2] /4 , result[3] /3)    # From the stoichiometry



def test_single_reaction_step_5():
    chem_data = ReactionData(names=["A", "B"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([3., 5.])

    # Reaction  2A <-> B , with 2nd-order kinetics in forward reaction, and 1st-order in reverse.
    # Based on experiment "1D/reactions/reaction7"
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.02, conc_array=conc_array)
    assert np.allclose(result, [-1.4 , 0.7])
    assert np.allclose(result[0] /2 , -result[1])      # From the stoichiometry



def test_single_reaction_step_6():
    chem_data = ReactionData(names=["A", "B", "C", "D", "E"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([3., 5., 1., 0.4, 0.1])

    # Coupled reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction8"
    chem_data.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)
    assert chem_data.number_of_reactions() == 2

    result = rxn._reaction_elemental_step(delta_time=0.02, conc_array=conc_array)
    assert np.allclose(result, [-1.46 , -1.46  , 1.404 , -0.056 ,  0.056])
    assert np.allclose(result[0] , result[1])                  # From the stoichiometry
    assert np.allclose(result[3] , -result[4])                 # From the stoichiometry
    assert np.allclose(result[0] + result[4], -result[2])      # From the stoichiometry
                                                                # The increase in [A] and [E] combined
                                                                # must match the decrease in [C]



def test_single_compartment_react_variable_steps_1():
    # Based on experiment "variable_steps_1"

    # Initialize the system
    chem_data = ReactionData(names=["U", "X", "S"])

    # Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
    chem_data.add_reaction(reactants=[(2, "S")], products="U",
                           forward_rate=8., reverse_rate=2.)

    # Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
    chem_data.add_reaction(reactants="S", products="X",
                           forward_rate=6., reverse_rate=3.)

    dynamics = ReactionDynamics(reaction_data=chem_data)
    dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})

    dynamics.single_compartment_react(time_step=0.01, stop_time=0.2,
                                      variable_steps=True, thresholds={"low": 0.25, "high": 0.64})

    df = dynamics.get_history()
    #print(df)
    assert len(df) == 23

    assert np.allclose(df.iloc[0][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.0000,  50.000000,  100.000000,  0.000000])

    assert np.allclose(df.iloc[1][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.0050,  49.500000,  98.500000,  2.500000], rtol=1e-03)     # Notice the halved step size

    assert np.allclose(df.iloc[2][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.0075,  49.302500,  97.798750,  3.596250], rtol=1e-03)

    assert np.allclose(df.iloc[22][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.2050,  55.600598,   69.127620,  19.671183], rtol=1e-03)

    #print(df.iloc[22][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'))



def test_single_compartment_correct_neg_conc_1():
    # Based on "Run 3" of experiment "negative_concentrations_1

    # Initialize the system
    chem_data = ReactionData(names=["U", "X", "S"])

    # Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
    chem_data.add_reaction(reactants=[(2, "S")], products="U",
                           forward_rate=8., reverse_rate=2.)

    # Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
    chem_data.add_reaction(reactants="S", products="X",
                           forward_rate=6., reverse_rate=3.)

    dynamics = ReactionDynamics(reaction_data=chem_data)
    dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})

    dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

    dynamics.single_compartment_react(time_step=0.1, stop_time=0.8, variable_steps=False)
    # Note: negative concentrations that would arise from the given step size, get automatically intercepted - and
    #       the step sizes get reduced as needed

    df = dynamics.get_history()
    #print(df)
    assert len(df) == 11

    assert np.allclose(df.iloc[0][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.0000,  50.000000,  100.000000,  0.000000])

    assert np.allclose(df.iloc[1][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.1,  40, 70, 50], rtol=1e-03)

    assert np.allclose(df.iloc[2][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.15,  56, 74.50, 13.5], rtol=1e-03)        # Notice the halved step size

    assert np.allclose(df.iloc[3][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.25,  55.6, 60.25, 28.55], rtol=1e-03)        # Back to the requested step size

    assert np.allclose(df.iloc[5][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.45,  58.7, 45.1465, 37.4535], rtol=1e-03)

    assert np.allclose(df.iloc[6][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.5,  67.8114, 49.610575, 14.766625], rtol=1e-03)        # Notice another instance of halved step size

    assert np.allclose(df.iloc[7][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.6,  66.062420, 43.587378, 24.287783], rtol=1e-03)        # Back to the requested step size

    assert np.allclose(df.iloc[10][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.9,  76.895206, 44.446655,	1.762933], rtol=1e-03)        # Final step




###########################  LOWER-LEVEL METHODS  ###########################


def test_step_determiner_A ():
    chem_data = ReactionData()
    rxn = ReactionDynamics(chem_data)
    rxn.variable_steps_threshold_low = 50
    rxn.variable_steps_threshold_high = 100
    rxn.variable_steps_threshold_abort = 200

    assert rxn.step_determiner_A (0) == 2
    assert rxn.step_determiner_A (49.9) == 2
    assert rxn.step_determiner_A (50.1) == 1
    assert rxn.step_determiner_A (99.0) == 1
    assert np.allclose(rxn.step_determiner_A (100.1), 0.5)
    assert np.allclose(rxn.step_determiner_A (199.9), 0.5)
    assert np.allclose(rxn.step_determiner_A (200.1), 0.5)
    assert np.allclose(rxn.step_determiner_A (10000), 0.5)



def test_compute_all_rate_deltas():
    chem_data = ReactionData(names=["A", "B", "C", "D"])
    rxn = ReactionDynamics(chem_data)

    # Start with reaction A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=20., reverse_rate=2.)  # Reaction 0
    conc_array = np.array([5., 8.])
    result = rxn.compute_all_reaction_deltas(conc_array=conc_array, delta_time=0.5) # {0: 42.0}
    assert len(result) == 1
    assert np.allclose(result[0], 42.0)

    # Add reaction 2B <-> 3C , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=[(2, "B")], products=[(3, "C")], forward_rate=10., reverse_rate=25.)   # Rxn 1
    conc_array = np.array([5., 8., 15.])
    result = rxn.compute_all_reaction_deltas(conc_array=conc_array, delta_time=0.5) # {0: 42.0, 1: -147.5}
    assert len(result) == 2
    assert np.allclose(result[0], 42.0)
    assert np.allclose(result[1], -147.5)

    # Add reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                     forward_rate=5., reverse_rate=2.)          # Rxn 2
    conc_array = np.array([5., 8., 15., 7.])
    result = rxn.compute_all_reaction_deltas(conc_array=conc_array, delta_time=0.5) # {0: 42.0, 1: -147.5, 2: -5.0}
    assert len(result) == 3
    assert np.allclose(result[0], 42.0)
    assert np.allclose(result[1], -147.5)
    assert np.allclose(result[2], -5)


    # Add reaction  2A <-> B , with 2nd-order kinetics in the forward direction
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=3., reverse_rate=2.)
    conc_array = np.array([5., 8., 15., 7.])
    result = rxn.compute_all_reaction_deltas(conc_array=conc_array, delta_time=0.5) # {0: 42.0, 1: -147.5, 2: -5.0, 3: 29.5}
    assert len(result) == 4
    assert np.allclose(result[0], 42.0)
    assert np.allclose(result[1], -147.5)
    assert np.allclose(result[2], -5)
    assert np.allclose(result[3], 29.5)

    # This time, only process reactions 0 and 2
    result_0_2 = rxn.compute_all_reaction_deltas(conc_array=conc_array, delta_time=0.5,
                                                 rxn_list=[0, 2])   # {0: 42.0, 2: -5.0}
    assert len(result_0_2) == 2
    assert np.allclose(result_0_2[0], 42.0)
    assert np.allclose(result_0_2[2], -5.0)


    # FLUSH OUT ALL REACTIONS (to start over)
    rxn.clear_reactions()
    # Start with reaction A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=20., reverse_rate=2.)   # Rxn 0
    conc_array = np.array([5., 8.])
    result = rxn.compute_all_reaction_deltas(conc_array=conc_array, delta_time=0.25)  #  {0: 21.0}
    assert len(result) == 1
    assert np.allclose(result[0], 21.0)



def test_compute_rate_delta():
    chem_data = ReactionData(names=["A", "B", "C", "D"])
    rxn = ReactionDynamics(chem_data)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=20., reverse_rate=2.)
    conc_array = np.array([5., 8.])
    result = rxn.compute_reaction_delta(rxn_index=0, conc_array=conc_array, delta_time=0.5)
    assert np.allclose(result, 42.0)

    # Reaction 2B <-> 3C , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=[(2, "B")], products=[(3, "C")], forward_rate=10., reverse_rate=25.)
    conc_array = np.array([0., 8., 15.])
    result = rxn.compute_reaction_delta(rxn_index=1, conc_array=conc_array, delta_time=1.5)
    assert np.allclose(result, -442.5)

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                    forward_rate=5., reverse_rate=2.)
    conc_array = np.array([3.5, 9., 11., 7.])
    result = rxn.compute_reaction_delta(rxn_index=2, conc_array=conc_array, delta_time=0.5)
    assert np.allclose(result, 1.75)

    # Reaction  2A <-> B , with 2nd-order kinetics in the forward direction
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    conc_array = np.array([4.5, 6.])
    result = rxn.compute_reaction_delta(rxn_index=3, conc_array=conc_array, delta_time=2.0)
    assert np.allclose(result, 178.5)

    # Reaction  B <-> 2C , with 2nd-order kinetics in the reverse direction
    chem_data.add_reaction(reactants=[("B")], products=[(2, "C", 2)], forward_rate=4., reverse_rate=2.)
    conc_array = np.array([0., 5., 4.])
    result = rxn.compute_reaction_delta(rxn_index=4, conc_array=conc_array, delta_time=0.5)
    assert np.allclose(result, -6.0)



def test_is_in_equilibrium():
    chem_data = ReactionData(names=["A", "B", "C", "D", "E", "F"])
    rxn = ReactionDynamics(chem_data)

    # Reaction 0 : A <-> B
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    c = {'A': 23.9931640625, 'B': 36.0068359375}
    assert rxn.is_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)
    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=1)      # Testing ALL reactions

    # Reaction 1 : A <-> F
    chem_data.add_reaction(reactants=["A"], products=["F"], forward_rate=20, reverse_rate=2.)
    c = {'A': 3, 'F': 32.999}
    assert rxn.is_in_equilibrium(rxn_index=1, conc=c, explain=False, tolerance=10)   # The deviation is just below the 10% tolerance

    with pytest.raises(Exception):
        assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=10) # Testing ALL reactions, but neglected to provide [B]

    c = {'A': 3, 'B': 4.6, 'F': 33.001}   # We're now including all the concentrations we need for both reactions
    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=10) \
           == {False: [1]}    # The deviation for reaction 1 is just above the 10% tolerance

    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=2.2) \
           == {False: [0, 1]}      # With a tolerance so low, not reaction meets the criteria

    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=2.3) \
           == {False: [1]}                # Reaction 0 barely passes with this tolerance (it deviates by 2.222 %)

    with pytest.raises(Exception):
        assert rxn.is_in_equilibrium(explain=False, tolerance=2.3)  # We're failing to provide concentrations

    rxn.set_conc(conc=[3, 4.6, 0, 0, 0, 33.001])    # The concentrations are in the same order as the declared chemicals
    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=2.3) \
           == {False: [1]}        # Now using the System concentrations



def test_reaction_in_equilibrium():
    chem_data = ReactionData(names=["A", "B", "C", "D", "E", "F"])
    rxn = ReactionDynamics(chem_data)

    # Reaction 0 : A <-> B
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    c = {'A': 23.9931640625, 'B': 36.0068359375}
    assert rxn.reaction_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)

    # Reaction 1 : A <-> F
    chem_data.add_reaction(reactants=["A"], products=["F"], forward_rate=20, reverse_rate=2.)
    c = {'A': 3, 'F': 32.999}
    assert rxn.reaction_in_equilibrium(rxn_index = 1, conc=c, explain=False, tolerance=10)   # Just below the 10% tolerance
    c = {'A': 3, 'F': 33.001}
    assert not rxn.reaction_in_equilibrium(rxn_index = 1, conc=c, explain=False, tolerance=10)  # Just above the 10% tolerance

    # Reaction 2 : A <-> 3B
    chem_data.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    c = {'A': 14.54545455, 'B': 36.36363636}
    assert rxn.reaction_in_equilibrium(rxn_index = 2, conc=c, explain=False, tolerance=1)

    # Reaction 3:  2A <-> 3B
    chem_data.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    c = {'A': 16.25, 'B': 40.625}
    assert rxn.reaction_in_equilibrium(rxn_index = 3, conc=c, explain=False, tolerance=1)

    # Reaction 4:  A + B <-> C , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 0.29487741, 'B': 40.29487741, 'C': 29.70512259}
    assert rxn.reaction_in_equilibrium(rxn_index = 4, conc=c, explain=False, tolerance=1)

    # Reaction 5:  A <-> 2C + D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[("A")], products=[(2, "C") , ("D")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 4.31058733, 'C': 6.37882534, 'D': 1.68941267}
    assert rxn.reaction_in_equilibrium(rxn_index = 5, conc=c, explain=False, tolerance=1)

    # Reaction 6:  2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 2.80284552, 'B': 4.00711381, 'C': 7.39430896, 'D': 3.79573172}
    assert rxn.reaction_in_equilibrium(rxn_index = 6, conc=c, explain=False, tolerance=1)


    rxn.clear_reactions()   # This will reset the reaction count to 0

    # Reaction 0:  2A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=[(2, "A")], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 2.16928427, 'B': 5.41535786}
    assert rxn.reaction_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)

    # Reaction 1:  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 1.51554944, 'B': 5.74222528}
    assert rxn.reaction_in_equilibrium(rxn_index = 1, conc=c, explain=False, tolerance=1)



def test_reaction_speeds():
    # To test that newly-added reactions are automatically tentatively marked "Fast",
    # and to test the methods managing the "Slow/Fast" data structure
    chem_data = ReactionData(names=["A", "B", "C"])
    rxn = ReactionDynamics(chem_data)
    assert rxn.slow_rxns() == []        # There are no reactions yet
    assert rxn.fast_rxns() == []        # There are no reactions yet
    assert rxn.are_all_slow_rxns()

    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    assert rxn.slow_rxns() == []        # The one reaction present is assumed to be fast
    assert rxn.fast_rxns() == [0]
    assert not rxn.are_all_slow_rxns()

    rxn.set_rxn_speed(0, "S")          # Mark the lone reaction as "Slow"
    assert rxn.slow_rxns() == [0]
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    chem_data.add_reaction(reactants=["B"], products=["C"], forward_rate=13., reverse_rate=12.)
    assert rxn.reaction_data.number_of_reactions() == 2
    assert rxn.slow_rxns() == [0]
    assert rxn.fast_rxns() == [1]       # The newly-added one
    rxn.set_rxn_speed(1, "S")
    assert rxn.slow_rxns() == [0, 1]
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    rxn.set_rxn_speed(0, "F")
    assert not rxn.are_all_slow_rxns()
    rxn.set_rxn_speed(1, "F")
    assert rxn.slow_rxns() == []
    assert rxn.fast_rxns() == [0, 1]
    assert not rxn.are_all_slow_rxns()



def test_validate_increment():
    chem_data = ReactionData(names=["A", "B", "C"])
    rxn = ReactionDynamics(chem_data)

    chem_data.add_reaction(reactants=["A"], products=["B"])

    rxn.validate_increment(delta_conc=50., baseline_conc=10., rxn_index=0, species_index=2, delta_time=0.02)
    rxn.validate_increment(delta_conc=-9.99, baseline_conc=10., rxn_index=0, species_index=2, delta_time=0.02)

    with pytest.raises(Exception):      # Would lead to a negative concentration
        rxn.validate_increment(delta_conc=-10.01, baseline_conc=10., rxn_index=0, species_index=2, delta_time=0.02)



def test_examine_increment_array():
    chem_data = ReactionData(names=["A", "B", "C"])
    rxn = ReactionDynamics(chem_data)

    rxn.fast_criterion_use_baseline = True   # TODO: This will probably be phased out

    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)   # 1st reaction: A <-> B
    rxn.set_rxn_speed(0, "S")          # Mark the lone reaction as "Slow"
    assert rxn.are_all_slow_rxns()


    # INITIAL THRESHOLD (fast_threshold_fraction=0.05)

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([4.98, 0, 0]), baseline_conc_array=np.array([100., 0, 0]),
                                time_subdivision=1, fast_threshold_fraction=0.05, use_baseline = True)
    assert rxn.are_all_slow_rxns()      # The small increment didn't trip the current "slow reaction" tag


    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 5.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.05, use_baseline = True)
    assert rxn.fast_rxns() == [0]       # The one reaction present got marked as "Fast" b/c of large delta_conc
    assert not rxn.are_all_slow_rxns()


    rxn.set_rxn_speed(0, "S")           # Reset the lone reaction to "Slow"
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -4.98, 0]), baseline_conc_array=np.array([ 0, 100.,0]),
                                time_subdivision=1, fast_threshold_fraction=0.05, use_baseline = True)
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([-5.01, 0, 0]), baseline_conc_array=np.array([100., 0, 0]),
                                time_subdivision=1, fast_threshold_fraction=0.05, use_baseline = True)
    assert not rxn.are_all_slow_rxns()  # The one reaction present got marked as "Fast" b/c of large abs(delta_conc)

    rxn.set_rxn_speed(0, "S")           # Reset the lone reaction to "Slow"
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -5.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.05, use_baseline = True)
    assert not rxn.are_all_slow_rxns()  # The one reaction present got marked as "Fast" b/c of large abs(delta_conc),
                                        # no matter which chemical species is being affected

    rxn.set_rxn_speed(0, "S")           # Reset the lone reaction to "Slow"
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([2.49, 0, 0]), baseline_conc_array=np.array([100., 0, 0]),
                                time_subdivision=2, fast_threshold_fraction=0.05, use_baseline = True)
    assert rxn.are_all_slow_rxns()      # The small increment didn't trip the current "slow reaction" tag

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 2.51, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=2, fast_threshold_fraction=0.05, use_baseline = True)
    assert rxn.fast_rxns() == [0]       # The reaction got marked as "Fast" b/c of large delta_conc for a time_subdivision of 2
    assert not rxn.are_all_slow_rxns()


    # NEW THRESHOLD (fast_threshold_fraction=0.1)

    chem_data.add_reaction(reactants=["B"], products=["C"], forward_rate=13., reverse_rate=12.)   # 2nd reaction: B <-> C
    assert rxn.reaction_data.number_of_reactions() == 2
    assert rxn.get_rxn_speed(1) == "F"  # Newly-added reactions are assumed "Fast" (until proven otherwise!)

    rxn.set_rxn_speed(0, "S")           # Reset reactions to "Slow"
    rxn.set_rxn_speed(1, "S")
    assert rxn.slow_rxns() == [0, 1]
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 0, 5.01]), baseline_conc_array=np.array([0, 0, 100.]),
                                time_subdivision=1, fast_threshold_fraction=0.1, use_baseline = True)
    assert rxn.slow_rxns() == [0, 1]    # No change, because the threshold is now higher
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 9.99, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1, use_baseline = True)
    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 0, -9.99]), baseline_conc_array=np.array([0, 0, 100.]),
                                time_subdivision=1, fast_threshold_fraction=0.1, use_baseline = True)
    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([9.99, 0, 0]), baseline_conc_array=np.array([100., 0, 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1, use_baseline = True)
    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -9.99, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1, use_baseline = True)
    assert rxn.slow_rxns() == [0, 1]    # Still no change, because all the abs(delta_conc) still below threshold
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=1, delta_conc_array=np.array([0, 10.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1, use_baseline = True)
    assert rxn.slow_rxns() == [0]       # The last function call was "the straw that broke the camel's back" for reaction 1 (no longer "slow")!
    assert rxn.fast_rxns() == [1]
    assert not rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -10.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1, use_baseline = True)
    assert rxn.slow_rxns() == []        # The last function call was "straw that broke the camel's back" for reaction 0 (no longer "slow")!
    assert rxn.fast_rxns() == [0, 1]
    assert not rxn.are_all_slow_rxns()

    rxn.set_rxn_speed(0, "S")           # Reset reactions to "Slow"
    rxn.set_rxn_speed(1, "S")
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 0, 3.32]), baseline_conc_array=np.array([0, 0, 100.]),
                                time_subdivision=3, fast_threshold_fraction=0.1, use_baseline = True)
    assert rxn.are_all_slow_rxns()      # The small increment didn't trip the current "slow reaction" tag for reaction 0

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 3.34, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=3, fast_threshold_fraction=0.1, use_baseline = True)
    assert rxn.fast_rxns() == [0]       # Reaction 0 got marked as "Fast" b/c of large delta_conc for a time_subdivision of 3

    rxn.examine_increment_array(rxn_index=1, delta_conc_array=np.array([0, 0, -0.99]), baseline_conc_array=np.array([0, 0, 100.]),
                                time_subdivision=10, fast_threshold_fraction=0.1, use_baseline = True)
    assert rxn.fast_rxns() == [0]       # The small increment didn't trip the current "slow reaction" tag for reaction 1

    rxn.examine_increment_array(rxn_index=1, delta_conc_array=np.array([0, -1.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=10, fast_threshold_fraction=0.1, use_baseline = True)
    assert rxn.fast_rxns() == [0, 1]    # Reaction 1 got marked as "Fast" b/c of large delta_conc for a time_subdivision of 3



def test_examine_increment_array_2():
    chem_data = ReactionData(names=["A", "B", "C"])
    rxn = ReactionDynamics(chem_data)


    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)   # 1st reaction: A <-> B
    rxn.set_rxn_speed(0, "S")          # Mark the lone reaction as "Slow"
    assert rxn.are_all_slow_rxns()


    # INITIAL THRESHOLD (fast_threshold_fraction=5.)

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([4.98, 0, 0]),
                                time_subdivision=1, fast_threshold_fraction=5.)
    assert rxn.are_all_slow_rxns()      # The small increment didn't trip the current "slow reaction" tag


    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 5.01, 0]),
                                time_subdivision=1, fast_threshold_fraction=5.)
    assert rxn.fast_rxns() == [0]       # The one reaction present got marked as "Fast" b/c of large delta_conc
    assert not rxn.are_all_slow_rxns()


    rxn.set_rxn_speed(0, "S")           # Reset the lone reaction to "Slow"
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -4.98, 0]),
                                time_subdivision=1, fast_threshold_fraction=5.)
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([-5.01, 0, 0]),
                                time_subdivision=1, fast_threshold_fraction=5.)
    assert not rxn.are_all_slow_rxns()  # The one reaction present got marked as "Fast" b/c of large abs(delta_conc)

    rxn.set_rxn_speed(0, "S")           # Reset the lone reaction to "Slow"
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -5.01, 0]),
                                time_subdivision=1, fast_threshold_fraction=5.)
    assert not rxn.are_all_slow_rxns()  # The one reaction present got marked as "Fast" b/c of large abs(delta_conc),
    # no matter which chemical species is being affected

    rxn.set_rxn_speed(0, "S")           # Reset the lone reaction to "Slow"
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([2.49, 0, 0]),
                                time_subdivision=2, fast_threshold_fraction=5.)
    assert rxn.are_all_slow_rxns()      # The small increment didn't trip the current "slow reaction" tag

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 2.51, 0]),
                                time_subdivision=2, fast_threshold_fraction=5.)
    assert rxn.fast_rxns() == [0]       # The reaction got marked as "Fast" b/c of large delta_conc for a time_subdivision of 2
    assert not rxn.are_all_slow_rxns()


    # NEW THRESHOLD (fast_threshold_fraction=10.)

    chem_data.add_reaction(reactants=["B"], products=["C"], forward_rate=13., reverse_rate=12.)   # 2nd reaction: B <-> C
    assert rxn.reaction_data.number_of_reactions() == 2
    assert rxn.get_rxn_speed(1) == "F"  # Newly-added reactions are assumed "Fast" (until proven otherwise!)

    rxn.set_rxn_speed(0, "S")           # Reset reactions to "Slow"
    rxn.set_rxn_speed(1, "S")
    assert rxn.slow_rxns() == [0, 1]
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 0, 5.01]),
                                time_subdivision=1, fast_threshold_fraction=10.)
    assert rxn.slow_rxns() == [0, 1]    # No change, because the threshold is now higher
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 9.99, 0]),
                                time_subdivision=1, fast_threshold_fraction=10.)
    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 0, -9.99]),
                                time_subdivision=1, fast_threshold_fraction=10.)
    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([9.99, 0, 0]),
                                time_subdivision=1, fast_threshold_fraction=10.)
    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -9.99, 0]),
                                time_subdivision=1, fast_threshold_fraction=10.)
    assert rxn.slow_rxns() == [0, 1]    # Still no change, because all the abs(delta_conc) still below threshold
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=1, delta_conc_array=np.array([0, 10.01, 0]),
                                time_subdivision=1, fast_threshold_fraction=10.)
    assert rxn.slow_rxns() == [0]       # The last function call was "the straw that broke the camel's back" for reaction 1 (no longer "slow")!
    assert rxn.fast_rxns() == [1]
    assert not rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -10.01, 0]),
                                time_subdivision=1, fast_threshold_fraction=10.)
    assert rxn.slow_rxns() == []        # The last function call was "straw that broke the camel's back" for reaction 0 (no longer "slow")!
    assert rxn.fast_rxns() == [0, 1]
    assert not rxn.are_all_slow_rxns()

    rxn.set_rxn_speed(0, "S")           # Reset reactions to "Slow"
    rxn.set_rxn_speed(1, "S")
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 0, 3.32]),
                                time_subdivision=3, fast_threshold_fraction=10.)
    assert rxn.are_all_slow_rxns()      # The small increment didn't trip the current "slow reaction" tag for reaction 0

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 3.34, 0]),
                                time_subdivision=3, fast_threshold_fraction=10.)
    assert rxn.fast_rxns() == [0]       # Reaction 0 got marked as "Fast" b/c of large delta_conc for a time_subdivision of 3

    rxn.examine_increment_array(rxn_index=1, delta_conc_array=np.array([0, 0, -0.99]),
                                time_subdivision=10, fast_threshold_fraction=10.)
    assert rxn.fast_rxns() == [0]       # The small increment didn't trip the current "slow reaction" tag for reaction 1

    rxn.examine_increment_array(rxn_index=1, delta_conc_array=np.array([0, -1.01, 0]),
                                time_subdivision=10, fast_threshold_fraction=10.)
    assert rxn.fast_rxns() == [0, 1]    # Reaction 1 got marked as "Fast" b/c of large delta_conc for a time_subdivision of 3




def test_stoichiometry_checker():
    chem = ReactionData(names=["A", "B", "C", "D"])
    rxn = ReactionDynamics(chem)

    chem.add_reaction(reactants=["A"], products=["B"])          # Reaction 0:   A <--> B

    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 50, 0, 0]), conc_arr_after=np.array([10, 40, 0, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 50, 0, 0]), conc_arr_after=np.array([10, 39.9, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 0, 0, 0]), conc_arr_after=np.array([90, 10, 0, 0]))


    chem.add_reaction(reactants=[(2, "A")], products=["B"])     # Reaction 1:   2A <--> B

    assert rxn.stoichiometry_checker(rxn_index=1, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=1, conc_arr_before=np.array([100, 0, 0, 0]), conc_arr_after=np.array([80, 10, 0, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=1, conc_arr_before=np.array([100, 0, 0, 0]), conc_arr_after=np.array([80, 10.1, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=1, conc_arr_before=np.array([0, 50, 0, 0]), conc_arr_after=np.array([10, 45, 0, 0]))


    chem.add_reaction(reactants=["A", "B"], products=["C"])     # Reaction 2:   A + B <--> C

    assert rxn.stoichiometry_checker(rxn_index=2, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=2, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 40, 10, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=2, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 40, 9.9, 0]))


    chem.add_reaction(reactants=["A", (3, "B")], products=["C"])     # Reaction 3:   A + 3B <--> C

    assert rxn.stoichiometry_checker(rxn_index=3, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=3, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 10, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=3, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 9.9, 0]))


    chem.add_reaction(reactants=["A", (3, "B")], products=[(4, "C")])                   # Reaction 4:   A + 3B <--> 4C

    assert rxn.stoichiometry_checker(rxn_index=4, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=4, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 40, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=4, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 39.9, 0]))


    chem.add_reaction(reactants=[(2, "A"), (3, "B")], products=[(4, "C"), (5, "D")])     # Reaction 5:   2A + 3B <--> 4C + 5D
    assert rxn.stoichiometry_checker(rxn_index=5, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=5, conc_arr_before=np.array([100, 100, 100, 100]), conc_arr_after=np.array([120, 130, 60, 50]))
    assert not rxn.stoichiometry_checker(rxn_index=5, conc_arr_before=np.array([100, 100, 100, 100]), conc_arr_after=np.array([120.1, 130, 60, 50]))
    assert rxn.stoichiometry_checker(rxn_index=5, conc_arr_before=np.array([100, 100, 100, 100]), conc_arr_after=np.array([80, 70, 140, 150]))
    assert not rxn.stoichiometry_checker(rxn_index=5, conc_arr_before=np.array([100, 100, 100, 100.1]), conc_arr_after=np.array([80, 70, 140, 150]))



def test_explain_time_advance():
    rxn = ReactionDynamics(None)

    with pytest.raises(Exception):
        rxn.explain_time_advance(return_times=True)      # Diagnostics weren't first enabled

    rxn.set_diagnostics()

    with pytest.raises(Exception):
        rxn.explain_time_advance(return_times=True)      # No diagnostic data yet present

    # Start out with uniform steps
    rxn.diagnostic_conc_data.store(par=20.,
                                   data_snapshot={"primary_timestep": 100.})

    assert rxn.explain_time_advance(return_times=True, silent=True) is None


    rxn.diagnostic_conc_data.store(par=30.,
                                   data_snapshot={"primary_timestep": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)    # TODO: also test the returned step sizes
    assert np.allclose(result, [20., 30.])

    rxn.diagnostic_conc_data.store(par=40.,
                                   data_snapshot={"primary_timestep": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40.])

    # Switching to smaller step
    rxn.diagnostic_conc_data.store(par=45.,
                                   data_snapshot={"primary_timestep": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 45.])

    rxn.diagnostic_conc_data.store(par=50.,
                                   data_snapshot={"primary_timestep": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50.])

    # Switching to larger step
    rxn.diagnostic_conc_data.store(par=70.,
                                   data_snapshot={"primary_timestep": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70.])

    # Yet larger
    rxn.diagnostic_conc_data.store(par=95.,
                                   data_snapshot={"primary_timestep": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95.])

    # Smaller again
    rxn.diagnostic_conc_data.store(par=96.,
                                   data_snapshot={"primary_timestep": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95., 96.])

    rxn.diagnostic_conc_data.store(par=97.,
                                   data_snapshot={"primary_timestep": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95., 97.])

    rxn.diagnostic_conc_data.store(par=98.,
                                   data_snapshot={"primary_timestep": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95., 98.])

    #print(rxn.diagnostic_data_baselines.get())
    #print(result)



def test__delta_names():
    chem_data = ReactionData(names=["A", "B", "X"])
    dyn = ReactionDynamics(chem_data)

    assert dyn._delta_names() == ["Delta A", "Delta B", "Delta X"]



def test__delta_conc_dict():
    chem_data = ReactionData(names=["A", "B", "X"])
    dyn = ReactionDynamics(chem_data)

    assert dyn._delta_conc_dict(np.array([10, 20, 30])) == \
           {"Delta A": 10, "Delta B": 20, "Delta X": 30}

    with pytest.raises(Exception):
        dyn._delta_conc_dict(np.array([10, 20, 30, 40]))    # One element too many



def test_save_diagnostic_rxn_data():
    chem_data = ReactionData(names=["A", "B", "C", "X"])
    # Add 3 reactions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A"], products=["X"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A", "B"], products=["X"], forward_rate=5., reverse_rate=2.)

    dyn = ReactionDynamics(chem_data)
    assert len(dyn.diagnostic_rxn_data) == 0

    dyn.system_time = 10.
    dyn.save_diagnostic_rxn_data(rxn_index=0, delta_time=0.5, increment_vector_single_rxn=np.array([2.4, -2.4, 0]))

    assert len(dyn.diagnostic_rxn_data) == 3    # Same as the number of reactions

    # TODO: verify that dyn.diagnostic_rxn_data got set correctly



def test_get_diagnostic_rxn_data():
    pass