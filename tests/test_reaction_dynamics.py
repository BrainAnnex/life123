import pytest
import numpy as np
from modules.reactions.reaction_data import ReactionData
from modules.reactions.reaction_dynamics import ReactionDynamics



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



def test_single_reaction_fixed_step():
    chem_data = ReactionData(names=["A", "B"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([10., 50.])

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn.reaction_step_FIXED_RESOLUTION(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [ 7. , -7.])
    assert result[0] == - result[1]         # From the stoichiometry


    rxn.clear_reactions()       # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction2"
    chem_data.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn.reaction_step_FIXED_RESOLUTION(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]      # From the stoichiometry


    rxn.clear_reactions()       # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction3"
    chem_data.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn.reaction_step_FIXED_RESOLUTION(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test_single_reaction_variable_step_1():
    chem_data = ReactionData(names=["A", "B"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([10., 50.])

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [ 7. , -7.])
    assert result[0] == - result[1]         # From the stoichiometry


    rxn.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction2"
    chem_data.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]      # From the stoichiometry


    rxn.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction3"
    chem_data.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)

    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test_single_reaction_step_2():
    chem_data = ReactionData(names=["A", "B", "C"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([10., 50., 20.])

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # # Based on experiment "reactions_single_compartment/react_1"
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.1, conc_array=conc_array)
    assert np.allclose(result, [ 7. , -7. , 0.])    # Chemical "C" not participating in this reaction; its delta conc. is 0
    assert result[0] == - result[1]         # From the stoichiometry


    chem_data.clear_reactions_data()   # Re-start with a blank slate of reactions
    # Reaction A + B <-> C , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction4"
    chem_data.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                     forward_rate=5., reverse_rate=2.)

    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.002, conc_array=conc_array)
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

    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.05, conc_array=conc_array)
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

    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.001, conc_array=conc_array)
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

    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.02, conc_array=conc_array)
    assert np.allclose(result, [-1.4 , 0.7])
    assert  np.allclose(result[0] /2 , -result[1])      # From the stoichiometry



def test_single_reaction_step_6():
    chem_data = ReactionData(names=["A", "B", "C", "D", "E"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([3., 5., 1., 0.4, 0.1])

    # Coupled reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction8"
    chem_data.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)
    assert chem_data.number_of_reactions() == 2

    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.02, conc_array=conc_array)
    assert np.allclose(result, [-1.46 , -1.46  , 1.404 , -0.056 ,  0.056])
    assert  np.allclose(result[0] , result[1])                  # From the stoichiometry
    assert  np.allclose(result[3] , -result[4])                 # From the stoichiometry
    assert  np.allclose(result[0] + result[4], -result[2])      # From the stoichiometry
                                                                # The increase in [A] and [E] combined
                                                                # must match the decrease in [C]



def test_adaptive_time_resolution_1():
    chem_data = ReactionData(names=["A", "B"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([10., 50.])

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_2"
    kF = 3.
    kR = 2.
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=kF, reverse_rate=kR)

    # Start testing the lowest-level function, and then proceed to testing progressively higher-level ones
    #rxn.verbose_list = [1]
    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=0.05, conc_array=conc_array,
                                                   time_subdivision=2, fast_threshold_fraction=0.05)
    # The above call is just 1 time step
    # Check the calculations, based on the forward Euler method
    fwd_delta = 0.05 * (-kF * conc_array[0] + kR * conc_array[1])   # 3.5
    rev_delta = -fwd_delta   # From the stoichiometry
    assert np.allclose(result, [fwd_delta ,rev_delta])    # [3.5 , -3.5]


    # Repeat at the next-higher level
    rxn.clear_reactions()   # IMPORTANT: because it'll reset all reaction in their default initial "fast" mode
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=kF, reverse_rate=kR)   # Re-add same reaction

    result = rxn.reaction_step_orchestrator(delta_time_full=0.1, conc_array=conc_array,
                                            dynamic_steps=2, fast_threshold=5)
    # The above call results in 2 time steps
    # Check the calculations, based on the forward Euler method
    half_step_conc = conc_array + [fwd_delta ,rev_delta]    # [13.5 46.5]   These are the conc's halfway thru delta_time_full
    new_fwd_delta = 0.05 * (-kF * half_step_conc[0] + kR * half_step_conc[1])   # 2.625
    new_rev_delta = -new_fwd_delta  # From the stoichiometry
    fwd_delta += new_fwd_delta      # 6.125
    rev_delta += new_rev_delta      # - 6.125
    assert np.allclose(result, [fwd_delta ,rev_delta])  # [6.125, -6.125]


    # Repeat at the yet-next-higher level
    rxn.clear_reactions()   # IMPORTANT: because it'll reset all reaction in their default initial "fast" mode
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=kF, reverse_rate=kR)   # Re-add same reaction
    rxn.system = conc_array.copy()      # The copy() is to avoid messing up conc_array
    rxn.system_time = 0.
    rxn.single_compartment_react(time_step=0.1, n_steps=1, dynamic_steps=2, fast_threshold=5)

    assert np.allclose(rxn.system_time, 0.1)
    assert np.allclose(rxn.system, conc_array + np.array([fwd_delta ,rev_delta]))     # [16.125, 43.875]
    # Note: in the previous scenario, we had computed the "delta's"; now, we're looking at the updated system state



def test_adaptive_time_resolution_2():
    chem_data = ReactionData(names=["A", "B", "C"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([10., 50., 20.])

    # Reaction A + B <-> C , with 1st-order kinetics for each species
    # Based on experiment "reactions_single_compartment/react_3"
    kF = 5.
    kR = 2.

    delta_time_full_interval = 0.004
    time_subdivision = 2
    delta_time_subinterval = delta_time_full_interval /time_subdivision

    chem_data.add_reaction(reactants=["A" , "B"], products=["C"],
                           forward_rate=kF, reverse_rate=kR)

    # Start testing the lower-level functions, and then proceed to testing progressively higher-level ones
    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=delta_time_subinterval, conc_array=conc_array,
                                                   time_subdivision=time_subdivision, fast_threshold_fraction=0.05)
    # Check the calculations, based on the forward Euler method
    delta_A = delta_time_subinterval * (-kF * conc_array[0] * conc_array[1] + kR * conc_array[2])   # -4.92
    delta_B = delta_A       # From the stoichiometry
    delta_C = -delta_A      # From the stoichiometry
    assert np.allclose(result, [delta_A ,delta_B, delta_C])    # [-4.92 -4.92  4.92]


    # Repeat at the next-higher level
    rxn.clear_reactions()   # IMPORTANT: because it'll reset all reaction in their default initial "fast" mode
    chem_data.add_reaction(reactants=["A" , "B"], products=["C"],
                           forward_rate=kF, reverse_rate=kR)        # Re-add the reaction
    result = rxn.reaction_step_orchestrator(delta_time_full=delta_time_full_interval, conc_array=conc_array,
                                            dynamic_steps=time_subdivision, fast_threshold=5)
    # Check the calculations, based on the forward Euler method
    half_step_conc = conc_array + [delta_A ,delta_B, delta_C]    # [5.08 45.08 24.92]  These are the conc's halfway thru delta_time_full

    new_delta_A = delta_time_subinterval * (-kF * half_step_conc[0] * half_step_conc[1] + kR * half_step_conc[2])   # -2.190384
    new_delta_B = new_delta_A       # From the stoichiometry
    new_delta_C = -new_delta_A      # From the stoichiometry

    delta_A += new_delta_A
    delta_B += new_delta_B
    delta_C += new_delta_C

    assert np.allclose(result, [delta_A ,delta_B, delta_C])  # [-7.110384 -7.110384 7.110384]


    # Repeat at the yet-next-higher level
    rxn.clear_reactions()   # IMPORTANT: because it'll reset all reaction in their default initial "fast" mode
    chem_data.add_reaction(reactants=["A" , "B"], products=["C"],
                           forward_rate=5., reverse_rate=2.)        # Re-add the reaction
    rxn.system = conc_array.copy()      # The copy() is to avoid messing up conc_array
    rxn.system_time = 0.
    rxn.single_compartment_react(time_step=delta_time_full_interval, n_steps=1, dynamic_steps=time_subdivision, fast_threshold=5)

    assert np.allclose(rxn.system_time, delta_time_full_interval)
    assert np.allclose(rxn.system,
                       conc_array + np.array([delta_A ,delta_B, delta_C]))     # [ 2.889616 42.889616 27.110384]


    # Do one more step at the high level
    rxn.single_compartment_react(time_step=delta_time_full_interval,
                                 n_steps=1, dynamic_steps=time_subdivision, fast_threshold=5)
    assert np.allclose(rxn.system_time, 2*delta_time_full_interval)     # system_time now is 0.008
    assert np.allclose(rxn.system, np.array([1.13726186, 41.13726186, 28.86273814]))


    # Do several (7) more steps at the high level
    rxn.single_compartment_react(time_step=delta_time_full_interval,
                                 n_steps=7, dynamic_steps=time_subdivision, fast_threshold=5)
    assert np.allclose(rxn.system_time, 0.036)  # Note that 0.036 is 0.008 + 7 * 0.004
    assert np.allclose(rxn.system, np.array([0.29501266, 40.29501266, 29.70498734]))



def test_adaptive_time_resolution_3():
    chem_data = ReactionData(names=["A", "C"])
    rxn = ReactionDynamics(chem_data)

    conc_array = np.array([200., 40.])

    # Reaction 2A <-> C,
    # with 2nd-order kinetics for A, and 1-st order kinetics for C
    # Based on experiment "reactions_single_compartment/react_4"
    kF = 3.
    kR = 2.

    delta_time_full_interval = 0.002
    time_subdivision = 4
    delta_time_subinterval = delta_time_full_interval /time_subdivision

    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["C"],
                           forward_rate=kF, reverse_rate=kR)
                           # Note: the first 2 in (2, "A", 2) is the stoichiometry coefficient,
                           #       while the other one is the order

    # Start testing the lower-level functions, and then proceed to testing progressively higher-level ones
    result = rxn.reaction_step_VARIABLE_RESOLUTION(delta_time=delta_time_subinterval, conc_array=conc_array,
                                                   time_subdivision=time_subdivision, fast_threshold_fraction=0.05)
    # Check the calculations, based on the forward Euler method
    delta_A = 2 * delta_time_subinterval * (-kF * conc_array[0] **2 + kR * conc_array[1])   # -119.92
    delta_C = -delta_A / 2.                           # From the stoichiometry
    assert np.allclose(result, [delta_A ,delta_C])    # [-119.92   59.96]


    # Repeat at the next-higher level
    rxn.clear_reactions()   # IMPORTANT: because it'll reset all reaction in their default initial "fast" mode
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["C"],
                           forward_rate=kF, reverse_rate=kR)        # Re-add the reaction
    result = rxn.reaction_step_orchestrator(delta_time_full=delta_time_full_interval, conc_array=conc_array,
                                            dynamic_steps=time_subdivision, fast_threshold=5)

    # Check the calculations, based on the forward Euler method

    # 1st substep, using the previously-computed [delta_A ,delta_C]
    substep_conc = conc_array + [delta_A ,delta_C]  # [80.08, 99.96]  These are the conc's 1/4 way thru delta_time_full
                                                    #                 i.e. at t = 0.0005
    # 2nd substep
    new_delta_A = 2 * delta_time_subinterval * (-kF * substep_conc[0] **2 + kR * substep_conc[1])
    new_delta_C = -new_delta_A / 2.       # From the stoichiometry

    delta_A += new_delta_A
    delta_C += new_delta_C

    substep_conc = conc_array + [delta_A ,delta_C]  # [61.0415008, 109.4792496]  These are the conc's 1/2 way thru delta_time_full
                                                    #                 i.e. at t = 0.0010
    # 3rd substep
    new_delta_A = 2 * delta_time_subinterval * (-kF * substep_conc[0] **2 + kR * substep_conc[1])   #
    new_delta_C = -new_delta_A / 2.       # From the stoichiometry

    delta_A += new_delta_A
    delta_C += new_delta_C

    substep_conc = conc_array + [delta_A ,delta_C]  # [50.08226484, 114.95886758]  These are the conc's 3/4 way thru delta_time_full
                                                    #                 i.e. at t = 0.0015

    # 4th substep (completing the full interval)
    new_delta_A = 2 * delta_time_subinterval * (-kF * substep_conc[0] **2 + kR * substep_conc[1])   #
    new_delta_C = -new_delta_A / 2.       # From the stoichiometry

    delta_A += new_delta_A
    delta_C += new_delta_C

    #substep_conc = conc_array + [delta_A ,delta_C]  # [42.78748282, 118.60625859]  These are the conc's 3/4 way thru delta_time_full
                                                    #                 i.e. at t = 0.0020

    assert np.allclose(result, [delta_A , delta_C])  # [-157.21251718   78.60625859]


    # Repeat at the yet-next-higher level
    rxn.clear_reactions()   # IMPORTANT: because it'll reset all reaction in their default initial "fast" mode
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["C"],
                           forward_rate=kF, reverse_rate=kR)        # Re-add the reaction
    rxn.system = conc_array.copy()      # The copy() is to avoid messing up conc_array
    rxn.system_time = 0.
    rxn.single_compartment_react(time_step=delta_time_full_interval, n_steps=1, dynamic_steps=time_subdivision, fast_threshold=5)

    assert np.allclose(rxn.system_time, delta_time_full_interval)
    assert np.allclose(rxn.system,
                       conc_array + np.array([delta_A, delta_C]))     # [42.78748282, 118.60625859]


    # Do one more step at the high level
    rxn.single_compartment_react(time_step=delta_time_full_interval,
                                 n_steps=1, dynamic_steps=time_subdivision, fast_threshold=5)
    assert np.allclose(rxn.system_time, 2 * delta_time_full_interval)     # system_time now is 0.004
    assert np.allclose(rxn.system, np.array([27.89238785, 126.05380607]))


    # Do several (3) more steps at the high level
    rxn.single_compartment_react(time_step=delta_time_full_interval,
                                 n_steps=3, dynamic_steps=time_subdivision, fast_threshold=5)
    assert np.allclose(rxn.system_time, 0.01)   # Note that 0.01 is 0.004 + 3 * 0.002
    assert np.allclose(rxn.system, np.array([15.34008717, 132.32995642]))



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
    assert rxn.is_in_equilibrium(rxn_index = 0, conc = c, explain=False)

    # Reaction 1 : A <-> F
    chem_data.add_reaction(reactants=["A"], products=["F"], forward_rate=20, reverse_rate=2.)
    c = {'A': 3, 'F': 32.999}
    assert rxn.is_in_equilibrium(rxn_index = 1, conc = c, explain=False, tolerance=10)   # Just below the 10% tolerance
    c = {'A': 3, 'F': 33.001}
    assert not rxn.is_in_equilibrium(rxn_index = 1, conc = c, explain=False, tolerance=10)  # Just above the 10% tolerance

    # Reaction 2 : A <-> 3B
    chem_data.add_reaction(reactants=["A"], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    c = {'A': 14.54545455, 'B': 36.36363636}
    assert rxn.is_in_equilibrium(rxn_index = 2, conc = c, explain=False)

    # Reaction 3:  2A <-> 3B
    chem_data.add_reaction(reactants=[(2,"A")], products=[(3,"B")], forward_rate=5., reverse_rate=2.)
    c = {'A': 16.25, 'B': 40.625}
    assert rxn.is_in_equilibrium(rxn_index = 3, conc = c, explain=False)

    # Reaction 4:  A + B <-> C , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 0.29487741, 'B': 40.29487741, 'C': 29.70512259}
    assert rxn.is_in_equilibrium(rxn_index = 4, conc = c, explain=False)

    # Reaction 5:  A <-> 2C + D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[("A")], products=[(2, "C") , ("D")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 4.31058733, 'C': 6.37882534, 'D': 1.68941267}
    assert rxn.is_in_equilibrium(rxn_index = 5, conc = c, explain=False)

    # Reaction 6:  2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 2.80284552, 'B': 4.00711381, 'C': 7.39430896, 'D': 3.79573172}
    assert rxn.is_in_equilibrium(rxn_index = 6, conc = c, explain=False)


    rxn.clear_reactions()   # This will reset the reaction count to 0

    # Reaction 0:  2A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=[(2, "A")], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 2.16928427, 'B': 5.41535786}
    assert rxn.is_in_equilibrium(rxn_index = 0, conc = c, explain=False)

    # Reaction 1:  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 1.51554944, 'B': 5.74222528}
    assert rxn.is_in_equilibrium(rxn_index = 1, conc = c, explain=False)



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

    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)   # 1st reaction: A <-> B
    rxn.set_rxn_speed(0, "S")          # Mark the lone reaction as "Slow"
    assert rxn.are_all_slow_rxns()


    # INITIAL THRESHOLD (fast_threshold_fraction=0.05)

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([4.98, 0, 0]), baseline_conc_array=np.array([100., 0, 0]),
                                time_subdivision=1, fast_threshold_fraction=0.05)
    assert rxn.are_all_slow_rxns()      # The small increment didn't trip the current "slow reaction" tag


    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 5.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.05)
    assert rxn.fast_rxns() == [0]       # The one reaction present got marked as "Fast" b/c of large delta_conc
    assert not rxn.are_all_slow_rxns()


    rxn.set_rxn_speed(0, "S")           # Reset the lone reaction to "Slow"
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -4.98, 0]), baseline_conc_array=np.array([ 0, 100.,0]),
                                time_subdivision=1, fast_threshold_fraction=0.05)
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([-5.01, 0, 0]), baseline_conc_array=np.array([100., 0, 0]),
                                time_subdivision=1, fast_threshold_fraction=0.05)
    assert not rxn.are_all_slow_rxns()  # The one reaction present got marked as "Fast" b/c of large abs(delta_conc)

    rxn.set_rxn_speed(0, "S")           # Reset the lone reaction to "Slow"
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -5.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.05)
    assert not rxn.are_all_slow_rxns()  # The one reaction present got marked as "Fast" b/c of large abs(delta_conc),
                                        # no matter which chemical species is being affected

    rxn.set_rxn_speed(0, "S")           # Reset the lone reaction to "Slow"
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([2.49, 0, 0]), baseline_conc_array=np.array([100., 0, 0]),
                                time_subdivision=2, fast_threshold_fraction=0.05)
    assert rxn.are_all_slow_rxns()      # The small increment didn't trip the current "slow reaction" tag

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 2.51, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=2, fast_threshold_fraction=0.05)
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
                                time_subdivision=1, fast_threshold_fraction=0.1)
    assert rxn.slow_rxns() == [0, 1]    # No change, because the threshold is now higher
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 9.99, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1)
    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 0, -9.99]), baseline_conc_array=np.array([0, 0, 100.]),
                                time_subdivision=1, fast_threshold_fraction=0.1)
    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([9.99, 0, 0]), baseline_conc_array=np.array([100., 0, 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1)
    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -9.99, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1)
    assert rxn.slow_rxns() == [0, 1]    # Still no change, because all the abs(delta_conc) still below threshold
    assert rxn.fast_rxns() == []
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=1, delta_conc_array=np.array([0, 10.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1)
    assert rxn.slow_rxns() == [0]       # The last function call was "the straw that broke the camel's back" for reaction 1 (no longer "slow")!
    assert rxn.fast_rxns() == [1]
    assert not rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, -10.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=1, fast_threshold_fraction=0.1)
    assert rxn.slow_rxns() == []        # The last function call was "straw that broke the camel's back" for reaction 0 (no longer "slow")!
    assert rxn.fast_rxns() == [0, 1]
    assert not rxn.are_all_slow_rxns()

    rxn.set_rxn_speed(0, "S")           # Reset reactions to "Slow"
    rxn.set_rxn_speed(1, "S")
    assert rxn.are_all_slow_rxns()

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 0, 3.32]), baseline_conc_array=np.array([0, 0, 100.]),
                                time_subdivision=3, fast_threshold_fraction=0.1)
    assert rxn.are_all_slow_rxns()      # The small increment didn't trip the current "slow reaction" tag for reaction 0

    rxn.examine_increment_array(rxn_index=0, delta_conc_array=np.array([0, 3.34, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=3, fast_threshold_fraction=0.1)
    assert rxn.fast_rxns() == [0]       # Reaction 0 got marked as "Fast" b/c of large delta_conc for a time_subdivision of 3

    rxn.examine_increment_array(rxn_index=1, delta_conc_array=np.array([0, 0, -0.99]), baseline_conc_array=np.array([0, 0, 100.]),
                                time_subdivision=10, fast_threshold_fraction=0.1)
    assert rxn.fast_rxns() == [0]       # The small increment didn't trip the current "slow reaction" tag for reaction 1

    rxn.examine_increment_array(rxn_index=1, delta_conc_array=np.array([0, -1.01, 0]), baseline_conc_array=np.array([0, 100., 0]),
                                time_subdivision=10, fast_threshold_fraction=0.1)
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

