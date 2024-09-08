import pytest
import numpy as np
from life123.reaction import Reaction
from life123.reaction_dynamics import ReactionDynamics, VariableTimeSteps




def test_solve_exactly():
    rxn = Reaction(reactants = "A", products="B", forward_rate=3., reverse_rate=2.)

    t = np.array([0, 0.005, 0.31739, 1.12, 100.])

    A_exact, B_exact = ReactionDynamics.solve_exactly(rxn=rxn, A0=80., B0=10., t_arr=t)

    assert np.allclose(A_exact, [80, 78.91363613, 45.0, 36.162706, 36])
    assert np.allclose(B_exact, [10, 11.08636387, 45.0, 53.837294, 54])



def test_approx_solution_combination_rxn():
    #rxn = Reaction(reactants = ["A", "B"], products="C", forward_rate=3., reverse_rate=2.)

    t = np.array([0, 0.000864, 0.001555, 0.009850, 0.067400])

    result = ReactionDynamics.approx_solution_combination_rxn(kF=5., kR=2., A0=10., B0=50., C0=20., t_arr=t)

    assert np.allclose(result[0], [10., 7.88742367,   6.53396759,  0.88583575,  0.29487746])
    assert np.allclose(result[1], [50., 47.88742367, 46.53396759, 40.88583575, 40.29487746])
    assert np.allclose(result[2], [20., 22.11257633, 23.46603241, 29.11416425, 29.70512254])



def test_exact_solution_combination_rxn():
    #rxn = Reaction(reactants = ["A", "B"], products="C", forward_rate=3., reverse_rate=2.)

    t = np.array([0, 0.000864, 0.001555, 0.009850, 0.067400])

    print()
    result = ReactionDynamics.exact_solution_combination_rxn(kF=5., kR=2., A0=10., B0=50., C0=20., t=0)
    print(result)
    print()
    result = ReactionDynamics.exact_solution_combination_rxn(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.03)
    print(result)



def test_exact_solution_combination_rxn_WRONG():
    #rxn = Reaction(reactants = ["A", "B"], products="C", forward_rate=3., reverse_rate=2.)

    t = np.array([0, 0.000864, 0.001555, 0.009850, 0.067400])

    result = ReactionDynamics.exact_solution_combination_rxn_WRONG(kF=5., kR=2., A0=10., B0=50., C0=20., t_arr=t)

    '''
    assert np.allclose(result[0], [10., 7.88742367,   6.53396759,  0.88583575,  0.29487746])
    assert np.allclose(result[1], [50., 47.88742367, 46.53396759, 40.88583575, 40.29487746])
    assert np.allclose(result[2], [20., 22.11257633, 23.46603241, 29.11416425, 29.70512254])
    '''




#########################################################################################################

def test_set_thresholds():
    rd = VariableTimeSteps()
    assert rd.thresholds == []

    with pytest.raises(Exception):
        rd.set_thresholds(norm=123)     # Bad norm name

    with pytest.raises(Exception):
        rd.set_thresholds(norm="")      # Bad norm name

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", low=5, high=1)     # Can't have low > high

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", low=5, abort=1)    # Can't have low > abort

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", high=8, abort=1)   # Can't have high > abort

    rd.set_thresholds(norm="norm_A", low=5)                 # Create a new rule, for norm_A
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 5}]

    rd.set_thresholds(norm="norm_A", low=6)                 # Update an existing value
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 6}]

    rd.set_thresholds(norm="norm_A", high=8)                # Add a new value to an existing rule
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 6, 'high': 8}]

    rd.set_thresholds(norm="norm_A", abort=10)               # Add a new value to an existing rule
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 6, 'high': 8, 'abort': 10}]

    # Bad values that violate low < high < abort
    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", low=8)

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", high=6)    # Too small

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", high=10)   # Too big

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", abort=8)

    assert rd.thresholds == [{'norm': 'norm_A', 'low': 6, 'high': 8, 'abort': 10}]  # Nothing got changed by the failed calls



def test_delete_thresholds():
    rd = VariableTimeSteps()
    assert rd.thresholds == []

    with pytest.raises(Exception):
        rd.delete_thresholds(norm="not found")  # No rule with that norm exists

    rd.set_thresholds(norm="norm_A", low=1, high=2, abort=3)
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 1, 'high': 2, 'abort': 3}]

    rd.delete_thresholds(norm="norm_A", high=True)
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 1, 'abort': 3}]

    with pytest.raises(Exception):
        rd.delete_thresholds(norm="norm_A", high=True)  # Trying to delete a non-existing threshold

    rd.delete_thresholds(norm="norm_A", low=True)
    assert rd.thresholds == [{'norm': 'norm_A', 'abort': 3}]

    rd.delete_thresholds(norm="norm_A", abort=True)
    assert rd.thresholds == []      # Nothing is left of that rule



def test_adjust_timestep():
    rd = VariableTimeSteps()

    prev =     np.array([1,   8, 8, 10, 10])
    baseline = np.array([2,   5, 5, 14, 14])
    delta =    np.array([0.5, 1, 4, -2, -5])

    n_chems = len(baseline)     # 5 chemicals (with indexes 0 thru 4)


    normA = rd.norm_A(delta_conc=delta)
    assert np.allclose(normA, 1.85)

    normB = rd.norm_B(baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normB, 0.8)


    rd.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.84)
    rd.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.79)
    rd.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4, error=0.25)

    indexes_of_active_chemicals = [0]     # To indicate that just the 0-th chemical is to be considered in the norms
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'step_factor': 1, 'norms': {'norm_A': 0.25, 'norm_B': 0.25}, 'applicable_norms': 'ALL'}

    indexes_of_active_chemicals = [0, 1]     # To indicate that just chemicals with indices 0 and 1 are to be considered in the norms
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'step_factor': 1, 'norms': {'norm_A': 0.3125, 'norm_B': 0.25}, 'applicable_norms': 'ALL'}

    indexes_of_active_chemicals = [0, 1, 2, 3, 4]       # All the chemicals are to be considered in the norms, from now on

    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,  delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85}, 'applicable_norms': ['norm_A']}

    rd.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.86)     # normA (1.85) no longer triggers abort, but normB (0.8) still does
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_B']}

    rd.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.81)    # normB (0.8) no longer triggers abort, but triggers a high.
                                                                        # normA (1.85) triggers a high, too
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'high', 'step_factor': 0.5, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_A', 'norm_B']}

    rd.set_thresholds(norm="norm_A", low=0.5, high=1.86, abort=1.87)    # normA (1.85) no longer triggers high nor abort
    rd.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.79)
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_B']}

    rd.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.81)    # normB (0.8) no longer triggers abort, but still triggers a high
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'high', 'step_factor': 0.5, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_B']}

    rd.set_thresholds(norm="norm_B", low=0.08, high=0.81, abort=0.82)    # normB (0.8) no longer triggers high nor abort
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'step_factor': 1, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'ALL'}

    rd.set_thresholds(norm="norm_A", low=1.86, high=1.87, abort=1.88)   # normA (1.85) will now trigger a "low"
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'step_factor': 1, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'ALL'}
    # We're still on the 'stay' action because we aren't below ALL the thresholds

    rd.set_thresholds(norm="norm_B", low=0.81, high=0.82, abort=0.83)   # normB (0.8) will now trigger a "low", too
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'low', 'step_factor': 1.2, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'ALL'}

    rd.set_thresholds(norm="norm_C", low=1.34, high=2.60, abort=2.61)   # normC (1.333) will still continue to trigger a "low"
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'low', 'step_factor': 1.2, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 4/3}, 'applicable_norms': 'ALL'}

    rd.set_thresholds(norm="norm_C", low=1.32, high=2.60, abort=2.61)   # normC (1.333) will no longer trigger a "low" - but not a "high" nor an "abort"
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'step_factor': 1, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 4/3}, 'applicable_norms': 'ALL'}

    rd.set_thresholds(norm="norm_C", low=1.3, high=1.32, abort=2.61)   # normC (1.333) will now trigger a "high"
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,  delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'high', 'step_factor': 0.5, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 4/3}, 'applicable_norms': ['norm_C']}

    rd.set_thresholds(norm="norm_C", low=1, high=1.31, abort=1.32)   # normC (1.333) will now trigger an "abort"
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 4/3}, 'applicable_norms': ['norm_C']}

    rd.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.79)    # normB (0.8) will now trigger an "abort" before we can even get to normC
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_B']}

    rd.set_thresholds(norm="norm_A", low=0.5, high=1.0, abort=1.84)    # normA (1.85) will now trigger an "abort" before we can even get to normB
    result = rd.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85}, 'applicable_norms': ['norm_A']}

    # TODO: also test with the other norms



def test_relative_significance():
    rd = VariableTimeSteps()
    
    # Assess the relative significance of various quantities
    #       relative to a baseline value of 10
    assert rd.relative_significance(1, 10) == "S"
    assert rd.relative_significance(4.9, 10) == "S"
    assert rd.relative_significance(5.1, 10) == "C"
    assert rd.relative_significance(10, 10) == "C"
    assert rd.relative_significance(19.9, 10) == "C"
    assert rd.relative_significance(20.1, 10) == "L"
    assert rd.relative_significance(137423, 10) == "L"






#################   INDIVIDUAL NORMS   #################

def test_norms():
    rd = VariableTimeSteps()
    
    prev =     np.array([1,   8, 8, 10, 10])
    baseline = np.array([2,   5, 5, 14, 14])
    delta =    np.array([0.5, 1, 4, -2, -5])

    normA = rd.norm_A(delta_conc=delta)
    assert np.allclose(normA, 1.85)

    normB = rd.norm_B(baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normB, 0.8)

    normC = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normC, 4/3)

    normD = rd.norm_D(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normD, 0.7833333333333333)



def test_norm_A():
    rd = VariableTimeSteps()

    delta_conc = np.array([1, 4])
    result = rd.norm_A(delta_conc)
    assert np.allclose(result, 4.25)

    delta_conc = np.array([.5, 2])
    result = rd.norm_A(delta_conc)
    assert np.allclose(result, 1.0625)

    delta_conc = np.array([.5, 2, 1])
    result = rd.norm_A(delta_conc)
    assert np.allclose(result, 0.5833333333)

    delta_conc = np.array([0.5, 1, 4, -2, -5])
    result = rd.norm_A(delta_conc)
    assert np.allclose(result, 1.85)



def test_norm_B():
    rd = VariableTimeSteps()

    delta = np.array([4, -1])
    base = np.array([10, 2])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.5)

    base = np.array([10, 0])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.4)     # The zero baseline concentration is disregarded

    delta = np.array([300,        6, 0, -1])
    base = np.array([0.00000001, 10, 0,  2])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.6)

    delta = np.array([300,   6, 0, -1])
    base = np.array([0.001, 10, 0,  2])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 300000)

    base = np.array([2,   5, 5, 14, 14])
    delta = np.array([0.5, 1, 4, -2, -5])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.8)

    with pytest.raises(Exception):
        rd.norm_B(baseline_conc=np.array([1, 2, 3]), delta_conc=delta) # Too many entries in array

    with pytest.raises(Exception):
        rd.norm_B(baseline_conc=base, delta_conc=np.array([1]))        # Too few entries in array



def test_norm_C():
    rd = VariableTimeSteps()

    result = rd.norm_C(prev_conc=np.array([1]), baseline_conc=np.array([2]), delta_conc=np.array([0.5]))
    assert result == 0

    result = rd.norm_C(prev_conc=np.array([8]), baseline_conc=np.array([5]), delta_conc=np.array([1]))
    assert result == 0

    result = rd.norm_C(prev_conc=np.array([8]), baseline_conc=np.array([5]), delta_conc=np.array([4]))
    assert np.isclose(result, 4/3)

    result = rd.norm_C(prev_conc=np.array([10]), baseline_conc=np.array([14]), delta_conc=np.array([-2]))
    assert result == 0

    result = rd.norm_C(prev_conc=np.array([10]), baseline_conc=np.array([14]), delta_conc=np.array([-5]))
    assert np.isclose(result, 5/4)

    prev =     np.array([1,   8, 8, 10, 10])
    baseline = np.array([2,   5, 5, 14, 14])
    delta =    np.array([0.5, 1, 4, -2, -5])
    result = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3)

    # A scenario where the 'prev' and 'baseline' values are almost identical
    prev = np.append(prev, 3)
    baseline = np.append(baseline, 2.999999999)
    delta = np.append(delta, 8)
    result = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3)

    # A scenario where the 'delta' dwarfs the change between 'prev' and 'baseline'
    prev = np.append(prev, 10)
    baseline = np.append(baseline, 10.05)
    delta = np.append(delta, -9)
    result = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3)

    prev = np.append(prev, 10)
    baseline = np.append(baseline, 10.2)
    delta = np.append(delta, -9)
    result = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 45)



def test_norm_D():
    rd = VariableTimeSteps()

    prev =     np.array([ 12.96672432,  31.10067726,  55.93259842,  44.72389482, 955.27610518])
    baseline = np.array([ 12.99244738,  31.04428765,  55.96326497,  43.91117372, 956.08882628])
    delta =    np.array([-2.56160549,   -0.03542113,   2.59702662,   1.36089356,  -1.36089356])
    result = rd.norm_D(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 37.64942285399873)
