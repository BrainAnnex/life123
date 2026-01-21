import pytest
import numpy as np
import math
from life123.reaction_kinetics import ReactionKinetics, VariableTimeSteps
from life123 import ThermoDynamics



def test_half_time_unimolecular_irreversible():
    # Reaction A -> B , with kF=5
    kF=5.
    half_time = ReactionKinetics.half_time_unimolecular_irreversible(kF=kF)
    assert np.allclose(half_time, math.log(2) / kF)
    a, _ = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=kF, A0=80., B0=10., t=half_time)
    assert np.allclose(a, 80./2)    # [A] has indeed dropped in half after half_time has elapsed



def test_half_time_relaxation_unimolecular_reversible():
    pass



def test_relaxation_time_unimolecular_reversible():
    relaxation_time = ReactionKinetics.relaxation_time_unimolecular_reversible(kF=8., kR=2.)
    assert np.allclose(relaxation_time, 0.1)    # 1 / (8.+2.)



def test_exact_advance_unimolecular_irreversible():
    # Reaction A -> B
    a, b = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., B0=10., t=0)
    assert np.allclose(a, 80.)
    assert np.allclose(b, 10.)

    incr = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., B0=10., t=0, incremental=True)
    assert np.allclose(incr, [0,0])


    a, b = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., B0=10., t=0.005)
    assert np.allclose(a, 78.808955168245)
    assert np.allclose(b, 11.191044831754994)
    assert np.allclose(a+b, 90)

    a, b = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., B0=10., t=0.31739)
    assert np.allclose(a, 30.87221642701797)
    assert np.allclose(b, 59.127783572982025)
    assert np.allclose(a+b, 90)

    incr = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., B0=10., t=0.31739, incremental=True)
    assert np.allclose(incr, [-49.127783572982025 , 49.127783572982025])
    assert np.allclose(incr[0] + incr[1], 0)


    a, b = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., B0=10., t=1.12)
    assert np.allclose(a, 2.778820715579084)
    assert np.allclose(b, 87.22117928442091)
    assert np.allclose(a+b, 90)

    a, b = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., B0=10., t=100.)
    assert np.allclose(a, 0)
    assert np.allclose(b, 90)

    incr = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., B0=10., t=100., incremental=True)
    assert np.allclose(incr, [-80, 80])
    assert np.allclose(incr[0] + incr[1], 0)



def test_exact_advance_unimolecular_reversible():
    # Reaction A <-> B
    a, b = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., B0=10., t=0)
    assert np.allclose(a, 80.)
    assert np.allclose(b, 10.)

    incr = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., B0=10., t=0, incremental=True)
    assert np.allclose(incr, [0,0])


    a, b = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., B0=10., t=0.005)
    assert np.allclose(a, 78.91363613)
    assert np.allclose(b, 11.08636387)
    assert np.allclose(a+b, 90)

    a, b = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., B0=10., t=0.31739)
    assert np.allclose(a, 45.)
    assert np.allclose(b, 45.)
    assert np.allclose(a+b, 90)

    a, b = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., B0=10., t=1.12)
    assert np.allclose(a, 36.162706)
    assert np.allclose(b, 53.837294)
    assert np.allclose(a+b, 90)

    incr = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., B0=10., t=1.12, incremental=True)
    assert np.allclose(incr, [-43.837294, 43.837294])
    assert np.allclose(incr[0] + incr[1], 0)


    a, b = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., B0=10., t=100.)
    assert np.allclose(a, 36.)
    assert np.allclose(b, 54.)
    assert np.allclose(a+b, 90)

    incr = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., B0=10., t=100., incremental=True)
    assert np.allclose(incr, [-44., 44.])
    assert np.allclose(incr[0] + incr[1], 0)



def test_approx_solution_synthesis_rxn():
    #rxn = Reaction(reactants = ["A", "B"], products="C", forward_rate=3., reverse_rate=2.)

    t = np.array([0, 0.000864, 0.001555, 0.009850, 0.067400])

    result = ReactionKinetics.approx_solution_synthesis_rxn(kF=5., kR=2., A0=10., B0=50., C0=20., t_arr=t)

    assert np.allclose(result[0], [10., 7.88742367,   6.53396759,  0.88583575,  0.29487746])
    assert np.allclose(result[1], [50., 47.88742367, 46.53396759, 40.88583575, 40.29487746])
    assert np.allclose(result[2], [20., 22.11257633, 23.46603241, 29.11416425, 29.70512254])


def test_exact_advance_synthesis_reversible():
    # Reaction A + B <-> C      TODO: Look into testing 2A -> C
    D_a, D_b, D_c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0, incremental=True)
    assert np.allclose(D_a, 0.)
    assert np.allclose(D_b, 0.)
    assert np.allclose(D_c, 0.)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0, incremental=False)
    assert np.allclose(a, 10.)
    assert np.allclose(b, 50.)
    assert np.allclose(c, 20.)


    D_a, D_b, D_c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.0001, incremental=True)
    assert np.allclose(D_a, -D_c)
    assert np.allclose(D_b, -D_c)
    assert np.allclose(D_c, 0.24233229989)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.0001, incremental=False)
    assert np.allclose(a, 10.-D_c)
    assert np.allclose(b, 50.-D_c)
    assert np.allclose(c, 20.24233229989)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.0003, incremental=False)
    assert np.allclose(a+c, 30.)    # From stoichiometry
    assert np.allclose(b+c, 70.)    # From stoichiometry
    assert np.allclose(c, 20.70580471094)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.0004, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 20.927461929486)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.0005, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 21.14272449633)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.0007, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 21.554973213)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.001, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 22.130796453845)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.0015, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 22.989395248289)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.002, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 23.73870897904)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.003, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 24.9720196494)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.005, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 26.681100361833)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.008, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 28.12354985934)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=0.01, incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 28.6688498538)

    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=1., incremental=False)
    assert np.allclose(a+c, 30.)
    assert np.allclose(b+c, 70.)
    assert np.allclose(c, 29.705122591242)

    # Verify reaching equilibrium at a large t
    a, b, c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., C0=20., t=10., incremental=False)
    eq = ReactionKinetics.compute_equilibrium_conc_first_order(kF=5., kR=2., a=1, A0=10., b=1, B0=50., c=1, C0=20.)
    assert np.allclose(a, eq["A"])
    assert np.allclose(b, eq["B"])
    assert np.allclose(c, eq["C"])



def test_estimate_rate_constants_simple():
    pass    # TODO


def test_estimate_rate_constants_synthesis():
    pass    # TODO



def test_compute_rate_elementary():
    # All reactions below are elementary reactions
    # that have 1st-order kinetics with respect to all the involved chemicals

    # Reaction A <-> B
    result = ReactionKinetics.compute_rate_elementary(reactants=["A"], products=["B"],
                                                      kF=20., kR=2., reversible=True,
                                                      conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5. - 2. * 8.)  # 84.0

    result = ReactionKinetics.compute_rate_elementary(reactants=["A"], products=["B"],
                                                      kF=20., kR=2., reversible=False,
                                                      conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 20. * 5)             # 100.0


    # Reaction A + B <-> C + D , with 1st-order kinetics for each species
    result = ReactionKinetics.compute_rate_elementary(reactants=["A", "B"], products=["C", "D"],
                                                      kF=5., kR=2., reversible=True,
                                                      conc_dict={"A": 3.5, "B": 9., "C": 11., "D": 7.})
    assert np.allclose(result,  5. * 3.5 * 9. - 2. * 11. * 7.)  # 3.5

    result = ReactionKinetics.compute_rate_elementary(reactants=["A", "B"], products=["C", "D"],
                                                      kF=5., kR=2., reversible=True,
                                                      conc_dict={"A": 5., "B": 8., "C": 15., "D": 7.})
    assert np.allclose(result,  -10.)




def test_compute_rate_mass_action_kinetics():

    # Reaction  2A <-> B , with 2nd-ORDER kinetics in the forward direction
    result = ReactionKinetics.compute_rate_mass_action_kinetics(reactant_terms=[(2, "A")], product_terms=[(1, "B")],
                                                    kF=5., kR=2.,
                                                    conc_dict={"A": 4.5, "B": 6.})
    assert np.allclose(result, 5. * 4.5 **2 - 2. * 6.)      # 89.25

    result = ReactionKinetics.compute_rate_mass_action_kinetics(reactant_terms=[(2, "A")], product_terms=[(1, "B")],
                                                    kF=5., kR=0,
                                                    conc_dict={"A": 4.5, "B": 6.})
    assert np.allclose(result, 5. * 4.5 **2)                # 101.25   (irreversible reaction)


    result = ReactionKinetics.compute_rate_mass_action_kinetics(reactant_terms=[(2, "A")], product_terms=[(1, "B")],
                                                    kF=3., kR=2.,
                                                    conc_dict={"A": 5., "B": 8.})
    assert np.allclose(result, 59.)


    # Reaction  B <-> 2C , with 2nd-ORDER kinetics in the reverse direction
    result = ReactionKinetics.compute_rate_mass_action_kinetics(reactant_terms=[(1, "B")], product_terms=[(2, "C")],
                                                kF=4., kR=2.,
                                                conc_dict={"B": 5., "C": 4})
    assert np.allclose(result, 4. * 5. - 2. * 4. **2)       # -12.0


    # Reaction A + B <-> C + D , with 1st-order kinetics for each species
    result = ReactionKinetics.compute_rate_mass_action_kinetics(reactant_terms=[(1, "A"), (1, "B")], product_terms=[(1, "C"), (1, "D")],
                                            kF=5., kR=2.,
                                            conc_dict={"A": 5., "B": 8., "C": 15., "D": 7.})
    assert np.allclose(result,  -10.)


def test_compute_rate_first_order():
    pass    #TODO



def test_compute_equilibrium_conc():
    # Unimolecular reaction A <-> C
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=3., kR=2., a=1, c=1, A0=80., C0=10.)
    assert np.allclose(result["A"], 36)
    assert np.allclose(result["C"], 54)
    # Further verify
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(1,"A")], product_data=[(1,"C")],
                                                 conc={"A": 36, "C": 54})
    assert np.allclose(K, 3./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # Only the forward reaction
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=3., kR=0, a=1, c=1, A0=80., C0=10.)
    assert np.allclose(result["A"], 0)
    assert np.allclose(result["C"], 90)

    # Only the reverse reaction
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=0, kR=2., a=1, c=1, A0=80., C0=10.)
    assert np.allclose(result["A"], 90)
    assert np.allclose(result["C"], 0)


    # A + B <-> C
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=5., kR=2., a=1, b=1, c=1, A0=10., B0=50., C0=20.)
    assert np.allclose(result["A"], 0.2948774087575341)
    assert np.allclose(result["B"], 40.294877408757536)
    assert np.allclose(result["C"], 29.705122591242464)
    # Further verify
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(1,"A"), (1,"B")], product_data=[(1,"C")],
                                                   conc={"A": 0.2948774087575341, "B": 40.294877408757536, "C": 29.705122591242464})
    assert np.allclose(K, 5./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # 2 A <-> C  (expressed as A + A <-> C)
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=3., kR=2., a=2, b=2, c=1, A0=200., B0=200., C0=40.)
    assert np.allclose(result["A"], 9.49568869375716)
    assert np.allclose(result["C"], 135.2521556531214)

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(2,"A")], product_data=[(1,"C")],
                                                 conc={"A": 9.49568869375716, "C": 135.2521556531214})
    assert np.allclose(K, 3./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # HYPOTHETICAL 1st order reaction A + B <-> C + D
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=10., kR=2., a=1, b=1, c=1, d=1, A0=1., B0=2., C0=3., D0=4.)
    expected = {'A': 1.0894541729001368, 'B': 2.089454172900137, 'C': 2.910545827099863, 'D': 3.910545827099863}
    for key in result.keys():
        assert np.allclose(result[key], expected[key])

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(1,"A"), (1,"B")], product_data=[(1,"C"), (1,"D")],
                                                   conc=expected)
    assert np.allclose(K, 10./2.)    # Non-elementary reaction, but still with kinetics following mass-action: the equilibrium constant equals to kF/kR


    # HYPOTHETICAL reaction A + 2 B <-> C + 3 D, 1st order in all chemicals
    # This is a scenario of generalized kinetics (non-mass-action with respect to the stoichiometry),
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=10., kR=2., a=1, b=2, c=1, d=3, A0=1., B0=2., C0=3., D0=4.)
    expected = {'A': 1.0598463308126616, 'B': 2.119692661625323, 'C': 2.9401536691873384, 'D': 3.8204610075620153}
    for key in result.keys():
        assert np.allclose(result[key], expected[key])

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(1,"A"), (2,"B")], product_data=[(1,"C"), (3,"D")],
                                                 conc=expected)
    corrective_factor = (result["B"] / result["D"]**2)      # Needed because we're dealing with generalized kinetics
                                                            # (non-mass-action with respect to the stoichiometry)
    assert np.allclose(K * corrective_factor, 10./2.)    # The "corrected" equilibrium constant equals to kF/kR






########    class VariableTimeSteps    ###########################################################################

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
