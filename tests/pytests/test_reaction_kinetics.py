import pytest
import numpy as np
import math
from life123.reaction_kinetics import ReactionKinetics
from life123.reaction_simulator import ReactionSimulator, AnalyticalReactionSolver
from life123.reactions import Stoichiometry




def test_half_time_unimolecular_irreversible():
    # Reaction A -> P , with kF=5
    kF=5.
    half_time = ReactionKinetics.half_time_unimolecular_irreversible(kF=kF)
    assert np.allclose(half_time, math.log(2) / kF)

    p_halftime = AnalyticalReactionSolver.exact_advance_unimolecular_irreversible(kF=kF, A0=80., P0=10., t=half_time)
    a_halftime = 80 - (p_halftime - 10)     # From the stoichiometry
    assert np.allclose(a_halftime, 80./2)   # [A] has indeed dropped in half after half_time has elapsed



def test_half_time_relaxation_unimolecular_reversible():
    # Reaction A <-> P , with kF=8 and kR=2
    kF=8.
    kR=2.
    half_time_relaxation = ReactionKinetics.half_time_relaxation_unimolecular_reversible(kF=kF, kR=kR)
    assert np.allclose(half_time_relaxation, math.log(2) / (kF + kR))

    # Simulate the reaction to t = half_time_relaxation
    p_halftime = AnalyticalReactionSolver.exact_advance_unimolecular_reversible(kF=kF, kR=kR, A0=80., P0=10., t=half_time_relaxation)
    a_halftime = 80 - (p_halftime - 10)      # From the stoichiometry

    # Determine the equilibrium concentrations (the other reference point for the halfway drop)
    equil_concs = ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=80., p=1, P0=10.)

    delta_a = 80. - equil_concs["A"]    # Change in [A] from initial state to equilibrium
    assert np.allclose(a_halftime, 80. - delta_a/2)    # [A] has indeed dropped by half of the overall delta_a after half_time has elapsed



def test_relaxation_time_unimolecular_reversible():
    relaxation_time = ReactionKinetics.relaxation_time_unimolecular_reversible(kF=8., kR=2.)
    assert np.allclose(relaxation_time, 0.1)    # 1 / (8.+2.)



def test_half_time_to_equilibrium_synthesis():
    # Reaction A + B -> P
    kF=8.

    half_time_relaxation = ReactionKinetics.half_time_to_equilibrium_irreversible_synthesis(kF=kF, A0=15, B0=5)
    assert np.allclose(half_time_relaxation, 0.006385320297074884)

    #equil_concs = ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=0, a=1, A0=15., b=1, B0=5, p=1, P0=0)
    #print(equil_concs)

    P_t = AnalyticalReactionSolver.exact_advance_synthesis_irreversible(kF=kF, A0=15., B0=5., P0=0, t=half_time_relaxation)
    assert np.allclose(P_t, 2.5)    # 2.5 is halfway between P0=0 and the final value of 5,
                                    # when all the limiting reagent (B) has been consumed


    half_time_relaxation = ReactionKinetics.half_time_to_equilibrium_irreversible_synthesis(kF=kF, A0=20, B0=20)
    print(half_time_relaxation)
    assert np.allclose(half_time_relaxation, 0.00625)
    P_t = AnalyticalReactionSolver.exact_advance_synthesis_irreversible(kF=kF, A0=20, B0=20, P0=0, t=half_time_relaxation)
    print(P_t)
    assert np.allclose(P_t, 10)     # 10 is halfway between P0=0 and the final value of 20,
                                    # when the reagents have been consumed



def test_estimate_rate_constants_simple():
    pass    # TODO


def test_estimate_rate_constants_synthesis():
    pass    # TODO




def test_compute_equilibrium_conc_elementary_decomposition():
    # A <-> 2 P
    result = ReactionKinetics.compute_equilibrium_conc_elementary_decomposition(kF=2, kR=3, A0=200., P0=40.)

    assert np.allclose(result["A"], 214.02745923365862)
    assert np.allclose(result["P"], 11.945081532682774)
    assert "B" not in result
    assert "Q" not in result

    # Verify from direct computation
    assert np.allclose(2 * 214.02745923365862 - 3 * 11.945081532682774 **2, 0)   # Net flow is zero at equilibrium
    st = Stoichiometry({"A": -1, "P": 2})
    st.consistency_checker(conc_before={"A": 200, "P": 40}, conc_after={'A': 214.02745923365862, 'P': 11.945081532682774})   # Stoichiometry check


    # A <-> 2 P
    result = ReactionKinetics.compute_equilibrium_conc_elementary_decomposition(kF=2, kR=3, A0=40., P0=200.)

    assert np.allclose(result["A"], 135.2521556531214)
    assert np.allclose(result["P"], 9.49568869375716)
    assert "B" not in result
    assert "Q" not in result

    # Verify from direct computation
    assert np.allclose(2 * 135.2521556531214 - 3 * 9.49568869375716 **2, 0)   # Net flow is zero at equilibrium
    st = Stoichiometry({"A": -1, "P": 2})
    st.consistency_checker(conc_before={"A": 40, "P": 200}, conc_after={'A': 135.2521556531214, 'P': 9.49568869375716})   # Stoichiometry check



def test_compute_equilibrium_conc_elementary_synthesis():

    # 2 A <-> P
    result = ReactionKinetics.compute_equilibrium_conc_elementary_synthesis(kF=3., kR=2., A0=200., P0=40.)
    assert np.allclose(result["A"], 9.49568869375716)
    assert np.allclose(result["P"], 135.2521556531214)
    assert "B" not in result
    assert "Q" not in result

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ReactionKinetics.compute_reaction_quotient(reactant_data=[(2,"A")], product_data=[(1,"P")],
                                                   conc={"A": 9.49568869375716, "P": 135.2521556531214})
    assert np.allclose(K, 3./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # Verify from direct computation
    assert np.allclose(3 * 9.49568869375716**2 - 2 * 135.2521556531214, 0)              # Net flow is zero at equilibrium
    st = Stoichiometry({"A": -2, "P": 1})
    st.consistency_checker(conc_before={"A": 200, "P": 40}, conc_after={"A": 9.49568869375716, "P": 135.2521556531214})   # Stoichiometry check



    # X <-> 2 Y  , with kF = 2, kR = 3, and X0 = 200, Y0 = 40
    # Reverse, to obtain 2 Y <-> X, with kF = 3, kR = 2 (reversed rates), and X0 = 200, Y0 = 40
    result = ReactionKinetics.compute_equilibrium_conc_elementary_synthesis(kF=3, kR=2, A0=200., P0=40.)    # TODO: ditch
    assert np.allclose(result["A"], 9.49568869375716)
    assert np.allclose(result["P"], 135.2521556531214)
    assert "B" not in result
    assert "Q" not in result

    # i.e. [X]eq = 135.2521556531214 , and [Y]eq = 9.49568869375716
    # Verify from direct computation
    assert np.allclose(2 * 135.2521556531214 - 3 * 9.49568869375716**2, 0)   # Net flow is zero at equilibrium



def test_compute_equilibrium_conc_mass_action():

    # Unimolecular reaction A <-> P
    result = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=3., kR=2., A0=80., P0=10.)
    assert np.allclose(result["A"], 36)
    assert np.allclose(result["P"], 54)
    assert "B" not in result
    assert "Q" not in result

    # Further verify
    K = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"A")], product_data=[(1,"P")],
                                                   conc={"A": 36, "P": 54})
    assert np.allclose(K, 3./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # Only the forward reaction
    result = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=3., kR=0, A0=80., P0=10.)
    assert np.allclose(result["A"], 0)
    assert np.allclose(result["P"], 90)
    assert "B" not in result
    assert "Q" not in result

    # Only the reverse reaction
    result = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=0, kR=2., A0=80., P0=10.)
    assert np.allclose(result["A"], 90)
    assert np.allclose(result["P"], 0)
    assert "B" not in result
    assert "Q" not in result

    # A + B <-> P
    result = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=5., kR=2., A0=10., B0=50., P0=20.)
    assert np.allclose(result["A"], 0.2948774087575341)
    assert np.allclose(result["B"], 40.294877408757536)
    assert np.allclose(result["P"], 29.705122591242464)
    assert "Q" not in result
    # Further verify
    K = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"A"), (1,"B")], product_data=[(1,"P")],
                                                   conc={"A": 0.2948774087575341, "B": 40.294877408757536, "P": 29.705122591242464})
    assert np.allclose(K, 5./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR



def test__compute_equilibrium_conc_first_order():
    # Unimolecular reaction A <-> P
    with pytest.raises(Exception):
        ReactionKinetics._compute_equilibrium_conc_first_order(kF=3., kR=2., a=1, p=1, A0=80., P0=10.,
                                                               b=0, B0=123., q=0, Q0=999.)      # Unexpected values for B0 and Q0

    result = ReactionKinetics._compute_equilibrium_conc_first_order(kF=3., kR=2., a=1, p=1, A0=80., P0=10.)
    assert np.allclose(result["A"], 36)
    assert np.allclose(result["P"], 54)
    assert "B" not in result
    assert "Q" not in result

    # Further verify
    K = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"A")], product_data=[(1,"P")],
                                                   conc={"A": 36, "P": 54})
    assert np.allclose(K, 3./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # Only the forward reaction
    result = ReactionKinetics._compute_equilibrium_conc_first_order(kF=3., kR=0, a=1, p=1, A0=80., P0=10.)
    assert np.allclose(result["A"], 0)
    assert np.allclose(result["P"], 90)
    assert "B" not in result
    assert "Q" not in result

    # Only the reverse reaction
    result = ReactionKinetics._compute_equilibrium_conc_first_order(kF=0, kR=2., a=1, p=1, A0=80., P0=10.)
    assert np.allclose(result["A"], 90)
    assert np.allclose(result["P"], 0)
    assert "B" not in result
    assert "Q" not in result

    # A + B <-> P
    result = ReactionKinetics._compute_equilibrium_conc_first_order(kF=5., kR=2., a=1, b=1, p=1, A0=10., B0=50., P0=20.)
    assert np.allclose(result["A"], 0.2948774087575341)
    assert np.allclose(result["B"], 40.294877408757536)
    assert np.allclose(result["P"], 29.705122591242464)
    assert "Q" not in result
    # Further verify
    K = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"A"), (1,"B")], product_data=[(1,"P")],
                                                   conc={"A": 0.2948774087575341, "B": 40.294877408757536, "P": 29.705122591242464})
    assert np.allclose(K, 5./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # 2 A <-> P  (expressed as A + A <-> P)
    result = ReactionKinetics._compute_equilibrium_conc_first_order(kF=3., kR=2., a=2, b=2, p=1, A0=200., B0=200., P0=40.)
    assert np.allclose(result["A"], 9.49568869375716)
    assert np.allclose(result["P"], 135.2521556531214)
    assert "Q" not in result

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ReactionKinetics.compute_reaction_quotient(reactant_data=[(2,"A")], product_data=[(1,"P")],
                                                 conc={"A": 9.49568869375716, "P": 135.2521556531214})
    assert np.allclose(K, 3./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # HYPOTHETICAL 1st order reaction A + B <-> P + Q
    result = ReactionKinetics._compute_equilibrium_conc_first_order(kF=10., kR=2., a=1, b=1, p=1, q=1, A0=1., B0=2., P0=3., Q0=4.)
    expected = {'A': 1.0894541729001368, 'B': 2.089454172900137, 'P': 2.910545827099863, 'Q': 3.910545827099863}
    for key in result.keys():
        assert np.allclose(result[key], expected[key])

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"A"), (1,"B")], product_data=[(1,"P"), (1,"Q")],
                                                   conc=expected)
    assert np.allclose(K, 10./2.)    # Non-elementary reaction, but still with kinetics following mass-action: the equilibrium constant equals to kF/kR


    # HYPOTHETICAL reaction A + 2 B <-> P + 3 Q, 1st order in all chemicals
    # This is a scenario of generalized kinetics (non-mass-action with respect to the stoichiometry),
    result = ReactionKinetics._compute_equilibrium_conc_first_order(kF=10., kR=2., a=1, b=2, p=1, q=3, A0=1., B0=2., P0=3., Q0=4.)
    expected = {'A': 1.0598463308126616, 'B': 2.119692661625323, 'P': 2.9401536691873384, 'Q': 3.8204610075620153}
    for key in result.keys():
        assert np.allclose(result[key], expected[key])

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"A"), (2,"B")], product_data=[(1,"P"), (3,"Q")],
                                                 conc=expected)
    corrective_factor = (result["B"] / result["Q"]**2)      # Needed because we're dealing with generalized kinetics
                                                            # (non-mass-action with respect to the stoichiometry)
    assert np.allclose(K * corrective_factor, 10./2.)    # The "corrected" equilibrium constant equals to kF/kR
