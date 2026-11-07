import numpy as np
import math
from life123.kinetics.reaction_kinetics import ReactionKinetics
from life123.reaction_simulator import AnalyticReactionSolver
from life123.kinetics.mass_action_analysis import MassAction_Analysis




def test_half_time_unimolecular_irreversible():
    # Reaction A -> P , with kF=5
    kF=5.
    half_time = ReactionKinetics.half_time_unimolecular_irreversible(kF=kF)
    assert np.allclose(half_time, math.log(2) / kF)

    p_halftime = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=kF, R0=80., P0=10., t=half_time)
    a_halftime = 80 - (p_halftime - 10)     # From the stoichiometry
    assert np.allclose(a_halftime, 80./2)   # [A] has indeed dropped in half after half_time has elapsed



def test_half_time_relaxation_unimolecular_reversible():
    # Reaction A <-> P , with kF=8 and kR=2
    kF=8.
    kR=2.
    half_time_relaxation = ReactionKinetics.half_time_relaxation_unimolecular_reversible(kF=kF, kR=kR)
    assert np.allclose(half_time_relaxation, math.log(2) / (kF + kR))

    # Simulate the reaction to t = half_time_relaxation
    p_halftime = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=kF, kR=kR, A0=80., P0=10., t=half_time_relaxation)
    a_halftime = 80 - (p_halftime - 10)      # From the stoichiometry

    # Determine the equilibrium concentrations (the other reference point for the halfway drop)
    equil_concs = MassAction_Analysis._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=80., p=1, P0=10.)

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

    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=kF, A0=15., B0=5., P0=0, t=half_time_relaxation)
    assert np.allclose(P_t, 2.5)    # 2.5 is halfway between P0=0 and the final value of 5,
                                    # when all the limiting reagent (B) has been consumed


    half_time_relaxation = ReactionKinetics.half_time_to_equilibrium_irreversible_synthesis(kF=kF, A0=20, B0=20)
    print(half_time_relaxation)
    assert np.allclose(half_time_relaxation, 0.00625)
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=kF, A0=20, B0=20, P0=0, t=half_time_relaxation)
    print(P_t)
    assert np.allclose(P_t, 10)     # 10 is halfway between P0=0 and the final value of 20,
                                    # when the reagents have been consumed
