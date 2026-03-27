import pytest
import numpy as np
import math
from life123.reaction_kinetics import ReactionKinetics, VariableTimeSteps
from life123 import ThermoDynamics, ReactionGeneric




def test_exact_advance_unimolecular_reversible():
    # Reaction A <-> P
    p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0)
    assert np.allclose(p, 10.)      # No change

    incr = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0, incremental=True)
    assert np.allclose(incr, 0)     # No change


    p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.005)
    assert np.allclose(p, 11.08636387)

    p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.2)
    assert np.allclose(p, 37.81330458845654)

    p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.31739)
    assert np.allclose(p, 45.)

    p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=1.12)
    assert np.allclose(p, 53.837294)

    incr = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=1.12, incremental=True)
    assert np.allclose(incr, 53.837294-10)      # P(t) - P0


    p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=100.)
    assert np.allclose(p, 54.)


def test_exact_advance_unimolecular_reversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law,
    # at the middle point of 3
    h = 0.0005
    t_start = 0.2
    p1 = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=t_start)
    p2 = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=t_start + h)
    a2 = 90 - p2    # From mass conservation
    p3 = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=t_start + 2 * h)
    #print(p2)
    gradient = np.gradient([p1, p2, p3], h)
    derivative_p2 = gradient[1]
    rate_at_mid_point = 3 * a2 - 2 * p2     # kF * a2 - kR * p2
    assert np.allclose(derivative_p2, rate_at_mid_point)


def test_exact_advance_unimolecular_reversible_3():
    # Compare the exact solution against a fine-grained forward Euler approximation
    rxn = ReactionGeneric(reactants="A", products="P", kF=3., kR=2.)

    t_final = 0.2
    n_steps = 50000
    t_step = t_final / n_steps
    a = 80
    p = 10
    for i in range(n_steps):
        increment_dict, _ = \
            rxn.step_simulation(delta_time=t_step, conc_dict={"A": a, "P": p})
        delta_p = increment_dict["P"]
        p += delta_p
        a -= delta_p

    exact_p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.2)   # 37.81330458845654
    assert np.allclose(p, exact_p)




def test_exact_advance_unimolecular_irreversible():
    # Reaction A -> P
    # Compare against the REVERSIBLE reaction solver with zero reverse rate constant
    p = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0)
    assert np.allclose(p, 10.)      # No change

    incr = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0, incremental=True)
    assert np.allclose(incr, 0)     # No change

    p = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0.005)
    assert np.allclose(p, 11.191044831754994)

    p = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0.31739)
    assert np.allclose(p, 59.127783572982025)

    incr = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0.31739, incremental=True)
    assert np.allclose(incr, p-10)      # P(t) - P0

    p = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=1.12)
    assert np.allclose(p, 87.22117928442091)

    p = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=100.)
    assert np.allclose(p, 90)

    incr = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=100., incremental=True)
    assert np.allclose(incr, 80)


def test_exact_advance_unimolecular_irreversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law, at the middle point
    h = 0.0005
    t_start = 0.3
    p1 = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=t_start)
    p2 = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=t_start + h)
    a2 = 90 - p2    # From mass conservation
    p3 = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=t_start + 2 * h)

    gradient = np.gradient([p1, p2, p3], h)
    derivative_p2 = gradient[1]
    rate_at_mid_point = 3 * a2   # kF * a2
    assert np.allclose(derivative_p2, rate_at_mid_point)




def test_approx_solution_synthesis_rxn():
    # Reaction A + B <-> P
    kF=5.
    kR=2.
    A0=10.
    B0=50.
    P0=20.

    with pytest.raises(Exception):
        ReactionKinetics.approx_solution_synthesis_rxn(kF=0, kR=kR, A0=A0, B0=B0, P0=P0, t=0)   # kF cannot be zero

    p = ReactionKinetics.approx_solution_synthesis_rxn(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=0)
    assert np.allclose(p, P0)      # No change from P0

    #equil = ReactionKinetics.compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0)
    #print(equil)       # {'A': 0.2948774087575341, 'B': 40.294877408757536, 'P': 29.705122591242464}

    p = ReactionKinetics.approx_solution_synthesis_rxn(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=0.1)
    assert np.allclose(p, 29.705122591242464)       # Reaching the equilibrium value after a sufficiently long time


    # Process a whole array of times in one call
    t = np.array([0, 0.000864, 0.001555, 0.009850, 0.067400])
    result = ReactionKinetics.approx_solution_synthesis_rxn(kF=5., kR=2., A0=10., B0=50., P0=20., t=t)
    assert np.allclose(result, [20., 22.11257633, 23.46603241, 29.11416425, 29.70512254])

    # Compare against exact solution
    for t in [0, 0.000864, 0.001555, 0.009850, 0.067400]:
        p_approx = ReactionKinetics.approx_solution_synthesis_rxn(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
        p_exact = ReactionKinetics.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
        assert abs(p_approx - p_exact) / p_exact < 0.02     # Less than 2% discrepancy



def test_exact_advance_synthesis_irreversible():
    # Reaction A + B -> P

    # We'll start with `A` as the limiting reagent
    # No change at time 0
    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0, incremental=False)
    assert np.allclose(P_t, 20.)

    Delta_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0, incremental=True)
    assert np.allclose(Delta_t, 0)

    # Compare against the REVERSIBLE reaction solver with a near-zero reverse rate constant : not exactly zero, because its
    # implementation reverts to the irreversible solver when the reverse rate constant is very close to zero
    eps = 0.00001
    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0.01, incremental=False)
    c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=10., B0=50., P0=20., t=0.01, incremental=False)
    assert np.allclose(P_t, c)  # 28.8871974442945

    Delta_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0.01, incremental=True)
    D_c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=10., B0=50., P0=20., t=0.01, incremental=True)
    assert np.allclose(Delta_t, D_c)

    for t  in [0.0003, 0.0004, 0.5, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.008, 0.01, 1.]:
        P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=t, incremental=False)
        c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=10., B0=50., P0=20., t=t, incremental=False)
        assert np.allclose(P_t, c)


    # A time so large that the reaction has gone to completion (`A` being the limiting reagent)
    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=3., incremental=False)
    assert np.allclose(P_t, 30.)    # All `A`(the limiting reagent) converted to `P`

    # A value of t so ridiculously large that an OverflowError is caused (but caught) in the internal math.exp usage
    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=200., incremental=False)
    assert np.allclose(P_t, 30.)    # All `A`(the limiting reagent) converted to `P`


    # Now, let `B` be the limiting reagent at large times
    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=5., P0=20., t=3., incremental=False)
    assert np.allclose(P_t, 25.)    # All `B`(the limiting reagent) converted to `P`

    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=5., P0=20., t=200., incremental=False)
    assert np.allclose(P_t, 25.)    # All `B`(the limiting reagent) converted to `P`


    # If either A0 or B0 is zero, the reaction doesn't proceed
    Delta_P = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=0, B0=50., P0=20., t=1., incremental=True)
    assert np.allclose(Delta_P, 0)

    Delta_P = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=0, P0=20., t=1., incremental=True)
    assert np.allclose(Delta_P, 0)


    # Case A0 = B0  (equivalently, 2 A -> P)

    # No change at time 0
    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30, P0=20, t=0, incremental=False)
    assert np.allclose(P_t, 20)     # No change

    for t  in [0.0003, 0.0004, 0.5, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.008, 0.01, 1.]:
        P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30, P0=20, t=t, incremental=False)
        P_t_approx = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30.000001, P0=20, t=t, incremental=False)
        c = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=30, B0=30, P0=20, t=t, incremental=False)
        assert np.allclose(P_t, c)
        assert np.allclose(P_t, P_t_approx)

    # A time so large that the reaction has gone to completion
    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30, P0=20, t=1000., incremental=False)
    assert np.allclose(P_t, 50.)    # The reagents fully converted to `P`


def test_exact_advance_synthesis_irreversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law,
    # at the middle point of 3
    h = 0.00001
    t_start = 0.003
    p1 = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50, P0=20., t=t_start)
    p2 = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50, P0=20., t=t_start + h)
    a2 = 10 - (p2 - 20)    # From mass conservation
    b2 = 50 - (p2 - 20)    # From mass conservation
    p3 = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50, P0=20., t=t_start + 2 * h)
    #print(a2, b2, p2)
    gradient = np.gradient([p1, p2, p3], h)
    derivative_p2 = gradient[1]
    rate_at_mid_point = 5 * a2 * b2   # kF * a2 * b2
    assert np.allclose(derivative_p2, rate_at_mid_point)


def test_exact_advance_synthesis_irreversible_3():
    # Compare the exact solution against a fine-grained forward Euler approximation
    # Reaction A + B -> P
    rxn = ReactionGeneric(reactants=["A", "B"], products="P", kF=5., kR=0)

    t_final = 0.01
    n_steps = 10000
    t_step = t_final / n_steps
    a = 10
    b = 50
    p = 20
    for i in range(n_steps):
        increment_dict, _ = \
            rxn.step_simulation(delta_time=t_step, conc_dict={"A": a, "B": b, "P": p})
        delta_p = increment_dict["P"]
        p += delta_p
        a -= delta_p
        b -= delta_p

    exact_p = ReactionKinetics.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0.01, incremental=False)   # 28.8871974442945
    assert np.allclose(p, exact_p)



def test_exact_advance_synthesis_reversible():
    # Reaction A + B <-> P

    Dp = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0, incremental=True)
    assert np.allclose(Dp, 0.)     # No change

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0, incremental=False)
    assert np.allclose(P_t, 20.)    # No change

    # The comparison values below are from Octave, v. 10.3.0  ;  for example:
    '''
    lsode_options("absolute tolerance", 1e-12);
    lsode_options("relative tolerance", 1e-10);
    
    kf = 5.0;
    kr = 2.0;
    a0 = 10.0;
    b0 = 50.0;
    p0 = 20.0;
    
    p_init = p0;
    
    t = [0, 0.0001, 0.0003, 0.0004, 0.0005, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.008, 0.01, 1.];
    
    function dpdt = bimol_ode(p, t, kf, kr, a0, b0, p0)
      dpdt = kf * (a0 - p + p0) .* (b0 - p + p0) - kr * p;
    end
    
    p = lsode(@(p,t) bimol_ode(p,t,kf,kr,a0,b0,p0), p_init, t);
    
    format long
    p
    '''
    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0001, incremental=False)
    assert np.allclose(P_t, 20.24233230049967)

    Dp = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0001, incremental=True)
    assert np.allclose(Dp, 0.24233230049967)    # This is a DELTA value

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0003, incremental=False)
    assert np.allclose(P_t, 20.70580470925962)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0004, incremental=False)
    assert np.allclose(P_t, 20.92746192779109)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0005, incremental=False)
    assert np.allclose(P_t, 21.14272449451510)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0007, incremental=False)
    assert np.allclose(P_t, 21.55497320836003)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.001, incremental=False)
    assert np.allclose(P_t, 22.13079644707523)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0015, incremental=False)
    assert np.allclose(P_t, 22.98939523899655)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.002, incremental=False)
    assert np.allclose(P_t, 23.73870896747882)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.003, incremental=False)
    assert np.allclose(P_t, 24.97201963453988)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.005, incremental=False)
    assert np.allclose(P_t, 26.68110034281166)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.008, incremental=False)
    assert np.allclose(P_t, 28.12354983857080)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.01, incremental=False)
    assert np.allclose(P_t, 28.66884983321557)

    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=1., incremental=False)
    assert np.allclose(P_t, 29.70512259124693)


    # Verify reaching equilibrium at a large t
    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=10., incremental=False)
    eq = ReactionKinetics.compute_equilibrium_conc_first_order(kF=5., kR=2., a=1, A0=10., b=1, B0=50., p=1, P0=20.)
    assert np.allclose(P_t, eq["P"])    # 29.705122591242464


    # Case A0 = B0
    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=35, B0=35, P0=20., t=0.002, incremental=False)
    assert np.allclose(P_t, 28.99659131142920)


    # No reverse reaction
    P_t = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=0, A0=30, B0=30, P0=20, t=0.0003, incremental=False)
    assert np.allclose(P_t, 21.29186602834257)


def test_exact_advance_synthesis_reversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law,
    # at the middle point of 3

    # Case where `A` is the limiting reagent
    A0=10
    B0=50
    P0=20
    kF=5
    kR=2

    #print(ReactionKinetics.compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0))
    #Equilibrium at {'A': 0.2948774087575341, 'B': 40.294877408757536, 'P': 29.705122591242464}

    h = 0.00001
    t_start = 0.002
    times = [t_start + i * h for i in range(3)]
    p = [ReactionKinetics.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
            for t in times
        ]
    p_middle = p[1]
    a_middle = A0 - (p_middle - P0)    # From mass conservation
    b_middle = B0 - (p_middle - P0)    # From mass conservation

    print(a_middle, b_middle, p)        # 6.247301538717998 46.247301538718 [23.738708979042038, 23.752698461282, 23.766650993990154]
    gradient = np.gradient(p, h)
    derivative_middle = gradient[1]
    rate_at_mid_point = kF * a_middle * b_middle - kR * p_middle  # kF [A] [B] - kR [P]
    assert np.allclose(derivative_middle, rate_at_mid_point)


    # Case where `B` is the limiting reagent
    A0=60
    B0=30
    P0=10
    kF=7
    kR=3

    #print(ReactionKinetics.compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0))
    #Equilibrium at {'A': 30.553318635734474, 'B': 0.5533186357344739, 'P': 39.44668136426553}

    h = 0.00001
    t_start = 0.002
    times = [t_start + i * h for i in range(3)]
    p = [ReactionKinetics.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
            for t in times
        ]
    p_middle = p[1]
    a_middle = A0 - (p_middle - P0)    # From mass conservation
    b_middle = B0 - (p_middle - P0)    # From mass conservation

    #print(a_middle, b_middle, p)   44.710779935913656 14.710779935913656 [25.243842573299553, 25.289220064086344, 25.334407842935416]
    gradient = np.gradient(p, h)
    derivative_middle = gradient[1]
    rate_at_mid_point = kF * a_middle * b_middle - kR * p_middle  # kF [A] [B] - kR [P]
    assert np.allclose(derivative_middle, rate_at_mid_point)


    # Case where A0 = B0 (and reaction mostly in reverse)
    A0=8
    B0=8
    P0=30
    kF=5
    kR=8

    #print(ReactionKinetics.compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0))
    #Equilibrium at {'A': 7.038367176906169, 'B': 7.038367176906169, 'P': 30.96163282309383}

    h = 0.00005
    t_start = 0.01
    times = [t_start + i * h for i in range(3)]
    #print(times)
    p = [ReactionKinetics.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
            for t in times
        ]
    p_middle = p[1]
    a_middle = A0 - (p_middle - P0)    # From mass conservation
    b_middle = B0 - (p_middle - P0)    # From mass conservation

    #print(a_middle, b_middle, p)    # 7.461626205565043 7.461626205565043 [30.536666653790505, 30.538373794434957, 30.54007389774443]
    gradient = np.gradient(p, h)
    derivative_middle = gradient[1]
    rate_at_mid_point = kF * a_middle * b_middle - kR * p_middle  # kF [A] [B] - kR [P]
    assert np.allclose(derivative_middle, rate_at_mid_point)


def test_exact_advance_synthesis_reversible_3():
    # Compare the exact solution against a fine-grained forward Euler approximation
    # Reaction A + B <-> P
    rxn = ReactionGeneric(reactants=["A", "B"], products="P", kF=5., kR=2.)

    t_final = 0.001
    n_steps = 2000
    t_step = t_final / n_steps
    a = 10
    b = 50
    p = 20
    for i in range(n_steps):
        increment_dict, _ = \
            rxn.step_simulation(delta_time=t_step, conc_dict={"A": a, "B": b, "P": p})
        delta_p = increment_dict["P"]
        p += delta_p
        a -= delta_p
        b -= delta_p

    exact_p = ReactionKinetics.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.001, incremental=False)   # 22.130796453845
    assert np.allclose(p, exact_p)



def test_half_time_unimolecular_irreversible():
    # Reaction A -> P , with kF=5
    kF=5.
    half_time = ReactionKinetics.half_time_unimolecular_irreversible(kF=kF)
    assert np.allclose(half_time, math.log(2) / kF)

    p_halftime = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=kF, A0=80., P0=10., t=half_time)
    a_halftime = 80 - (p_halftime - 10)     # From the stoichiometry
    assert np.allclose(a_halftime, 80./2)   # [A] has indeed dropped in half after half_time has elapsed



def test_half_time_relaxation_unimolecular_reversible():
    # Reaction A <-> P , with kF=8 and kR=2
    kF=8.
    kR=2.
    half_time_relaxation = ReactionKinetics.half_time_relaxation_unimolecular_reversible(kF=kF, kR=kR)
    assert np.allclose(half_time_relaxation, math.log(2) / (kF + kR))

    # Simulate the reaction to t = half_time_relaxation
    p_halftime = ReactionKinetics.exact_advance_unimolecular_reversible(kF=kF, kR=kR, A0=80., P0=10., t=half_time_relaxation)
    a_halftime = 80 - (p_halftime - 10)      # From the stoichiometry

    # Determine the equilibrium concentrations (the other reference point for the halfway drop)
    equil_concs = ReactionKinetics.compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=80., p=1, P0=10.)

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

    #equil_concs = ReactionKinetics.compute_equilibrium_conc_first_order(kF=kF, kR=0, a=1, A0=15., b=1, B0=5, p=1, P0=0)
    #print(equil_concs)

    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=kF, A0=15., B0=5., P0=0, t=half_time_relaxation)
    assert np.allclose(P_t, 2.5)    # 2.5 is halfway between P0=0 and the final value of 5,
                                    # when all the limiting reagent (B) has been consumed


    half_time_relaxation = ReactionKinetics.half_time_to_equilibrium_irreversible_synthesis(kF=kF, A0=20, B0=20)
    print(half_time_relaxation)
    assert np.allclose(half_time_relaxation, 0.00625)
    P_t = ReactionKinetics.exact_advance_synthesis_irreversible(kF=kF, A0=20, B0=20, P0=0, t=half_time_relaxation)
    print(P_t)
    assert np.allclose(P_t, 10)     # 10 is halfway between P0=0 and the final value of 20,
                                    # when the reagents have been consumed



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



def test_compute_equilibrium_conc_first_order():
    # Unimolecular reaction A <-> P
    with pytest.raises(Exception):
        ReactionKinetics.compute_equilibrium_conc_first_order(kF=3., kR=2., a=1, p=1, A0=80., P0=10.,
                                                              b=0, B0=123., q=0, Q0=999.)      # Bogus values for B0 and D0

    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=3., kR=2., a=1, p=1, A0=80., P0=10.)
    assert np.allclose(result["A"], 36)
    assert np.allclose(result["P"], 54)
    assert "B" not in result
    assert "Q" not in result

    # Further verify
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(1,"A")], product_data=[(1,"P")],
                                                 conc={"A": 36, "P": 54})
    assert np.allclose(K, 3./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # Only the forward reaction
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=3., kR=0, a=1, p=1, A0=80., P0=10.)
    assert np.allclose(result["A"], 0)
    assert np.allclose(result["P"], 90)
    assert "B" not in result
    assert "Q" not in result

    # Only the reverse reaction
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=0, kR=2., a=1, p=1, A0=80., P0=10.)
    assert np.allclose(result["A"], 90)
    assert np.allclose(result["P"], 0)
    assert "B" not in result
    assert "Q" not in result

    # A + B <-> P
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=5., kR=2., a=1, b=1, p=1, A0=10., B0=50., P0=20.)
    assert np.allclose(result["A"], 0.2948774087575341)
    assert np.allclose(result["B"], 40.294877408757536)
    assert np.allclose(result["P"], 29.705122591242464)
    assert "Q" not in result
    # Further verify
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(1,"A"), (1,"B")], product_data=[(1,"P")],
                                                   conc={"A": 0.2948774087575341, "B": 40.294877408757536, "P": 29.705122591242464})
    assert np.allclose(K, 5./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # 2 A <-> P  (expressed as A + A <-> P)
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=3., kR=2., a=2, b=2, p=1, A0=200., B0=200., P0=40.)
    assert np.allclose(result["A"], 9.49568869375716)
    assert np.allclose(result["P"], 135.2521556531214)
    assert "Q" not in result

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(2,"A")], product_data=[(1,"P")],
                                                 conc={"A": 9.49568869375716, "P": 135.2521556531214})
    assert np.allclose(K, 3./2.)    # At the equilibrium concentrations, the reaction quotient of elementary equations equals kF/kR

    # HYPOTHETICAL 1st order reaction A + B <-> P + Q
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=10., kR=2., a=1, b=1, p=1, q=1, A0=1., B0=2., P0=3., Q0=4.)
    expected = {'A': 1.0894541729001368, 'B': 2.089454172900137, 'P': 2.910545827099863, 'Q': 3.910545827099863}
    for key in result.keys():
        assert np.allclose(result[key], expected[key])

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(1,"A"), (1,"B")], product_data=[(1,"P"), (1,"Q")],
                                                   conc=expected)
    assert np.allclose(K, 10./2.)    # Non-elementary reaction, but still with kinetics following mass-action: the equilibrium constant equals to kF/kR


    # HYPOTHETICAL reaction A + 2 B <-> P + 3 Q, 1st order in all chemicals
    # This is a scenario of generalized kinetics (non-mass-action with respect to the stoichiometry),
    result = ReactionKinetics.compute_equilibrium_conc_first_order(kF=10., kR=2., a=1, b=2, p=1, q=3, A0=1., B0=2., P0=3., Q0=4.)
    expected = {'A': 1.0598463308126616, 'B': 2.119692661625323, 'P': 2.9401536691873384, 'Q': 3.8204610075620153}
    for key in result.keys():
        assert np.allclose(result[key], expected[key])

    # Further verify.  Compute the reaction quotient with the equilibrium concentrations, i.e. the thermodynamic equilibrium constant
    K = ThermoDynamics.compute_reaction_quotient(reactant_data=[(1,"A"), (2,"B")], product_data=[(1,"P"), (3,"Q")],
                                                 conc=expected)
    corrective_factor = (result["B"] / result["Q"]**2)      # Needed because we're dealing with generalized kinetics
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
