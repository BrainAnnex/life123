import pytest
import numpy as np
from life123.analytic_reaction_solver import AnalyticReactionSolver
from life123.reaction_kinetics import ReactionKinetics
from life123.species_registry import SpeciesRegistry
from life123.reactions import ReactionDefinition



def test_exact_advance_unimolecular_reversible():
    # Reaction A <-> P
    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0)
    assert np.allclose(p, 10.)      # No change

    incr = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0, incremental=True)
    assert np.allclose(incr, 0)     # No change


    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.005)
    assert np.allclose(p, 11.08636387)

    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.2)
    assert np.allclose(p, 37.81330458845654)

    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.31739)
    assert np.allclose(p, 45.)

    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=1.12)
    assert np.allclose(p, 53.837294)

    incr = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=1.12, incremental=True)
    assert np.allclose(incr, 53.837294-10)      # P(t) - P0


    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=100.)
    assert np.allclose(p, 54.)

    equil = ReactionKinetics._compute_equilibrium_conc_first_order(kF=3., kR=2., a=1, A0=80., p=1, P0=10.)
    assert np.allclose(p, equil["P"])


def test_exact_advance_unimolecular_reversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law,
    # at the middle point of 3 sampled points
    h = 0.0005
    t_start = 0.2
    p1 = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=t_start)
    p2 = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=t_start + h)
    a2 = 90 - p2    # From mass conservation
    p3 = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=t_start + 2 * h)
    #print(p2)
    gradient = np.gradient([p1, p2, p3], h)
    derivative_p2 = gradient[1]
    rate_at_mid_point = 3 * a2 - 2 * p2     # kF * a2 - kR * p2
    assert np.allclose(derivative_p2, rate_at_mid_point)


def test_exact_advance_unimolecular_reversible_3():
    # Compare the exact solution against a fine-grained forward-Euler approximation
    sr = SpeciesRegistry()
    rxn_defn = ReactionDefinition(reactants="A", products="P",
                                  species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF": 3., "kR": 2.})
    rxn = rxn_defn.sim_reactions[0]

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

    exact_p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.2)   # 37.81330458845654
    assert np.allclose(p, exact_p)



def test_exact_advance_unimolecular_irreversible():
    # Reaction A -> P
    # Compare against the REVERSIBLE reaction solver with zero reverse rate constant
    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0)
    assert np.allclose(p, 10.)      # No change

    incr = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0, incremental=True)
    assert np.allclose(incr, 0)     # No change

    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0.005)
    assert np.allclose(p, 11.191044831754994)

    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0.31739)
    assert np.allclose(p, 59.127783572982025)

    incr = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0.31739, incremental=True)
    assert np.allclose(incr, p-10)      # P(t) - P0

    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=1.12)
    assert np.allclose(p, 87.22117928442091)

    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=100.)
    assert np.allclose(p, 90)

    incr = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=100., incremental=True)
    assert np.allclose(incr, 80)


def test_exact_advance_unimolecular_irreversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law, at the middle point
    h = 0.0005
    t_start = 0.3
    p1 = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=t_start)
    p2 = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=t_start + h)
    a2 = 90 - p2    # From mass conservation
    p3 = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=t_start + 2 * h)

    gradient = np.gradient([p1, p2, p3], h)
    derivative_p2 = gradient[1]
    rate_at_mid_point = 3 * a2   # kF * a2
    assert np.allclose(derivative_p2, rate_at_mid_point)



def test_exact_advance_synthesis_reversible():
    # Reaction A + B <-> P

    # General case A0 != B0

    Dp = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0, incremental=True)
    assert np.allclose(Dp, 0.)     # No change

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0, incremental=False)
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
    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0001, incremental=False)
    assert np.allclose(P_t, 20.24233230049967)

    Dp = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0001, incremental=True)
    assert np.allclose(Dp, 0.24233229989273042)    # This is a DELTA value

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0003, incremental=False)
    assert np.allclose(P_t, 20.70580470925962)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0004, incremental=False)
    assert np.allclose(P_t, 20.92746192779109)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0005, incremental=False)
    assert np.allclose(P_t, 21.14272449451510)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0007, incremental=False)
    assert np.allclose(P_t, 21.55497320836003)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.001, incremental=False)
    assert np.allclose(P_t, 22.13079644707523)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0015, incremental=False)
    assert np.allclose(P_t, 22.98939523899655)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.002, incremental=False)
    assert np.allclose(P_t, 23.73870896747882)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.003, incremental=False)
    assert np.allclose(P_t, 24.97201963453988)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.005, incremental=False)
    assert np.allclose(P_t, 26.68110034281166)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.008, incremental=False)
    assert np.allclose(P_t, 28.12354983857080)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.01, incremental=False)
    assert np.allclose(P_t, 28.66884983321557)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=1., incremental=False)
    assert np.allclose(P_t, 29.70512259124693)


    # Verify reaching equilibrium at a large t
    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=10., incremental=False)
    eq = ReactionKinetics._compute_equilibrium_conc_first_order(kF=5., kR=2., a=1, A0=10., b=1, B0=50., p=1, P0=20.)
    assert np.allclose(P_t, eq["P"])    # 29.705122591242464


    # Special case A0 = B0
    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=35, B0=35, P0=20., t=0.002, incremental=False)
    assert np.allclose(P_t, 28.99659131142920)



def test_exact_advance_synthesis_reversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law,
    # at the middle point of 3

    # Case where `A` is the limiting reagent
    A0=10
    B0=50
    P0=20
    kF=5
    kR=2

    #print(ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0))
    #Equilibrium at {'A': 0.2948774087575341, 'B': 40.294877408757536, 'P': 29.705122591242464}

    h = 0.00001
    t_start = 0.002
    times = t_start + h * np.arange(3)

    # Sample at 3 closely-spaced points
    p = [AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
         for t in times
         ]
    #print(p)    # [23.738708979042038, 23.752698461282, 23.766650993990154]

    p_middle = p[1]
    a_middle = A0 - (p_middle - P0)    # From mass conservation
    b_middle = B0 - (p_middle - P0)    # From mass conservation

    #print(a_middle, b_middle)        # 6.247301538717998 46.247301538718
    gradient = np.gradient(p, h)
    derivative_middle = gradient[1]
    rate_at_mid_point = kF * a_middle * b_middle - kR * p_middle    # kF [A] [B] - kR [P]
    assert np.allclose(derivative_middle, rate_at_mid_point)


    # Case where `B` is the limiting reagent
    A0=60
    B0=30
    P0=10
    kF=7
    kR=3

    #print(ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0))
    #Equilibrium at {'A': 30.553318635734474, 'B': 0.5533186357344739, 'P': 39.44668136426553}

    h = 0.00001
    t_start = 0.002
    times = t_start + h * np.arange(3)

    # Sample at 3 closely-spaced points
    p = [AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
         for t in times
         ]
    #print(p)    # [25.243842573299553, 25.289220064086344, 25.334407842935416]
    p_middle = p[1]
    a_middle = A0 - (p_middle - P0)    # From mass conservation
    b_middle = B0 - (p_middle - P0)    # From mass conservation

    #print(a_middle, b_middle)   44.710779935913656 14.710779935913656
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

    #print(ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0))
    #Equilibrium at {'A': 7.038367176906169, 'B': 7.038367176906169, 'P': 30.96163282309383}

    h = 0.00005
    t_start = 0.01
    times = t_start + h * np.arange(3)

    # Sample at 3 closely-spaced points
    p = [AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
         for t in times
         ]
    #print(p)    # [30.536666653790505, 30.538373794434957, 30.54007389774443]
    p_middle = p[1]
    a_middle = A0 - (p_middle - P0)    # From mass conservation
    b_middle = B0 - (p_middle - P0)    # From mass conservation

    #print(a_middle, b_middle, p)    # 7.461626205565043 7.461626205565043
    gradient = np.gradient(p, h)
    derivative_middle = gradient[1]
    rate_at_mid_point = kF * a_middle * b_middle - kR * p_middle  # kF [A] [B] - kR [P]
    assert np.allclose(derivative_middle, rate_at_mid_point)



def test_exact_advance_synthesis_reversible_3():
    # Compare the exact solution against a fine-grained forward Euler approximation
    # Reaction A + B <-> P
    sr = SpeciesRegistry()
    rxn_defn = ReactionDefinition(reactants=["A", "B"], products="P",
                                  species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF": 5., "kR": 2.})
    rxn = rxn_defn.sim_reactions[0]

    t_final = 0.001
    n_steps = 1800
    t_step = t_final / n_steps
    a = 10
    b = 50
    p = 20

    #print(ReactionKinetics._compute_equilibrium_conc_first_order(kF=5., kR=2., a=1, A0=a, b=1, B0=b, p=1, P0=p))
    #Equilibrium at {'A': 0.2948774087575341, 'B': 40.294877408757536, 'P': 29.705122591242464}

    for i in range(n_steps):
        increment_dict, _ = \
            rxn.step_simulation(delta_time=t_step, conc_dict={"A": a, "B": b, "P": p})
        delta_p = increment_dict["P"]
        p += delta_p
        a -= delta_p
        b -= delta_p

    #print(p)    # 22.130930189693114   Value at t_final, from the fine-grained forward Euler approximation
    exact_p = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=t_final, incremental=False)   # 22.130796453845
    assert np.allclose(p, exact_p)


    # Let's do a second identical round, to reach t_final * 2
    for i in range(n_steps):
        increment_dict, _ = \
            rxn.step_simulation(delta_time=t_step, conc_dict={"A": a, "B": b, "P": p})
        delta_p = increment_dict["P"]
        p += delta_p
        a -= delta_p
        b -= delta_p

    #print(p)    # 23.738906198053474  Value at t_final * 2, from the fine-grained forward Euler approximation
    exact_p = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=t_final * 2, incremental=False)
    assert np.allclose(p, exact_p)



def test_exact_advance_synthesis_irreversible():
    # Reaction A + B -> P

    # We'll start with `A` as the limiting reagent
    # No change at time 0
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0, incremental=False)
    assert np.allclose(P_t, 20.)

    Delta_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0, incremental=True)
    assert np.allclose(Delta_t, 0)

    # Compare against the REVERSIBLE reaction solver with a near-zero reverse rate constant : not exactly zero, because its
    # implementation reverts to the irreversible solver when the reverse rate constant is very close to zero
    eps = 0.00001
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0.01, incremental=False)
    c = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=10., B0=50., P0=20., t=0.01, incremental=False)
    assert np.allclose(P_t, c)  # 28.8871974442945

    Delta_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0.01, incremental=True)
    D_c = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=10., B0=50., P0=20., t=0.01, incremental=True)
    assert np.allclose(Delta_t, D_c)

    for t  in [0.0003, 0.0004, 0.5, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.008, 0.01, 1.]:
        P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=t, incremental=False)
        c = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=10., B0=50., P0=20., t=t, incremental=False)
        assert np.allclose(P_t, c)


    # A time so large that the reaction has gone to completion (`A` being the limiting reagent)
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=3., incremental=False)
    assert np.allclose(P_t, 30.)    # All `A`(the limiting reagent) converted to `P`

    # A value of t so ridiculously large that an OverflowError is caused (but caught) in the internal math.exp usage
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=200., incremental=False)
    assert np.allclose(P_t, 30.)    # All `A`(the limiting reagent) converted to `P`


    # Now, let `B` be the limiting reagent at large times
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=5., P0=20., t=3., incremental=False)
    assert np.allclose(P_t, 25.)    # All `B`(the limiting reagent) converted to `P`

    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=5., P0=20., t=200., incremental=False)
    assert np.allclose(P_t, 25.)    # All `B`(the limiting reagent) converted to `P`


    # If either A0 or B0 is zero, the reaction doesn't proceed
    Delta_P = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=0, B0=50., P0=20., t=1., incremental=True)
    assert np.allclose(Delta_P, 0)

    Delta_P = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=0, P0=20., t=1., incremental=True)
    assert np.allclose(Delta_P, 0)


    # Case A0 = B0  (equivalently, 2 A -> P)

    # No change at time 0
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30, P0=20, t=0, incremental=False)
    assert np.allclose(P_t, 20)     # No change

    for t  in [0.0003, 0.0004, 0.5, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.008, 0.01, 1.]:
        P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30, P0=20, t=t, incremental=False)
        P_t_approx = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30.000001, P0=20, t=t, incremental=False)
        c = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=30, B0=30, P0=20, t=t, incremental=False)
        assert np.allclose(P_t, c)
        assert np.allclose(P_t, P_t_approx)

    # A time so large that the reaction has gone to completion
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30, P0=20, t=1000., incremental=False)
    assert np.allclose(P_t, 50.)    # The reagents fully converted to `P`


def test_exact_advance_synthesis_irreversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law,
    # at the middle point of 3
    h = 0.00001
    t_start = 0.003
    p1 = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50, P0=20., t=t_start)
    p2 = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50, P0=20., t=t_start + h)
    a2 = 10 - (p2 - 20)    # From mass conservation
    b2 = 50 - (p2 - 20)    # From mass conservation
    p3 = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50, P0=20., t=t_start + 2 * h)
    #print(a2, b2, p2)
    gradient = np.gradient([p1, p2, p3], h)
    derivative_p2 = gradient[1]
    rate_at_mid_point = 5 * a2 * b2   # kF * a2 * b2
    assert np.allclose(derivative_p2, rate_at_mid_point)


def test_exact_advance_synthesis_irreversible_3():
    # Compare the exact solution against a fine-grained forward Euler approximation
    # Reaction A + B -> P
    sr = SpeciesRegistry()
    rxn_defn = ReactionDefinition(reactants=["A", "B"], products="P",
                                  species_registry= sr, autoregister_species=True,
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF": 5., "kR": 0})
    rxn = rxn_defn.sim_reactions[0]

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

    exact_p = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0.01, incremental=False)   # 28.8871974442945
    assert np.allclose(p, exact_p)



def test_approx_solution_synthesis_rxn():
    # Reaction A + B <-> P
    kF=5.
    kR=2.
    A0=10.
    B0=50.
    P0=20.

    with pytest.raises(Exception):
        AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=0, kR=kR, A0=A0, B0=B0, P0=P0, t=0)   # kF cannot be zero

    p = AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=0)
    assert np.allclose(p, P0)      # No change from P0

    #equil = ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0)
    #print(equil)       # {'A': 0.2948774087575341, 'B': 40.294877408757536, 'P': 29.705122591242464}

    p = AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=0.1)
    assert np.allclose(p, 29.705122591242464)       # Reaching the equilibrium value after a sufficiently long time


    # Process a whole array of times in one call
    t = np.array([0, 0.000864, 0.001555, 0.009850, 0.067400])
    result = AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=5., kR=2., A0=10., B0=50., P0=20., t=t)
    assert np.allclose(result, [20., 22.11257633, 23.46603241, 29.11416425, 29.70512254])

    # Compare against exact solution
    for t in [0, 0.000864, 0.001555, 0.009850, 0.067400]:
        p_approx = AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
        p_exact = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
        assert abs(p_approx - p_exact) / p_exact < 0.02     # Less than 2% discrepancy



def test_approx_solution_synthesis_rxn_ALT():
    pass    # TODO
