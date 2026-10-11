import math
import cmath
import numpy as np



class AnalyticReactionSolver:
    """
    For reactions that have known analytical solutions (exact or approximate)
    """
    @staticmethod
    def exact_advance_unimolecular_reversible(kF, kR, A0, P0, t, incremental=False) -> float:
        """
        Exactly advance the concentrations
        in the reversible elementary Reaction A <-> P,
        from time 0 to time t,
        with the specified parameters.

        :param kF:  Forward reaction rate constant
        :param kR:  Reverse reaction rate constant
        :param A0:  Initial concentration of the reactant A
        :param P0:  Initial concentration of the product P
        :param t:   The end time of the reaction that started at time zero
        :param incremental: [OPTIONAL] If True, the changes in concentrations is returned,
                                rather than the final one.  Default: False

        :return:    The concentrations of P at time t (if `incremental` is False)
                        or its concentration change during the time interval (if `incremental` is True)
        """
        TOT = A0 + P0

        # Formula is:  A(t) = (A0 - [(kR TOT) / (kF + kR)]) Exp[-(kF + kR) t] + kR TOT / (kF + kR)
        # In python:          (A0 - (kR * TOT) / (kF + kR)) * math.exp(-(kF + kR) * t) + kR * TOT / (kF + kR)

        sum_rates = kF + kR
        ratio = (kR * TOT) / sum_rates

        A_t = (A0 - ratio) * math.exp(-sum_rates * t) + ratio
        P_t = TOT - A_t     # From mass conservation

        if incremental:
            return P_t - P0

        return P_t


    @staticmethod
    def exact_advance_unimolecular_irreversible(kF, A0, P0, t, incremental=False) -> float:
        """
        Exactly advance the concentrations
        in the ir-reversible 1st Order Reaction A -> P,
        from time 0 to time t,
        with the specified parameters.

        :param kF:  Forward reaction rate constant (the reverse one is taken to be zero)
        :param A0:  Initial concentration of the reactant A
        :param P0:  Initial concentration of the product P
        :param t:   The end time of the reaction that started at time zero
        :param incremental: [OPTIONAL] If True, the change in concentration is returned,
                                rather than the final ones.  Default: False

        :return:    The concentrations of P at time t (if `incremental` is False)
                        or its concentration change during the time interval (if `incremental` is True)
        """

        # Formula is:  A(t) = A0 Exp(-kF t)

        if incremental:
            A_incr = A0 * (np.exp(-kF * t) - 1)
            P_incr = - A_incr   # From mass conservation
            return P_incr


        TOT = A0 + P0
        A_t = A0 * np.exp(-kF * t)
        P_t = TOT - A_t     # From mass conservation

        return P_t



    @staticmethod
    def exact_advance_synthesis_reversible(kF, kR, A0, B0, P0, t, incremental=False) -> float:
        """
        Exactly advance the concentrations
        in the reversible elementary Reaction A + B <-> P,
        from time 0 to time t,
        with the specified parameters.

        :param kF:  Forward reaction rate constant
        :param kR:  Reverse reaction rate constant
        :param A0:  Initial concentration of the 1st reactant A
        :param B0:  Initial concentration of the 2nd reactant B
        :param P0:  Initial concentration of the product P
        :param t:   The end time of the reaction that started at time zero
        :param incremental: [OPTIONAL] If True, the changes in concentration of P is returned,
                                rather than the final ones.  Default: False

        :return:   The concentration of P at time t (if `incremental` is False)
                        or its concentration change during the time interval (if `incremental` is True)
        """
        if np.allclose(kR, 0):
            '''
            Note - the case when KR=0 and A0=B0 cannot be passed to the general solver below,
                   because in such as case:
                        T = A0 + P0     # A collapse of the general cases of AP_tot and BP_tot
                        alpha = kF
                        beta = 2 * kF * T
                        gamma = kF * T**2 
            In this scenario:
                - beta**2 + 4 * alpha * gamma = - [2 kF T]**2 + 4 kF * (kF * T**2) 
                                              = - 4 kF**2 T**2 + 4 kF**2 T**2 
                                              = 0  
            and the general solution cannot be used, because in it we divide by the square root of that term!       
            '''
            #print("****** Switching to IR-reversible version")
            return  AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=kF, A0=A0, B0=B0, P0=P0, t=t, incremental=incremental)


        AP_tot = A0 + P0        # Quantity conserved thru the rxn, from the stoichiometry
        BP_tot = B0 + P0        # Quantity conserved thru the rxn, from the stoichiometry

        alpha = kF
        beta = kF * (AP_tot + BP_tot) + kR
        gamma = kF * AP_tot * BP_tot
        '''
        The ODE to solve for P(t) is:
            P'(t) = kF A(t) B(t) - kR P(t)      # Forward rate minus reverse rate
            
        The stoichiometry dictates that the drop in each reagent equals the increase in the product:
            A(t) = A0 - (P(t) - P0)
            B(t) = B0 - (P(t) - P0)
            
        which can be re-arranged to:
            P'(t) = kF [A0 - P(t) + P0] [B0 - P(t) + P0] - kR P(t)
            P'(t) = kF [AC_tot - P(t)] [BC_tot - P(t)] - kR P(t)
            P'(t) = kF [AC_tot * BC_tot - AC_tot * P(t) - P(t) * BC_tot + P(t)**2] - kR P(t)
            P'(t) = kF * P(t)**2 
                    - kF * AC_tot * P(t) - kF * BC_tot * P(t) - kR P(t) 
                    + kF * AC_tot * BC_tot 
            P'(t) = kF * P(t)**2 
                    - (kF * AC_tot + kF * BC_tot + kR) * P(t) 
                    +  kF * AC_tot * BC_tot 
            P'(t) = alpha [P(t)]**2 - beta P(t) + gamma
        
        The following is based on the answer by Mathematica v.7 to:
            DSolve[{p'[t] == alpha (p[t])**2 - beta p[t] + gamma, p[0] == P0}, p[t], t]
            Simplify[%]
        '''
        term = - beta**2 + 4 * alpha * gamma
        #print(f"term = {term}")

        if np.allclose(term, 0):
            # The term = 0 must be intercepted because we later divide by the square root of that term
            # TODO: can this scenario ever occur when kR > 0  ?
            raise Exception("exact_advance_synthesis_reversible(): No provision for this case yet, when beta**2 = 4 * alpha * gamma")
        else:
            sqrt_term = cmath.sqrt(term)            # Note: `term` could be negative; hence, the complex square root

            arctan_arg = (beta - 2 * alpha * P0) / sqrt_term    # This will be the argument to pass to the arctan() function, below
            arctan_value = np.arctan(arctan_arg)
            #print(f"arctan_value = {arctan_value}")

            # Note that every variable so far is a constant, NOT dependent on time

            tan_arg = (0.5 * sqrt_term * t) - arctan_value     # This will be the argument to pass to the tan() function, below

            result = (beta + sqrt_term * np.tan(tan_arg) ) / (2*alpha)

            P_t = result.real             # Drop the imaginary components (which ought to be very close to zero)


        if incremental:
            return P_t - P0

        return P_t



    @staticmethod
    def exact_advance_synthesis_irreversible(kF, A0, B0, P0, t, incremental=False) -> float:
        """
        Exactly advance the concentrations
        in the irreversible elementary Reaction A + B -> P,
        from time 0 to time t,
        with the specified parameters.

        :param kF:  Forward reaction rate constant
        :param A0:  Initial concentration of the 1st reactant A
        :param B0:  Initial concentration of the 2nd reactant B
        :param P0:  Initial concentration of the product P
        :param t:   The end time of the reaction that started at time zero
        :param incremental: [OPTIONAL] If True, the changes in concentrations are returned,
                                rather than the final ones.  Default: False
        :return:    The concentration of P at time t (if `incremental` is False)
                        or its concentration change during the time interval (if `incremental` is True)
        """
        # Adapted from "Athel Cornish-Bowden, Fundamentals of Enzyme Kinetics, 4th edn", section 1.2.3
        # Note: the above author assumes P0 = 0, and overlooks the case A0 = B0 (which requires special treatment)
        '''
        Start with the ODE:
            d P(t) / dt = kF * A(t) * B(t) = kF * (A0 - (P(t)-P0)) * (B0 - (P(t)-P0))
                     
        Special case - when A0 = B0, the ODE becomes:
            d P(t) / dt = kF * A(t) * A(t) = kF * ( A0 - (P(t)-P0) )**2
            which can, for example, be plugged into Mathematica as:
            DSolve[{ p'[t] == kF (A0 + P0 - p[t])^2 , p[0] == P0 }, p[t], t]
        '''
        alpha =  P0 + A0

        if np.allclose(A0, B0):
            # When A0 = B0, the ODE has a special solution:
            #   P(t) = (P0 + alpha A0 kF t) / (1 + A0 kF t)
            prod = A0 * kF * t
            P_t = (P0 + alpha * prod) / (1 + prod)  # Concentration of `P` at time t

        else:
            # When A0 != B0, the solution P(t) of the ODE must satisfy:
            #   (A0 * (P0+B0 - P(t))) / (B0 * (P0+A0 - P(t))) = e ^ (B0-A0)*kF*t
            #   [see "Athel Cornish-Bowden, Fundamentals of Enzyme Kinetics, 4th edn", section 1.2.3]
            #   which gets solved for P(t) below
            #   Note: when A0 = B0, the above equation doesn't hold (at, at any rate, turns into a trivial 1 = 1)
            beta =   P0 + B0
            try:
                # Evaluate the exponential on the right-hand size
                eps = math.exp((B0-A0)* kF * t)
            except OverflowError:
                # Overflow occurs if B0 > A0 and t is very large;
                # in such a case, `A` is the limiting reagent, and fully converts to `P`
                if incremental:
                    return A0
                return alpha    # P0 + A0

            ''' 
            With the above variable assignments, and writing P for P(t), the earlier equation can be re-stated as:
                (A0 * (beta - P)) / (B0 * (alpha - P)) = eps 
                           
            i.e.:
                (A0 * (beta - P)) = eps * (B0 * (alpha - P))
                A0 * beta - A0 * P = eps * (B0 * alpha - B0 * P)
                A0 * beta - A0 * P = B0 * alpha * eps - B0 * P * eps
                B0 * P * eps - A0 * P = B0 * alpha * eps - A0 * beta
                (B0 * eps - A0) * P   = B0 * alpha * eps - A0 * beta
            '''
            #print(f"alpha: {alpha} , beta: {beta} , eps: {eps}")
            P_t = (B0 * alpha * eps - A0 * beta) / (B0 * eps - A0)    # Concentration of `P` at time t


        if incremental:
            return P_t-P0

        return P_t



    @staticmethod
    def approx_solution_synthesis_rxn(kF, kR, A0, B0, P0, t :float|np.ndarray) -> float|np.ndarray:
        """
        Return the APPROXIMATE analytical solution, by way of exponentials,
        of the reversible 2nd Order Reaction A + B <-> P,
        with the specified parameters, sampled at the given time(s).

        This approximation is relatively coarse; in the pytests, errors of up to about 2% were seen,
        with the worst errors in the mid-portion between the starting state and the final equilibrium

        :param kF:  Forward reaction rate constant (cannot be zero)
        :param kR:  Reverse reaction rate constant
        :param A0:  Initial concentration of the 1st reactant A
        :param B0:  Initial concentration of the 2nd reactant B
        :param P0:  Initial concentration of the product P
        :param t:   The end time of the reaction that started at time zero
                        OR a Numpy array with the desired times at which the solutions are desired

        :return:    A concentration value, or Numpy arrays with the concentrations,
                        of P at the time - or times - given by the argument `t`
        """
        assert not np.allclose(kF, 0), \
            "approx_solution_synthesis_rxn(): this approximation cannot be used when kF is zero"

        # Calculate the equilibrium concentrations
        K_inv = kR / kF     # Inverse of the equilibrium constant.  Note that we're assuming kF isn't zero

        # The following derive from solving the quadratic equation  (A0 - m) * (B0 - m) / (C0 + m) == K_inv , for m
        # where m is the Product concentration change ("moles/liter of forward reaction")
        TOT_reactants = A0+B0
        r = (K_inv**2 + TOT_reactants**2) / 4 + (K_inv * TOT_reactants / 2) + P0 * K_inv - A0 * B0
        m = (K_inv + TOT_reactants)/2  - math.sqrt(r)         # Product concentration change

        # The reactants get consumed by m, while the product increases by m
        A_eq = A0 - m
        B_eq = B0 - m

        #print(f"\nProduct concentration change: {m} | A_eq: {A_eq}  |  B_eq: {B_eq}")

        l = (kF + kR) * (A_eq + B_eq)  # Relaxation rate constant λ

        AC_tot = A0 + P0        # Quantity conserved thru the rxn, from the stoichiometry

        A_arr = A_eq + (A0 - A_eq) * np.exp(-l * t)     # Approximate analytical solution
        P_arr = AC_tot - A_arr

        return P_arr



    @staticmethod
    def approx_solution_synthesis_rxn_ALT(kF, kR, A0, B0, P0, t, incremental) -> float:
        """
        Provided by ChatGPT.  Not fully tested.

        Closed-form solution for:
        dp/dt = kf (a0 - p + p0)(b0 - p + p0) - kr p
        with p(0) = p0

        Parameters
        ----------
        t : float or array_like
            Time(s) at which to evaluate p(t)
        P0, kF, kR, A0, B0 : float
            Reaction parameters

        Returns
        -------
        p : float or ndarray
            Value(s) of p(t)
        """

        # Coefficients of Riccati equation
        alpha = kF
        beta  = kF * (A0 + B0 + 2 * P0) + kR
        gamma = kF * (A0 + P0) * (B0 + P0)

        # Discriminant
        discr = beta**2 - 4*alpha*gamma
        if np.any(discr < 0):
            raise ValueError("Negative discriminant: parameters give complex roots.")

        sqrt_disc = np.sqrt(discr)

        # Equilibria
        p_plus  = (beta + sqrt_disc) / (2*alpha)
        p_minus = (beta - sqrt_disc) / (2*alpha)

        # Relaxation rate
        lam = alpha * (p_plus - p_minus)

        #t = np.asarray(t, dtype=float)

        exp_term = np.exp(-lam * t)

        numerator = (
                p_plus * (P0 - p_minus)
                - p_minus * (P0 - p_plus) * exp_term
        )

        denominator = (
                (P0 - p_minus)
                - (P0 - p_plus) * exp_term
        )

        P_t = numerator / denominator

        if incremental:
            return P_t - P0

        return P_t
