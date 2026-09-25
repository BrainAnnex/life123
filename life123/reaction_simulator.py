# Classes ReactionSimulator, AnalyticalReactionSolver, and VariableTimeSteps:

import math
import cmath
import numpy as np




class ReactionSimulator:
    """

    """
    pass
    





####################################################################################################

class AnalyticalReactionSolver:
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
            return  AnalyticalReactionSolver.exact_advance_synthesis_irreversible(kF=kF, A0=A0, B0=B0, P0=P0, t=t, incremental=incremental)


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






####################################################################################################

class VariableTimeSteps:
    """
    Methods for managing variable time steps during reactions
    """


    def __init__(self, uc=None):
        # ***  PARAMETERS FOR AUTOMATED ADAPTIVE TIME STEP SIZES  ***
        # Note: The "aborts" below are "elective" aborts - i.e. not aborts from hard errors (further below)
        #       The default values get packed into a "preset", specified by a name, and optionally passed when
        #       instantiating this object.
        #       Users in need to change these values post-instantiation will generally use use_adaptive_preset(), below;
        #       or, for more control, set_thresholds() and set_step_factors()
        self.thresholds = []    # A list of "rules"
                                # EXAMPLE:                  [
                                #                             {"norm": "norm_A", "low": 0.5, "high": 0.8, "abort": 1.44},
                                #                             {"norm": "norm_B", "low": 0.08, "high": 0.5, "abort": 1.5}
                                #                           ]
                                #       low:    The value below which the norm is considered to be good (and the variable steps are being made larger);
                                #                if None, it doesn't get used
                                #       high:   The value above which the norm is considered to be excessive (and the variable steps get reduced);
                                #               if None, it doesn't get used
                                #       abort:  The value above which the norm is considered to be dangerously large (and the last variable step gets
                                #               discarded and back-tracked with a smaller size)
                                #               if None, it doesn't get used

        self.step_factors = {}
                            # EXAMPLE: {"upshift": 1.2, "downshift": 0.5, "abort": 0.4, "error": 0.2}
                            # "upshift" must be > 1 ; all the other values must be < 1
                            # Generally, error <= abort <= downshift
                            # "Error" value: Factor by which to multiply the time step
                            #   in case of negative-concentration error from excessive step size
                            #   NOTE: this is from ERROR aborts,
                            #   not to be confused with general aborts based on reaching high threshold


        # Zero the number of times that each norm got involved in step-size decision
        self.norm_usage = {}
        self.reset_norm_usage_stats()

        self.uc = uc    # Object of type "UniformCompartment"



    def set_thresholds(self, norm :str, low=None, high=None, abort=None) -> None:
        """
        Create or update a rule based on the given norm
        (simulation parameters that affect the adaptive variable step sizes.)
        If no rule with the specified norm was previously set, it will be added.

        All None values in the arguments are ignored.

        :param norm:    The name of the "norm" (criterion) being used - either a previously-set one or not
        :param low:     The value below which the norm is considered to be good (and the variable steps are being made larger);
                            if None, it doesn't get changed
        :param high:    The value above which the norm is considered to be excessive (and the variable steps get reduced);
                            if None, it doesn't get changed
        :param abort:   The value above which the norm is considered to be dangerously large (and the last variable step gets
                            discarded and back-tracked with a smaller size) ;
                            if None, it doesn't get changed
        :return:        None
        """
        assert type(norm) == str and norm != "", \
            "set_thresholds(): the `norm` argument must be a non-empty string"

        if (abort is not None) and (high is not None):
            assert abort > high, \
                f"set_thresholds(): `abort` value ({abort}) must be > `high' value ({high})"

        if (high is not None) and (low is not None):
            assert high > low, \
                f"set_thresholds(): `high` value ({high}) must be > `low' value ({low})"

        if (abort is not None) and (low is not None):
            assert abort > low, \
                f"add_thresholds(): `abort` value ({high}) must be > `low' value ({low})"

        for i, t in enumerate(self.thresholds):
            if t.get("norm") == norm:
                # Found a rule using the requested norm

                t_original = t.copy()  # Create a backup copy in case of error

                if low is not None:
                    t["low"] = low

                if high is not None:
                    t["high"] = high

                if abort is not None:
                    t["abort"] = abort

                if len(t) == 1:
                    del self.thresholds[i]      # Completely eliminate this un-used norm

                low, high, abort = t.get("low"), t.get("high"), t.get("abort")

                if (abort is not None) and (high is not None):
                    if abort <= high:
                        self.thresholds[i] = t_original     # Restore original values
                        raise Exception(f"set_thresholds(): `abort` value ({abort}) must be > `high' value ({high})")

                if (high is not None) and (low is not None):
                    if high <= low:
                        self.thresholds[i] = t_original     # Restore original values
                        raise Exception(f"set_thresholds(): `high` value ({high}) must be > `low' value ({low})")

                if (abort is not None) and (low is not None):
                    if abort <= low:
                        self.thresholds[i] = t_original     # Restore original values
                        raise Exception(f"add_thresholds(): `abort` value ({high}) must be > `low' value ({low})")

                return

        # If we get here, it means that we're handling a norm
        # not currently present in the list self.thresholds
        new_t = {"norm": norm}
        if low is not None:
            new_t["low"] = low
        if high is not None:
            new_t["high"] = high
        if abort is not None:
            new_t["abort"] = abort

        if self.thresholds is None:
            self.thresholds = [new_t]
        else:
            self.thresholds.append(new_t)



    def delete_thresholds(self, norm :str, low=False, high=False, abort=False) -> None:
        """
        Delete one or more of the threshold values associated to a rule using the specified norm.
        If none of the threshold values remains in place, the whole rule gets eliminated altogether.
        Attempting to delete something not present, will raise an Exception

        :param norm:    The name of the "norm" (criterion) being used - for a previously-set rule
        :param low:     If True, this threshold value will get deleted
        :param high:    If True, this threshold value will get deleted
        :param abort:   If True, this threshold value will get deleted
        :return:        None
        """
        for i, t in enumerate(self.thresholds):
            if t.get("norm") == norm:
                # Found a rule using the requested norm

                if low:
                    del t["low"]

                if high:
                    del t["high"]

                if abort:
                    del t["abort"]

                if len(t) == 1:
                    del self.thresholds[i]      # Completely eliminate this un-used norm
                return

        # If we get here, it means that we're handling a norm
        # not currently present in the list self.thresholds
        raise Exception(f"delete_thresholds(): no norm named '{norm}' was found")



    def display_value_against_thresholds(self, all_norms) -> None:
        """

        :param all_norms:
        :return:            None
        """
        for rule in self.thresholds:
            value = all_norms.get(rule['norm'])
            print(self.display_value_against_thresholds_single_rule(rule, value))


    def display_value_against_thresholds_single_rule(self, rule :dict, value) -> str:
        """
        Examine how the specified value fits
        relatively to the 'low', 'high', and 'abort' stored in the given rule

        :param rule:    A dict that must contain the key 'norm',
                            and may contain the keys: 'low', 'high', and 'abort'
                            (referring to 3 increasingly-high thresholds)
        :param value:   Either None, or a number to compare to the thresholds
        :return:        A string that visually highlights the relative position of the value
                            relatively to the given thresholds
        """
        s = f"                   {rule['norm']} : "     # Name of the norm being used

        if value is None:
            return s + " (skipped; not needed)"

        # Extract the 3 thresholds (some might be missing)
        low = rule.get('low')
        high = rule.get('high')
        abort = rule.get('abort')

        if low is not None and value <= low:
            # The value is below the `low` threshold
            return f"{s}(VALUE {value:.5g}) | low {low} | high {high} | abort {abort}"

        if high is not None and value < high:
            # The value is between the `low` and `high` thresholds
            return f"{s}low {low} | (VALUE {value:.5g}) | high {high} | abort {abort}"

        if abort is not None and value < abort:
            # The value is above the `high` threshold
            return f"{s}low {low} | high {high} | (VALUE {value:.5g}) | abort {abort}"

        # If we get thus far, the value is above the `abort` threshold
        return f"{s}low {low} | high {high} | abort {abort} | (VALUE {value:.5g})"



    def set_step_factors(self, upshift=None, downshift=None, abort=None, error=None) -> None:
        """
        Over-ride current values for simulation parameters that affect the adaptive variable step sizes.
        Values not explicitly passed will remain the same.

        :param upshift:     Fraction by which to increase the variable step size
        :param downshift:   Fraction by which to decrease the variable step size
        :param abort:       Fraction by which to decrease the variable step size in cases of step re-do
        :param error:       Fraction by which to decrease the variable step size in cases of error
        :return:            None
        """
        if upshift is not None:
            assert upshift > 1, "set_step_factors(): `upshift` value must be > 1"
            self.step_factors["upshift"] = upshift

        if downshift is not None:
            assert 0 < downshift < 1, "set_step_factors(): `downshift` value must be a positive number < 1"
            self.step_factors["downshift"] = downshift

        if abort is not None:
            assert 0 < abort < 1, "set_step_factors(): `abort` value must be a positive number < 1"
            self.step_factors["abort"] = abort

        if error is not None:
            assert 0 < error < 1, "set_step_factors(): `error` value must be a positive number < 1"
            self.step_factors["error"] = error



    def show_adaptive_parameters(self) -> None:
        """
        Print out the current values for the adaptive time-step parameters

        :return:    None
        """
        print("Parameters used for the automated adaptive time step sizes -")
        print("    THRESHOLDS: ", self.thresholds)
        print("    STEP FACTORS: ", self.step_factors)



    def use_adaptive_preset(self, preset :str) -> None:
        """
        Lets the user choose a preset to use from now on, unless later explicitly changed,
        for use in ALL reaction simulations involving adaptive time steps.
        The preset will affect the degree to which the simulation will be "risk-taker" vs. "risk-averse" about
        taking larger steps.

        Note:   for more control, use set_thresholds() and set_step_factors()

                For example, using the "mid" preset is the same as issuing:
                    ReactionSimulator.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
                    ReactionSimulator.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
                    ReactionSimulator.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4, error=0.25)

        :param preset:  String with one of the available preset names;
                            allowed values are (in generally-increasing speed):
                            'heavy_brakes', 'slower', 'slow', 'mid', 'fast'
        :return:        None
        """
        if preset == "heavy_brakes":   # It slams on the "brakes" hard in case of abort or errors
            self.thresholds = [{"norm": "norm_A", "low": 0.02, "high": 0.025, "abort": 0.03},
                               {"norm": "norm_B", "low": 0.05, "high": 1.0, "abort": 2.0}]
            self.step_factors = {"upshift": 1.6, "downshift": 0.15, "abort": 0.08, "error": 0.05}

        elif preset == "small_rel_change":
            self.thresholds = [{"norm": "norm_A", "low": 2., "high": 5., "abort": 10.},
                               {"norm": "norm_B", "low": 0.008, "high": 0.5, "abort": 2.0}]     # The "low" value of "norm_B" is very strict
            self.step_factors = {"upshift": 1.5, "downshift": 0.25, "abort": 0.25, "error": 0.2}

        elif preset == "slower":   # Very conservative about taking larger steps
            self.thresholds = [{"norm": "norm_A", "low": 0.2, "high": 0.5, "abort": 0.8},
                               {"norm": "norm_B", "low": 0.03, "high": 0.05, "abort": 0.5}]
            self.step_factors = {"upshift": 1.01, "downshift": 0.5, "abort": 0.1, "error": 0.1}

        elif preset == "slow":   # Conservative about taking larger steps
            self.thresholds = [{"norm": "norm_A", "low": 0.2, "high": 0.5, "abort": 0.8},
                               {"norm": "norm_B", "low": 0.05, "high": 0.4, "abort": 1.3}]
            self.step_factors = {"upshift": 1.1, "downshift": 0.3, "abort": 0.2, "error": 0.1}

        elif preset == "mid":     # A "middle-of-the road" heuristic: somewhat "conservative" but not overly so
            self.thresholds = [{"norm": "norm_A", "low": 0.5, "high": 0.8, "abort": 1.44},
                               {"norm": "norm_B", "low": 0.08, "high": 0.5, "abort": 1.5}]
            self.step_factors = {"upshift": 1.2, "downshift": 0.5, "abort": 0.4, "error": 0.25}

        elif preset == "fast":   # Less conservative (more "risk-taker") about taking larger steps
            self.thresholds = [{"norm": "norm_A", "low": 0.8, "high": 1.2, "abort": 1.7},
                               {"norm": "norm_B", "low": 0.15, "high": 0.8, "abort": 1.8}]
            self.step_factors = {"upshift": 1.5, "downshift": 0.8, "abort": 0.6, "error": 0.5}

        elif preset == "mid_inclusive": # A "middle-of-the road" heuristic that makes use of more norms
            self.thresholds = [{'norm': 'norm_A', 'low': 0.2, 'high': 0.8, 'abort': 1.44},
                               {'norm': 'norm_B', 'low': 0.08, 'high': 0.5, 'abort': 1.5},
                               {'norm': 'norm_C', 'low': 0.5, 'high': 1.2, 'abort': 1.6},
                               {'norm': 'norm_D', 'low': 1.3, 'high': 1.7, 'abort': 1.8}]
            self.step_factors = {'upshift': 1.1, 'downshift': 0.5, 'abort': 0.4, 'error': 0.25}

        elif preset == "mid_inclusive_slow":
            self.thresholds = [{'norm': 'norm_A', 'low': 0.15, 'high': 0.8, 'abort': 1.44},
                               {'norm': 'norm_B', 'low': 0.05, 'high': 0.5, 'abort': 1.5},
                               {'norm': 'norm_C', 'low': 0.5, 'high': 1.2, 'abort': 1.6},
                               {'norm': 'norm_D', 'low': 1.1, 'high': 1.7, 'abort': 1.8}]
            self.step_factors = {'upshift': 1.1, 'downshift': 0.5, 'abort': 0.4, 'error': 0.25}

        else:
            raise Exception(f"set_adaptive_parameters(): unknown value for the `preset` argument ({preset}); "
                            f"allowed values are 'heavy_brakes', 'slower', 'slow', 'mid', 'fast', 'mid_inclusive'")



    def adjust_timestep(self, n_chems: int, indexes_of_active_chemicals :[int],
                        delta_conc: np.array, baseline_conc=None, prev_conc=None,
                        ) -> dict:
        """
        Computes some measures of the change of concentrations, from the last step, in the context of the
        baseline initial concentrations of that same step, and the concentrations in the step before that.
        Based on the magnitude of the measures, propose a course of action about what to do for the next step.

        :param n_chems:         The total number of registered species - exclusive of water and of macro-molecules
        :param indexes_of_active_chemicals: The ordered list (numerically sorted) of the INDEX numbers of all the chemicals
                                                involved in ANY of the registered reactions,
                                                but NOT counting chemicals that always appear in a catalytic role in all the reactions they
                                                participate in

        :param delta_conc:      A numpy array of changes in concentrations for the chemicals of interest,
                                    across a simulation time step (typically, the current step a run in progress)
        :param baseline_conc:   A numpy array of baseline concentration values for those same chemicals,
                                    prior to the above change, at the start of a simulation time step
        :param prev_conc:       A numpy array of concentration values for those same chemicals,
                                    in the step prior to the current one (i.e. an "archive" value)
        :return:                A dict:
                                    "action"           - String with the name of the computed recommended action:
                                                            either "low", "stay", "high" or "abort"
                                    "step_factor"      - A factor by which to multiply the time step at the next iteration round;
                                                            if no change is deemed necessary, 1
                                    "norms"            - A dict of all the computed norm name/values (any of the norms, except norm_A,
                                                            may be missing)
                                    "applicable_norms" - The name of the norm that triggered the decision; if all norms were involved,
                                                            it will be "ALL"
        """
        if baseline_conc is not None:
            assert len(baseline_conc) == len(delta_conc), \
                f"adjust_timestep(): the number of entries in the passed array `delta_conc` ({len(delta_conc)}) " \
                f"does not match the number of entries in the passed array `baseline_conc` ({len(baseline_conc)})"

        assert n_chems == len(delta_conc), \
            f"adjust_timestep(): the number of entries in the passed array `delta_conc` ({len(delta_conc)}) " \
            f"does not match the number of registered chemicals ({n_chems})"


        # If some chemicals are not dynamically involved in the reactions
        # (i.e. if they don't occur in any reaction, or occur as enzyme),
        # restrict our consideration to only the dynamically involved ones
        # CAUTION: the concept of "active chemical" might change in future versions, where only SOME of
        #          the reactions are simulated
        #if self.species_data.number_of_active_chemicals() < n_chems:
        if len(indexes_of_active_chemicals) < n_chems:
            delta_conc = delta_conc[indexes_of_active_chemicals]
            #print(f"\nadjust_timestep(): restricting adaptive time step analysis to {n_chems} chemicals only; their delta_conc is {delta_conc}")
            if baseline_conc is not None:
                baseline_conc = baseline_conc[indexes_of_active_chemicals]
            if prev_conc is not None:
                prev_conc = prev_conc[indexes_of_active_chemicals]
            # Note: setting delta_conc, etc, only affects local variables, and won't mess up the arrays passed as arguments

        all_norms = {}

        all_small = True            # Provisional answer to the question: "do ALL the rule yield a 'low'?"
        # Any rule failing to yield a 'low' will flip this status

        high_seen_at = []           # List of rule names at which a "high" is encountered, if applicable

        for rule in self.thresholds:
            norm_name = rule["norm"]

            if norm_name == "norm_A":
                result = self.norm_A(delta_conc)
            elif norm_name == "norm_B":
                result = self.norm_B(baseline_conc, delta_conc)
            elif norm_name == "norm_C":
                result = self.norm_C(prev_conc, baseline_conc, delta_conc)
            else:
                result = self.norm_D(prev_conc, baseline_conc, delta_conc)

            all_norms[norm_name] = result

            if ("abort" in rule) and (result > rule["abort"]):
                # If any rules declares an abort, no need to proceed further: it's an abort
                #self.norm_usage[norm_name] += 1
                self.increase_norm_count(norm_name)
                return {"action": "abort", "step_factor": self.step_factors["abort"], "norms": all_norms, "applicable_norms": [norm_name]}

            if ("high" in rule) and (result > rule["high"]):
                # If any rules declares a "high", still need to consider the other rules - in case any of them over-rides
                # the "high" with an "abort"
                high_seen_at.append(norm_name)
                all_small = False           # No longer the case of all 'low` yields

            if all_small and ("low" in rule) and (result > rule["low"]):
                all_small = False           # No longer the case of all 'low` yields
        # END for


        if high_seen_at:
            for n in high_seen_at:
                #self.norm_usage[n] += 1
                self.increase_norm_count(n)
            return {"action": "high", "step_factor": self.step_factors["downshift"], "norms": all_norms, "applicable_norms": high_seen_at}


        if all_small:
            for i in self.norm_usage:
                # All the norms were used
                self.increase_norm_count(i)
                #self.norm_usage[i] += 1

            return {"action": "low", "step_factor": self.step_factors["upshift"], "norms": all_norms, "applicable_norms": "ALL"}


        # If we get thus far, none of the rules were found applicable
        return {"action": "stay", "step_factor": 1, "norms": all_norms, "applicable_norms": "ALL"}



    def relative_significance(self, value :float, baseline :float) -> str:
        """
        Estimate, in a loose categorical fashion, the magnitude of the quantity "value"
        in proportion to the quantity "baseline".
        Both are assumed non-negative (NOT checked.)
        Return one of:
            "S" ("Small" ; up to 1/2 the size)
            "C" ("Comparable" ; from 1/2 to double)
            "L" ("Large" ; over double the size)
        This method is meant for large-scale computations, and on purpose avoids doing divisions.

        NOT IN CURRENT ACTIVE USAGE (in former use for the discontinued substep implementation)

        :param value:
        :param baseline:
        :return:        An assessment of relative significance, as one of
                        "S" ("Small"), "C" ("Comparable"), "L" ("Large")
        """
        #TODO: a Numpy array version

        if value < baseline:
            if value + value < baseline:
                return "S"
            else:
                return "C"

        else:
            if baseline + baseline < value:
                return "L"
            else:
                return "C"





    #####################################################################################################

    '''                                         ~  NORMS  ~                                           '''

    def ________NORMS________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def reset_norm_usage_stats(self):
        """
        Reset the count of the number of times that each norm got involved in step-size decision

        :return:    None
        """
        self.norm_usage = {"norm_A": 0, "norm_B": 0, "norm_C": 0, "norm_D": 0}



    def increase_norm_count(self, norm_name :str) -> None:
        """

        :param norm_name:
        :return:            None
        """
        assert norm_name in self.norm_usage, \
            f"increase_norm_count(): unknown norm named `{norm_name}`"

        self.norm_usage[norm_name] += 1



    def norm_A(self, delta_conc :np.array) -> float:
        """
        Return a measure of system change, based on the average concentration changes
        of ALL the specified chemicals across a time step, adjusted for the number of chemicals.
        A square-of-sums computation (the square of an L2 norm) is used.

        :param delta_conc:  A Numpy array with the concentration changes
                                of the chemicals of interest across a time step
        :return:            A measure of change in the concentrations across the simulation step
        """
        n_active_chemicals = len(delta_conc)

        assert n_active_chemicals > 0, \
            "norm_A(): zero-sized array was passed as argument"

        # The following are normalized by the number of chemicals
        #L2_rate = np.linalg.norm(delta_concentrations) / n_chems
        #L2_rate = np.sqrt(np.sum(delta_concentrations * delta_concentrations)) / n_chems
        #print("    L_inf norm:   ", np.linalg.norm(delta_concentrations, ord=np.inf) / delta_time)
        #print("    Adjusted L1 norm:   ", np.linalg.norm(delta_concentrations, ord=1) / n_chems)

        adjusted_L2_rate = np.sum(delta_conc * delta_conc) / (n_active_chemicals * n_active_chemicals)
        return adjusted_L2_rate



    def norm_B(self, baseline_conc: np.array, delta_conc: np.array) -> float:
        """
        Return a measure of system change, based on the max absolute relative concentration
        change of all the chemicals across a time step (based on an L infinity norm - but disregarding
        any baseline concentration that is very close to zero)

        :param baseline_conc:   A Numpy array with the concentration of the chemicals of interest
                                    at the start of a simulation time step
        :param delta_conc:      A Numpy array with the concentration changes
                                    of the chemicals of interest across a time step
        :return:                A measure of change in the concentrations across the simulation step
        """
        arr_size = len(baseline_conc)
        assert len(delta_conc) == arr_size, "norm_B(): mismatch in the sizes of the 2 passed array arguments"

        to_keep = ~ np.isclose(baseline_conc, 0)    # Element-wise negation; this will be an array of Booleans
        # with True for all the elements of baseline_conc that aren't too close to 0

        ratios = delta_conc[to_keep] / baseline_conc[to_keep]   # Using boolean indexing to only select some of the elements :
        # the non-zero denominators, and their corresponding numerators
        if len(ratios) == 0:
            return 0.

        return max(abs(ratios))



    def norm_C(self, prev_conc: np.array, baseline_conc: np.array, delta_conc: np.array) -> float:
        """
        Return a measure of system short-period oscillation; larger values might be heralding
        onset of simulation instability

        :param prev_conc:       A numpy array with the concentration of the chemicals of interest,
                                    in the step prior to the current one (i.e. an "archive" value)
        :param baseline_conc:   A Numpy array with the concentration of the chemicals of interest
                                    at the start of a simulation time step
        :param delta_conc:      A Numpy array with the concentration changes
                                    of the chemicals of interest across a time step
        :return:
        """
        if prev_conc is None:
            return 0            # Unable to compute a norm; 0 represents "perfect"

        D1 = baseline_conc - prev_conc
        D2 = delta_conc

        #print("****** D1: ", D1)
        #print("****** D2: ", D2)

        sign_flip = ((D1 >= 0) & (D2 < 0)) | ((D1 < 0) & (D2 >= 0))

        criterion_met = sign_flip & (abs(D2) > abs(D1)) & ~np.isclose(D1, 0) & (abs(D2) < 50 * abs(D1))
        #criterion_met = sign_flip & ~np.isclose(D1, 0) & (abs(D2) < 50 * abs(D1))
        # Note: values with D1 very close to zero are ignored
        #       likewise, values where |D2| dwarfs |D1| are ignored
        #print("criterion_met: ", criterion_met)

        # Use boolean indexing to select elements from delta1 where the criterion is True
        D1_selected = D1[criterion_met]
        D2_selected = D2[criterion_met]

        #print("**** delta1_selected: ", D1_selected)
        #print("*** delta2_selected: ", D2_selected)

        ratios = abs(D2_selected / D1_selected)
        #print(ratios)

        if len(ratios) == 0:
            return 0
        else:
            return np.max(ratios)    # An argument might be made for taking the avg SUM instead



    def norm_D(self, prev_conc: np.array, baseline_conc: np.array, delta_conc: np.array) -> float:
        """
        Return a measure of curvature in the concentration vs. time curves; larger values might be heralding
        onset of simulation instability

        :param prev_conc:       A numpy array with the concentration of the chemicals of interest,
                                    in the step prior to the current one (i.e. an "archive" value)
        :param baseline_conc:   A Numpy array with the concentration of the chemicals of interest
                                    at the start of a simulation time step
        :param delta_conc:      A Numpy array with the concentration changes
                                    of the chemicals of interest across a time step
        :return:
        """
        if prev_conc is None:
            return 0            # Unable to compute a norm; 0 represents "perfect"

        D1 = baseline_conc - prev_conc  # Change from prev state to current one
        D2 = delta_conc                 # Change from current state to next one

        #print("\nnorm_D ****** D1: ", D1)
        #print("norm_D ****** D2: ", D2)

        criterion_met = ~np.isclose(D1, 0) & (abs(D2) < 100 * abs(D1))
        # Note: values with D1 very close to zero are ignored
        #       likewise, values where |D2| dwarfs |D1| are ignored
        #print("criterion_met: ", criterion_met)

        # Use boolean indexing to select elements from delta1 where the criterion is True
        D1_selected = D1[criterion_met]
        D2_selected = D2[criterion_met]

        #print("norm_D **** delta1_selected: ", D1_selected)
        #print("norm_D *** delta2_selected: ", D2_selected)

        ratios = abs(D2_selected / D1_selected) # How big are the next changes relative to the previous ones
        #print(ratios)

        res = np.sum(ratios)    # An argument might be made for taking the MAX instead
        # TODO: turn into a separate norm
        arr_size = len(baseline_conc)
        normalized_res = res / arr_size
        #print("norm_D *** : ", normalized_res)

        return normalized_res
