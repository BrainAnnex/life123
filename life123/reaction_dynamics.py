# Classes ReactionDynamics and VariableTimeSteps:

import numpy as np
from life123.reaction import Reaction



class ReactionDynamics:
    """
    Static methods about the dynamic evolution of reactions
    """

    @classmethod
    def solve_exactly(cls, rxn :Reaction, A0 :float, B0 :float, t_arr) -> (np.array, np.array):
        """
        Return the exact solution of the given reaction,
        PROVIDED that it is a 1st Order Reaction of the type A <=> B.

        Use the given initial conditions,
        and return the solutions sampled at the specified times.

        For details, see https://life123.science/reactions

        :param rxn:     Object of type "Reaction", containing data for the reaction of interest
        :param A0:      The initial concentration of the reactant A
        :param B0:      The initial concentration of the product B
        :param t_arr:   A Numpy array with the desired times at which the solutions are desired
        :return:        A pair of Numpy arrays with, respectively, the concentrations of A and B
        """

        reactants, products, kF, kR = rxn.unpack_for_dynamics()

        assert len(reactants) == 1, "Currently only works for `A <-> B` reactions"
        assert len(products) == 1, "Currently only works for `A <-> B` reactions"
        assert rxn.extract_stoichiometry(reactants[0]) == 1, \
            "Currently only works for `A <-> B` reactions"
        assert rxn.extract_stoichiometry(products[0]) == 1, \
            "Currently only works for `A <-> B` reactions"
        # TODO: should also verify the reaction orders to be 1

        return cls._exact_solution(kF, kR, A0, B0, t_arr)



    @classmethod
    def _exact_solution(cls, kF, kR, A0, B0, t_arr :np.ndarray) -> (np.ndarray, np.ndarray):
        """
        Return the exact solution of the 1st Order Reaction A <=> B,
        with the specified parameters,
        sampled at the given times.

        For details, see https://life123.science/reactions   (TODO: make B(t) more symmetric, as suggested by ChatGPT)

        :param kF:
        :param kR:
        :param A0:
        :param B0:
        :param t_arr:   A Numpy array with the desired times at which the solutions are desired
        :return:        A pair of Numpy arrays
        """
        TOT = A0 + B0
        # Formula is:  A = (A0 - (kR TOT) / (kF + kR)) Exp[-(kF + kR) t] + kR TOT / (kF + kR)

        sum_rates = kF + kR
        A_arr = (A0 - (kR * TOT) / sum_rates) * np.exp(-sum_rates * t_arr) + (kR * TOT / sum_rates)
        B_arr = TOT - A_arr
        return (A_arr, B_arr)







###########################################################################################

class VariableTimeSteps:
    """
    Methods for managing variable time steps during reactions
    """
    # TODO: maybe rename "VariableTimeSteps"

    def __init__(self):
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



    def display_value_against_thresholds(self, all_norms):
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
                    dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
                    dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
                    dynamics.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4, error=0.25)

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

        :param n_chems:         The total number of registered chemicals - exclusive of water and of macro-molecules
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
        #if self.chem_data.number_of_active_chemicals() < n_chems:
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
