# Classes:
#   1) ReactionSimulator
#   2) AnalyticalReactionSolver
#   3) VariableTimeSteps:

import math
import cmath
import numpy as np
import time
from dataclasses import dataclass, field
from typing import Any
from life123.species_index_map import SpeciesIndexMap
from life123.diagnostics import Diagnostics
from life123.history import HistoryUniformConcentration, HistoryReactionRate




class ExcessiveTimeStepHard(Exception):
    """
    Used to raise Exceptions arising from excessively large time steps
    (that lead to negative concentration values, i.e. "HARD" errors)

    The single argument passed in the call to:
        raise ExcessiveTimeStepHard(my_exception_details)

    is made conveniently available as a `details` attribute in the caught Exception:
        except ExcessiveTimeStepHard as ex:
            ex.details      # Contains exactly what you passed in my_exception_details
    """
    def __init__(self, details):
        super().__init__(details)
        self.details = details


class ExcessiveTimeStepSoft(Exception):
    """
    Used to raise Exceptions arising from excessively large time steps
    (that lead to norms regarded as excessive because of user-specified values,
     or to other signs of overshoots, i.e. "SOFT" errors)
    """
    def __init__(self, details):
        super().__init__(details)
        self.details = details




#############################################################################################

@dataclass(slots=True)      # Note: (slots=True) has the effect of prohibiting non-listed fields,
                            #       and of making the class more efficient
class ReactionStep:
    """
    Diagnostic data for the last reaction step.
    Note that if an intended step, in variable-step mode, lead to multiple steps,
    this data gets overwritten at each step
    """
    delta_time : float                          # REQUIRED
                                                # Note: this is the size of the step eventually ACTUALLY taken,
                                                #       not necessarily what first attempted (in the case of variable steps)

    system_rxn_rates: dict[int, float] | None = field(default_factory=dict)
                                        # Reaction rates at the start of the last (current) step, for each of the involved reactions
                                        # Keys are the reaction indexes
                                        # EXAMPLE: {0: 0.42, 1: 6.2, 2: -3.7}


    # All the attributes below are only applicable to variable steps

    number_neg_concs : int = 0                  # Zero default upon starting a new sim step
    number_soft_aborts : int = 0                # Zero default upon starting a new sim step

    norms: dict[str, Any] | None = None         #
                                                # EXAMPLE: {'norm_A': 0.98, 'norm_B': 0.14}
    decision_data: dict[str, Any] | None = None #
                                                # EXAMPLE: {'action': 'stay', 'operation': 'stay',
                                                #           'step_factor': 1,
                                                #           'applicable_norms': 'ALL'}



#######################################################################

class ReactionSimulator:
    """
    Used to simulate the dynamics of reactions (in a single compartment.)
    This might be thought of as a "zero-dimensional system".
    """
    
    def __init__(self, system=None, species_index_map=None, reaction_registry=None, exact=True,
                 diagnostics=None, diagnostics_enabled=False, method="forward_euler", preset="mid",
                 uniform_compartment=None):
        """

        :param system:
        :param species_index_map:
        :param reaction_registry:
        :param exact:
        :param diagnostics:
        :param diagnostics_enabled:
        :param method:
        :param preset:
        :param uniform_compartment:
        """

        if uniform_compartment is not None:
            self.uniform_compartment = uniform_compartment  # Use "UniformCompartment" object to manage the species index
        else:
            self.uniform_compartment = self         # Self-manage the species index

        self.system :np.ndarray = system

        self.species_index_map :SpeciesIndexMap = species_index_map 
        self.reaction_registry = reaction_registry
        self.exact :bool = exact

        self.diagnostics_enabled = diagnostics_enabled  # Flag indicating whether using diagnostics
        self.diagnostics = diagnostics  # Object of class "Diagnostics"
        self.diagnostic_data_snapshot = {}      # TODO: unclear if still needed, now that we have
                                                #       `self.reaction_step_diagnostics.system_rxn_rates` (maybe merge?)

        self.reaction_step_diagnostics = None   # Internal data about a single reaction step
                                                #   (taken in full OR aborted!)
                                                #   Object of dataclass "ReactionStep"


        self.method :str = method

        self.system_time = 0.       # Global time of the system, from initialization
        self.previous_system = None # Concentration data of all the species at the previous simulation step


        self.rate_history = HistoryReactionRate(active=True)        # Object used to store user-requested snapshots
                                                                    # of (some of) the chemical reaction rates:
                                                                    # 'SYSTEM TIME', 'rxn0_rate', 'rxn1_rate', ...

        self.conc_history = HistoryUniformConcentration(active=True)    # Object used to store user-requested snapshots
                                                                        # of (some of) the species concentrations:
                                                                        # 'SYSTEM TIME', 'A', 'B', ..., 'comments'

        if (species_index_map is None) and (uniform_compartment is None) and (reaction_registry is not None):
            # Assign an array index to all the species we're dealing with    TODO: use set_system_conc()
            all_species = reaction_registry.get_species_in_any_reaction(sort=True)
            self.species_index_map = SpeciesIndexMap(all_species)


        # FOR AUTOMATED ADAPTIVE TIME STEP SIZES
        self.adaptive_steps = VariableTimeSteps()

        if (self.diagnostics_enabled) and (not self.diagnostics):
            self.diagnostics = Diagnostics(reactions=self.reaction_registry, species_to_index=self.uniform_compartment.species_index_map.species_to_index,
                                           species_index_map=self.uniform_compartment.species_index_map)

        if preset:
            self.adaptive_steps.use_adaptive_preset(preset)




    #########  MISC. UTILITIES  #########

    def number_of_system_species(self) -> int:
        """
        Number of species being simulated (and kept in the system state)
        
        :return:
        """
        # TODO: alternatively, use the size of self.system
        return self.uniform_compartment.species_index_map.number_of_system_species()
 
 
    def get_reactions(self):
        """
        Return all the reactions associated to this Uniform Compartment

        :return:    Object ot type "ReactionRegistry" (with data about all the reactions)
        """
        return self.reaction_registry   



    def indexes_of_active_chemicals(self) -> list[int]:
        """
        Return the ordered list (numerically SORTED) of the INDEX numbers of all the species
        involved in ANY of the registered reactions,
        but NOT counting species that always appear in a catalytic role in all the reactions they
        participate in
        (if a species participates in a non-catalytic role in ANY reaction, it'll appear here.)

        EXAMPLE: [2, 7, 8]  if only those 3 chemicals (with indexes of, respectively, 2, 7 and 8)
                            are actively involved in ANY of the registered reactions

        CAUTION: the concept of "active species" might change in future versions, where only SOME of
                 the reactions are simulated
        """
        set_active_species = self.get_reactions().active_chemicals

        index_list = list(
                            map(lambda species_id: self.uniform_compartment.species_index_map.index_of(species_id), set_active_species)
                         )
        return sorted(index_list)        



    def capture_rate_snapshot(self, step_count=None, caption=None,
                             force=False, system_time=None) -> None:
        """
        Take the reaction rates for the last (current) step of all reactions,
        stored as a dict in the object property `self.reaction_step_diagnostics.system_rxn_rates`,
        and save them in the ongoing table storing the "rate history".

        Note: reaction rate history is regarded as akin to concentration history - something to periodically
              sample and record, separate from diagnostics

        :param step_count:  [OPTIONAL] Step count in the simulation; used to possibly pare down the frequency of snapshot saving
        :param caption:     [OPTIONAL] String to save as a caption field, alongside the rate fields
        :param force:       [OPTIONAL] If True, take a snapshot regardless of step_count; default is False
        :param system_time: [OPTIONAL] If not specified, use the current SYSTEM TIME
        :return:            None
        """
        if not force and not self.rate_history.to_capture(step_count):
            return

        if system_time is None:
            system_time = self.system_time

        data_snapshot = {}     # rxn_rates_snapshot
        for k, v in self.reaction_step_diagnostics.system_rxn_rates.items():
            if type(v) == tuple:
                # Only pairs are currently used (for sub-reactions of enzymatic reactions) -> TODO: this is being phased out
                data_snapshot[f"rxn{k}_rate_1"] = v[0]      # EXAMPLE:  "rxn4_rate_1" = 18.2
                data_snapshot[f"rxn{k}_rate_2"] = v[1]      # EXAMPLE:  "rxn4_rate_2" = 5.52
            else:
                data_snapshot[f"rxn{k}_rate"] = v      # EXAMPLE:  "rxn4_rate" = 18.2
        '''
           EXAMPLE of data_snapshot:
                {"rxn1_rate": 6.3, "rxn2_rate": 14.3}        
        '''
        self.rate_history.save_snapshot(step_count=step_count, system_time=system_time,
                                        data_snapshot=data_snapshot,
                                        caption=caption)



    def capture_conc_snapshot(self, step_count=None, caption="", extra=False) -> None:
        """

        :param step_count:  [OPTIONAL] Step count in the simulation
        :param caption:     [OPTIONAL]
        :param extra:       [OPTIONAL] If True, it means that this is a special extra capture;
                                the capture frequency will NOT considered in the
                                decision about saving it, but a check will be performed
                                to make sure it's not a duplicate of the earlier capture

        :return:                    None
        """
        if not self.conc_history.to_capture(step_count, extra=extra):
            return

        data_snapshot = self.get_conc_dict(chem_labels=self.conc_history.restrict_chemicals)
        '''
           EXAMPLE of data_snapshot:
                {"A": 1.3, "B": 4.9}        
        '''
        self.conc_history.save_snapshot(step_count=step_count, system_time=self.system_time,
                                        data_snapshot=data_snapshot,
                                        caption=caption)


    def set_system_conc(self, conc_array :np.ndarray):
        """

        :param conc_array:
        :return:
        """
        self.system = conc_array
        if (self.uniform_compartment == self) and (self.reaction_registry is not None):
            # Assign an array index to all the species we're dealing with
            all_species = self.reaction_registry.get_species_in_any_reaction(sort=True)
            self.species_index_map = SpeciesIndexMap(all_species)



    def get_species_conc(self, label :str) -> float:
        """
        Return the current system concentration of the given species, specified by its id.
        If no species by that name exists, an Exception is raised

        :param label:   The id of a species
        :return:        The current system concentration of the above chemical
        """
        # TODO: this ought to be returned to UniformCompartment
        species_index = self.uniform_compartment.species_index_map.index_of(label)
        return self.system[species_index]



    def get_conc_dict(self, chem_labels=None, system_data=None) -> dict|None:
        """
        Retrieve the concentrations of the requested species (by default all),
        as a dictionary indexed by the species id

        :param chem_labels: [OPTIONAL] List, tuple or set of the id's of the species;
                                by default, return all
        :param system_data: [OPTIONAL] A Numpy array of concentration values, in the same order as the
                                index of the chemical species; by default, use the SYSTEM DATA

        :return:            A dictionary, indexed by the chemical labels, of the concentration values;
                                or None if no data available
                                EXAMPLE: {"A": 1.2, "D": 4.67}
        """
        # TODO: this ought to be returned to UniformCompartment
        # TODO: probably change the None returns to empty dict's
        if system_data is None:
            system_data = self.system
        else:
            assert system_data.size == self.number_of_system_species(), \
                f"UniformCompartment.get_conc_dict(): the argument `system_data` must be a 1-D Numpy array with as many entries " \
                f"as the declared number of chemicals ({self.number_of_system_species()})"


        if chem_labels is None:
            if system_data is None:
                return {}
            else:
                return {self.uniform_compartment.species_index_map.species_at(index): system_data[index]
                        for index, conc in enumerate(system_data)}
        else:
            assert isinstance(chem_labels, (list, tuple, set)), \
                f"UniformCompartment.get_conc_dict(): the argument `species` must be a list or tuple" \
                f" (it was of type {type(chem_labels)})"

            conc_dict = {}
            for name in chem_labels:
                species_index = self.uniform_compartment.species_index_map.index_of(name)
                conc_dict[name] = system_data[species_index]

            return conc_dict



    def _save_soft_abort_diagnostics(self, delta_concentrations):
        # TODO: maybe this call could wait until the interception of ExcessiveTimeStepSoft?

        if self.diagnostics_enabled:
            # Define the dict `self.diagnostic_data_snapshot`
            all_norms = self.reaction_step_diagnostics.norms
            self.diagnostic_data_snapshot['norm_A'] = all_norms.get('norm_A')
            self.diagnostic_data_snapshot['norm_B'] = all_norms.get('norm_B')
            self.diagnostic_data_snapshot['norm_C'] = all_norms.get('norm_C')
            self.diagnostic_data_snapshot['norm_D'] = all_norms.get('norm_D')
            self.diagnostic_data_snapshot['action'] = "ABORT"
            self.diagnostic_data_snapshot['step_factor'] = self.reaction_step_diagnostics.decision_data["step_factor"]
            self.diagnostic_data_snapshot['time_step'] = self.reaction_step_diagnostics.delta_time
            self.diagnostics.save_diagnostic_decisions_data(system_time=self.system_time,
                                                data=self.diagnostic_data_snapshot, delta_conc_arr=delta_concentrations,
                                                caption="excessive norm value(s)")

            # Make a note of the abort action in all the reaction-specific diagnostics
            self.diagnostics.annotate_abort_rxn_data("aborted: excessive norm value(s)")



    def _save_hard_abort_diagnostics_single_rxn(self, exception_data :dict):
        """
        Save 1 diagnostic entry under "decision data" diagnostics
        and 1 under "rxn_data" diagnostics

        :param exception_data:
        :return:
        """
        # Save under "decision data"
        diagnostics_data = {  "action": "ABORT",
                              "caption": f"neg. conc. in {exception_data['species_id']} from rxn # {exception_data['rxn_index']}",
                              "time_step": exception_data["delta_time"]}
        if self.adaptive_steps:     # Add more diagnostics data if available
            diagnostics_data["step_factor"] = self.adaptive_steps.step_factors.get('error')

        self.diagnostics.save_diagnostic_decisions_data(system_time=self.system_time,
                                                        data=diagnostics_data,
                                                        delta_conc_arr=None)

        # Save under "rxn_data" diagnostics
        self.diagnostics.save_rxn_data(rxn_index=exception_data["rxn_index"], system_time=self.system_time,
                                       time_step=exception_data["delta_time"],
                                       increment_dict_single_rxn=None,
                                       aborted=True,
                                       rate=exception_data["rate"], caption=f"aborted: neg. conc. in `{exception_data['species_id']}`")



    def _save_hard_abort_diagnostics_all_rxns(self, explain_variable_steps, delta_time):
        """

        :param explain_variable_steps:
        :param delta_time:
        :return:
        """
        # TODO: maybe this call could wait until the interception of ExcessiveTimeStepHard ?
        # All this can way until the interception of ExcessiveTimeStepHard
        if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
            print(f"*** CAUTION: negative concentration resulting from the combined effect of all reactions, "
                  f"upon advancing reactions from system time t={self.system_time:,.5g}\n"
                  f"         It'll be AUTOMATICALLY CORRECTED with a reduction in time step size")

        # A type of HARD ABORT is detected and raised (a negative concentration resulting from the combined effect of all reactions)
        if self.diagnostics_enabled:
            self.diagnostics.save_diagnostic_decisions_data(system_time=self.system_time,
                                                            data={"action": "ABORT",
                                                                  "step_factor": self.adaptive_steps.step_factors["error"],
                                                                  "caption": "neg. conc. from combined effect of all rxns",
                                                                  "time_step": delta_time})
            # Save up diagnostic data for ALL reactions
            self.diagnostics.save_diagnostic_aborted_rxns(system_time=self.system_time, time_step=delta_time,
                                                         caption=f"aborted: neg. conc. from combined multiple rxns")
    


    def specify_steps(self, duration=None, initial_step=None, n_steps=None) -> (float, int):
        """
        If either the `initial_step` or `n_steps` is not provided (but at least 1 of them must be present),
        determine the other one from `duration`

        Their desired relationship is: total_duration = time_step * n_steps

        :param duration:    Float with the overall time advance (i.e. time_step * n_steps)
        :param initial_step:Float with the size of each time step
        :param n_steps:     Integer with the desired number of steps
        :return:            The pair (time_step, n_steps)
        """
        DEFAULT_N_STEPS = 10

        assert (not duration or not initial_step or not n_steps), \
            "UniformCompartment.specify_steps(): cannot specify all 3 arguments: " \
            "`duration`, `initial_step`, `n_steps` (specify at most 2 of them)"

        assert (duration or initial_step), \
            "UniformCompartment.specify_steps(): at least 1 of the arguments `duration` and `initial_step` must be provided"

        if (not initial_step) and (not n_steps):    # If only `duration` is provided
            n_steps = DEFAULT_N_STEPS
            initial_step = duration
            return (initial_step, n_steps)

        if (not duration) and (not n_steps):        # If only `initial_step` is provided
            n_steps = DEFAULT_N_STEPS
            return (initial_step, n_steps)

        # If we get this far, both `duration` and `initial_step` must be present

        if not initial_step:
            initial_step = duration / n_steps

        if not n_steps:
            n_steps = math.ceil(duration / initial_step)

        return (initial_step, n_steps)     # Note: could opt to also return total_duration if there's a need for it





    ###########   THE MAIN SIMULATION SEQUENCE STARTS HERE   ###########


    def single_compartment_react(self, duration=None, target_end_time=None, stop=None,
                                 initial_step=None, n_steps=None, max_steps=None,
                                 variable_steps=True, explain_variable_steps=None,
                                 silent=False, report_interval=0.5) -> dict:
        """
        Simulate ALL the (previously-registered) reactions in the single uniform ("well-stirred") compartment -
        based on the INITIAL concentrations stored in self.system

        Update the system state and the system time accordingly
        (object attributes self.system and self.system_time)

        :param duration:        [OPTIONAL] The overall time advance for the reactions (it may be exceeded in case of variable steps)
        :param target_end_time: [OPTIONAL] The final time at which to stop the reaction; it may be exceeded in case of variable steps
                                    If both `target_end_time` and `duration` are specified, an error will result

        :param initial_step:    [OPTIONAL] The suggested size of the first step (it might be reduced automatically,
                                    in case of "hard" errors resulting from overly-large steps)

        :param stop:            [OPTIONAL] Pair of the form (termination_keyword, termination_parameter), to indicate
                                    the criterion to use to stop the reaction
                                    EXAMPLES:
                                        ("conc_below", (chem_name, conc))  Stop when conc first dips below
                                        ("conc_above", (chem_name, conc))  Stop when conc first rises above

                                        TODO: add more options, such as
                                        ("before_time", t)                  Stop just before the given target time
                                        ("after_time", t)                   Stop just after the given target time
                                        ("equilibrium", tolerance)          Stop when equilibrium reached

        :param n_steps:         [OPTIONAL] The desired number of steps

        :param max_steps:       [OPTIONAL] Max numbers of steps; if reached, it'll terminate regardless of any other criteria

        :param variable_steps:          [OPTIONAL] If True (default), the steps sizes will get automatically adjusted, based on thresholds
        :param explain_variable_steps:  [OPTIONAL] If not None, a brief explanation is printed about how the variable step sizes were chosen,
                                            when the System time inside that range;
                                            only applicable if variable_steps is True

        :param silent:                  [OPTIONAL] If True, less output is generated
        :param report_interval:         [OPTIONAL] How frequently, in terms of elapsed running time, in minutes,
                                            to inform the user of the current status

        :return:                        A dictionary containing the key "recommended_next_step"
                                        Note: the object attributes self.system and self.system_time get updated
        """

        # Default values
        initial_step_caption = "1st reaction step"
        final_step_caption = "last reaction step"

        # Validation
        assert self.system is not None, "UniformCompartment.single_compartment_react(): " \
                                        "the concentration values of the various chemicals must be set first"

        if variable_steps and (n_steps is not None):
            raise Exception("UniformCompartment.single_compartment_react(): if `variable_steps` is True, cannot specify `n_steps` "
                            "(because the number of steps will vary); specify `duration` or `target_end_time` instead")

        if stop is not None:
            assert type(stop) == tuple and len(stop) == 2, \
                f"UniformCompartment.single_compartment_react(): the argument `stop`, if passed, " \
                f"must be a pair of values, of the form (termination_keyword, termination_parameter)"

        assert self.reaction_registry.number_of_reactions() > 0, \
            f"UniformCompartment.single_compartment_react(): no reactions are present.  Make sure to first add them with add_reaction()"


        self.conc_history.initial_caption = initial_step_caption     # TODO: turn into method


        """
        Determine all the various time parameters that were not explicitly provided
        """

        if stop is not None:
            assert initial_step > 0, \
                "single_compartment_react(): when using the `stop` argument, an `initial_step` argument must be provided"
            assert max_steps is not None, \
                "single_compartment_react(): when using the `stop` argument, a `max_steps` argument must be provided"
            time_step = initial_step

        else:
            if target_end_time is not None:
                if duration is not None:
                    raise Exception("single_compartment_react(): cannot provide values for BOTH `target_end_time` and `duration`")
                else:
                    assert target_end_time > self.system_time, \
                        f"single_compartment_react(): `target_end_time` must be larger than the current System Time ({self.system_time})"
                    duration = target_end_time - self.system_time

            # Determine the time step,
            # as well as the required number of such steps
            # TODO: if the following call results in an Exception, the reported arguments are confusing because
            #       the names don't match
            time_step, n_steps = self.specify_steps(duration=duration,
                                                    initial_step=initial_step,
                                                    n_steps=n_steps)
            #print(f"time_step: {time_step} , n_steps: {n_steps}")
            # Note: if variable steps are requested then `n_steps` stops being particularly meaningful; it becomes a
            #       hypothetical value, in the (unlikely) event that the step sizes were never changed - and is only
            #       used to detect a very excessive number of actual attempted steps

            if target_end_time is None:
                if variable_steps:
                    target_end_time = self.system_time + duration
                else:
                    target_end_time = self.system_time + time_step * n_steps


        step_count = 0

        # Reset some diagnostic variables
        #self.reaction_step_diagnostics.number_neg_concs = 0
        #self.reaction_step_diagnostics.number_soft_aborts = 0
        self.adaptive_steps.reset_norm_usage_stats()

        # Time-related
        t_start = time.perf_counter()
        t_report = t_start
        report_interval *= 60.      # Convert to seconds


        try:
            while True:     # Loop until one of various criteria becomes applicable

                # Check various criteria for termination
                if (max_steps is not None) and (step_count >= max_steps):
                    print(f"single_compartment_react(): computation stopped because max # of steps ({max_steps}) reached")
                    break       # We have reached the max allowable number of steps

                if (target_end_time is not None) and (self.system_time >= target_end_time):
                    break       # The system time has reached the target endtime

                if (stop is not None):
                    (termination_keyword, termination_parameter) = stop
                    if termination_keyword == "conc_below":
                        chem_name, conc_threshold = termination_parameter
                        if self.get_species_conc(chem_name) < conc_threshold:
                            break   # The concentration of the specified chemical has dropped the requested threshold
                    elif termination_keyword == "conc_above":
                        chem_name, conc_threshold = termination_parameter
                        if self.get_species_conc(chem_name) > conc_threshold:
                            break   # The concentration of the specified chemical has risen above the requested threshold

                if (not variable_steps) and (step_count == n_steps)\
                        and (target_end_time is not None) and np.allclose(self.system_time, target_end_time):
                    break       # When dealing with fixed steps, catch scenarios where after performing n_steps,
                                #   the System Time is below the target_end_time because of roundoff error


                # ---  CORE OPERATION OF MAIN LOOP  ---
                step_count, recommended_next_step = self._single_compartment_react_main_loop(step_count=step_count, n_steps=n_steps,
                                                                    variable_steps=variable_steps, time_step=time_step,
                                                                    explain_variable_steps=explain_variable_steps)


                if self.diagnostics_enabled:
                    # Save up the current time and System State as "diagnostic 'concentration' data"
                    system_data = self.get_conc_dict(system_data=self.system)   # The current System State, as a dict
                    self.diagnostics.save_diagnostic_conc_data(system_data=system_data, system_time=self.system_time)

                t_now = time.perf_counter()
                t_elapsed = t_now - t_report    # Time elapsed since the last report
                if (not silent) and (t_elapsed > report_interval):
                    if variable_steps:
                        info_on_step = f"(doing step size {time_step:,.2g})"
                    else:
                        info_on_step = ""
                    print(f"... running : currently at System Time {self.system_time:,.4g} {info_on_step} after running for {(t_now - t_start)/60:.1f} min")
                    t_report = t_now            # Reset

                if variable_steps:
                    time_step = recommended_next_step   # Follow the recommendation of the ODE solver for the next time step to take

        # --- END while ---

        except KeyboardInterrupt:
            print("\n*** KeyboardInterrupt exception caught")



        # We're now at the end of the computation
        # Report whether extra steps were automatically added
        n_steps_taken = step_count


        if (not variable_steps) and (n_steps is not None):
            extra_steps = n_steps_taken - n_steps
            if extra_steps > 0:
                print(f"The computation took {extra_steps} extra step(s) - "
                      f"automatically added to prevent negative concentrations")


        if not silent:
            # Print out a summary, at the termination of the run
            t_now = time.perf_counter()
            step_type_str = "variable " if variable_steps else "fixed "
            time_taken = t_now - t_start
            if time_taken < 60:
                display_time_taken = f"{time_taken:.3f} sec"
            else:
                display_time_taken = f"{time_taken/60.:.2f} min"
            print(f"{n_steps_taken} total {step_type_str}step(s) taken in {display_time_taken}")
            if variable_steps:
                if self.reaction_step_diagnostics.number_neg_concs:
                    print(f"Number of step re-do's because of negative concentrations: {self.reaction_step_diagnostics.number_neg_concs}")
                if self.reaction_step_diagnostics.number_soft_aborts:
                    print(f"Number of step re-do's because of elective soft aborts: {self.reaction_step_diagnostics.number_soft_aborts}")

                print("Norm usage:", self.adaptive_steps.norm_usage)
                print(f"System Time is now: {self.system_time:,.5g}")


        # One final snapshot, unless already taken for the last step done
        self.capture_conc_snapshot(step_count=step_count, caption=final_step_caption, extra=True)

        # Add a caption to the very last entry in the system history
        self.conc_history.set_caption_last_snapshot(final_step_caption)

        return {"recommended_next_step": recommended_next_step}



    def _single_compartment_react_main_loop(self, time_step, variable_steps :bool,
                                            step_count, n_steps :int,
                                            explain_variable_steps=None) -> (int, float):
        """
        Helper function to single_compartment_react(), for its main loop.
        Perform the reaction step (either fixed or variable step, as appropriate),
        then UPDATE THE SYSTEM STATE, and preserve the RATES data.
        If the number of steps taken so far is getting excessive, raise an Exception

        :param time_step:               The requested time duration of the reaction step;
                                            in case of adaptive variable steps, a smaller step might actually get taken
        :param variable_steps:          If True, the steps sizes will get automatically adjusted, based on thresholds
        :param step_count:              Step count in the simulation (the initial value at start of simulation should be zero)
        :param n_steps:                 The desired number of steps;
                                            the actual number might be considerably higher in case of adaptive variable steps.
                                            Here it's used purely as a safeguard against infinite loops
        :param explain_variable_steps:  [OPTIONAL] If provided, it must be a pair of numbers of the form [t_start, t_end];
                                            a brief explanation will be printed about how the variable step sizes were chosen,
                                            when the System time inside that range.
                                            Only applicable if variable_steps is True
        :return:                        The pair (new step count, recommended next step);
                                            the latter is only meaningful in case of adaptive variable steps
        """
        # ----------  CORE OPERATION OF MAIN LOOP  ----------
        if variable_steps:
            delta_concentrations, step_actually_taken, recommended_next_step = \
                self.reaction_step_common_variable_step(delta_time=time_step,
                                                        explain_variable_steps=explain_variable_steps,
                                                        step_counter=step_count)
        else:
            # Fixed steps
            delta_concentrations = \
                self.reaction_step_common_fixed_step(delta_time=time_step, step_counter=step_count)
            step_actually_taken = time_step
            recommended_next_step = time_step


        # Update the System State
        self.previous_system = self.system.copy()   # Clone the earlier state (note: this is used for adaptive time steps. TODO: maybe jettison for fixed steps)
        self.system += delta_concentrations
        if min(self.system) < 0:    # Check for negative concentrations. TODO: redundant, since reaction_step_common() now does that
            print(f"***********  SYSTEM STATE ERROR: FAILED TO CATCH negative concentration "
                  f"upon advancing reactions from system time t={self.system_time:,.5g}")


        # Preserve the RATES data, as requested ("part1", BEFORE updating the System Time, because reaction rates are
        # based on the *start* time of the simulation step)
        if step_count == 0:
            self.capture_rate_snapshot(force=True, step_count=0)    # Always save the initial rate
        else:
            self.capture_rate_snapshot(step_count=step_count)       # Save historical rate values (if enabled)


        # UPDATE THE SYSTEM TIME and step count (now we're at the END of the current time step)
        # Note: we've taken exactly 1 more step - though not necessarily of the requested size, in case of adaptive variable steps
        self.system_time += step_actually_taken
        step_count += 1


        # Preserve the CONCENTRATION data, as requested ("part2", AFTER updating the System Time, because current concentrations
        # refer to the System Time at the end of the simulation step)
        self.capture_conc_snapshot(step_count=step_count) # Save historical concentration values (if enabled)
                                                          # We use the updated step count because we save the conc. values at the END of the step


        # The following is only applicable to variable steps (TODO: maybe move to calling module)
        if variable_steps:
            if (n_steps is not None) and (step_count > 1000 * n_steps):  # Another approach to catch infinite loops
                raise Exception(f"_single_compartment_react_main_loop(): "
                                f"the computation is taking a very large number of steps, probably from automatically trying to correct instability;"
                                f" try reducing the time_step")   # TODO: is the explanation correctly phrased?


        return step_count, recommended_next_step



    def reaction_step_common_fixed_step(self, delta_time: float, conc_array=None,
                                        step_counter=1) -> np.array:
        """
        This is the common entry point for FIXED-step simulations,
        both for single-compartment reactions,
        and for the reaction component of reaction-diffusions in 1D, 2D and 3D.

        "Compartments" may or may not correspond to the "bins" of the higher layers;
        the calling code might have opted to merge some bins into a single "compartment".

        Using the given concentration data for all the applicable species in a single compartment,
        do a single reaction time step for ALL the reactions -
        based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        Return the increment vector for all the chemical species concentrations in the compartment

        NOTES:  * the actual system concentrations are NOT changed
                * this method doesn't modify the step sizes: in case of any error caused by an excessively-large
                  step size, an Exception is raised.

        :param delta_time:      The requested time duration of the reaction step
        :param conc_array:      [OPTIONAL]All initial concentrations at the start of the reaction step,
                                    as a Numpy array for ALL the chemical species, in their index order.
                                    If not provided, self.system is used instead
        :param step_counter:    [OPTIONAL] Integer currently only used for diagnostics

        :return:                The increment vector for the concentrations of ALL the chemical species,
                                    in their array index order, as a Numpy array
                                    EXAMPLE (for a single-reaction reactant and product with a 1:3 stoichiometry):
                                        array([7. , -21.])
        """
        # TODO: no longer pass conc_array .  Use the object variable self.system instead
        #       Determine whether 1 or multiple UC objects are to be used by Bio1D, etc.

        if conc_array is not None:
            self.system = conc_array    # For historical reasons, as a convenience to Bio1D, etc.
                                        # TODO: maybe it ought to be kept separate in separate instances of the UC object
            if self.uniform_compartment:
                self.uniform_compartment.system = self.system   # Re-align

        # Validate the setup
        assert self.system is not None, "ReactionSimulator.reaction_step_common_fixed_step(): " \
                                        "the concentration values of the various species must be set first"

        #print(f"************ At SYSTEM TIME: {self.system_time:,.4g}, calling reaction_step_common_fixed_step() with:")
        #print(f"             delta_time={delta_time}, system={self.system}, ")

        try:
            delta_concentrations, _  =  \
                self.attempt_reaction_step(delta_time, variable_steps=False,
                                           explain_variable_steps=None, step_counter=step_counter)

        # CATCH any 'ExcessiveTimeStepHard' exception raised in the loop  (i.e. a HARD ABORT),
        #       in order to re-raise an Exception with expanded error message and additional error data
        except ExcessiveTimeStepHard as ex:
            # Single reactions steps can fail with this error condition if the attempted time step was too large,
            # under the following scenarios:
            #       1. negative concentrations from any one reaction
            #       2. negative concentration from the combined effect of multiple reactions
            #print("*** CAUGHT a HARD ABORT in reaction_step_common_fixed_step()")
            exception_data = ex.details     # Start with the details of the caught exception...
            # ... and add/edit some
            exception_data["previous_function"] = exception_data["function"]
            exception_data["function"] = "reaction_step_common_fixed_step"   # Overwrite
            old_message = exception_data.get("message")
            exception_data["message"] = f"reaction_step_common_fixed_step(): unable to complete the reaction step.  " \
                                        f"Try REDUCING the time step, or switching to variable time steps. \n" \
                                        f"DETAILS: \n{old_message}"

            raise ExcessiveTimeStepHard(exception_data)


        return  delta_concentrations    # TODO: consider returning tentative_updated_system , since we already computed it



    def reaction_step_common_variable_step(self, delta_time: float, conc_array=None,
                                           explain_variable_steps=None, step_counter=1) -> (np.array, float, float):
        """
        TODO: corresponds to reaction_step_common() in UniformCompartment -> to phase out

        This is the common entry point for VARIABLE-step simulations,
        both for single-compartment reactions,
        and for the reaction component of reaction-diffusions in 1D, 2D and 3D.

        "Compartments" may or may not correspond to the "bins" of the higher layers;
        the calling code might have opted to merge some bins into a single "compartment".

        Using the given concentration data for all the applicable species in a single compartment,
        do a single reaction time step for ALL the reactions -
        based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        Return the increment vector for all the chemical species concentrations in the compartment

        NOTES:  * the actual system concentrations are NOT changed
                * this method doesn't decide on step sizes - except in case of ("hard" or "soft") aborts, which are
                    followed by repeats with a smaller step.  Also, it makes suggestions
                    to the calling module about the next step to best take (whether as a result of an abort,
                    or for other considerations)

        :param delta_time:      The requested time duration of the reaction step
        :param conc_array:      [OPTIONAL] All initial concentrations at the start of the reaction step,
                                    as a Numpy array for ALL the species, in their array index order.
                                    If not provided, self.system is used instead
        :param explain_variable_steps:  [OPTIONAL] If provided, it must be a pair of numbers of the form (t_start, t_end);
                                            a brief explanation will be printed about how the variable step sizes were chosen,
                                            when the System time inside that range
        :param step_counter:    [OPTIONAL] Integer currently only used for diagnostics

        :return:                The triplet:
                                    1) increment vector for the concentrations of ALL the species,
                                        in their array index order, as a Numpy array
                                        EXAMPLE (for a single-reaction reactant and product with a 1:3 stoichiometry):
                                            array([7. , -21.])
                                    2) time step size actually taken - which might be smaller than the requested one
                                        because of reducing the step to avoid negative-concentration errors
                                    3) recommended_next_step : a suggestions to the calling module
                                       about the next step to best take
        """
        # TODO: no longer pass conc_array .  Use the object variable self.system instead
        #       Determine whether 1 or multiple UC objects are to be used by Bio1D, etc.

        if conc_array is not None:
            self.system = conc_array    # For historical reasons, as a convenience to Bio1D, etc.
                                        # TODO: maybe it ought to be kept separate in separate instances of the UC object
            if self.uniform_compartment:
                self.uniform_compartment.system = self.system   # Re-align

        # Validate arguments
        assert self.system is not None, "ReactionSimulator.reaction_step_common_variable_step(): " \
                                        "the concentration values of the various chemicals must be set first"


        if explain_variable_steps:
            assert isinstance(explain_variable_steps, (list, tuple))  and  (len(explain_variable_steps) == 2), \
                "reaction_step_common_variable_step(): the argument `explain_variable_steps`, " \
                "if provided, must be a pair of numbers [t_start, t_end]"


        #print(f"************ At SYSTEM TIME: {self.system_time:,.4g}, calling reaction_step_common() with:")
        #print(f"             delta_time={delta_time}, system={self.system}, ")


        recommended_next_step = delta_time     # Baseline value; no reason yet to suggest a change in step size


        delta_concentrations = None

        SMALLEST_VALUE_TO_TRY = delta_time / 2000.       # Used to prevent infinite loops

        normal_exit = False

        while delta_time > SMALLEST_VALUE_TO_TRY:       # TODO: consider moving the inside of the WHILE loop into a separate function
            try:    # We want to catch Exceptions that can arise from excessively large time steps
                    #   that lead to negative concentrations or violation of user-set thresholds ("HARD" or "SOFT" aborts)
                (delta_concentrations, recommended_next_step) = \
                        self.attempt_reaction_step(delta_time=delta_time, variable_steps=True,
                                                   explain_variable_steps=explain_variable_steps, step_counter=step_counter)
                normal_exit = True

                break       # IMPORTANT: this is needed because, in the absence of errors, we need to go thru the WHILE loop only once!


            # CATCH any 'ExcessiveTimeStepHard' exception raised in the loop  (i.e. a HARD ABORT)
            except ExcessiveTimeStepHard as ex:
                # Single reactions steps can fail with this error condition if the attempted time step was too large,
                # under the following scenarios:
                #       1. negative concentrations from any one reaction - caught by  validate_increment()
                #       2. negative concentration from the combined effect of multiple reactions - caught in this function
                #print("*** CAUGHT a HARD ABORT")
                self.reaction_step_diagnostics.number_neg_concs += 1
                if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                    explanation = ex.details.get("message")
                    explanation += f"\n      -> will backtrack, and re-do step with a SMALLER delta time, " \
                                    f"multiplied by {self.adaptive_steps.step_factors['error']} " \
                                    f"(set to {delta_time * self.adaptive_steps.step_factors['error']:.5g}) " \
                                    f"\n      [Step started at t={self.system_time:.5g}, and will rewind there]"
                    print(explanation)


                delta_time *= self.adaptive_steps.step_factors["error"]       # Reduce the excessive time step by a pre-set factor
                recommended_next_step = delta_time
                # At this point, the loop will generally try the simulation again, with a smaller step (a revised delta_time)


            # CATCH any 'ExcessiveTimeStepSoft' exception raised in the loop  (i.e. a SOFT ABORT)
            except ExcessiveTimeStepSoft as ex:
                # Single reactions steps, in the variable step scenario,
                # can fail with this error condition if the attempted time step was too large,
                # under the following scenario:
                #       * excessive norm(s) measures in the overall step - caught in this function
                #print("*** CAUGHT a soft ABORT")
                self.reaction_step_diagnostics.number_soft_aborts += 1
                if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                    msg = ex.details.get("message")
                    print(f"       {msg}")
                delta_time *= self.adaptive_steps.step_factors["abort"]       # Reduce the excessive time step by a pre-set factor
                recommended_next_step = delta_time
                # At this point, the loop will generally try the simulation again, with a smaller step (a revised delta_time)

        # END while


        if not normal_exit:         # i.e., if no reaction simulation took place in the WHILE loop, above
            raise Exception(f"reaction_step_common_variable_step(): unable to complete the reaction step.  "
                            f"In spite of numerous automated reductions of the time step, "
                            f"it continues to lead to concentration changes that are considered excessive; "
                            f"try reducing the original time step, and/or increasing the 'abort' thresholds with set_thresholds(). "
                            f"Current values: {self.adaptive_steps.thresholds}")


        # If we get thus far, it's the normal exit of the reaction step

        return  (delta_concentrations, delta_time, recommended_next_step)     # TODO: consider returning tentative_updated_system , since we already computed it



    def attempt_reaction_step(self, delta_time :float, variable_steps :bool, explain_variable_steps=None, step_counter=None) -> (np.array, float):
        """
        Attempt to perform a single simulation step for ALL reactions - and then raise an Exception if it needs to be aborted,
        based on various criteria.
        If variable_steps is True, determine a new value for the "recommended next step"

        :param delta_time:              The requested time duration of the reaction step
        :param variable_steps:          If True, the step sizes will get automatically adjusted with an adaptive algorithm
        :param explain_variable_steps:  [OPTIONAL] If provided, it must be a pair of numbers of the form [t_start, t_end];
                                            a brief explanation will be printed about how the variable step sizes were chosen,
                                            when the System time inside that range.
                                            Only applicable if variable_steps is True
        :param step_counter:            A pair with a time range to show in the explanations about the variable step sizes;
                                            only applicable if explain_variable_steps is True

        :return:                The pair (delta_concentrations, recommended_next_step)
                                    If variable_steps is False, recommended_next_step is trivially the same as `delta_time`
        """
        # TODO: explain_variable_steps should be a boolean

        # *****  CORE OPERATION  *****
        delta_concentrations = self.single_step_all_rxns(delta_time=delta_time, rxn_list=None) # ALL reactions


        if self.diagnostics_enabled:
            self.diagnostic_data_snapshot = {}  # Reset


        recommended_next_step = delta_time      # (Only applicable if variable_steps is True)
                                                # Baseline value; no reason yet to suggest a change in step size

        if variable_steps:
            decision_data = self.adaptive_steps.adjust_timestep(n_chems=self.number_of_system_species(),
                                                                indexes_of_active_chemicals= self.indexes_of_active_chemicals(),
                                                                delta_conc=delta_concentrations,
                                                                baseline_conc=self.system, prev_conc=self.previous_system)
            #print("decision_data: ", decision_data)
            # EXAMPLE: {'action': 'stay', 'step_factor': 1, 'norms': {'norm_A': 0.98, 'norm_B': 0.14}, 'applicable_norms': 'ALL'}

            saved_decision_data = {k : v  for k, v in decision_data.items() if k != "norms"}   # Drop the "norms" key
            self.reaction_step_diagnostics.decision_data = saved_decision_data
            self.reaction_step_diagnostics.norms = decision_data["norms"]

            step_factor = decision_data['step_factor']
            action = decision_data['action']
            applicable_norms = decision_data['applicable_norms']

            self._explain_variable_timestep_prelude(explain_variable_steps=explain_variable_steps,
                                                    step_counter=step_counter,
                                                    delta_concentrations=delta_concentrations)


            # Abort the current step if some rate of change is deemed excessive.
            # TODO: maybe ALWAYS check this, regardless of variable-steps option
            if action == "abort":       # NOTE: this is a "strategic" abort, not a hard one from error

                # TODO: maybe this call could wait until the interception of ExcessiveTimeStepSoft?
                self._save_soft_abort_diagnostics(delta_concentrations=delta_concentrations)


                exception_data = {
                    "message": f"* INFO: the tentative time step ({delta_time:.5g}) "
                        f"leads to a value of {applicable_norms} > its ABORT threshold:\n"
                        f"       -> will backtrack, and re-do step with a SMALLER Δt, x{step_factor:.5g} (now set to {delta_time * step_factor:.5g}) "
                        f"[Step started at t={self.system_time:.5g}, and will rewind there]",
                    "delta_time": delta_time
                }
                raise ExcessiveTimeStepSoft(exception_data)    # ABORT THE CURRENT STEP


            # Put together a recommendation to the higher-level functions, about the next best step size
            recommended_next_step = delta_time * step_factor

            # Append data to self.diagnostic_data_snapshot
            self._gather_diagnostics_for_var_step(step_factor=step_factor, delta_time=delta_time, action=action)
            #print(self.diagnostic_data_snapshot)

            self._explain_variable_timestep_upon_success(step_factor=step_factor, delta_time=delta_time,
                                                         explain_variable_steps=explain_variable_steps,
                                                         recommended_next_step=recommended_next_step, applicable_norms=applicable_norms)

        # END if variable_steps


        self._save_diagnostics_end_of_step(delta_concentrations=delta_concentrations)    # "diagnostic_decisions_data"


        # Check whether the COMBINED delta_concentrations will make any conc negative;
        # if so, raise an "ExcessiveTimeStepHard" exception (a custom exception)
        tentative_updated_system = self.system + delta_concentrations
        if min(tentative_updated_system) < 0:
            # TODO: maybe this call could wait until the interception of ExcessiveTimeStepHard ?
            self._save_hard_abort_diagnostics_all_rxns(explain_variable_steps=explain_variable_steps, delta_time=delta_time)

            neg_indices = np.where(tentative_updated_system < 0)[0]
            first_neg_index = neg_indices[0]
            chem_name = self.uniform_compartment.species_index_map.species_at(int(first_neg_index))  # The int() is to convert the NumPy integer type
            exception_data = {
                "message": f"      The tentative time step ({delta_time:.6g}) "
                                        f"would lead to a NEGATIVE concentration "
                                        f"\n      in one or more of the chemicals (for instance `{chem_name}`, of index {first_neg_index}), from the combined reactions."
                                        f"\n      Baseline concentration values: {self.system} at system time {self.system_time:.5g}; requested changes (NOT carried out): {delta_concentrations}",
                "function": "attempt_reaction_step",
                "delta_time": delta_time
            }
            raise ExcessiveTimeStepHard(exception_data)


        return  (delta_concentrations, recommended_next_step)       # Maybe also return tentative_updated_system



    def _explain_variable_timestep_prelude(self, explain_variable_steps, step_counter,
                                           delta_concentrations) -> None:
        """
        If requested, print out an explanation for the user about the current variable time step,
        and the decision made

        :param explain_variable_steps:  [OPTIONAL] If provided, it must be a pair of numbers of the form [t_start, t_end];
                                            a brief explanation will be printed about how the variable step sizes were chosen,
                                            when the System time inside that range.
        :param step_counter:
        :param delta_concentrations:
        :return:                        None
        """
        if explain_variable_steps \
                and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
            step_factor = self.reaction_step_diagnostics.decision_data['step_factor']
            action      = self.reaction_step_diagnostics.decision_data['action']
            operation   = self.reaction_step_diagnostics.decision_data['operation']
            all_norms  = self.reaction_step_diagnostics.norms
            delta_time = self.reaction_step_diagnostics.delta_time

            if action == "abort":
                step_status = "aborted"
            else:
                step_status = "completed"

            print(f"\n(STEP {step_counter} {step_status}) SYSTEM TIME {self.system_time:.5g} : Examining Conc. changes "
                  f"due to tentative Δt={delta_time:.5g} ...")
            print("    Previous: ", self.previous_system)
            print("    Baseline: ", self.system)
            print("    Deltas:   ", delta_concentrations)

            if len(self.reaction_registry.active_chemicals) < self.number_of_system_species():
                print(f"    Restricting adaptive time step analysis to {len(self.reaction_registry.active_chemicals)} "
                f"species only: {self.reaction_registry.get_species_in_any_reaction()} , with indexes: {self.indexes_of_active_chemicals()}")

            self.adaptive_steps.display_overview(action=action, operation=operation, step_factor=step_factor, all_norms=all_norms)

            #print(f"    => Action: '{action.upper()}'  (with step size factor of {step_factor})")


    def _explain_variable_timestep_upon_success(self, step_factor, delta_time, explain_variable_steps, recommended_next_step, applicable_norms) -> None:
        """
        Print out, if requested, a 2-line wrap-up informational message
        at the end of a successful variable step of simulation

        :param step_factor:
        :param delta_time:              The requested time duration of the reaction step
        :param explain_variable_steps:  [OPTIONAL] If provided, it must be a pair of numbers of the form [t_start, t_end];
                                            a brief explanation will be printed about how the variable step sizes were chosen,
                                            when the System time inside that range.
        :param recommended_next_step:
        :param applicable_norms:
        :return:                        None
        """
        if explain_variable_steps \
                and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
            msg = "       "

            if step_factor > 1:         # "INCREASE
                msg +=  f"INFO: COMPLETED STEP NORMALLY and MADE INTERVAL LARGER, " \
                        f"multiplied by {step_factor} (set to {recommended_next_step:.5g}) at the next round, because all norms are low"
            elif step_factor < 1:       # "DECREASE"
                msg +=  f"INFO: COMPLETED STEP NORMALLY and MADE INTERVAL SMALLER, " \
                        f"multiplied by {step_factor} (set to {recommended_next_step:.5g}) at the next round, because {applicable_norms} is high"
            else:                       # "STAY THE COURSE"
                msg +=  f"INFO: COMPLETED STEP NORMALLY - we're inside the target range of all norms.  No change to step size."

            msg += f"\n    [The current step started at System Time: {self.system_time:.5g} , and will continue to {self.system_time + delta_time:.5g}]"
            print(msg)



    def _gather_diagnostics_for_var_step(self, step_factor, delta_time, action) -> None:
        """
        Gather diagnostic data, if enabled, about the variable step.
        Add this data to the dict self.diagnostic_data_snapshot

        :param step_factor:
        :param delta_time:
        :param action:
        :return:            None
        """
        if self.diagnostics_enabled:
            all_norms = self.reaction_step_diagnostics.norms
            # Populate the dict self.diagnostic_data_snapshot
            self.diagnostic_data_snapshot['norm_A'] = all_norms.get('norm_A')    # TODO: combine all norms in 1 step
            self.diagnostic_data_snapshot['norm_B'] = all_norms.get('norm_B')
            self.diagnostic_data_snapshot['norm_C'] = all_norms.get('norm_C')
            self.diagnostic_data_snapshot['norm_D'] = all_norms.get('norm_D')
            self.diagnostic_data_snapshot['action'] = f"OK ({action})"
            self.diagnostic_data_snapshot['step_factor'] = step_factor
            self.diagnostic_data_snapshot['time_step'] = delta_time


    def _save_diagnostics_end_of_step(self, delta_concentrations):
        """
        THIS DATA GETS SAVED INTO  "diagnostic_decisions_data".

        Save diagnostic data, if enabled, at the conclusion of the simulation step
        (which might be variable or fixed).

        Used to save the diagnostic concentration values, and concentration changes,
        indexed by the given System Time.
        Note: if an interval run is aborted, by our convention an entry is STILL created here

        :param delta_concentrations:
        :return:
        """
        # NOTE: self.diagnostic_data_snapshot is only applicable to VARIABLE steps
        if self.diagnostics_enabled:
            self.diagnostics.save_diagnostic_decisions_data(system_time=self.system_time,
                                                            delta_conc_arr=delta_concentrations,
                                                            data=self.diagnostic_data_snapshot)



    def single_step_all_rxns(self, delta_time: float, rxn_list=None) -> np.array:
        """
        TODO: this corresponds to UniformCompartment._reaction_elemental_step() , to eventually ditch

        Using the system concentration data,
        do the specified SINGLE TIME STEP
        for ONLY the requested reactions (by default all).

        Return the Numpy increment vector for ALL the species concentrations, in their index order
        (whether involved in these reactions or not)

        The object variable `self.reaction_step_diagnostics.system_rxn_rates` is also set.

        NOTES:  - the actual System Concentrations
                    and the System Time (stored in object variables) are NOT changed
                - if any of the concentrations go negative, an Exception is raised

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span.
        :param rxn_list:    OPTIONAL list of reactions (specified by their indices) to include in this simulation step ;
                                EXAMPLE: [1, 3, 7]
                                If None, do all the reactions

        :return:            The increment vector caused by all the specified reactions
                                for the concentrations of ALL the chemical species,
                                (whether involved in the reactions or not),
                                as a Numpy array for all the chemical species, in their index order
                            EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):
                                array([7. , -21.])
        """
        assert self.reaction_registry is not None, \
            "ReactionSimulator.single_step_all_rxns(): Cannot perform simulations because the 'ReactionRegistry' object isn't available;\n" \
            "    did you pass the argument `reaction_registry` at instantiation time, or set it later?"

        self.reaction_step_diagnostics = ReactionStep(delta_time=delta_time)     # RESET, ahead of this reaction step

        # The increment vector is cumulative for ALL the requested reactions.  Initialize it to all zeros
        number_species = self.uniform_compartment.species_index_map.number_of_system_species()
        increment_vector = np.zeros(number_species, dtype=float)       # One element per species


        if rxn_list is None:    # Meaning ALL reactions
            # A list of the reaction indices of all the active reactions
            rxn_list = self.reaction_registry.active_reaction_indices()


        # For each applicable reaction, find the needed adjustments ("deltas")
        #   to the concentrations of the reactants and products,
        #   based on the forward and reverse rates of the reaction
        for rxn_index in rxn_list:      # Consider each reaction in turn
            rxn = self.reaction_registry.get_reaction(rxn_index)
            rxn_rate = self.single_step_single_rxn(rxn=rxn, rxn_index=rxn_index, delta_time=delta_time,
                                                   increment_vector=increment_vector)
            # Compute and save up the rates ("velocities") of all the reactions we're looking into, as a dict;
            # the keys are the reaction indexes
            self.reaction_step_diagnostics.system_rxn_rates[rxn_index] = rxn_rate       # Save the value

        return increment_vector



    def single_step_single_rxn(self, increment_vector, delta_time,
                                rxn, rxn_index) -> float:
        """
        Update the Numpy array passed as the argument `increment_vector`.

        If the object variable `self.diagnostics_enabled` is True,
        then a variety of routine reaction diagnostics are logged;
        and addition diagnostics are logged in case of Exceptions
        resulting from the proposed concentration increments

        :param increment_vector:Numpy array that gets updated in place
        :param delta_time:      The duration of the single time step to take
        :param rxn:             Object of type "SimulationReaction"

        (The remaining arguments are just for diagnostics and error printing
        :param rxn_index:       The index (0-based) of the above reaction in self.reaction_registry  (ONLY USED for diagnostics and error printing)

        :return:                The initial rate of the reaction
        """
        conc_init = self._fetch_concs_for_rnx(rxn=rxn)
        # For the species in this rxn only.  EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}

        increment_dict_single_rxn, rxn_rate = self.dispatcher_single_rxn(rxn=rxn, conc_init=conc_init,
                                                                         delta_time=delta_time)
        # EXAMPLE of increment_dict_single_rxn: {"B": -1.3, "F": 2.9, "D": -1.6}


        for (species_id, delta_conc) in increment_dict_single_rxn.items():
            species_index = self.uniform_compartment.species_index_map.index_of(species_id)
            # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
            # for any of the proposed concentration changes for this reaction.
            # Note: it's not enough to detect conc going negative from combined changes from multiple reactions!
            #       Further testing done upstream
            self._validate_increment(delta_conc=delta_conc, baseline_conc=self.system[species_index],
                                     rxn=rxn, rxn_index=rxn_index, species_id=species_id, rxn_rate=rxn_rate,
                                     delta_time=delta_time)

            # Accumulate the increment vector from the species in this reaction
            increment_vector[species_index] += delta_conc  # Accumulate  all the increments from this reaction


        if self.diagnostics_enabled:
            self.diagnostics.save_rxn_data(rxn_index=rxn_index,
                                           system_time=self.system_time, time_step=delta_time,
                                           increment_dict_single_rxn=increment_dict_single_rxn,
                                           rate=rxn_rate)   # Save to "diagnostic_rxn_data"

        return rxn_rate



    def _validate_increment(self, delta_conc :float, baseline_conc :float,
                            rxn, rxn_index :int, species_id: str, rxn_rate, delta_time) -> None:
        """
        Examine the single requested concentration change `delta_conc`
        (typically, as computed by an ODE solver),
        relative to the baseline (pre-reaction) value `baseline_conc`,
        for the given SINGLE species and SINGLE reaction.

        If the requested concentration change would render the concentration negative,
        save diagnostic data if diagnostics are enabled, and then
        raise an Exception of the custom type "ExcessiveTimeStepHard"

        :param delta_conc:      The change in concentration that we're considering
                                    for the specified species, in the given reaction
        :param baseline_conc:   The initial concentration value for that species

        [The remaining arguments are ONLY USED for diagnostics and error printing]
        :param rxn_index:       The index (0-based) to identify the reaction of interest (ONLY USED for diagnostics and error message)
        :param species_id:      The id of the species under consideration (ONLY USED for diagnostics and error message)
        :param delta_time:      The time duration of the reaction step (ONLY USED for diagnostics and error message)

        :return:                None.  An Exception is raised if a negative new concentration would result
                                    from the requested concentration change
        """
        if (baseline_conc + delta_conc) < 0:
            # If the requested concentration change would lead to a negative concentration

            # A type of HARD ABORT is detected (a single reaction that, by itself, would lead to a negative concentration;
            #   while it's possible that other coupled reactions might counterbalance this - nonetheless,
            #   it's taken as a sign of excessive step size)

            # After having saved the appropriate diagnostic data, raise the custom Exception "ExcessiveTimeStepHard"
            exception_data = {
                "message": f"      The tentative time step ({delta_time:.6g}) "
                           f"would lead to a NEGATIVE concentration in the species `{species_id}` "
                           f"from the reaction `{rxn.describe(concise=True)}` (rxn # {rxn_index})\n"
                           f"      Baseline concentration value of `{species_id}` : {baseline_conc:.6g} at system time {self.system_time:.5g}; requested change (NOT carried out): {delta_conc:.6g}",
                "function": "_validate_increment",
                "delta_time": delta_time,
                #"action": "ABORT",
                "caption": f"aborted: neg. conc. in `{species_id}` from rxn # {rxn_index}",
                "system_time":  self.system_time,
                #"increment_dict_single_rxn": None,
                #"aborted": True,
                "rate": rxn_rate,
                "rxn_index": rxn_index,
                "species_id": species_id
            }

            #if self.adaptive_steps:     # Add more diagnostics data if available
            #    exception_data["step_factor"] = self.adaptive_steps.step_factors.get('error')


            if self.diagnostics_enabled:
                self._save_hard_abort_diagnostics_single_rxn(exception_data)

            # This Exception will get caught by either reaction_step_common_variable_step()
            # or by reaction_step_common_fixed_step(), leading to different remediation
            raise ExcessiveTimeStepHard(exception_data)



    def _fetch_concs_for_rnx(self, rxn):
        """
        Extract, out of the Numpy array of the system concentrations,
        just the concentrations of relevance for the specified reaction

        :param rxn:         An object of type "SimulationReaction"
        :return:            A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        """
        # Get the SET of the id's of ALL the species appearing in this reaction
        species_ids = rxn.stoichiometry.get_all_species_ids()   # EXAMPLE: {"B", "F", "D"}

        return self.uniform_compartment.get_conc_dict(chem_labels=species_ids)

        

    def dispatcher_single_rxn(self, rxn, conc_init :dict, delta_time :float):
        """
        Dispatch to the appropriate ODE solver for this single reaction,
        based on the object variable self.method

        If the requested method is "analytic", but no analytic solver is available
        for this reaction type, then fall back to "forward_euler" for this reaction

        :param rxn:         Object of type "SimulationReaction"
        :param conc_init:
        :param delta_time:  The duration of the single time step to take
        :return:            EXAMPLE: (  {"B": -1.3, "F": 2.9, "D": -1.6} ,
                                        35.
                                     )
        """
        method = self.method

        if method == "analytic" and not rxn.analytic_solution_family:
            method = "forward_euler"    # Default when no analytic solver is available

        if self.method == "forward_euler":
            increment_dict_single_rxn, rxn_rate = ReactionSimulator.forward_euler_single_rxn(rxn=rxn, conc_init=conc_init, delta_time=delta_time)
        elif self.method == "heun":
            increment_dict_single_rxn, rxn_rate = ReactionSimulator.heun_single_rxn(rxn=rxn, conc_init=conc_init, delta_time=delta_time)
        elif method == "analytic":
            rxn_rate = rxn.model.rate(conc_dict = conc_init)    # Rate at start of time step
            increment_dict_single_rxn = ReactionSimulator.analytic_solver_single_rxn(rxn=rxn, conc_init=conc_init, delta_time=delta_time)
        else:
            raise Exception(f"dispatcher_single_rxn(): Unknown reaction-solver method: '{self.method}'")


        return (increment_dict_single_rxn, rxn_rate)



    @staticmethod
    def analytic_solver_single_rxn(rxn, conc_init, delta_time) -> dict:
        """
        Unpack the reaction's parameters, and dispatch to the appropriate analytic solver for its type

        :param rxn:         Object of type "SimulationReaction"
        :param conc_init:   A dictionary of initial concentrations for all the species involved
        :param delta_time:  The duration of the single time step to take
        :return:            A dictionary of concentration increments for all the species involved
        """
        reactants = rxn.stoichiometry.get_reactant_list()     # A list of pairs of the form (stoichiometry coefficient, species id))
        products = rxn.stoichiometry.get_product_list()       # A list of pairs of the form (stoichiometry coefficient, species id))

        if rxn.analytic_solution_family == "ONE_TO_ONE":
            r = reactants[0][1]           # EXAMPLE: "R"
            p = products[0][1]            # EXAMPLE: "P"

            R0 = conc_init[r]
            P0 = conc_init[p]
            # Compute the respective increments of R0 and P0
            if rxn.model.reversible:
                delta_p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=rxn.model.kF, kR=rxn.model.kR,
                                                                                       A0=R0, P0=P0, t=delta_time, incremental=True)
            else:
                delta_p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=rxn.model.kF,
                                                                                         A0=R0, P0=P0, t=delta_time, incremental=True)

            # Work out the stoichiometry for all the species
            increment_dict_single_rxn = {r: -delta_p, p: delta_p}
            return increment_dict_single_rxn
        else:       # TODO: implement "TWO_TO_ONE" and "ONE_TO_TWO"
            raise Exception(f"analytic_solver_single_rxn(): no exact analytic solution "
                            f"is currently implemented for reactions of type {rxn.analytic_solution_family}")




    @staticmethod
    def forward_euler_single_rxn(rxn, conc_init :dict, delta_time :float) -> tuple[dict, float]:
        """
        Simulate the given reaction, over the specified single time step,
        using the "forward Euler" method

        :param rxn:         The "SimulationReaction" object for this reaction
        :param conc_init:   A dictionary of initial concentrations for all the species involved
        :param delta_time:  The duration of the single time step to take
        :return:            The pair (increment_dict_single_rxn, rxn_rate)
        """
        # Compute the reaction rate ("velocity"), at the current system concentrations, for this reaction
        rate_initial = rxn.model.rate(conc_dict = conc_init)    # Rate at start of time step

        delta_rxn = rate_initial * delta_time      # forward reaction - reverse reaction

        #Note: conc_final = conc_init + v_i * rate_initial * delta_time
        #      We return conc_final - conc_init , which we call delta_conc

        stoich = rxn.stoichiometry
        all_species = stoich.get_all_species_ids(exclude_catalysts=True)
        #print(delta_rxn)
        #print(rate_initial)

        # Determine the concentration adjustments as a result of this reaction step,
        #       for this individual reaction being considered
        # Note: the SIGNED stoichiometry coefficient ensure that,
        #       if delta_rxn is positive,
        #       the reactants decrease in concentration, and the products increase
        delta_conc = {species: (stoich.vector[species] * delta_rxn)
                                    for species in all_species}

        # TODO: maybe raise an "ExcessiveTimeStepHard" Exception, if any of the delta_conc
        #       components would make its final concentration negative

        return (delta_conc, rate_initial)




    @staticmethod
    def heun_single_rxn(rxn, conc_init :dict, delta_time :float) -> tuple[dict, float]:
        """
        Simulate the given reaction, over the specified single time step,
        using the "Heun" method (aka "explicit trapezoidal rule")

        :param rxn:         The "SimulationReaction" object for this reaction
        :param conc_init:   A dictionary of initial concentrations for all the species involved
        :param delta_time:  The duration of the single time step to take
        :return:            The pair (increment_dict_single_rxn, rxn_rate)
        """
        # The first part is just like the forward Euler method: we'll compute the
        # final concentrations after the given single time step

        # Compute the reaction rate ("velocity"), at the current system concentrations, for this reaction
        rate_initial = rxn.model.rate(conc_dict = conc_init)    # Rate at start of time step

        delta_rxn_prelim = rate_initial * delta_time      # forward reaction - reverse reaction (early pass)

        stoich = rxn.stoichiometry
        all_species = stoich.get_all_species_ids(exclude_catalysts=True)

        final_conc = {species: (stoich.vector[species] * delta_rxn_prelim) + conc_init[species]
                                    for species in all_species}
        # So far, same as the "Forward Euler" method
        #print("Euler final_conc: ", final_conc)

        min_conc = min(final_conc.values())
        if min_conc < 0:
            exception_data = {
                "message": f"heun_single_rxn(): excessive time step ({delta_time:.6g}), "
                           f"leading to negative concentrations",
                "function": "heun_single_rxn",
                "delta_time": delta_time,
                "rate": rate_initial
            }
            raise ExcessiveTimeStepHard(exception_data)

        # Now compute the rate at END of time step, UNDER THE ASSUMPTION that the reaction
        # proceeded as predicted by the "Forward Euler" approximation
        rate_final = rxn.model.rate(conc_dict = final_conc)
        #print(f"final rate: {rate_final}")

        # We want to intercept scenarios where the reaction rate changes sign,
        # and the Heun update would reverse the reaction's direction or result in zero changes
        if (np.sign(rate_final) != np.sign(rate_initial)) \
            and (np.abs(rate_final) >= np.abs(rate_initial)):
                exception_data = {
                    "message":  f"heun_single_rxn(): excessive time step ({delta_time:.6g}), "
                                f"leading to final rate of {rate_final:.5g} which, when averaged with the "
                                f"initial rate of {rate_initial} would flip the reaction's direction",
                    "function": "heun_single_rxn",
                    "delta_time": delta_time,
                    "rate": rate_initial,
                    "rate_final": rate_final
                }
                raise ExcessiveTimeStepHard(exception_data)

        # Finally, average the two rate, and use the average the advance the reaction
        rate_heun = (rate_initial + rate_final) / 2
        #print(f"rate_heun: {rate_heun}")

        delta_rxn = rate_heun * delta_time      # forward reaction - reverse reaction (corrected)

        delta_conc = {species: (stoich.vector[species] * delta_rxn)
                                    for species in all_species}

        return (delta_conc, rate_initial)
    



    ############  DIAGNOSTICS  ############

    def enable_diagnostics(self):
        """
        Turn on the diagnostics mode.

        CAUTION: if not done early enough in the simulation, some desired diagnostic data may be missing.
                 You might consider using the argument  enable_diagnostics=True
                 when first instantiating UniformCompartment.

        :return: None
        """
        #TODO: maybe automatically capture a snapshot of the current data
        if self.diagnostics_enabled:
            print("*** INFO: enable_diagnostics() - diagnostics were ALREADY enabled")

        self.diagnostics_enabled = True
        if not self.diagnostics:
            self.diagnostics = Diagnostics(reactions=self.reaction_registry,
                                           species_to_index=self.uniform_compartment.species_index_map.species_to_index)



    def pause_diagnostics(self):
        """
        Turn off the overall diagnostics mode; existing diagnostics data, if any, is left untouched

        :return:    None
        """
        if not self.diagnostics_enabled:
            print("*** INFO: pause_diagnostics() - diagnostics were ALREADY off")

        self.diagnostics_enabled = False



    def get_diagnostics(self):
        """

        :return:    Object of type life123.diagnostics.Diagnostics
        """
        assert self.diagnostics is not None, \
            "get_diagnostics(): no diagnostics data is available.  " \
            "Did you call enable_diagnostics() prior to running the simulation?"

        return self.diagnostics





####################################################################################################

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

        #self.uc = uc    # Object of type "UniformCompartment"  TODO: NOT USED



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



    def display_overview(self, action :str, operation :str, step_factor :dict, all_norms :dict) -> None:
        """

        :param action:
        :param all_norms:
        :return:            None
        """
        # Round off all the norm values before printing them
        print("    Norms:     {",
              ", ".join(f"'{key}': {value:.5g}" for key, value in all_norms.items()),
              "}")      # EXAMPLE:     Norms:    { 'norm_A': 0.98, 'norm_B': 0.14 }

        print("    Thresholds:    ")
        self.display_value_against_thresholds(all_norms)



        if action != "stay":    # The step is trivially 1 when the action is "stay"
            print("    Step Factors:    ", self.step_factors)

        if action == 'stay':
            print(f"    => Action: '{action.upper()}'  (with step size factor of {step_factor})")
        else:
            print(f"    => Action: '{action.upper()}'  ('{operation}' with step size factor of {step_factor})")



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
        :param indexes_of_active_chemicals: The ordered list (numerically sorted) of the INDEX numbers of all the species
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
                                                            either "low", "stay", "high" or "abort"     TODO: possibly rename "action" to "determination"
                                    "operation"        - String: either 'upshift', 'stay', 'downshift', 'abort'
                                    "step_factor"      - A factor by which to multiply the time step at the next iteration round;
                                                            if no change is deemed necessary, 1
                                    "norms"            - A dict of all the computed norm name/values (any of the norms, except norm_A,
                                                            may be missing)
                                    "applicable_norms" - The name of the norm(s) that triggered the decision; if all norms were involved,
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
                return {"action": "abort", "operation": "abort", "step_factor": self.step_factors["abort"], "norms": all_norms, "applicable_norms": [norm_name]}

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
            return {"action": "high", "operation": "downshift", "step_factor": self.step_factors["downshift"], "norms": all_norms, "applicable_norms": high_seen_at}


        if all_small:
            for i in self.norm_usage:
                # All the norms were used
                self.increase_norm_count(i)
                #self.norm_usage[i] += 1

            return {"action": "low", "operation": "upshift", "step_factor": self.step_factors["upshift"], "norms": all_norms, "applicable_norms": "ALL"}


        # If we get thus far, none of the rules were found applicable
        return {"action": "stay", "operation": "stay", "step_factor": 1, "norms": all_norms, "applicable_norms": "ALL"}



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
