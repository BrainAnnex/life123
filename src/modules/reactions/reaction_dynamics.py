import math
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from typing import Union
from src.modules.movies.movies import MovieTabular
from src.modules.numerical.numerical import Numerical as num



class ExcessiveTimeStepHard(Exception):
    """
    Used to raise Exceptions arising from excessively large time steps
    (that lead to negative concentration values, i.e. "HARD" errors)
    """
    pass

class ExcessiveTimeStepSoft(Exception):
    """
    Used to raise Exceptions arising from excessively large time steps
    (that lead to norms regarded as excessive because of user-specified values, i.e. "SOFT" errors)
    """
    pass

#############################################################################################

class ReactionDynamics:
    """
    Used to simulate the dynamics of reactions (in a single compartment.)
    In the context of Life123, this may be thought of as a "zero-dimensional system"
    """

    def __init__(self, chem_data):
        """
        :param chem_data:   Object of type "ChemData" (with data about the chemicals and their reactions)
                                    It's acceptable to pass None,
                                    and take care of it later (though probably a bad idea!)
                                    TODO: maybe offer an option to let the constructor instantiate that object?
        """
        self.chem_data = chem_data  # Object of type "ChemData" (with data about the chemicals and their reactions)

        self.system_time = 0.   # Global time of the system, from initialization on

        self.system = None  # Concentration data in the single compartment we're simulating, for all the (bulk) chemicals
                            # A Numpy array of the concentrations of all the chemical species, in their index order
                            # 1-dimensional NumPy array of floats, whose size is the number of chemical species.
                            # Each entry is the concentration of the species with that index (in the "ChemData" object)
                            # Note that this is the counterpart - with 1 less dimension - of the array by the same name
                            #       in the class BioSim1D

        self.macro_system = None    # The counterpart of the system data for macro-molecules, if present
                                    # Binding fractions of the applicable transcription factors, for all the macro-molecules
                                    # EXAMPLE:   {"M1": {"A": 0.2, "F": 0.93, "P": 0.},
                                    #             "M2": {"C": 0.5, "F": 0.1}}
                                    # For background, see: https://www.annualreviews.org/doi/10.1146/annurev-cellbio-100617-062719

        self.history = MovieTabular()   # To store user-selected snapshots of (some of) the chemical concentrations,
                                        #   whenever requested by the user.


        # ***  FOR AUTOMATED ADAPTIVE STEP SIZES  ***
        # Note: the "aborts" below are "elective" aborts - not aborts from hard errors (further below)
        self.thresholds = [{"norm": "norm_A", "low": 0.5, "high": 0.8, "abort": 1.44},
                           {"norm": "norm_B", "low": 0.08, "high": 0.5, "abort": 1.5}]
        self.step_factors = {"upshift": 1.5, "downshift": 0.5, "abort": 0.5}
                            # TODO: consider more conservative defaults upshift=1.2, downshift=0.5, abort=0.4

        self.error_abort_step_factor = 0.5      # MUST BE < 1.  Factor by which to multiply the time step
                                                #   in case of negative-concentration error from excessive step size
                                                #   NOTE: this is from ERROR aborts, not to be confused with high-threshold aborts
        # TODO: consider more conservative default 0.25


        self.reaction_speeds = {}       # TODO: not in active use; possibly obsolete
        # A dictionary, where the keys are reaction indices,
        #   and the values are either "S" (Slow) or "F" (Fast)
        #   EXAMPLE : { 1: "F", 4: "S", 5: "F" , 8: "S" }
        #   Any reaction with a missing entry is regarded as "F" (Fast)


        # ***  FOR DIAGNOSTICS  ***     TODO: maybe turn all diagnostic data/methods into a separate object

        self.verbose_list = []          # A list of integers or strings with the codes of the desired verbose checkpoints
                                        #   EXAMPLE: [1, "my_ad_hoc_tag"] to invoke sections of code marked as 1 or 3
                                        #   Those sections will have entry points such as:  if "my_ad_hoc_tag" in self.verbose_list


        self.diagnostics = False        # Overall flag about whether using diagnostics

        self.diagnostic_rxn_data = {}   # "Diagnostic reaction data", PER REACTION: a dict with as many entries as reactions.
                                        #   The keys are the reaction indices; the values are objects of type "MovieTabular",
                                        #   which contain Pandas dataframes with the following columns
                                        #   (referring to 1 reaction):
                                        #           'START TIME' , 'Delta A' , 'Delta B' ...'time_step' , 'caption'
                                        #
                                        #   Notes:  - entries are always added, even if an interval run is aborted
                                        #           - the various 'Delta concentrations' are for ALL the chemicals (in the reaction or not),
                                        #                   over the time interval that *STARTS* at the value in the "TIME" column

        self.diagnostic_conc_data = MovieTabular(parameter_name="TIME")
                                        # An expanded version of the normal System History.
                                        #   Columns of the dataframes:
                                        #       'TIME' 	'A' 'B' ...  'caption'
                                        #
                                        #   Note: if an interval run is aborted, NO entry is created here
                                        #         (this approach DIFFERS from that of other diagnostic data)

        self.diagnostic_decisions_data = MovieTabular(parameter_name="START_TIME")
                                        #   Columns of the dataframes:
                                        #       'START_TIME' 	'Delta A' 'Delta B' ...
                                        #               [plus, if applicable, other fields such as 'action', 'norm_A', 'norm_B', 'step_factors']
                                        #
                                        #   Note: entries are always added, even if an interval run is aborted



    #####################################################################################################

    '''                             ~   TO SET/READ DATA   ~                                          '''

    def ________TO_SET_AND_READ_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def set_conc(self, conc: Union[list, tuple, dict], snapshot=True) -> None:     # TODO: maybe rename set_all_conc()
        """
        Set the concentrations of ALL the chemicals at once   TODO: maybe a dict indicates a selection of chems, while a list, etc means "ALL"

        :param conc:        A list or tuple of concentration values for ALL the chemicals, in their index order;
                                alternatively, a dict indexed by the chemical names, again for ALL the chemicals
                                EXAMPLE of the latter:  {"A": 12.4, "B": 0.23, "C": 2.6}
                                                        (assuming that "A", "B", "C" are ALL the chemicals)

                                Note: any previous values will get over-written)
                                TODO: [maybe the dict version should be the responsibility of set_chem_conc() instead;] OR
                                      allow the dict to be a subset of all chemicals
                                TODO: also allow a Numpy array; make sure to do a copy() to it!
                                TODO: pytest for the dict option
        :param snapshot:    (OPTIONAL) boolean: if True, add to the history
                                a snapshot of this state being set.  Default: True
        :return:            None
        """
        # TODO: more validations, incl. of individual values being wrong data type

        assert len(conc) == self.chem_data.number_of_chemicals(), \
            f"ReactionDynamics.set_conc(): the number of concentration values passed in the argument 'conc' ({len(conc)}) " \
            f"must match the number of declared chemicals ({self.chem_data.number_of_chemicals()})"

        if type(conc) == dict:
            conc_list = []
            for species_index in range(self.chem_data.number_of_chemicals()):
                name = self.chem_data.get_name(species_index)
                if name is None:
                    raise Exception(f"ReactionDynamics.set_conc(): to use a dictionary for the `conc` arguments, "
                                    f"all the chemicals must first be given names - to be used as keys for the dictionary "
                                    f"(Chemical in index position {species_index} lacks a name")
                conc_list.append(conc[name])
            # END for
            conc = conc_list


        assert min(conc) >= 0, \
            f"ReactionDynamics.set_conc(): values meant to be chemical concentrations cannot be negative " \
            f"(such as the passed value {min(conc)})"


        self.system = np.array(conc)

        self.set_rxn_speed_all_fast()   # Reset all the reaction speeds to "Fast"    TODO: obsolete

        if snapshot:
            self.add_snapshot(caption="Initial state")



    def set_chem_conc(self, conc, species_index=None, species_name=None, snapshot=True) -> None:
        # TODO Maybe rename set_single_conc()
        #      [OR call it set_conc(), and also handle dict's]
        """
        Set the concentrations of 1 chemical
        Note: if both species_index and species_name are provided, species_name is used     TODO: pytest this part

        :param conc:            A non-negative number with the desired concentration value for the above value.
                                    (Any previous value will get over-written)
        :param species_index:   (OPTIONAL) An integer that indexes the chemical of interest (numbering starts at 0)
        :param species_name:    (OPTIONAL) A name for the chemical of interest.
                                    If both species_index and species_name are provided, species_name is used
        :param snapshot:        (OPTIONAL) boolean: if True, add to the history
                                    a snapshot of this state being set.  Default: True
        :return:                None
        """
        # Validate the arguments
        assert conc >= 0, \
            f"ReactionDynamics.set_chem_conc(): chemical concentrations cannot be negative (value passed: {conc})"

        if species_name is not None:
            species_index = self.chem_data.get_index(species_name)
        elif species_index is not None:
            self.chem_data.assert_valid_species_index(species_index)
        else:
            raise Exception("ReactionDynamics.set_chem_conc(): at least one "
                            "of the arguments `species_index` or `species_name` must be provided")


        # TODO: if setting concentrations of a newly-added chemical, needs to first expand self.system
        self.system[species_index] = conc

        self.set_rxn_speed_all_fast()   # Reset all the reaction speeds to "Fast"      TODO: obsolete
                                        # TODO: this is overkill; ought to only reset the affected reactions -
                                        #       as returned by get_reactions_participating_in()

        if snapshot:
            self.add_snapshot(caption=f"Set concentration of `{self.chem_data.get_name(species_index)}`")



    def get_system_conc(self) -> np.array:
        """
        Retrieve the concentrations of ALL the chemicals as a Numpy array

        :return:        Either a Numpy array with the concentrations of ALL the chemicals, in their index order
                            EXAMPLE:  array([12.3, 4.56, 0.12])    The 0-th chemical has concentration 12.3, and so on...
        """
        return self.system



    def get_conc_dict(self, species=None, system_data=None) -> dict:
        """
        Retrieve the concentrations of the requested chemicals (by default all),
        as a dictionary

        :param species:     (OPTIONAL) list or tuple of names of the chemical species; by default, return all
        :param system_data: (OPTIONAL) a Numpy array of concentration values, in the same order as the
                                index of the chemical species; by default, use the SYSTEM DATA
                                (which is set and managed by various functions)
        :return:            A dictionary, indexed by chemical name, of the concentration values;
                                EXAMPLE: {"A": 1.2, "D": 4.67}
        """
        if system_data is None:
            system_data = self.system
        else:
            assert system_data.size == self.chem_data.number_of_chemicals(), \
                f"ReactionDynamics.get_conc_dict(): the argument `system_data` must be a 1-D Numpy array with as many entries " \
                f"as the declared number of chemicals ({self.chem_data.number_of_chemicals()})"


        if species is None:
            if system_data is None:
                return {}
            else:
                return {self.chem_data.get_name(index): system_data[index]
                        for index, conc in enumerate(system_data)}
        else:
            assert type(species) == list or  type(species) == tuple, \
                f"ReactionDynamics.get_conc_dict(): the argument `species` must be a list or tuple" \
                f" (it was of type {type(species)})"

            conc_dict = {}
            for name in species:
                species_index = self.chem_data.get_index(name)
                conc_dict[name] = system_data[species_index]

            return conc_dict




    '''
    Management of reactions
    '''

    def slow_rxns(self) -> [int]:
        """
        Return a list of all the reactions that are marked as "slow"
        :return:
        """
        return [k  for k, v in self.reaction_speeds.items()  if v == "S"]
        # Alternate approach:
        # return list(filter(lambda k:self.reaction_speeds[k] == "S", d ))


    def fast_rxns(self) -> [int]:
        """
        Return a list of all the reactions that are marked as "fast"
        :return:
        """
        return [i for i in range(self.chem_data.number_of_reactions())
                if i not in self.slow_rxns()]
        # Alternate way:
        # return list(set(range(self.reaction_data.number_of_reactions()).difference(self.slow_rxns()))



    def are_all_slow_rxns(self) -> bool:
        """
        Return True iff all the reactions are marked as "slow"
        :return:
        """
        return len(self.slow_rxns()) == self.chem_data.number_of_reactions()



    def get_rxn_speed(self, rxn_index: int) -> str:     # TODO: obsolete
        """
        For the requested reaction, get the string code that it was marked with
        to classify its speed.
        If the reaction has not been previously classified, regard it as "F" (Fast)

        :param rxn_index:   The index (0-based) to identify the reaction of interest
        :return:            A 1-letter string with the code "F" (for Fast) or "S" (Slow)
        """
        return self.reaction_speeds.get(rxn_index, "F")     # Any missing entry is regarded to be "F" (Fast)


    def set_rxn_speed(self, rxn_index: int, speed: str) -> None:     # TODO: obsolete
        """
        Set a code value that classifies the reaction speed to tag the given reaction to

        :param rxn_index:   The index (0-based) to identify the reaction of interest
        :param speed:       A 1-letter string with the code "F" (for Fast) or "S" (Slow)
        :return:            None
        """
        assert speed in ["S", "F"], "set_rxn_speed(): `speed` argument must be either 'S' or 'F'"
        self.chem_data.assert_valid_rxn_index(rxn_index)
        self.reaction_speeds[rxn_index] = speed


    def set_rxn_speed_all_fast(self) -> None:     # TODO: obsolete
        """
        Reset all the reaction speeds to "Fast"

        :return:    None
        """
        self.reaction_speeds = {}   # Now they'll all be assumed to be "Fast"
                                    #   (which is the default for missing values)



    def clear_reactions(self) -> None:
        """
        Get rid of all reactions; start again with "an empty slate" (but still with reference
        to the same data object about the chemicals)

        # TODO: maybe offer an option to clear just one reaction, or a list of them
        # TODO: provide support for "inactivating" reactions

        :return:    None
        """
        self.chem_data.clear_reactions_data()
        self.set_rxn_speed_all_fast()   # Reset all the reaction speeds to "Fast".      TODO: obsolete




    #############################################################################################
    #                                                                                           #
    #                                   TO VISUALIZE SYSTEM                                     #
    #                                                                                           #
    #############################################################################################
    def ________TO_VISUALIZE_SYSTEM________(DIVIDER):  # Used to get a better structure view in IDEs such asPycharm
        pass


    def describe_state(self) -> None:
        """
        A simple printout of the state of the system
        :return:        None
        """
        print(f"SYSTEM STATE at Time t = {self.system_time:,.8g}:")

        n_species = self.chem_data.number_of_chemicals()
        print(f"{n_species} species:")

        # Show a line of line of data for each chemical species in turn
        for species_index, name in enumerate(self.chem_data.get_all_names()):
            if name:    # If a name was provided, show it
                name = f" ({name})"
            else:
                name = ""

            print(f"  Species {species_index}{name}. Conc: {self.system[species_index]}")




    #############################################################################################
    #                                                                                           #
    #                                TO PERFORM THE REACTIONS                                   #
    #                                                                                           #
    #############################################################################################
    def ________TO_PERFORM_THE_REACTIONS________(DIVIDER):  # Used to get a better structure view in IDEs such asPycharm
        pass


    def specify_steps(self, total_duration=None, time_step=None, n_steps=None) -> (float, int):
        """
        If either the time_step or n_steps is not provided (but at least 1 of them must be present),
        determine the other one from total_duration

        Their desired relationship is: total_duration = time_step * n_steps

        :param total_duration:  Float with the overall time advance (i.e. time_step * n_steps)
        :param time_step:       Float with the size of each time step
        :param n_steps:         Integer with the desired number of steps
        :return:                The pair (time_step, n_steps)
        """
        assert (not total_duration or not time_step or not n_steps), \
            "ReactionDynamics.specify_steps(): cannot specify all 3 arguments: `total_duration`, `time_step`, `n_steps` (specify only any 2 of them)"

        assert (total_duration and time_step) or (total_duration and n_steps) or (time_step and n_steps), \
            "ReactionDynamics.specify_steps(): must provide exactly 2 arguments from:  `total_duration`, `time_step`, `n_steps`"

        if not time_step:
            time_step = total_duration / n_steps

        if not n_steps:
            n_steps = math.ceil(total_duration / time_step)

        return (time_step, n_steps)     # Note: could opt to also return total_duration if there's a need for it



    def set_thresholds(self, norm="norm_A", low=None, high=None, abort=None) -> None:
        """
        Over-ride default values for simulation parameters.
        The default None values can be used to eliminate some of the threshold rules, for the specified norm

        :param norm:
        :param low:
        :param high:
        :param abort:
        :return:            None
        """
        for i, t in enumerate(self.thresholds):
            if t.get("norm") == norm:
                if low is None:
                    del t["low"]
                else:
                    t["low"] = low

                if high is None:
                    del t["high"]
                else:
                    t["high"] = high

                if abort is None:
                    del t["abort"]
                else:
                    t["abort"] = abort

                if len(t) == 1:
                    del self.thresholds[i]      # Completely eliminate this un-used norm
                return

        # If we get here, it means that we're handling a norm not currently present in the list self.thresholds
        new_t = {}
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



    def set_step_factors(self, upshift, downshift, abort) -> None:
        """
        Over-ride default values for simulation parameters

        :param upshift:
        :param downshift:
        :param abort:
        :return:            None
        """
        if self.step_factors is None:
            self.step_factors = {}

        if abort is not None:
            self.step_factors["abort"] = abort
        if downshift is not None:
            self.step_factors["downshift"] = downshift
        if upshift is not None:
            self.step_factors["upshift"] = upshift


    def set_error_step_factor(self, value) -> None:
        """
        Over-ride the default value for the simulation parameter error_abort_step_factor

        :param value:
        :return:        None
        """
        assert value < 1, "set_error_step_factor(): the argument must strictly be < 1"
        assert value > 0, "set_error_step_factor(): the argument must be a non-zero positive number"

        self.error_abort_step_factor = value




    def single_compartment_react(self, reaction_duration=None, target_end_time=None,
                                 initial_step=None, n_steps=None,
                                 snapshots=None, silent=False,
                                 variable_steps=False, explain_variable_steps=False) -> None:
        """
        Perform ALL the reactions in the single compartment -
        based on the INITIAL concentrations,
        which are used as the basis for all the reactions.

        Update the system state and the system time accordingly
        (object attributes self.system and self.system_time)

        :param reaction_duration:  The overall time advance for the reactions (it might be exceeded in case of variable steps)
        :param target_end_time: The final time at which to stop the reaction
                                    If both target_end_time and reaction_duration are specified, an error will result

        :param initial_step:    The suggested size of the first step (it might be reduced automatically,
                                    in case of "hard" errors from large steps)

        :param n_steps:         The desired number of steps

        :param snapshots:       OPTIONAL dict that may contain any the following keys:
                                        -"frequency" (default 1)
                                        -"show_intermediates" (default True)
                                        -"species" (default None, meaning all species)
                                        -"initial_caption" (default blank)
                                        -"final_caption" (default blank)
                                    If provided, take a system snapshot after running a multiple
                                    of "frequency" reaction steps (default 1, i.e. at every step.)
                                    EXAMPLE: snapshots={"frequency": 2, "species": ["A", "H"]}

        :param silent:              If True, less output is generated

        :param variable_steps:      If True, the steps sizes will get automatically adjusted, based on thresholds
        :param explain_variable_steps:

        :return:                None.   The object attributes self.system and self.system_time get updated
        """

        # Validation
        assert self.system is not None, "ReactionDynamics.single_compartment_react(): " \
                                        "the concentration values of the various chemicals must be set first"

        if variable_steps and n_steps is not None:
            raise Exception("ReactionDynamics.single_compartment_react(): if `variable_steps` is True, cannot specify `n_steps` "
                            "(because the number of steps will vary); specify `reaction_duration` or `target_end_time` instead")


        """
        Determine all the various time parameters that were not explicitly provided
        """

        if target_end_time is not None:
            if reaction_duration is not None:
                raise Exception("single_compartment_react(): cannot provide values for BOTH `target_end_time` and `reaction_duration`")
            else:
                assert target_end_time > self.system_time, \
                    f"single_compartment_react(): `target_end_time` must be larger than the current System Time ({self.system_time})"
                reaction_duration = target_end_time - self.system_time

        # Determine the time step,
        # as well as the required number of such steps
        time_step, n_steps = self.specify_steps(total_duration=reaction_duration,
                                                time_step=initial_step,
                                                n_steps=n_steps)
        # Note: if variable steps are requested then n_steps stops being particularly meaningful; it becomes a
        #       hypothetical value, in the (unlikely) event that the step size were never changed

        if target_end_time is None:
            if variable_steps:
                target_end_time = self.system_time + reaction_duration
            else:
                target_end_time = self.system_time + time_step * n_steps


        if snapshots:
            frequency = snapshots.get("frequency", 1)   # If not present, it will be 1
            species = snapshots.get("species")          # If not present, it will be None (meaning show all)
            first_snapshot = True
            if not snapshots.get("show_intermediates"):
                snapshots["show_intermediates"] = True
        else:
            snapshots = {"frequency": 1, "show_intermediates": True, "species": None, "first_snapshot": True}
            frequency = 1
            species = None
            first_snapshot = True


        if self.diagnostics:
            # Save up the current System State, with some extra info, as "diagnostic 'concentration' data"
            system_data = self.get_conc_dict(system_data=self.system)   # The current System State, as a dict
            self.save_diagnostic_conc_data(system_data)


        i = 0
        while self.system_time < target_end_time:
            if (not variable_steps) and (i == n_steps) and np.allclose(self.system_time, target_end_time):
                break       # When dealing with fixed steps, catch scenarios where after performing n_steps,
                            #   the System Time is below the target_end_time because of roundoff error


            # ----------  CORE OPERATION OF MAIN LOOP  ----------
            delta_concentrations, step_actually_taken, recommended_next_step = \
                    self.reaction_step_common(delta_time=time_step, conc_array=self.system,
                                              variable_steps=variable_steps, explain_variable_steps=explain_variable_steps,
                                              step_counter=i)

            # Update the System State
            self.system += delta_concentrations
            if min(self.system) < 0:    # Check for negative concentrations. TODO: redundant, since reaction_step_common() now does that
                print(f"+++++++++++ SYSTEM STATE ERROR: FAILED TO CATCH negative concentration upon advancing reactions from system time t={self.system_time:,.4g}")

            self.system_time += step_actually_taken

            # Preserve some of the data, as requested
            if snapshots and ((i+1)%frequency == 0):
                if first_snapshot and "initial_caption" in snapshots:
                    self.add_snapshot(species=species, caption=snapshots["initial_caption"])
                    first_snapshot = False
                else:
                    self.add_snapshot(species=species)

            i += 1

            if i > 1000 * n_steps:  # To catch infinite loops
                raise Exception("single_compartment_react(): "
                                "the computation is taking a very large number of steps, probably from automatically trying to correct instability;"
                                " trying reducing the time_step")   # TODO: is the explanation correctly phrased?

            if self.diagnostics:
                # Save up the current System State, with some extra info, as "diagnostic 'concentration' data"
                system_data = self.get_conc_dict(system_data=self.system)   # The current System State, as a dict
                self.save_diagnostic_conc_data(system_data)

            if variable_steps:
                time_step = recommended_next_step   # Follow the recommendation of the ODE solver for the next time step to take
        # END while

        # Report as to whether extra steps were automatically added, as well as the total # taken
        n_steps_taken = i

        extra_steps = n_steps_taken - n_steps
        if extra_steps > 0:
            if variable_steps:
                print(f"Some steps were backtracked and re-done, "
                      f"to prevent negative concentrations or excessively large concentration changes")
            else:
                print(f"The computation took {extra_steps} extra step(s) - "
                      f"automatically added to prevent negative concentrations")


        if not silent:
            print(f"{n_steps_taken} total step(s) taken")


        if snapshots and "final_caption" in snapshots:
            self.history.set_caption_last_snapshot(snapshots["final_caption"])



    def reaction_step_common(self, delta_time: float, conc_array,
                             variable_steps=False, explain_variable_steps=False, step_counter=1) -> (np.array, float, float):
        """
        This is the common entry point for both single-compartment reactions,
        and the reaction part of reaction-diffusions in 1D, 2D and 3D.

        "Compartments" may or may not correspond to the "bins" of the higher layers;
        the calling code might have opted to merge some bins into a single "compartment".

        Using the given concentration data for all the applicable species in a single compartment,
        do a single reaction time step for ALL the reactions -
        based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        Return the increment vector for all the chemical species concentrations in the compartment

        NOTES:  - the actual system concentrations are NOT changed
                - this method doesn't decide on step sizes - except in case of aborts and repeats with smaller step - but
                  makes suggestions to the calling module about the next step to best take

        :param delta_time:      The requested time duration of the overall reaction step
        :param conc_array:      All initial concentrations at the start of the reaction step,
                                    as a Numpy array for ALL the chemical species, in their index order;
                                    this can be thought of as the "SYSTEM STATE"
        :param variable_steps:  If True, the step sizes will get automatically adjusted
        :param explain_variable_steps:
        :param step_counter:

        :return:                The triplet:
                                    1) increment vector for the concentrations of ALL the chemical species,
                                        in their index order, as a Numpy array
                                        EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):
                                            array([7. , -21.])
                                    2) step actually taken - which might be smaller than the requested one
                                        because of reducing the step to avoid negative-concentration errors
                                    3) recommended_next_step
        """
        recommended_next_step = delta_time     # Baseline; no reason yet to suggest a change in step size

        # Validate arguments
        assert conc_array is not None, "reaction_step_common(): the argument 'conc_array' must be a Numpy array"


        if 2 in self.verbose_list:
            print(f"************ At SYSTEM TIME: {self.system_time:,.4g}, calling reaction_step_common() with:")
            print(f"             delta_time={delta_time}, conc_array={conc_array}, ")


        delta_concentrations = None

        SMALLEST_VALUE_TO_TRY = delta_time / 1000.       # Used to prevent infinite loops

        normal_exit = False

        while delta_time > SMALLEST_VALUE_TO_TRY:       # TODO: consider moving the inside of the WHILE loop into a separate function
            try:    # We want to catch Exceptions that can arise from excessively large time steps
                    #   that lead to negative concentrations or violation of user-set thresholds ("HARD" or "SOFT" aborts)

                (delta_concentrations, recommended_next_step) = \
                        self.attempt_reaction_step(recommended_next_step, delta_time, conc_array, variable_steps, explain_variable_steps, step_counter)
                normal_exit = True

                break       # IMPORTANT: this is needed because, in the absence of errors, we need to go thru the WHILE loop only once!


            # CATCH any 'ExcessiveTimeStepHard' exception raised in the loop  (i.e. a HARD ABORT)
            except ExcessiveTimeStepHard as ex:
                # Single reactions steps can fail with this error condition if the attempted time step was too large,
                # under the following scenarios:
                #       1. negative concentrations from any one reaction - caught by  validate_increment()
                #       2. negative concentration from the combined effect of multiple reactions - caught in this function
                #print("*** CAUGHT a HARD ABORT")
                print(ex)
                delta_time *= self.error_abort_step_factor       # Reduce the excessive time step by a pre-set factor
                recommended_next_step = delta_time


            # CATCH any 'ExcessiveTimeStepSoft' exception raised in the loop  (i.e. a SOFT ABORT)
            except ExcessiveTimeStepSoft as ex:
                # Single reactions steps can fail with this error condition if the attempted time step was too large,
                #   under the following scenario:
                #       3. excessive norm(s) measures in the overall step - caught in this function (currently only checked in case of the variable-steps option)
                #print("*** CAUGHT a soft ABORT")
                print(ex)
                delta_time *= self.step_factors["abort"]       # Reduce the excessive time step by a pre-set factor
                recommended_next_step = delta_time

        # END while


        if not normal_exit:         # i.e., if no reaction simulation took place in the WHILE loop, above
            raise Exception(f"reaction_step_common(): unable to complete the reaction step.  "
                            f"In spite of numerous automated reductions of the time step, it continues to lead to concentration changes that are considered excessive; "
                            f"consider reducing the original time step, and/or increasing the 'abort' thresholds with set_thresholds(). Current values: {self.thresholds}")

        return  (delta_concentrations, delta_time, recommended_next_step)     # TODO: consider returning tentative_updated_system , since we already computed it




    def attempt_reaction_step(self, recommended_next_step, delta_time, conc_array, variable_steps, explain_variable_steps, step_counter):
        """

        :param recommended_next_step:
        :param delta_time:
        :param conc_array:
        :param variable_steps:
        :param explain_variable_steps:
        :param step_counter:

        :return:                        The pair (delta_concentrations, recommended_next_step)
        """

        # *****  CORE OPERATION  *****
        delta_concentrations = self._reaction_elemental_step(delta_time=delta_time, conc_array=conc_array, rxn_list=None)


        if self.diagnostics:
            diagnostic_data_snapshot = self._delta_conc_dict(delta_concentrations)      # A dict


        if variable_steps:
            decision_data = self.adjust_speed(delta_conc=delta_concentrations, baseline_conc=conc_array)
            step_factor = decision_data[1]
            action = decision_data[0]
            all_norms = decision_data[2]

            if explain_variable_steps:
                print(f"\n(STEP {step_counter}) ANALYSIS: Examining Conc. Changes from System Time {self.system_time:.5g} "
                      f"due to tentative single step of {delta_time:.5g}:")
                print("    Baseline: ", conc_array)
                print("    Deltas:   ", delta_concentrations)
                with np.errstate(divide='ignore'):  # Suppress warning about divisions by zero,
                                                    # which are expected, and no issue here
                    print("    Relative Deltas:   ", delta_concentrations / conc_array)
                print("    Norms:    ", all_norms)

                print("    Thresholds:    ")
                for rule in  self.thresholds:
                    print(self.display_thresholds(rule, all_norms.get(rule['norm'])))

                print("    Step Factors:    ", self.step_factors)
                print(f"    => Action:    {action}  (with step size factor of {step_factor})")


            # Abort the current step if the rate of change is deemed excessive.
            # TODO: maybe ALWAYS check this, regardless of variable-steps option
            if action == "abort":       # NOTE: this is a "strategic" abort, not a hard one from error
                msg =   f"* INFO: the tentative time step ({delta_time:.5g}) " \
                        f"leads to a least one norm value > its ABORT threshold:\n" \
                        f"      -> will backtrack, and re-do step with a SMALLER delta time, multiplied by {step_factor} (set to {delta_time * step_factor:.5g}) " \
                        f"[Step started at t={self.system_time:.5g}, and will rewind there]"
                #print("WARNING: ", msg)
                if self.diagnostics:
                    # Expand the dict diagnostic_data_snapshot
                    diagnostic_data_snapshot['norm_A'] = all_norms.get('norm_A')
                    diagnostic_data_snapshot['norm_B'] = all_norms.get('norm_B')
                    diagnostic_data_snapshot['action'] = "ABORT"
                    diagnostic_data_snapshot['step_factor'] = step_factor
                    diagnostic_data_snapshot['time_step'] = delta_time
                    self.save_diagnostic_decisions_data(data=diagnostic_data_snapshot, caption="excessive norm value(s)")

                    # Make a note of the abort action in all the reaction-specific diagnostics
                    self.comment_diagnostic_rxn_data("aborted: excessive norm value(s)")

                raise ExcessiveTimeStepSoft(msg)    # ABORT THE CURRENT STEP


            recommended_next_step = delta_time * step_factor

            if self.diagnostics:
                # Expand the dict diagnostic_data_snapshot
                diagnostic_data_snapshot['norm_A'] = all_norms.get('norm_A')
                diagnostic_data_snapshot['norm_B'] = all_norms.get('norm_B')
                diagnostic_data_snapshot['action'] = f"OK ({action})"
                diagnostic_data_snapshot['step_factor'] = step_factor
                diagnostic_data_snapshot['time_step'] = delta_time


            if explain_variable_steps:
                msg =   f"the tentative time step ({delta_time:.5g}) results in norm values that leads to the following:\n"

                if step_factor > 1:         # "INCREASE
                    msg +=  f"ACTION: COMPLETE STEP NORMALLY and MAKE THE INTERVAL LARGER, " \
                            f"multiplied by {step_factor} (set to {recommended_next_step:.5g}) at the next round, because all norms are low"
                elif step_factor < 1:     # "DECREASE"
                    msg +=  f"ACTION: COMPLETE STEP NORMALLY and MAKE THE INTERVAL SMALLER, " \
                            f"multiplied by {step_factor} (set to {recommended_next_step:.5g}) at the next round, because at least one norm is high"
                else:   # "STAY THE COURSE"
                    msg +=  f"ACTION: COMPLETE NORMALLY - we're inside the target range.  No change to step size."

                msg += f"\n        [The current step started at System Time: {self.system_time:.5g}, and will continue to {self.system_time + delta_time:.5g}]"
                print("NOTICE:", msg)


        if self.diagnostics:
            self.save_diagnostic_decisions_data(data=diagnostic_data_snapshot)


        # Check whether the COMBINED delta_concentrations will make any conc negative
        if self.system is not None:     # TODO - IMPORTANT: usage of this function doesn't always involve self.system
            tentative_updated_system = self.system + delta_concentrations
            if min(tentative_updated_system) < 0:
                print(f"*** CAUTION: negative concentration resulting from the combined effect of all reactions, "
                      f"upon advancing reactions from system time t={self.system_time:,.5g}\n"
                      f"         It'll be AUTOMATICALLY CORRECTED with a reduction in time step size")

                # A type of HARD ABORT is detected (a negative concentration resulting from the combined effect of all reactions)
                if self.diagnostics:
                    self.save_diagnostic_decisions_data(data={"action": "ABORT",
                                                              "step_factor": self.error_abort_step_factor,
                                                              "caption": "neg. conc. from combined effect of all rxns",
                                                              "time_step": delta_time})
                    for rxn_index in range(self.chem_data.number_of_reactions()):
                        self.save_diagnostic_rxn_data(rxn_index=rxn_index, time_step=delta_time,
                                                      increment_vector_single_rxn=None,
                                                      caption=f"aborted: neg. conc. from combined multiple rxns")

                raise ExcessiveTimeStepHard(f"INFO: the tentative time step ({delta_time:.5g}) "
                                            f"leads to a NEGATIVE concentration of one of the chemicals: "
                                            f"\n      -> will backtrack, and re-do step with a SMALLER delta time, "
                                            f"multiplied by {self.error_abort_step_factor} (set to {delta_time * self.error_abort_step_factor:.5g}) "
                                            f"[Step started at t={self.system_time:.5g}, and will rewind there.  Baseline values: {self.system} ; delta conc's: {delta_concentrations}]")


        return  (delta_concentrations, recommended_next_step)





    def norm_A(self, delta_conc: np.array) -> float:
        """
        Return "version A" of a measure of change, based on the average concentration changes
        of ALL chemicals across a time step, adjusted for the number of chemicals

        :param delta_conc:  A Numpy array with the concentration changes
                                of all the chemicals across a time step
        :return:            A measure of change in the concentrations across a simulation step
        """
        n_chems = self.chem_data.number_of_chemicals()

        assert n_chems == len(delta_conc), \
            f"norm_A(): the number of entries in the passed array ({len(delta_conc)}) " \
            f"does not match the number of registered chemicals ({n_chems})"


        # The following are normalized by the number of chemicals
        #L2_rate = np.linalg.norm(delta_concentrations) / n_chems
        #L2_rate = np.sqrt(np.sum(delta_concentrations * delta_concentrations)) / n_chems
        #print("    L_inf norm:   ", np.linalg.norm(delta_concentrations, ord=np.inf) / delta_time)
        #print("    Adjusted L1 norm:   ", np.linalg.norm(delta_concentrations, ord=1) / n_chems)

        adjusted_L2_rate = np.sum(delta_conc * delta_conc) / (n_chems * n_chems)   # The square of the rate above
        return adjusted_L2_rate


    def norm_B(self, baseline_conc: np.array, delta_conc: np.array) -> float:
        """
        Return "version B" of a measure of change, based on the max absolute relative concentration
        change of all the chemicals across a time step (based on an L infinity norm - but disregarding
        any baseline concentration that is very close to zero)

        :param baseline_conc:   A Numpy array with the concentration of all the chemicals
                                    at the start of a simulation time step
        :param delta_conc:      A Numpy array with the concentration changes
                                    of all the chemicals across a time step
        :return:                A measure of change in the concentrations across a simulation step
        """
        n_chems = self.chem_data.number_of_chemicals()

        assert n_chems == len(baseline_conc), \
            f"norm_B(): the number of entries in the array passed in the arg `baseline_conc` ({len(baseline_conc)}) " \
            f"does not match the number of registered chemicals ({n_chems})"

        assert n_chems == len(delta_conc), \
            f"norm_B(): the number of entries in the array passed in the arg `delta_conc` ({len(delta_conc)}) " \
            f"does not match the number of registered chemicals ({n_chems})"

        largest = 0.
        for i in range(len(baseline_conc)):
            if np.allclose(baseline_conc[i], 0):
                continue            # We ignore any baseline concentrations that are very close to zero

            abs_ratio = abs(delta_conc[i] / baseline_conc[i])
            if abs_ratio > largest:
                largest = abs_ratio

        return largest



    def adjust_speed(self, delta_conc: np.array, baseline_conc=None) -> (str, Union[float, int], dict):
        """

        :param delta_conc:
        :param baseline_conc:
        :return:                A triplet:
                                    1) String with the name of the action to take: either "low", "stay", "high" or "abort"
                                    2) A factor by which to multiple the time step at the next iteration round;
                                       if no change is deemed necessary, 1 is returned
                                    3) A dict of all the computed norms (any of the last ones, except the first one, may be missing),
                                       indexed by their names
        """
        all_norms = {}

        all_small = True

        for i, rule in enumerate(self.thresholds):
            norm_name = rule["norm"]

            if norm_name == "norm_A":
                result = self.norm_A(delta_conc)
            else:
                result = self.norm_B(baseline_conc, delta_conc)

            all_norms[norm_name] = result

            if ("abort" in rule) and (result > rule["abort"]):
                return ("abort", self.step_factors["abort"], all_norms)

            if ("high" in rule) and (result > rule["high"]):
                return ("high", self.step_factors["downshift"], all_norms)

            if all_small and ("low" in rule):
                if result >= rule["low"]:
                    all_small = False


        if all_small:
            return ("low", self.step_factors["upshift"], all_norms)

        return ("stay", 1, all_norms)




    def display_thresholds(self, rule, value):
        """

        :param rule:
        :param value:
        :return:
        """
        s = f"                   {rule['norm']} : "

        if value is None:
            return s + " (skipped; not needed)"

        low = rule.get('low')
        high = rule.get('high')
        abort = rule.get('abort')
        if low is not None and value < low:
            return f"{s}(VALUE {value:.5g}) | low {low} | high {high} | abort {abort}"

        if high is not None and value < high:
            return f"{s}low {low} | (VALUE {value:.5g}) | high {high} | abort {abort}"

        if abort is not None and value < abort:
            return f"{s}low {low} | high {high} | (VALUE {value:.5g}) | abort {abort}"

        return f"{s}low {low} | high {high} | abort {abort} | (VALUE {value:.5g})"



    def _reaction_elemental_step(self, delta_time: float, conc_array: np.array, rxn_list=None) -> np.array:
        """
        Using the given concentration data of ALL the chemical species,
        do the specified SINGLE TIME STEP for ONLY the requested reactions (by default all).

        All computations are based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions (in "forward Euler" approach.)

        Return the Numpy increment vector for ALL the chemical species concentrations, in their index order
        (whether involved in these reactions or not)

        NOTES:  - the actual System Concentrations and the System Time (stored in object variables) are NOT changed
                - if any of the concentrations go negative, an Exception is raised

        :param delta_time:      The time duration of this individual reaction step - assumed to be small enough that the
                                    concentration won't vary significantly during this span.
        :param conc_array:      All initial concentrations at the start of the reaction step,
                                    as a Numpy array for all the chemical species, in their index order;
                                    this can be thought of (as a copy of) the "SYSTEM STATE"
        :param rxn_list:        OPTIONAL list of reactions (specified by their indices) to include in this simulation step ;
                                    EXAMPLE: [1, 3, 7]
                                    If None, do all the reactions

        :return:            The increment vector for the concentrations of ALL the chemical species
                                (whether involved in the reactions or not),
                                as a Numpy array for all the chemical species, in their index order
                            EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):   array([7. , -21.])
        """
        # The increment vector is cumulative for ALL the requested reactions
        increment_vector = np.zeros(self.chem_data.number_of_chemicals(), dtype=float)       # One element per chemical species

        # Compute the forward and reverse "conversions" of all the applicable reactions
        delta_dict = self.compute_all_reaction_deltas(conc_array=conc_array, delta_time=delta_time, rxn_list=rxn_list)
        if 3 in self.verbose_list:
            print(f"      delta_list: {delta_dict}")


        if rxn_list is None:    # Meaning ALL reactions
            rxn_list = range(self.chem_data.number_of_reactions())  # This will be a list of all the reaction index numbers

        number_chemicals = self.chem_data.number_of_chemicals()

        # For each applicable reaction, find the needed adjustments ("deltas")
        #   to the concentrations of the reactants and products,
        #   based on the forward and reverse rates of the reaction
        for rxn_index in rxn_list:      # Consider each reaction in turn
            # TODO: maybe switch to a call to the experimental _reaction_elemental_step_SINGLE_REACTION()

            # One element per chemical species; notice that this array is being RESET for EACH REACTION
            # TODO: instead of using a potentially very large array (mostly of zeros) for each rxn, consider a dict instead
            #       then combine then at end
            increment_vector_single_rxn = np.zeros(number_chemicals, dtype=float)

            rxn = self.chem_data.get_reaction(rxn_index)
            reactants = rxn.extract_reactants()
            products = rxn.extract_products()

            """
            Determine the concentration adjustments as a result of this reaction step, 
            for the reaction being considered
            """

            # The reactants DECREASE based on the quantity (forward reaction - reverse reaction)
            for r in reactants:
                # Unpack data from the reactant r
                species_index = rxn.extract_species_index(r)
                if species_index == rxn.enzyme:
                    #print(f"*** SKIPPING reactant ENZYME {species_index} in reaction {rxn_index}")
                    continue    # Skip if r is an enzyme for this reaction

                stoichiometry = rxn.extract_stoichiometry(r)

                delta_conc = stoichiometry * (- delta_dict[rxn_index])  # Increment to this reactant from the reaction being considered
                # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
                # Note: not enough to detect conc going negative from combined changes from multiple reactions!
                #       Further testing done upstream
                self.validate_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                increment_vector_single_rxn[species_index] += delta_conc


            # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
            for p in products:
                # Unpack data from the reactant r
                species_index = rxn.extract_species_index(p)
                if species_index == rxn.enzyme:
                    #print(f"*** SKIPPING product ENZYME {species_index} in reaction {rxn_index}")
                    continue    # Skip if p is an enzyme for this reaction

                stoichiometry = rxn.extract_stoichiometry(p)

                delta_conc = stoichiometry * delta_dict[rxn_index]  # Increment to this reaction product from the reaction being considered
                # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
                # Note: not enough to detect conc going negative from combined changes from multiple reactions!
                #       Further testing done upstream
                self.validate_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                increment_vector_single_rxn[species_index] += delta_conc


            increment_vector += increment_vector_single_rxn     # Accumulate the increment vector across ALL reactions
                                                                # TODO: consider using a small dict in lieu of increment_vector_single_rxn

            if self.diagnostics:
                self.save_diagnostic_rxn_data(rxn_index=rxn_index, time_step=delta_time,
                                              increment_vector_single_rxn=increment_vector_single_rxn)
        # END for (over rxn_list)

        return increment_vector



    def _reaction_elemental_step_SINGLE_REACTION(self, delta_time: float, conc_array: np.array, increment_vector,
                                                 rxn_index: int,
                                                 delta_dict):     # TODO: EXPERIMENTAL; not yet in use

        if 2 in self.verbose_list:
            print(f"      Determining the conc.'s changes as a result of rxn # {rxn_index}")

        # NOTE: instead of using a potentially very large array (mostly of zeros) for each rxn,
        #       we use a dict instead; then combine all at end
        #increment_vector_single_rxn = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)
        increment_dict_single_rxn = {}      # The key will be the elements in the reaction

        # TODO: turn into a more efficient single step, as as:
        #(reactants, products) = cls.all_reactions.unpack_terms(rxn_index)
        reactants = self.chem_data.get_reactants(rxn_index)
        products = self.chem_data.get_products(rxn_index)


        """
        Determine the concentration adjustments as a result of this reaction step, 
        for the reaction being considered
        """

        # The reactants DECREASE based on the quantity (forward reaction - reverse reaction)
        for r in reactants:
            stoichiometry, species_index, order = r                 # Unpack
            delta_conc = stoichiometry * (- delta_dict[rxn_index])  # Increment to this reactant from the reaction being considered
            # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
            # Note: not enough to detect conc going negative from combined changes from multiple reactions!  Further testing done upstream
            # TODO: it might be more efficient to check this in bulk than with many function calls
            self.validate_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                    rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

            increment_dict_single_rxn[species_index] = increment_dict_single_rxn.get(species_index, 0) + delta_conc
            #increment_vector_single_rxn[species_index] += delta_conc


        # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
        for p in products:
            stoichiometry, species_index, order = p             # Unpack
            delta_conc = stoichiometry * delta_dict[rxn_index]  # Increment to this reaction product from the reaction being considered
            # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
            # Note: not enough to detect conc going negative from combined changes from multiple reactions!  Further testing done upstream
            # TODO: it might be more effecient to check this in bulk than with many function calls
            self.validate_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                    rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

            increment_dict_single_rxn[species_index] = increment_dict_single_rxn.get(species_index, 0) + delta_conc
            #increment_vector_single_rxn[species_index] += delta_conc


        for k, v in increment_dict_single_rxn.items():
            increment_vector[k] += v                # Accumulate the increment vector across all the increments from this reaction
        #increment_vector += increment_vector_single_rxn     # Accumulate the increment vector across all reactions


        # If requested to evaluate the reaction "speeds"
        # TODO: need a version of examine_increment_array() and save_diagnostic_data() that accept increment_dict_single_rxn
        #       instead of increment_vector_single_rxn [also, maybe move one or both of them to calling function]
        '''  
        if self.diagnostics:
            self.save_diagnostic_data(delta_time, increment_vector_single_rxn, rxn_index)
        '''



    def criterion_fast_reaction(self, delta_conc, fast_threshold_fraction,
                                baseline_conc=None, use_baseline=False) -> bool:
        """
        Apply a criterion to determine, from the given data,
        whether the originating reaction (the source of the data) needs to be classified as "Fast".
        All the passed data is for the concentration changes in 1 chemical from 1 reaction

        :param delta_conc:
        :param fast_threshold_fraction:
        :param baseline_conc:           # TODO: probably phase out
        :param use_baseline:            # TODO: gave poor results when formerly used for substeps

        :return:                        True if the concentration change is so large (based on some criteria)
                                            that the reaction that caused it, ought to be regarded as "fast"
        """
        if use_baseline:    # TODO: this criterion gives poor results, and will probably be eliminated
            # if abs(delta_conc) / baseline_conc > fast_threshold_fraction / time_subdivision
            # Perhaps more intuitive as shown above;
            # but, to avoid time-consuming divisions (and potential divisions by zero), re-formulated as below:
            return abs(delta_conc) > fast_threshold_fraction * baseline_conc

        else:
            # Perhaps more intuitive written as:  if abs(delta_conc) > fast_threshold_fraction / time_subdivision
            return abs(delta_conc) > fast_threshold_fraction



    def validate_increment(self,  delta_conc, baseline_conc: float,
                           rxn_index: int, species_index: int, delta_time) -> None:
        """
        Examine the requested concentration change given by delta_conc
        (typically, as computed by an ODE solver),
        relative to the baseline (pre-reaction) value baseline_conc,
        for the given SINGLE chemical species and SINGLE reaction.

        If the concentration change would render the concentration negative,
        raise an Exception (of custom type "ExcessiveTimeStepHard")

        :param delta_conc:              The change in concentration computed by the ode solver
                                            (for the specified chemical, in the given reaction)
        :param baseline_conc:           The initial concentration

        [The remaining arguments are ONLY USED for error printing]
        :param rxn_index:               The index (0-based) to identify the reaction of interest (ONLY USED for error printing)
        :param species_index:           The index (0-based) to identify the chemical species of interest (ONLY USED for error printing)
        :param delta_time:              The time duration of the reaction step (ONLY USED for error printing)

        :return:                        None (an Exception is raised if a negative concentration is detected)
        """
        if (baseline_conc + delta_conc) < 0:
            print(f"\n*** CAUTION: negative concentration in chemical `{self.chem_data.get_name(species_index)}` in step starting at t={self.system_time:.4g}.  "
                  f"It will be AUTOMATICALLY CORRECTED with a reduction in time step size, as follows:")

            # A type of HARD ABORT is detected (a single reaction that, by itself, would lead to a negative concentration;
            #   while it's possible that other coupled reactions might counterbalance this - nonetheless, it's taken as a sign of excessive step size)
            if self.diagnostics:
                self.save_diagnostic_decisions_data(data={"action": "ABORT",
                                                          "step_factor": self.error_abort_step_factor,
                                                          "caption": f"neg. conc. in {self.chem_data.get_name(species_index)} from rxn # {rxn_index}",
                                                          "time_step": delta_time})
                self.save_diagnostic_rxn_data(rxn_index=rxn_index, time_step=delta_time,
                                              increment_vector_single_rxn=None,
                                              caption=f"aborted: neg. conc. in `{self.chem_data.get_name(species_index)}`")

            raise ExcessiveTimeStepHard(f"    INFO: the tentative time step ({delta_time:.5g}) "
                                    f"leads to a NEGATIVE concentration of `{self.chem_data.get_name(species_index)}` "
                                    f"from reaction {self.chem_data.single_reaction_describe(rxn_index=rxn_index, concise=True)} (rxn # {rxn_index}): "
                                    f"\n      Baseline value: {baseline_conc:.5g} ; delta conc: {delta_conc:.5g}"
                                    f"\n      -> will backtrack, and re-do step with a SMALLER delta time, "
                                    f"multiplied by {self.error_abort_step_factor} (set to {delta_time * self.error_abort_step_factor:.5g}) "
                                    f"[Step started at t={self.system_time:.5g}, and will rewind there]")



    def examine_increment_array_OBSOLETE(self, rxn_index: int,
                                         delta_conc_array: np.array,
                                         fast_threshold_fraction,
                                         baseline_conc_array=None,
                                         use_baseline=False
                                         ) -> None:
        """
        TODO: unclear if there's any use for this anymore.       TODO: obsolete

        Examine the requested concentration changes given by delta_conc_array
        (typically, as computed by an ODE solver),
        relative to their baseline (pre-reaction) values baseline_conc_array,
        for the given SINGLE reaction,
        across all the chemicals affected by it.

        If the concentration change is large, relative to its baseline value,
        for ANY of the  chemicals involved in the reaction,
        then mark the given reaction as "Fast" (in its data structure);
        this will OVER-RIDE any previous annotation about the speed of that reaction.

        Note: if reaction is already marked as "Fast" then no action is taken

        :param rxn_index:               An integer with the (zero-based) index to identify the reaction of interest
        :param delta_conc_array:        The change in concentrations computed by the ODE solver
                                            (for ALL the chemicals,
                                            though only the ones involved in the above reaction are looked at)
        :param baseline_conc_array:     The initial concentration (for ALL the chemicals),
                                            before the last simulated reaction step  TODO: phase out
        :param fast_threshold_fraction: The minimum relative size of the concentration baseline over its change, AS A FRACTION,
                                            for a reaction to be regarded as "Slow".
                                            IMPORTANT: this refers to the FULL step size, and will be scaled for time_subdivision
                                            (e.g. a time_subdivision of 2 signifies 1/2 a time step,
                                            and will mean that just 1/2 of the change will constitute a threshold)
        :param use_baseline:            # TODO: gave poor results when formerly used for substeps

        :return:                        None (the equation is marked as "Fast", if appropriate, in its data structure)
        """

        if self.get_rxn_speed(rxn_index) == "F":
            # If the reaction was already marked as "Fast" then no action is taken
            return


        # If we get this far, the reaction came tagged as "SLOW"

        # If the reaction was (tentatively) labeled as "Slow", decide whether to flip it to "Fast"
        #   Note: the status will be used for the current simulation sub-cycle


        # Loop over just the chemicals in the given reaction (not over ALL the chemicals!)
        for i in self.chem_data.get_chemicals_in_reaction(rxn_index):

            delta_conc = delta_conc_array[i]
            if baseline_conc_array is None:
                baseline_conc = None
            else:
                baseline_conc = baseline_conc_array[i]

            # Determine whether the reaction needs to be classified as "Fast"
            if self.criterion_fast_reaction(delta_conc=delta_conc, baseline_conc=baseline_conc,
                                            fast_threshold_fraction=fast_threshold_fraction,
                                            use_baseline=use_baseline):
                self.set_rxn_speed(rxn_index, "F")
                return

        # END for
        # If we get thus far, the reaction is still regarded as "Slow", and its status is left unchanged



    def compute_all_reaction_deltas(self, delta_time: float, conc_array=None, rxn_list=None) -> dict:
        """
        For an explanation of the "reaction delta", see compute_reaction_delta().
        Compute the "reaction delta" for all the specified reaction (by default, all.)
        Return a list with an entry for each reaction, in their index order.

        For background info: https://life123.science/reactions

        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param conc_array:  All initial concentrations at the start of the reaction step,
                                as a Numpy array for all the chemical species, in their index order;
                                this can be thought of as the "SYSTEM STATE"
        :param rxn_list:    OPTIONAL list of reactions (specified by their integer index);
                                if None, do all the reactions.  EXAMPLE: [1, 3, 7]

        :return:            A dict of the differences between forward and reverse "conversions" -
                                for explanation, see compute_reaction_delta().
                                The dict is indexed by the reaction number, and contains as many entries as the
                                number of reactions being investigated
        """
        delta_dict = {}

        if rxn_list is None:    # Meaning ALL reactions
            rxn_list = range(self.chem_data.number_of_reactions())

        # Process the requested reactions
        for i in rxn_list:      # Consider each reaction in turn
            rxn = self.chem_data.get_reaction(i)
            delta = self.compute_reaction_delta(rxn=rxn, conc_array=conc_array, delta_time=delta_time)
            delta_dict[i] = delta

        return delta_dict



    def compute_reaction_delta(self, rxn, delta_time: float, conc_array: np.array) -> float:
        """
        For the SINGLE specified reaction, the given time interval, and the specified concentrations of chemicals,
        compute the difference of the reaction's forward and back "conversions",
        a non-standard term we're using here to refer to delta_time * (Forward_Rate − Reverse_Rate)

        TODO: maybe take the multiplication with delta_time to the calling function

        For background info: https://life123.science/reactions
        What we're computing here, is referred to as:  (Δt)∗delta_forward(n)

        :param rxn:         An object of type "Reaction"
        :param conc_array:  ALL initial concentrations at the start of the reaction step,
                                as a Numpy array for ALL the chemical species, in their index order
                                (regardless of their involvement in the reaction of interest);
                                this can be thought of as the "SYSTEM STATE"
        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                            concentration won't vary significantly during this span

        :return:            The differences between the reaction's forward and reverse "conversions",
                                a non-standard term we're using here to refer to delta_time * (Forward_Rate − Reverse_Rate),
                                for the given reaction during the specified time span
        """

        # TODO: maybe turn into a more efficient single step, as as:
        #(reactants, products, fwd_rate_constant, rev_rate_constant) = cls.all_reactions.unpack_reaction(i)
        reactants = rxn.extract_reactants()
        products = rxn.extract_products()
        fwd_rate_constant = rxn.extract_forward_rate()
        rev_rate_constant = rxn.extract_reverse_rate()

        forward_rate = fwd_rate_constant
        for r in reactants:
            # Unpack data from the reactant r
            species_index = rxn.extract_species_index(r)
            order = rxn.extract_rxn_order(r)
            conc = conc_array[species_index]
            #assert conc is not None, \
            #   f"ReactionDynamics.compute_reaction_delta(): lacking the value for the concentration of the chemical species `{self.reaction_data.get_name(species_index)}`"
            forward_rate *= conc ** order      # Raise to power

        reverse_rate = rev_rate_constant
        for p in products:
            #stoichiometry, species_index, order = p     # Unpack the data of the reaction product p
            species_index = rxn.extract_species_index(p)
            order = rxn.extract_rxn_order(p)
            conc = conc_array[species_index]
            #assert conc is not None, \
            #   f"ReactionDynamics.compute_reaction_delta(): lacking the concentration value for the species `{self.reaction_data.get_name(species_index)}`"
            reverse_rate *= conc ** order     # Raise to power

        return delta_time * (forward_rate - reverse_rate)   # TODO: maybe take the multiplication with delta_time to the calling function




    #############################################################################################
    #                                                                                           #
    #                                      FOR DIAGNOSTICS                                      #
    #                                                                                           #
    #############################################################################################
    def ________FOR_DIAGNOSTICS________(DIVIDER):  # Used to get a better structure view in IDEs such asPycharm
        pass


    def stoichiometry_checker(self, rxn_index: int, conc_arr_before: np.array, conc_arr_after: np.array,
                              suppress_warning=False) -> bool:
        """
        For the indicated reaction, investigate the change in the concentration of the involved chemicals,
        to ascertain whether the change is consistent with the reaction's stoichiometry.
        See https://life123.science/reactions

        IMPORTANT: this function is currently meant for simulations involving only 1 reaction (TODO: generalize)

        NOTE: the concentration changes in chemicals not involved in the specified reaction are ignored

        :param rxn_index:       Integer to identify the reaction of interest
        :param conc_arr_before: Numpy array with the concentrations of ALL the chemicals (whether involved
                                    in the reaction or not), in their index order, before the reaction
        :param conc_arr_after:  Same as above, but after the reaction
        :param suppress_warning:
        :return:                True if the change in reactant/product concentrations is consistent with the
                                    reaction's stoichiometry, or False otherwise
        """
        return self.stoichiometry_checker_from_deltas(rxn_index = rxn_index,
                                                      delta_arr = conc_arr_after-conc_arr_before,
                                                      suppress_warning = suppress_warning)



    def stoichiometry_checker_from_deltas(self, rxn_index: int, delta_arr: np.array, suppress_warning=False) -> bool:
        """
        For the indicated reaction, investigate the change in the concentration of the involved chemicals,
        to ascertain whether the change is consistent with the reaction's stoichiometry.
        See https://life123.science/reactions

        IMPORTANT: this function is currently meant for simulations involving only 1 reaction (TODO: generalize)

        NOTE: the concentration changes in chemicals not involved in the specified reaction are ignored

        :param rxn_index:   Integer to identify the reaction of interest
        :param delta_arr:   Numpy array of numbers, with the concentrations changes of ALL the chemicals (whether involved
                                in the reaction or not), in their index order,
                                as a result of JUST the reaction of interest
                            TODO: maybe also accept a Panda's data frame row
        :param suppress_warning:
        :return:            True if the change in reactant/product concentrations is consistent with the
                                reaction's stoichiometry, or False otherwise
                                Note: if any of the elements of the passed Numpy array is NaN, then True is returned
                                      (because NaN values are indicative of aborted steps; can't invalidate the stoichiometry
                                      check because of that)
        """
        if (not suppress_warning) and (self.chem_data.number_of_reactions() > 1):
            print(f"*** WARNING: {self.chem_data.number_of_reactions()} reactions are present.  "
                  f"stoichiometry_checker() currently only works for 1-reaction simulations")

        self.chem_data.assert_valid_rxn_index(rxn_index)

        if np.isnan(delta_arr).any():
            return True         # The presence of a NaN, anywhere in delta_arr, is indicative of an aborted step

        rxn = self.chem_data.get_reaction(rxn_index)
        reactants = self.chem_data.get_reactants(rxn_index)
        products = self.chem_data.get_products(rxn_index)

        # Pick (arbitrarily) the first reactant,
        # to establish a baseline change in concentration relative to its stoichiometric coefficient
        baseline_term = reactants[0]
        baseline_species = rxn.extract_species_index(baseline_term)
        baseline_stoichiometry = rxn.extract_stoichiometry(baseline_term)
        baseline_ratio =  (delta_arr[baseline_species]) / baseline_stoichiometry
        #print("baseline_ratio: ", baseline_ratio)

        for i, term in enumerate(reactants):
            if i != 0:
                species = rxn.extract_species_index(term)
                stoichiometry = rxn.extract_stoichiometry(term)
                ratio =  (delta_arr[species]) / stoichiometry
                #print(f"ratio for `{self.reaction_data.get_name(species)}`: {ratio}")
                if not np.allclose(ratio, baseline_ratio):
                    return False

        for term in products:
            species = rxn.extract_species_index(term)
            stoichiometry = rxn.extract_stoichiometry(term)
            ratio =  - (delta_arr[species]) / stoichiometry     # The minus in front is b/c we're on the other side of the eqn
            #print(f"ratio for `{self.reaction_data.get_name(species)}`: {ratio}")
            if not np.allclose(ratio, baseline_ratio):
                return False

        return True



    def stoichiometry_checker_entire_run(self) -> bool:
        """
        Verify that the stoichiometry is satisfied in all the reaction (sub)steps,
        using the diagnostic data from an earlier run

        IMPORTANT: this function is currently meant for simulations involving only 1 reaction (TODO: generalize)

        :return:    True if everything checks out, or False otherwise
        """
        if self.diagnostic_rxn_data == {}:
            print("WARNING *** In order to run stoichiometry_checker_entire_run(), "
                  "the diagnostics must be turned on, with set_diagnostics(), prior to the simulation run!")
            return False

        number_rxns = self.chem_data.number_of_reactions()
        assert number_rxns == 1, \
                f"stoichiometry_checker_entire_run(): this function is currently designed for just 1 reaction " \
                f"(whereas {number_rxns} are present)"

        for rxn_index in range(number_rxns):
            diagnostic_data = self.get_diagnostic_rxn_data(rxn_index=rxn_index, print_reaction=False)
            for row_index in range(len(diagnostic_data)):
                df_row = diagnostic_data.loc[row_index]     # Row in the Panda's data frame of diagnostic data
                chemical_delta_list = self._delta_names()           # EXAMPLE: ["Delta A", "Delta B", "Delta C"]
                delta = df_row[chemical_delta_list].to_numpy(dtype='float32')      # Extract select columns from the data frame row, and turn into Numpy array
                status = self.stoichiometry_checker_from_deltas(rxn_index=rxn_index, delta_arr=delta, suppress_warning=True)
                if not status:
                    print(f"Stoichiometry NOT satisfied by reaction # {rxn_index}: "
                          f"see row # {row_index} in the diagnostic data for that reaction, from get_diagnostic_rxn_data()")
                    print(df_row )
                    return False

        return True



    def _delta_names(self) -> [str]:
        """
        Return a list of strings, with the names of ALL the registered chemicals,
        in their index order, each prefixed by the string "Delta "
        EXAMPLE: ["Delta A", "Delta B", "Delta X"]

        :return:    A list of strings
        """
        chemical_list = self.chem_data.get_all_names()      # EXAMPLE: ["A", "B", "X"]
        chemical_delta_list = ["Delta " + name
                               for name in chemical_list]
        return chemical_delta_list



    def _delta_conc_dict(self, delta_conc_arr: np.ndarray) -> dict:
        """
        Convert a Numpy array into a dict, based on all the registered chemicals.
        The keys are the chemical names, prefixed by "Delta "

        :param delta_conc_arr:  A Numpy array of "delta concentrations".  EXAMPLE: array[1.23, 52.2]
        :return:                A dictionary such as {"Delta A": 1.23, "Delta X": 52.2}
        """
        chemical_delta_list = self._delta_names()    # EXAMPLE: ["Delta A", "Delta X"]
        assert len(chemical_delta_list) == len(delta_conc_arr), \
            f"_delta_conc_dict(): mismatch in number of chemicals ({len(chemical_delta_list)} vs. {len(delta_conc_arr)})"

        delta_conc_dict = {}
        for i, delta_name in enumerate(chemical_delta_list):
            delta_conc_dict[delta_name] = delta_conc_arr[i]

        return delta_conc_dict




    #####  DIAGNOSTICS  #####

    def set_diagnostics(self):
        # Turn on the overall diagnostics mode
        self.diagnostics = True

    def unset_diagnostics(self):
        # Turn off the overall diagnostics mode; existing diagnostics data, if any, is left untouched
        self.diagnostics = False



    #####  1. diagnostic_rxn_data  #####

    def save_diagnostic_rxn_data(self, time_step, increment_vector_single_rxn: Union[np.array, None],
                                 rxn_index: int, caption="") -> None:
        """
        Save up diagnostic data for 1 reaction, for a simulation step
        (by convention, regardless of whether the step is completed or aborted)

        :param time_step:                   The duration of the current simulation step
        :param increment_vector_single_rxn: A Numpy array of size equal to the total number of chemical species,
                                                containing the "delta concentrations" for
                                                ALL the chemicals (whether involved in the reaction or not)
        :param rxn_index:                   An integer that indexes the reaction of interest (numbering starts at 0)
        :param caption:                     OPTIONAL string to describe the snapshot
        :return:                            None
        """
        # Validate the reaction index
        self.chem_data.assert_valid_rxn_index(rxn_index)

        # Validate increment_vector_single_rxn
        if increment_vector_single_rxn is None:     # If the values aren't available (as is the case in aborts)
            increment_vector_single_rxn = np.full(self.chem_data.number_of_chemicals(), np.nan)
        else:
            assert len(increment_vector_single_rxn) == self.chem_data.number_of_chemicals(), \
                f"save_diagnostic_rxn_data(): the length of the Numpy array increment_vector_single_rxn " \
                f"({len(increment_vector_single_rxn)}) doesn't match the total numbers of chemicals ({self.chem_data.number_of_chemicals()})"


        if self.diagnostic_rxn_data == {}:    # INITIALIZE the dictionary self.diagnostic_rxn_data, if needed
            for i in range(self.chem_data.number_of_reactions()):
                self.diagnostic_rxn_data[i] = MovieTabular(parameter_name="START_TIME")       # One per reaction

        data_snapshot = {}
        # Add entries to the above dictionary, starting with the Delta conc. for all the chemicals
        for index, conc in enumerate(increment_vector_single_rxn):
            data_snapshot["Delta " + self.chem_data.get_name(index)] = conc

        #data_snapshot["reaction"] = rxn_index           # TODO: now redundant because factored out into separate data frames
        data_snapshot["time_step"] = time_step
        self.diagnostic_rxn_data[rxn_index].store(par=self.system_time,
                                                  data_snapshot=data_snapshot, caption=caption)



    def comment_diagnostic_rxn_data(self, msg: str) -> None:
        """
        Set the comment field of the last record for EACH of the reaction-specific dataframe

        :param msg: Value to set the comment field to
        :return:    None
        """
        for rxn_index in range(self.chem_data.number_of_reactions()):
            self.diagnostic_rxn_data[rxn_index].set_caption_last_snapshot(msg)



    def get_diagnostic_rxn_data(self, rxn_index: int, head=None, tail=None,
                                t=None, print_reaction=True) -> pd.DataFrame:
        """
        Return a Pandas dataframe with the diagnostic run data of the requested SINGLE reaction,
        from the time that the diagnostics were activated by a call to set_diagnostics().

        In particular, the dataframe contains the "Delta" values for each of the chemicals
        involved in the reaction - i.e. the change in their concentrations
        over the time interval that *STARTS* at the value in the "TIME" column.
        (So, there'll be no row with the final current System Time)

        Note: entries are always added, even if an interval run is aborted, and automatically re-done.

        Optionally, print out a brief description of the reaction.

        Optionally, limit the dataframe to a specified numbers of rows at the end,
        or just return one entry corresponding to a specific time
        (the row with the CLOSEST time to the requested one, which will appear in an extra column
        called "search_value")

        :param rxn_index:       An integer that indexes the reaction of interest (numbering starts at 0)
                                TODO: if not specified, show all reactions in turn

        :param head:            (OPTIONAL) Number of records to return,
                                    from the start of the diagnostic dataframe.
        :param tail:            (OPTIONAL) Number of records to return,
                                    from the end of the diagnostic dataframe.
                                    If either the "head" arguments is passed, this argument will get ignored

        :param t:               (OPTIONAL) Individual time to pluck out from the dataframe;
                                    the row with closest time will be returned.
                                    If this parameter is specified, an extra column - called "search_value" -
                                    is inserted at the beginning of the dataframe
                                    If either the "head" or the "tail" arguments are passed, this argument will get ignored

        :param print_reaction:  (OPTIONAL) If True (default), concisely print out the requested reaction

        :return:                A Pandas data frame with (all or some of)
                                    the diagnostic data of the specified reaction.
                                    Columns of the dataframes:
                                    'START_TIME' 'Delta A' 'Delta B'... 'time_step' 'caption'
        """
        # Validate the reaction index
        self.chem_data.assert_valid_rxn_index(rxn_index)

        if print_reaction:
            print("Reaction: ", self.chem_data.single_reaction_describe(rxn_index=rxn_index, concise=True))

        movie_obj = self.diagnostic_rxn_data.get(rxn_index)    # Object of type "MovieTabular"

        if movie_obj is None:
            raise Exception(f"get_diagnostic_rxn_data(): no diagnostics data exists for reaction index {rxn_index} ;"
                            f" did you call set_diagnostics() prior to the simulation run?")

        return movie_obj.get_dataframe(head=head, tail=tail, search_col="START_TIME", search_val=t)



    #####  2. diagnostic_conc_data  #####

    def save_diagnostic_conc_data(self, system_data) -> None:
        """
        To save the diagnostic concentration data during the run, indexed by the current System Time.
        Note: if an interval run is aborted, by convention NO entry is created here

        :return: None
        """
        self.diagnostic_conc_data.store(par=self.system_time,
                                        data_snapshot=system_data)



    def get_diagnostic_conc_data(self) -> pd.DataFrame:
        """
        Return the diagnostic concentration data saved during the run.
        This will be a complete set of simulation steps,
        even if we only saved part of the history during the run

        Note: if an interval run is aborted, by convention NO entry is created here

        :return: A Pandas dataframe, with the columns:
                    'TIME' 	'A' 'B' ... 'caption'
                 where 'A', 'B', ... are all the chemicals
        """
        return self.diagnostic_conc_data.get_dataframe()



    #####  3. save_diagnostic_decisions_data  #####

    def save_diagnostic_decisions_data(self, data, caption="") -> None:
        """
        Used to save the diagnostic concentration data during the run, indexed by the current System Time.
        Note: if an interval run is aborted, by convention an entry is STILL created here

        :return: None
        """
        self.diagnostic_decisions_data.store(par=self.system_time,
                                             data_snapshot=data, caption=caption)


    def get_diagnostic_decisions_data(self) -> pd.DataFrame:
        """
        Determine and return the diagnostic data about concentration changes at every step - EVEN aborted ones

        :return:    A Pandas dataframe with a "TIME" column, and columns for all the "Delta concentration" values
        """
        return self.diagnostic_decisions_data.get_dataframe()



    def get_diagnostic_decisions_data_ALT(self) -> pd.DataFrame:
        """
        Determine and return the diagnostic data about concentration changes at every step - EVEN aborted ones
        TODO: OBSOLETE - BEING PHASED OUT

        :return:    A Pandas dataframe with a "TIME" column, and columns for all the "Delta concentration" values
        """
        fields = self._delta_names()    # EXAMPLE: ["Delta A", "Delta B", "Delta C"]

        rxn_0 = self.get_diagnostic_rxn_data(rxn_index=0, print_reaction=False)

        df_sum = rxn_0[fields]

        time_col = rxn_0["START_TIME"]

        n_rows = len(df_sum)

        for rxn_index in range(1, self.chem_data.number_of_reactions()):
            rxn = self.get_diagnostic_rxn_data(rxn_index=rxn_index, print_reaction=False)
            assert len(rxn)  == n_rows, "get_diagnostic_delta_data(): this function cannot be used because different dataframes" \
                                        "of diagnostic reaction data have different numbers of rows"
            df_sum += rxn[fields]

        df_sum.insert(0, "START_TIME", time_col)

        return df_sum



    def explain_reactions(self) -> bool:
        """
        Provide a detailed explanation of all the steps of the reactions,
        from the saved diagnostic data

        WARNING: Currently designed only for exactly 2 reactions!  TODO: generalize to any number of reactions

        TODO: test and validate usefulness, now that substeps got eliminated
        TODO: allow arguments to specify the min and max reaction times during which to display the explanatory data

        :return:    True if the diagnostic data is consistent for all the steps of all the reactions,
                    or False otherwise
        """
        if not self.diagnostics:
            print("No diagnostic information is available:\nIn order to run explain_reactions(), "
                  "call set_diagnostics() prior to running single_compartment_react()")
            return False

        number_reactions = self.chem_data.number_of_reactions()

        assert number_reactions == 2, \
            "explain_reactions() currently ONLY works when exactly 2 reactions are present. " \
            "Future versions will lift this restriction"

        row_baseline = 0
        row_list = [0, 0]       # TODO: generalize
        active_list = [0, 1]    # ALL the reactions.  TODO: generalize

        self._explain_reactions_helper(active_list=active_list,
                                       row_baseline=row_baseline, row_list=row_list)


        while row_baseline < len(self.get_diagnostic_conc_data()) - 2:
            row_baseline += 1

            print("Advance a single-step across all tables")

            for i in active_list:
                row_list[i] += 1

            status = self._explain_reactions_helper(active_list=active_list, row_baseline=row_baseline, row_list=row_list) # TODO: generalize


            if not status:
                return False    # Error termination

        # END while

        return True             # Successful termination



    def _explain_reactions_helper(self, active_list, row_baseline, row_list) -> bool:
        """
        Helper function for explain_reactions()

        :param active_list:
        :param row_baseline:
        :param row_list:
        :return:            True is the diagnostic data is consistent for this step, or False otherwise
        """
        print("-----------")
        print("ROW of baseline data: ", row_baseline)

        current_time = self.get_diagnostic_conc_data().loc[row_baseline]['TIME']
        print(f"TIME = {current_time:.5g}")

        print("row_list: ", row_list)
        print("active_list: ", active_list)

        chemical_list = self.chem_data.get_all_names()
        chemical_delta_list = self._delta_names()

        conc_arr_before = self.get_diagnostic_conc_data().loc[row_baseline][chemical_list].to_numpy(dtype='float16')
        print("baseline concentrations: ", conc_arr_before)

        delta_cumulative = np.zeros(self.chem_data.number_of_chemicals(),
                                    dtype=float)  # One element per chemical species

        # For each reaction
        for rxn_index in range(self.chem_data.number_of_reactions()):
            if (rxn_index in active_list):  # If reaction is tagged as "fast"
                row = row_list[rxn_index]   # Row in the data frame for the diagnostic data on this reaction
                delta_rxn = self.get_diagnostic_rxn_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy(dtype='float16')
                print(f"From fast rxn {rxn_index}: delta_rxn = {delta_rxn}")
            else:                           # If reaction is tagged as "slow"
                row = row_list[rxn_index]   # Row in the data frame for the diagnostic data on this reaction
                delta_rxn = self.get_diagnostic_rxn_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy(dtype='float16')
                #delta_rxn = np.zeros(self.reaction_data.number_of_chemicals(),
                #                     dtype=float)   # This is the way it was in release beta 19
                print(f"From slow rxn {rxn_index}: delta_rxn = {delta_rxn}")


            delta_cumulative += delta_rxn
        # END for

        print("delta_cumulative: ", delta_cumulative)

        conc_after = conc_arr_before + delta_cumulative
        print("updated concentrations: ", conc_after)

        next_system_state = self.get_diagnostic_conc_data().loc[row_baseline + 1][chemical_list].to_numpy(dtype='float16')
        print(f"concentrations from the next row ({row_baseline + 1}) of the system state: ", next_system_state)

        status = np.allclose(conc_after.astype(float), next_system_state.astype(float))
        if status:
            print("Match OK")
        else:
            print("****************************************   MISMATCH!!!  ****************************************")

        print("-----------")

        return status



    def explain_time_advance(self, return_times=False, silent=False, use_history=False) \
                                                                    -> Union[None, tuple]:
        """
        Use the saved-up diagnostic data, to print out details of the timescales of the reaction run

        If diagnostics weren't enabled ahead of calling this function, an Exception is raised

        EXAMPLE of output:
            From time 0 to 0.0304, in 17 FULL steps of 0.0008
            (for a grand total of 38 FULL steps)

        :param return_times:    If True, all the critical times (times where the interval steps change)
                                    are saved and returned as a list
        :param silent:          If True, nothing gets printed out
        :param use_history:     If True, use the system history in lieu of the diagnostic data;
                                    to keep in mind is the fact that the user might only have asked
                                    for PART of the history to be saved
        :return:                Depending on the argument return_times, either None, or a pair with 2 lists:
                                        1 - list of time values
                                        2 - list of step sizes  (will have one less element than the first list)
        """
        if use_history:
            df = self.get_history()
            assert not df.empty , \
                "ReactionDynamics.explain_time_advance(): no history data found.  " \
                "Did you run the reaction simulation prior to calling this function?"
        else:
            df = self.get_diagnostic_conc_data()
            assert not df.empty , \
                "ReactionDynamics.explain_time_advance(): no diagnostic data found.  " \
                "Did you call set_diagnostics() prior to the reaction run?"

        n_entries = len(df)
        if use_history:
            t = list(df["SYSTEM TIME"])    # List of times (simulation step points)
        else:
            t = list(df["TIME"])    # List of times (simulation step points)

        if n_entries == 1:
            print(f"ReactionDynamics.explain_time_advance(): no time advance found. "
                  f"Diagnostics only contain an initial System Time of {t[0]:.3g} ;  "
                  f"did you run the reaction simulation?")
            return


        # If we get thus far, we have at least 2 entries in the list "t" of times
        grad = np.diff(t)
        grad_shifted = np.insert(grad, 0, 0)  # Insert a zero at the top (shifting down the array)
        #print(grad)
        #print(grad_shifted)
        start_interval = t[0]

        critical_times = [start_interval]  # List of times to report on (start time, plus times where the steps change)
        step_sizes = []

        total_number_full_steps = 0
        for i in range(1, n_entries-1):
            #print(i)
            if not np.allclose(grad[i], grad_shifted[i]):
                #print (f"   Detection at element {i} : time={t[i]}")
                #print (f"   From time {start_interval} to {t[i]}, in steps of {grad_shifted[i]}")
                n_full_steps_taken = self._explain_time_advance_helper(t_start=start_interval, t_end=t[i],
                                                                       delta_baseline=grad_shifted[i], silent=silent)

                start_interval = t[i]
                if n_full_steps_taken > 0:
                    total_number_full_steps += n_full_steps_taken
                    critical_times.append(t[i])
                    step_sizes.append(grad_shifted[i])


        # Final wrap-up of the interval's endpoint
        n_full_steps_taken = self._explain_time_advance_helper(t_start=start_interval, t_end=t[-1],
                                                               delta_baseline=grad_shifted[-1], silent=silent)

        if n_full_steps_taken > 0:
            total_number_full_steps += n_full_steps_taken
            critical_times.append(t[-1])
            step_sizes.append(grad_shifted[-1])

        if not silent:
            print(f"({total_number_full_steps} steps total)")

        if return_times:
            assert len(critical_times) == len(step_sizes) + 1, \
                "explain_time_advance(): validation error in the values to return"
            return (critical_times, step_sizes)



    def _explain_time_advance_helper(self, t_start, t_end, delta_baseline, silent: bool) -> Union[int, float]:
        """
        Using the provided data, about a group of same-size steps, create and print a description of it for the user

        :param t_start:
        :param t_end:
        :param delta_baseline:
        :param silent:          If True, nothing gets printed; otherwise, a line is printed out
        :return:                The corresponding number of FULL steps taken
        """
        if np.allclose(t_start, t_end):
            #print(f"   [Ignoring interval starting and ending at same time {t_start:.3g}]")
            return 0

        n_steps = round((t_end - t_start) / delta_baseline)
        step_s = "step" if n_steps == 1 else "steps"  # singular vs. plural
        if not silent:
            print(f"From time {t_start:.4g} to {t_end:.4g}, in {n_steps} {step_s} of {delta_baseline:.3g}")
        return n_steps




    #####################################################################################################

    '''                                    ~   GRAPHICS   ~                                           '''

    def ________GRAPHICS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def plot_curves(self, chemicals=None, colors=None, title=None, title_prefix=None,
                    vertical_lines=None, show_intervals=False, suppress=False) -> Union[None, go.Figure]:
        """
        Using plotly, draw the plots of concentration values over time, based on the saved history data.
        TODO: offer an option to just display a part of the timeline (e.g. a t_start and t_end)
        TODO: allow alternate labels for x-axis

        EXAMPLE - to combine plots:

            import plotly.graph_objects as go
            fig0 = plot_curves(chemicals=["A", "B", "C"], suppress=True)
            fig1 = px.line(x=[2,2], y=[0,100], color_discrete_sequence = ['gray'])
            all_fig = go.Figure(data=fig0.data + fig1.data, layout = fig0.layout)    # Note that the + is concatenating lists
            all_fig.update_layout(title="My title")
            all_fig.show()

        :param chemicals:   (OPTIONAL) List of the names of the chemicals to plot;
                                if None, then display all
        :param colors:      (OPTIONAL) List of the colors names to use;
                                if None, then use the hardwired defaults
        :param title:       (OPTIONAL) Title for the plot;
                                if None, use default titles that will vary based on the # of reactions; EXAMPLES:
                                    "Changes in concentrations for 5 reactions"
                                    "Reaction `A <-> 2 B` .  Changes in concentrations with time"
                                    "Changes in concentration for `2 S <-> U` and `S <-> X`"
        :param vertical_lines: (OPTIONAL) List or tuple or Numpy array or Pandas series
                                    of x-coordinates at which to draw thin vertical dotted gray lines
        :param show_intervals: (OPTIONAL) If True, it over-rides any value in vertical_lines, and draws
                                    thin vertical dotted gray lines at all the x-coords of the data points in the saved history data;
                                    also, it adds a comment to the title
        :param suppress:    If True, nothing gets shown - and a plotly "Figure" object gets returned instead;
                                this is useful to combine multiple plots (see example above)

        :return:            None or a plotly "Figure" object, depending on the "suppress" flag
        """
        default_colors = ['blue', 'green', 'brown', 'red', 'gray', 'orange', 'purple', 'cyan', 'darkorange', 'navy', 'darkred', 'black']

        df = self.get_history()     # The expected columns are "SYSTEM TIME", followed by the various chemicals

        if chemicals is None:
            chemicals = self.chem_data.get_all_names()      # List of the chemical names.  EXAMPLE: ["A", "B", "H"]

        number_of_curves = len(chemicals)

        if colors is None:
            colors = default_colors[:number_of_curves]      # Pick the first default colors; TODO: rotate if needing more


        if title is None:   # If no title was specified, create one based on how many reactions are present
            number_of_rxns = self.chem_data.number_of_reactions()
            if number_of_rxns > 2:
                title = f"Changes in concentrations for {number_of_rxns} reactions"
            elif number_of_rxns == 1:
                rxn_text = self.chem_data.single_reaction_describe(rxn_index=0, concise=True)   # The only reaction
                title = f"Reaction `{rxn_text}` .  Changes in concentrations with time"
            else:   # Exactly 2 reactions
                rxn_text_0 = self.chem_data.single_reaction_describe(rxn_index=0, concise=True)
                rxn_text_1 = self.chem_data.single_reaction_describe(rxn_index=1, concise=True)
                title = f"Changes in concentration for `{rxn_text_0}` and `{rxn_text_1}`"

        if title_prefix is not None:
            title = f"{title_prefix}.  {title}"

        if show_intervals:
            vertical_lines = df["SYSTEM TIME"]
            title += " (time steps shown in dashed lines)"

        # Create the main plot
        fig = px.line(data_frame=df, x="SYSTEM TIME", y=chemicals,
                      title=title,
                      color_discrete_sequence = colors,
                      labels={"value":"concentration", "variable":"Chemical"})


        if vertical_lines is not None:
            assert (type(vertical_lines) == list) or (type(vertical_lines) == tuple) or (type(vertical_lines) == np.ndarray) or (type(vertical_lines) == pd.core.series.Series), \
                "plot_curves(): the argument `vertical_lines`, if not None, must be a list or tuple or Numpy array or Pandas series of numbers (x-axis coords)"
            for xi in vertical_lines:
                fig.add_vline(x=xi, line_width=1, line_dash="dot", line_color="gray")


        if not suppress:
            fig.show()
        else:
            return fig



    def plot_step_sizes(self, show_intervals=False) -> None:
        """
        Using plotly, draw the plot of the step sizes vs. time
        (only meaningful when the variable-step option was used).
        The same scale as plot_curves() will be used.
        This function requires the diagnostics option to be turned on, prior to running the simulation

        :param show_intervals:  If True, will add to the plot thin vertical dotted gray lines
                                    at the time steps
        :return:                None
        """
        (transition_times, step_sizes) = self.explain_time_advance(return_times=True, silent=True)

        x=transition_times
        y=step_sizes

        # Create a step plot (TODO: there might be a way to directly do this in plotly)
        new_x = [x[0]]
        new_y = [y[0]]
        for i, xi in enumerate(x[1:-1]) :   # Drop the first and last elements
            new_x.append(xi)
            new_y.append(y[i])

            new_x.append(xi)
            new_y.append(y[i+1])

        new_x.append(x[-1])
        new_y.append(y[-1])


        df = self.get_diagnostic_conc_data()    # Pandas dataframe with a column called "TIME"

        # Note: the step size at the final end time isn't a defined quantity - so, we'll just repeat
        #       the last value, to maintain the full x-axis size
        #fig = px.line(x=transition_times, y=step_sizes+[step_sizes[-1]])

        fig = px.line(x=new_x, y=new_y)

        if show_intervals:
            for xi in df["TIME"]:
                fig.add_vline(x=xi, line_width=1, line_dash="dot", line_color="gray")


        fig.update_layout(title='Simulation step sizes',
                          xaxis_title='SYSTEM TIME',
                          yaxis_title='Step size')
        '''
        if show_intervals:
            for xi in transition_times:
                fig.add_vline(x=xi, line_width=1, line_dash="dot", line_color="gray")
        '''
        fig.show()




    #####################################################################################################

    '''                                      ~   HISTORY   ~                                          '''

    def ________HISTORY________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def add_snapshot(self, species=None, caption="", time=None, system_data=None) -> None:
        """
        Preserve some or all the chemical concentrations into the history,
        linked to the passed time (by default the current System Time),
        with an optional caption.

        EXAMPLES:  add_snapshot()
                    add_snapshot(species=['A', 'B']), caption="Just prior to infusion")

        :param species:     (OPTIONAL) list of name of the chemical species whose concentrations we want to preserve for later use.
                                If not specified, save all
        :param caption:     (OPTIONAL) caption to attach to this preserved data
        :param time:        (OPTIONAL) time value to attach to the snapshot (default: current System Time)
        :param system_data: (OPTIONAL) a Numpy array of concentration values, in the same order as the
                                index of the chemical species; by default, use the SYSTEM DATA
                                (which is set and managed by various functions)
        :return:            None
        """

        data_snapshot = self.get_conc_dict(species=species, system_data=system_data)    # This will be a dict

        if time is None:
            time = self.system_time     # By default, use the system time

        self.history.store(par=time,
                           data_snapshot=data_snapshot, caption=caption)



    def get_history(self, t_start=None, t_end=None, head=None, tail=None, t=None) -> pd.DataFrame:
        """
        Retrieve and return a Pandas dataframe with the system history that had been saved
        using add_snapshot()
        Optionally, restrict the result with a start and/or end times,
        or by limiting to a specified numbers of rows at the end

        :param t_start: (OPTIONAL) Start time in the "SYSTEM TIME" column
        :param t_end:   (OPTIONAL) End time
        :param head:    (OPTIONAL) Number of records to return,
                                   from the start of the diagnostic dataframe.
        :param tail:    (OPTIONAL) Number of records to consider, from the end of the dataframe
        :param t:       (OPTIONAL) Individual time to pluck out from the dataframe;
                                   the row with closest time will be returned.
                                   If this parameter is specified, an extra column - called "search_value" -
                                   is inserted at the beginning of the dataframe.
                                   If either the "head" or the "tail" arguments are passed, this argument will get ignored
        :return:        A Pandas dataframe
        """
        df = self.history.get_dataframe(head=head, tail=tail, search_val=t,
                                        search_col="SYSTEM TIME", val_start=t_start, val_end=t_end)

        return df



    def get_historical_concentrations(self, row: int, df=None) -> np.array:
        """
        Return a Numpy array with ALL the chemical concentrations (in their index order)
        from the specified row number of given Pandas data frame (by default, the system history)

        :param row: Integer with the zero-based row number of the system history (which is a Pandas data frame)
        :param df:  (OPTIONAL) A Pandas data frame with concentration information in columns that have
                        the names of the chemicals (if None, the system history is used)
        :return:    A Numpy array of floats.  EXAMPLE: array([200., 40.5])
        """
        if df is None:
            df = self.get_history()

        chem_list = self.chem_data.get_all_names()  # List of all the chemicals' names
        arr = df.loc[row][chem_list].to_numpy(dtype='float32')
        return arr




    #####################################################################################################

    '''                                ~   RESULT ANALYSIS   ~                                        '''

    def ________RESULT_ANALYSIS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def is_in_equilibrium(self, rxn_index=None, conc=None, tolerance=1, explain=True) -> Union[bool, dict]:
        """
        Ascertain whether the given concentrations are in equilibrium for the specified reactions
        (by default, for all reactions)
        TODO: optionally display last lines in diagnostic data

        :param rxn_index:   The index (0-based integer) to identify the reaction of interest;
                                if None, then check all the reactions
        :param conc:        Dict with the concentrations of the species involved in the reaction(s).
                            The keys are the chemical names
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
                            If None, then use the current System concentrations
        :param tolerance:   Allowable relative tolerance, as a PERCENTAGE, to establish satisfactory equality
        :param explain:     If True, print out details about the analysis,
                                incl. the formula(s) being used to check the equilibrium
                                EXAMPLES:   "([C][D]) / ([A][B])"
                                            "[B] / [A]^2"

        :return:            Return True if ALL the reactions are close enough to an equilibrium,
                                as allowed by the requested tolerance;
                                otherwise, return a dict of the form {False: [list of reaction index]}
                                for all the reactions that failed the criterion
                                (EXAMPLE:  {False:  [3, 6:]})
        """
        if conc is None:
            conc=self.get_conc_dict()   # Use the current System concentrations

        failures_dict = {False: []}     # 1-element dict whose value is
        # a list of reactions that fail to meet the criterion for equilibrium

        if rxn_index is not None:
            # Check the 1 given reaction
            if explain:
                description = self.chem_data.single_reaction_describe(rxn_index=rxn_index, concise=True)
                print(description)

            status = self.reaction_in_equilibrium(rxn_index=rxn_index, conc=conc, tolerance=tolerance, explain=explain)
            if not status:
                failures_dict = {False: [rxn_index]}

        else:
            # Check all the reactions
            status = True   # Overall status
            description_list = self.chem_data.multiple_reactions_describe(concise=True)
            for i in range(self.chem_data.number_of_reactions()):
                # For each reaction
                if explain:
                    print(description_list[i])

                single_status = self.reaction_in_equilibrium(rxn_index=i, conc=conc, tolerance=tolerance, explain=explain)

                if not single_status:
                    status = False
                    failures_dict[False].append(i)

        if status:
            return True
        else:
            return failures_dict



    def reaction_in_equilibrium(self, rxn_index, conc, tolerance, explain: bool) -> bool:
        """
        Ascertain whether the given concentrations are in equilibrium for the specified SINGLE reaction;
        return True or False, accordingly.

        :param rxn_index:   The index (0-based integer) to identify the reaction of interest
        :param conc:        Dict with the concentrations of the species involved in the reaction.
                            The keys are the chemical names
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
        :param tolerance:   Allowable relative tolerance, as a PERCENTAGE, to establish satisfactory equality
        :param explain:     If True, print out the formula being used, as well as some other info.
                                EXAMPLES of formulas:   "([C][D]) / ([A][B])"
                                                        "[B] / [A]^2"
        :return:            True if the given reaction is close enough to an equilibrium,
                            as allowed by the requested tolerance
        """
        rxn = self.chem_data.get_reaction(rxn_index)

        reactants = rxn.extract_reactants()    # A list of triplets
        products = rxn.extract_products()      # A list of triplets
        kF = rxn.extract_forward_rate()
        kR = rxn.extract_reverse_rate()

        rate_ratio = kF / kR                    # Ratio of forward/reverse reaction rates

        conc_ratio = 1.
        numerator = ""
        denominator = ""
        all_concs = []      # List of strings

        for p in products:
            # Loop over the reaction products
            species_index =  rxn.extract_species_index(p)
            rxn_order =  rxn.extract_rxn_order(p)

            species_name =  self.chem_data.get_name(species_index)
            species_conc = conc.get(species_name)
            assert species_conc is not None, f"reaction_in_equilibrium(): unable to proceed because the " \
                                             f"concentration of `{species_name}` was not provided"
            conc_ratio *= (species_conc ** rxn_order)
            if explain:
                all_concs.append(f"[{species_name}] = {species_conc:,.4g}")
                numerator += f"[{species_name}]"
                if rxn_order > 1:
                    numerator += f"^{rxn_order} "

        if explain and len(products) > 1:
            numerator = f"({numerator})"

        for r in reactants:
            # Loop over the return
            species_index =  rxn.extract_species_index(r)
            rxn_order =  rxn.extract_rxn_order(r)

            species_name = self.chem_data.get_name(species_index)
            species_conc = conc.get(species_name)
            assert species_conc is not None, f"reaction_in_equilibrium(): unable to proceed because the " \
                                             f"concentration of `{species_name}` was not provided"
            conc_ratio /= (species_conc ** rxn_order)
            if explain:
                all_concs.append(f"[{species_name}] = {species_conc:,.4g}")
                denominator += f"[{species_name}]"
                if rxn_order > 1:
                    denominator += f"^{rxn_order} "


        status = np.allclose(conc_ratio, rate_ratio, rtol=tolerance/100., atol=0)

        if explain:
            if len(reactants) > 1:
                denominator = f"({denominator})"

            print(f"Final concentrations: ", " ; ".join(all_concs))
            print(f"1. Ratio of reactant/product concentrations, adjusted for reaction orders: {conc_ratio:,.6g}")
            print(f"    Formula used:  {numerator} / {denominator}")
            print(f"2. Ratio of forward/reverse reaction rates: {rate_ratio}")
            print(f"Discrepancy between the two values: {100 * abs(conc_ratio - rate_ratio)/rate_ratio :,.4g} %")
            if status:
                print(f"Reaction IS in equilibrium (within {tolerance:.2g}% tolerance)\n")
            else:
                print(f"Reaction is NOT in equilibrium (not within {tolerance:.2g}% tolerance)\n")

        return status



    def extract_delta_concentrations(self, df, row_from: int, row_to: int, chem_list: [str]) -> np.array:
        """
        Extract the concentration changes of the specified chemical species from a Pandas dataframe
        of concentration values

        EXAMPLE:  extract_delta_concentrations(my_dataframe, 7, 8, ['A', 'B'])

        :param df:
        :param row_from:
        :param row_to:
        :param chem_list:
        :return:            A Numpy array of floats
        """
        from_values = df.loc[row_from][chem_list]
        to_values = df.loc[row_to][chem_list]
        return (to_values - from_values).astype(float).to_numpy(dtype='float32')



    def curve_intersection(self, chem1, chem2, t_start, t_end) -> (float, float):
        """
        Find and return the intersection of the 2 curves in the columns var1 and var2,
        in the time interval [t_start, t_end]
        If there's more than one intersection, only one - in an unpredictable choice - is returned
        TODO: the current implementation fails in cases where the 2 curves stay within some distance of each other,
              and then one curve jumps on the opposite side of the other curve, at at BIGGER distance.
              See the missed intersection at the end of experiment "reactions_single_compartment/up_regulate_1"

        :param chem1:   The name of the 1st chemical of interest
        :param chem2:   The name of the 2nd chemical of interest
        :param t_start: The start of the time interval being considered
        :param t_end:   The end of the time interval being considered
        :return:        The pair (time of intersection, common value)
        """
        df = self.get_history(t_start=t_start, t_end=t_end) # A Pandas dataframe

        row_index = abs(df[chem1] - df[chem2]).idxmin()     # The index of the Pandas dataframe row
                                                            #   with the smallest absolute value
                                                            #   of difference in concentrations
        print(f"Min abs distance found at data row: {row_index}")

        return num.curve_intersect_interpolate(df, row_index,
                                               x="SYSTEM TIME", var1=chem1, var2=chem2)
