import math
import numpy as np
import pandas as pd
import plotly.express as px
from typing import Union
from src.modules.movies.movies import MovieTabular
from src.modules.numerical.numerical import Numerical as num


class ExcessiveTimeStep(Exception):
    """
    Used to raise Exceptions arising from excessively large time steps
    (that lead to negative concentration values)
    """
    pass


#############################################################################################

class ReactionDynamics:
    """
    Used to simulate the dynamics of reactions (in a single compartment.)
    In the context of Life123, this may be thought of as a "zero-dimensional system"
    """

    def __init__(self, reaction_data):
        """
        :param reaction_data:   Object of type "ReactionData" (with data about the chemicals and their reactions)
                                    It's acceptable to pass None,
                                    and take care of it later (though probably a bad idea!)
                                    TODO: maybe offer an option to let the constructor instantiate that object?
        """
        self.reaction_data = reaction_data

        self.system = None  # Concentration data in the single compartment we're simulating, for all the chemicals
                            # A Numpy array of the concentrations of all the chemical species, in their index order
                            # 1-dimensional NumPy array of floats, whose size is the number of chemical species.
                            # Each entry is the concentration of the species with that index (in the "ReactionData" object)
                            # Note that this is the counterpart - with 1 less dimension - of the array by the same name
                            #       in the class BioSim1D

        self.system_time = 0.   # Global time of the system, from initialization on

        self.history = MovieTabular()   # To store user-selected snapshots of (some of) the chemical concentrations,
                                        #   whenever requested by the user.

        self.reaction_speeds = {}       # A dictionary, where the keys are reaction indices,
                                        #   and the values are either "S" (Slow) or "F" (Fast)
                                        #   EXAMPLE : { 1: "F", 4: "S", 5: "F" , 8: "S" }
                                        #   Any reaction with a missing entry is regarded as "F" (Fast)

        self.last_abs_fast_threshold = None # The last-used value of the threshold for fast reactions

        self.variable_steps = True      # EXPERIMENTAL

        self.verbose_list = []          # A list of integers with the codes of the desired verbose checkpoints
                                        #   EXAMPLE: [1, 3] to invoke sections of code marked as 1 or 3
                                        #   Those sections will have entry points such as "if 1 in self.verbose_list"
                                        #   MAX code currently used: 5

        self.diagnostic_data = {}       # This will be a dict with as many entries as reactions
        self.diagnostic_data_baselines = MovieTabular(parameter_name="TIME")    # An expanded version of the normal System History

        self.diagnostics = False

        self.fast_criterion_use_baseline = False    # TODO: EXPERIMENTAL - Probably to be phased out,
                                                    #       because True gives poor results




    #############################################################################################
    #                                                                                           #
    #                                ~  TO SET/READ DATA  ~                                     #
    #                                                                                           #
    #############################################################################################
    def ___________TO_SET_AND_READ_DATA___________(DIVIDER):  # Used to get a better structure listing in IDEs such asPycharm
        pass


    def set_conc(self, conc: Union[list, tuple, dict], snapshot=True) -> None:     # TODO: maybe rename set_all_conc()
        """
        Set the concentrations of ALL the chemicals at once   TODO: maybe a dict indicates a selection of chems, while a list, etc means "ALL"

        :param conc:        A list or tuple of concentration values for ALL the chemicals, in their index order;
                                alternatively, a dict indexed by the chemical names, again for ALL the chemicals
                                EXAMPLE of the latter: {"A": 12.4, "B": 0.23, "C": 2.6}  (assuming that "A", "B", "C" are ALL the chemicals)
                                (Any previous values will get over-written)
                                TODO: [maybe the dict version should be the responsibility of set_chem_conc() instead;] OR
                                      allow the dict to be a subset of all chemicals
                                TODO: also allow a Numpy array; make sure to do a copy() to it!
                                TODO: pytest for the dict option
        :param snapshot:    (OPTIONAL) boolean: if True, add to the history
                                a snapshot of this state being set.  Default: True
        :return:            None
        """
        # TODO: more validations, incl. of individual values being wrong data type

        assert len(conc) == self.reaction_data.number_of_chemicals(), \
            f"ReactionDynamics.set_conc(): the number of concentration values passed in the argument 'conc' ({len(conc)}) " \
            f"must match the number of declared chemicals ({self.reaction_data.number_of_chemicals()})"

        if type(conc) == dict:
            conc_list = []
            for species_index in range(self.reaction_data.number_of_chemicals()):
                name = self.reaction_data.get_name(species_index)
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

        self.set_rxn_speed_all_fast()   # Reset all the reaction speeds to "Fast"

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
            species_index = self.reaction_data.get_index(species_name)
        elif species_index is not None:
            self.reaction_data.assert_valid_species_index(species_index)
        else:
            raise Exception("ReactionDynamics.set_chem_conc(): at least one "
                            "of the arguments `species_index` or `species_name` must be provided")


        # TODO: if setting concentrations of a newly-added chemical, needs to first expand self.system
        self.system[species_index] = conc

        self.set_rxn_speed_all_fast()   # Reset all the reaction speeds to "Fast"
                                        # TODO: this is overkill; ought to only reset the affected reactions -
                                        #       as returned by get_reactions_participating_in()

        if snapshot:
            self.add_snapshot(caption=f"Set concentration of `{self.reaction_data.get_name(species_index)}`")



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
            assert system_data.size == self.reaction_data.number_of_chemicals(), \
                f"ReactionDynamics.get_conc_dict(): the argument `system_data` must be a 1-D Numpy array with as many entries " \
                f"as the declared number of chemicals ({self.reaction_data.number_of_chemicals()})"


        if species is None:
            if system_data is None:
                return {}
            else:
                return {self.reaction_data.get_name(index): system_data[index]
                        for index, conc in enumerate(system_data)}
        else:
            assert type(species) == list or  type(species) == tuple, \
                f"ReactionDynamics.get_conc_dict(): the argument `species` must be a list or tuple" \
                f" (it was of type {type(species)})"

            conc_dict = {}
            for name in species:
                species_index = self.reaction_data.get_index(name)
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
        return [i for i in range(self.reaction_data.number_of_reactions())
                if i not in self.slow_rxns()]
        # Alternate way:
        # return list(set(range(self.reaction_data.number_of_reactions()).difference(self.slow_rxns()))



    def are_all_slow_rxns(self) -> bool:
        """
        Return True iff all the reactions are marked as "slow"
        :return:
        """
        return len(self.slow_rxns()) == self.reaction_data.number_of_reactions()



    def get_rxn_speed(self, rxn_index: int) -> str:
        """
        For the requested reaction, get the string code that it was marked with
        to classify its speed.
        If the reaction has not been previously classified, regard it as "F" (Fast)

        :param rxn_index:   The index (0-based) to identify the reaction of interest
        :return:            A 1-letter string with the code "F" (for Fast) or "S" (Slow)
        """
        return self.reaction_speeds.get(rxn_index, "F")     # Any missing entry is regarded to be "F" (Fast)


    def set_rxn_speed(self, rxn_index: int, speed: str) -> None:
        """
        Set a code value that classifies the reaction speed to tag the given reaction to

        :param rxn_index:   The index (0-based) to identify the reaction of interest
        :param speed:       A 1-letter string with the code "F" (for Fast) or "S" (Slow)
        :return:            None
        """
        assert speed in ["S", "F"], "set_rxn_speed(): `speed` argument must be either 'S' or 'F'"
        self.reaction_data.assert_valid_rxn_index(rxn_index)
        self.reaction_speeds[rxn_index] = speed


    def set_rxn_speed_all_fast(self) -> None:
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
        self.reaction_data.clear_reactions_data()
        self.set_rxn_speed_all_fast()   # Reset all the reaction speeds to "Fast"




    #############################################################################################
    #                                                                                           #
    #                                   TO VISUALIZE SYSTEM                                     #
    #                                                                                           #
    #############################################################################################
    def ___________TO_VISUALIZE_SYSTEM___________(DIVIDER):  # Used to get a better structure listing in IDEs such asPycharm
        pass


    def describe_state(self) -> None:
        """
        A simple printout of the state of the system
        :return:        None
        """
        print(f"SYSTEM STATE at Time t = {self.system_time:,.8g}:")

        n_species = self.reaction_data.number_of_chemicals()
        print(f"{n_species} species:")

        # Show a line of line of data for each chemical species in turn
        for species_index, name in enumerate(self.reaction_data.get_all_names()):
            if name:    # If a name was provided, show it
                name = f" ({name})"
            else:
                name = ""

            print(f"  Species {species_index}{name}. Conc: {self.system[species_index]}")




    #############################################################################################
    #                                                                                           #
    #                                       HISTORY                                             #
    #                                                                                           #
    #############################################################################################
    def ___________HISTORY___________(DIVIDER):  # Used to get a better structure listing in IDEs such asPycharm
        pass


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



    def get_history(self, t_start=None, t_end=None, tail=None) -> pd.DataFrame:
        """
        Retrieve and return a Pandas dataframe with the system history that had been saved
        using add_snapshot()
        Optionally, restrict the result with a start and/or end times,
        or by limiting to a specified numbers of rows at the end

        :param t_start: (OPTIONAL) Start time in the "SYSTEM TIME" column
        :param t_end:   (OPTIONAL) End time
        :param tail:    (OPTIONAL) Number of records to consider, from the end of the dataframe
        :return:        A Pandas dataframe
        """
        df = self.history.get(tail=tail)

        if (t_start is not None) and (t_end is not None):
            return df[df["SYSTEM TIME"].between(t_start, t_end)]

        if t_start is not None:
            return df[df["SYSTEM TIME"] >= t_start]

        if t_end is not None:
            return df[df["SYSTEM TIME"] <= t_end]

        return df



    def get_historical_concentrations(self, row: int, df=None) -> np.array:
        """
        Return a Numpy array with ALL the chemical concentrations (in their index order)
        from the specified row number of given Pandas data frame (by default, the system history)

        :param row: Integer with the zero-based row number of the system history (which is a Pandas data frame)
        :param df:  (OPTIONAL) A Pandas data frame with concentration information in columns that have
                        the names of the chemicals (if None, the system history is used)
        :return:    A Numpy array.  EXAMPLE: array([200., 40.5])
        """
        if df is None:
            df = self.get_history()

        chem_list = self.reaction_data.get_all_names()  # List of all the chemicals' names
        arr = df.loc[row][chem_list].to_numpy()
        return arr




    #############################################################################################
    #                                                                                           #
    #                                TO PERFORM THE REACTIONS                                   #
    #                                                                                           #
    #############################################################################################
    def ___________TO_PERFORM_THE_REACTIONS___________(DIVIDER):  # Used to get a better structure listing in IDEs such asPycharm
        pass


    def specify_steps(self, total_duration=None, time_step=None, n_steps=None) -> (float, int):
        """
        Receive 2 out of 3 possible parameters, and determine the 3rd one:
                        total_duration, time_step, n_steps

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

        return (time_step, n_steps)     # Note: could opt to also return total_duration is there's a need for it



    def single_compartment_react(self, reaction_duration=None, stop_time=None, time_step=None, n_steps=None,
                                 snapshots=None, silent=False,
                                 dynamic_substeps=1, rel_fast_threshold=None, abs_fast_threshold=None) -> None:
        """
        Perform ALL the reactions in the single compartment -
        based on the INITIAL concentrations,
        which are used as the basis for all the reactions.

        Update the system state and the system time accordingly
        (object attributes self.system and self.system_time)

        :param reaction_duration:  The overall time advance for the reactions (i.e. time_step * n_steps)
        :param stop_time:       The final time at which to stop the reaction
                                    If both stop_time and reaction_duration are specified, an error will result
        :param time_step:       The size of each time step.
                                Note: if a dynamic_substeps > 1 is passed, then the time step will get internally
                                subdivided for the "fast" reactions
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

        :param dynamic_substeps:    An integer >= 1.  If > 1, individual steps may get divided by that factor,
                                    on a reaction-by-reaction basis,
                                    if that reaction has fast dynamics,
                                    or multiplied by that factor, if that reaction has slow dynamics

        :param rel_fast_threshold:  Only used when dynamic_substeps > 1
                                    The minimum relative size of the concentration [baseline over its] change, AS A PERCENTAGE,
                                    for a reaction to be regarded as "Slow".
                                    IMPORTANT: this refers to the FULL step size : it's RELative to it
        :param abs_fast_threshold:  Similar to fast_threshold, but in terms of time units - NOT relative to the requested step size
                                    TODO: maybe rename fast_threshold

        :param silent:              If True, less output is generated

        :return:                None.   The object attributes self.system and self.system_time get updated
        """
        # Validation
        assert self.system is not None, "ReactionDynamics.single_compartment_react(): " \
                                        "the concentration values of the various chemicals must be set first"


        """
        Determine all the various time parameters that are not explicitly provided
        """

        if stop_time is not None:
            if reaction_duration is not None:
                raise Exception("single_compartment_react(): cannot provide values for BOTH `stop_time` and `reaction_duration`")
            else:
                assert stop_time > self.system_time, \
                    f"single_compartment_react(): `stop_time` must be larger than the current System Time ({self.system_time})"
                reaction_duration = stop_time - self.system_time

        # Determine the time step (the full step, regardless of optional substeps),
        # as well as the required number of such steps
        time_step, n_steps = self.specify_steps(total_duration=reaction_duration,
                                                time_step=time_step,
                                                n_steps=n_steps)

        if stop_time is None:
            stop_time = self.system_time + time_step * n_steps


        if dynamic_substeps > 1:   # If variable substep size was requested
            if abs_fast_threshold is not None:
                if rel_fast_threshold is not None:
                    raise Exception("single_compartment_react(): "
                                    "cannot provide values for BOTH `rel_fast_threshold` and `abs_fast_threshold`")
                else:
                    rel_fast_threshold = abs_fast_threshold * time_step * 100.  # The *100. is b/c rel_fast_threshold is expressed as %
                    print(f"single_compartment_react(): setting rel_fast_threshold to {rel_fast_threshold}")
            elif rel_fast_threshold is None:
                raise Exception("single_compartment_react(): at least one of `rel_fast_threshold` and `abs_fast_threshold` must be provided")
            else:
                abs_fast_threshold = rel_fast_threshold / (time_step * 100.)
                print(f"single_compartment_react(): setting abs_fast_threshold to {abs_fast_threshold}")

            if self.last_abs_fast_threshold and (abs_fast_threshold < self.last_abs_fast_threshold):
                # If we're changing (from the previous run) our standards of what constitutes a "fast" reaction,
                #   in the direction of possibly making more reaction fast,
                #   then play it safe and initially regard all reactions as fast
                self.set_rxn_speed_all_fast()   # Reset all the reaction speeds to "Fast"

            self.last_abs_fast_threshold = abs_fast_threshold
        # END if dynamic_substeps


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
            # Save up the current System State, with some extra info, as "diagnostic 'baseline' data"
            system_data = self.get_conc_dict(system_data=self.system)   # The current System State, as a dict
            system_data["is_primary"] = True
            system_data["primary_timestep"] = time_step
            system_data["n_substeps"] = dynamic_substeps
            system_data["substep_number"] = 0
            self.diagnostic_data_baselines.store(par=self.system_time,
                                                 data_snapshot=system_data)


        i = 0
        extra_steps = 0     # To record whether any extra steps, not requested by the user,
                            #   were automatically taken to avoid negative-concentration errors
                            # NOTE: dynamic substeps, if requested, have no bearing on the number of (full) steps that were taken
        while self.system_time < stop_time:
            delta_concentrations, step_actually_taken = self.reaction_step_orchestrator(delta_time_full=time_step, conc_array=self.system,
                                                                                        snapshots=snapshots,
                                                                                        dynamic_substeps=dynamic_substeps, rel_fast_threshold=rel_fast_threshold)
            # Update the System State
            self.system += delta_concentrations
            if min(self.system) < 0:    # Check for negative concentrations. TODO: redundant, since reaction_step_orchestrator() now does that
                print(f"+++++++++++ SYSTEM STATE ERROR: FAILED TO CATCH negative concentration upon advancing reactions from system time t={self.system_time:,.4g}")

            self.system_time += step_actually_taken
            if step_actually_taken < time_step:
                extra_steps += 1

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
                                " trying reducing the delta_time")

            if self.diagnostics:
                # Save up the current System State, with some extra info, as "diagnostic 'baseline' data"
                system_data = self.get_conc_dict(system_data=self.system)   # The current System State, as a dict
                system_data["is_primary"] = True
                system_data["primary_timestep"] = time_step
                system_data["n_substeps"] = dynamic_substeps
                system_data["substep_number"] = 0
                self.diagnostic_data_baselines.store(par=self.system_time,
                                                     data_snapshot=system_data)
        # END while

        # Report as to whether extra steps were automatically added, as well as the total # taken
        n_steps_taken = i

        if extra_steps > 0:
            #if dynamic_substeps > 1:
                #print(f"The computation took {extra_steps} extra step(s) - "
                      #f"automatically added to prevent negative concentrations and/or as a result of the requested dynamic steps")
            #else:
            print(f"The computation took {extra_steps} extra step(s) - "
                  f"automatically added to prevent negative concentrations")
            if dynamic_substeps > 1:
                print("(an extra dynamic substeps that might have been taken, are NOT part of the above count)")

        if not silent:
            print(f"{n_steps_taken} total step(s) taken")


        if snapshots and "final_caption" in snapshots:
            self.history.set_caption_last_snapshot(snapshots["final_caption"])



    def reaction_step_orchestrator(self, delta_time_full: float, conc_array, snapshots=None,
                                   dynamic_substeps=1, rel_fast_threshold=5) -> (np.array, float):
        """     TODO: the word "orchestrator" may no longer be a good descriptor
        This is the common entry point for both single-compartment reactions,
        and the reaction part of reaction-diffusions in 1D, 2D and 3D.

        "Compartments" may or may not correspond to the "bins" of the higher layers;
        the calling code might have opted to merge some bins into a single "compartment".

        Using the given concentration data for all the applicable species in a single compartment,
        do a single reaction time step for ALL the reactions -
        based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        Return the increment vector for all the chemical species concentrations in the compartment

        NOTE:  the actual system concentrations are NOT changed

        :param delta_time_full: The requested time duration of the overall reaction step,
                                    which will be subdivided for the "fast" reactions, if dynamic_substeps is greater than 1
        :param conc_array:      All initial concentrations at the start of the reaction step,
                                    as a Numpy array for ALL the chemical species, in their index order;
                                    this can be thought of as the "SYSTEM STATE"
        :param snapshots:       See explanation under single_compartment_react()
        :param dynamic_substeps: An integer by which to subdivide the time intervals for "fast" reactions;
                                    default 1, i.e. no subdivisions
        :param rel_fast_threshold:  The minimum relative size of the concentration baseline over its change, AS A PERCENTAGE,
                                for a reaction to be regarded as "Slow".  IMPORTANT: this refers to the FULL step size

        :return:                The pair:
                                    1) increment vector for the concentrations of ALL the chemical species,
                                        as a Numpy array for all the chemical species, in their index order
                                        EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):
                                            array([7. , -21.])
                                    2) step actually taken  - which might be smaller than the requested one
                                        because of reducing the step to avoid negative-concentration errors
                                        NOTE: dynamic substeps, if requested, have no bearing on the (full) step size being returned
        """
        if 5 in self.verbose_list:
            print(f"Calling reaction_step_orchestrator() with delta_time={delta_time_full}, "
                  f"conc_array={conc_array}, dynamic_substeps={dynamic_substeps}")

        assert conc_array is not None, "reaction_step_orchestrator(): the argument 'conc_array' must be set to a Numpy array"
        assert type(dynamic_substeps) == int, "reaction_step_orchestrator(): the argument 'dynamic_substeps' must be an integer"
        assert dynamic_substeps >= 1, "reaction_step_orchestrator(): the argument 'dynamic_substeps' must be an integer greater or equal than 1"

        if dynamic_substeps > 1:   # If the variable time step option was requested
            assert rel_fast_threshold > 0, "reaction_step_orchestrator(): the argument 'fast_threshold' must be greater than 0"
            fast_threshold_fraction = rel_fast_threshold / 100.     # Here we switch over from percentages to fractions


        if 2 in self.verbose_list:
            print(f"************ SYSTEM TIME: {self.system_time:,.4g}")


        delta_concentrations = None

        while delta_time_full > 0.000001:   # TODO: tweak this number (used to prevent infinite loops)
            try:
                # We want to catch errors that can arise from excessively large time steps that lead to negative concentrations
                if 1 in self.verbose_list:
                    print(f"reaction_step_orchestrator(): entering WHILE loop at System Time={self.system_time:.5g} "
                          f"with delta_time_full={delta_time_full:.5g}")

                if dynamic_substeps == 1:   # If the variable time step option was NOT requested
                    if 2 in self.verbose_list:
                        print("    NO adaptive variable time resolution used")
                        print(f"    Processing ALL the {self.reaction_data.number_of_reactions()} reaction(s) with a single step")

                    delta_concentrations = self._reaction_elemental_step(delta_time=delta_time_full, conc_array=conc_array, rxn_list=None,
                                                                         tag_reactions=False)
                    if self.variable_steps:
                        n_chems = self.reaction_data.number_of_chemicals()
                        print(f"EXAMINING CONCENTRATION CHANGES at System Time {self.system_time} from the upcoming single step (for all rxns):")
                        print("    Baseline: ", conc_array)
                        print("    Deltas:   ", delta_concentrations)
                        #print("    Relative Deltas:   ", delta_concentrations / conc_array)
                        #print("    L_inf norm:   ", np.linalg.norm(delta_concentrations, ord=np.inf) / delta_time_full)
                        # The following are normalized by the number of steps and the number of chemicals
                        print("    L2 norm:   ", np.linalg.norm(delta_concentrations) / (delta_time_full * n_chems))
                        print("    L1 norm:   ", np.linalg.norm(delta_concentrations, ord=1) / (delta_time_full * n_chems))

                else:                       # Using variable time steps
                    delta_concentrations = self._advance_variable_substeps(delta_time_full=delta_time_full, conc_array= conc_array,
                                                                           snapshots=snapshots,
                                                                           dynamic_substeps=dynamic_substeps, fast_threshold_fraction=fast_threshold_fraction)

                # Check whether the COMBINED delta_concentrations will make any conc negative
                if self.system is not None:     # IMPORTANT: usage of this function doesn't always involve self.system
                    tentative_updated_system = self.system + delta_concentrations
                    if min(tentative_updated_system) < 0:
                        print(f"******** CAUTION: negative concentration resulting from the combined effect of multiple reactions, "
                              f"upon advancing reactions from system time t={self.system_time:,.5g}\n"
                              f"         It will be AUTOMATICALLY CORRECTED with a reduction in the time step size")
                        raise ExcessiveTimeStep(f"The chosen time interval ({delta_time_full}) "
                                                f"leads to a NEGATIVE concentration of one of the chemicals: "
                                                f"MUST MAKE THE INTERVAL SMALLER!\n"
                                                f"[System Time: {self.system_time:.5g} : Baseline values: {self.system} ; delta conc's: {delta_concentrations}]")


                break       # IMPORTANT: this is needed because, in the absence of errors, we need to go thru the WHILE loop only once!

            except ExcessiveTimeStep as ex:
                # Single reactions steps can fail if the attempted time step was too large (leading to negative concentrations)
                if 1 in self.verbose_list:
                    print(ex)

                delta_time_full /= 2.       # Reduce the excessive time step in 1/2
                if 1 in self.verbose_list:
                    print(f"reaction_step_orchestrator(): RE-DOING THE LAST REACTION STEP with the smaller time interval of {delta_time_full}\n")


        if delta_concentrations is None:
            raise Exception("reaction_step_orchestrator(): unable to complete the reaction step")

        return  (delta_concentrations, delta_time_full)     # TODO: consider returning tentative_updated_system , since we already computed it



    def _advance_variable_substeps(self, delta_time_full: float, conc_array,
                                   snapshots, dynamic_substeps: int, fast_threshold_fraction) -> np.array:
        """
        This version is for when using adaptive VARIABLE TIME resolution.

        Using the given concentration data of ALL the chemical species,
        do the specified FULL TIME STEP for ALL the reactions - both those marked "slow" and those marked "fast" -
        but for the fast reaction, break down the full time interval into dynamic_substeps parts.

        All computations are based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions (in "forward Euler" approach.)

        Return the Numpy increment vector for ALL the chemical species concentrations, in their index order
        (whether involved in these reactions or not)

        NOTE: the actual system concentrations are NOT changed

        :param delta_time_full: The requested time duration of the overall reaction step,
                                    which will be subdivided for the "fast" reactions, if dynamic_substeps is greater than 1
        :param conc_array:      All initial concentrations at the start of the reaction step,
                                    as a Numpy array for all the chemical species, in their index order;
                                    this can be thought of as the "SYSTEM STATE"
        :param snapshots:       See explanation under single_compartment_react()
        :param dynamic_substeps:   An integer > 1 by which to subdivide the time intervals for "fast" reactions
        :param fast_threshold_fraction:  The minimum relative size of the concentration baseline over its change, AS A FRACTION,
                                            for a reaction to be regarded as "Slow".  IMPORTANT: this refers to the FULL step size

        :return:                The increment vector for the concentrations of ALL the chemical species,
                                    as a Numpy array for all the chemical species, in their index order
                                    EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):
                                        array([7. , -21.])
        """
        assert dynamic_substeps > 1, \
            "_advance_variable_substeps(): the function is being called in a scenario where FIXED time steps should be used"

        if self.are_all_slow_rxns():
            # If all the reactions are labeled as "slow"
            if 1 in self.verbose_list:
                print("    All the reactions are SLOW")
                print(f"    Processing ALL the {self.reaction_data.number_of_reactions()} reaction(s)")

            delta_concentrations = self._reaction_elemental_step(delta_time=delta_time_full, conc_array=conc_array, rxn_list=None,
                                                                 tag_reactions=True,
                                                                 time_subdivision=1, substep_number=0, fast_threshold_fraction=fast_threshold_fraction,
                                                                 )
            return delta_concentrations


        # If we get thus far, not all reactions are slow (i.e., some are fast)

        '''
            Process all the SLOW reactions first
        '''
        slow_rxns = self.slow_rxns()
        if 1 in self.verbose_list:
            if slow_rxns == []:
                print(f"    There are NO SLOW reactions")
            else:
                print(f"    * SLOW REACTIONS: {slow_rxns}")
                print(f"    Processing SLOW reactions")

        if slow_rxns == []:
            # If there are no slow reactions
            delta_concentrations_slow = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)  # One element per chemical species
        else:
            # Time-advance the simulation for all the slow reactions
            delta_concentrations_slow = self._reaction_elemental_step(delta_time=delta_time_full, conc_array=conc_array, rxn_list=slow_rxns,
                                                                      tag_reactions=True,
                                                                      time_subdivision=1, fast_threshold_fraction=fast_threshold_fraction)

        '''
            Next, process all the FAST reactions
        '''
        fast_rxns = self.fast_rxns()    # Look up which reactions are currently tagged as "fast"
        assert fast_rxns != [], "_advance_variable_substeps(): INTERNAL ERROR.  No fast reactions found"

        if 1 in self.verbose_list:
            print(f"    * FAST REACTIONS: {fast_rxns}")

        local_conc_array = conc_array.copy()    # Duplicate the array conc_array, to avoid altering it
                                                # local_conc_array is the best estimate of the System state at the start of each substep

        # TODO: maybe a better name for local_conc_array could be "system_state_estimate", and "conc_array" could be renamed "initial_system_state"

        delta_concentrations_fast = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)     # One element per chemical species
        reduced_time_step = delta_time_full / dynamic_substeps

        for substep in range(dynamic_substeps):
            # For each of the substeps into which delta_time_full has been divided up into
            if 1 in self.verbose_list:
                local_system_time = self.system_time + substep * reduced_time_step
                print(f"    - SUBSTEP: {substep} (in processing of FAST reactions, i.e. {fast_rxns})."
                      f"  'Local' system time: {local_system_time:,.4g}")

            if substep == dynamic_substeps-1:
                tag_reactions = True
            else:
                tag_reactions = False
                if 1 in self.verbose_list:
                    print(f"        Skipping evaluation of rxn speed for all reactions because NOT at last substep "
                          f"(substep_number = {substep}, time_subdivision = {dynamic_substeps})")

            incr_vector = self._reaction_elemental_step(delta_time=reduced_time_step, time_subdivision=dynamic_substeps,
                                                        fast_threshold_fraction=fast_threshold_fraction,
                                                        conc_array=local_conc_array, rxn_list=fast_rxns,
                                                        substep_number=substep, tag_reactions=tag_reactions)
            delta_concentrations_fast += incr_vector

            # TODO: the next 2 lines don't need to be run, if at the last step
            local_conc_array += incr_vector     # This is the contribution to the advance of the system state from the fast reactions
            # Also incorporate a scaled-down fraction of delta_concentrations_slow (IF there are any slow rxn's)
            if slow_rxns != []:
                local_conc_array += delta_concentrations_slow / dynamic_substeps


            # Preserve the intermediate-state data (the "secondary staging points"), if requested
            # NOTE: skipping the last substep, which will be handled by the caller function
            if snapshots and snapshots.get("show_intermediates") and substep < dynamic_substeps-1:
                species_to_show = snapshots.get("species")          # If not present, it will be None (meaning show all)
                time = self.system_time + (substep+1) * reduced_time_step
                self.add_snapshot(time=time, species=species_to_show,
                                  system_data=local_conc_array,
                                  caption=f"Interm. step, due to the fast rxns: {fast_rxns}")
                if self.diagnostics:
                    system_data = self.get_conc_dict(system_data=local_conc_array)
                    system_data["is_primary"] = False
                    system_data["primary_timestep"] = delta_time_full
                    system_data["n_substeps"] = dynamic_substeps   # This must be set, or the data type of the whole column changes to accomodate NaN's!
                    system_data["substep_number"] = substep+1
                    self.diagnostic_data_baselines.store(par=time,
                                                         data_snapshot=system_data)
        # END for

        # At this point, all the slow AND all the fast reactions have been processed

        # Combine the results of all the slow reactions and all the fast reactions
        if 2 in self.verbose_list:
            print("      delta_time: ", delta_time_full)
            print("          delta_concentrations_slow: ", delta_concentrations_slow)
            print("          delta_concentrations_fast: ", delta_concentrations_fast)

        delta_concentrations = delta_concentrations_slow + delta_concentrations_fast

        assert np.allclose(conc_array + delta_concentrations, local_conc_array), \
                    f"_advance_variable_substeps(): FAILED VALIDATION"         # TODO: eventually ditch this check

        return  delta_concentrations



    def _reaction_elemental_step(self, delta_time: float, conc_array: np.array, rxn_list=None, tag_reactions=False,
                                 time_subdivision=1, substep_number=0, fast_threshold_fraction=0.
                                 ) -> np.array:
        """
        Using the given concentration data of ALL the chemical species,
        do the specified SINGLE TIME STEP for ONLY the requested reactions (by default all).

        NOTE: this time step might be
              either the "full" time step of the high-level calling function,
              or a substep (if adaptive variable time resolution is being used)

        All computations are based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions (in "forward Euler" approach.)

        Return the Numpy increment vector for ALL the chemical species concentrations, in their index order
        (whether involved in these reactions or not)

        NOTES:  - the actual System Concentrations and the System Time are NOT changed
                - if any of the concentrations go negative, an Exception is raised

        :param delta_time:      The time duration of this individual reaction step - assumed to be small enough that the
                                    concentration won't vary significantly during this span.
                                    NOTE: this may be a "full step" or a "substep", depending on the adaptive time scale
                                          being used by the caller function
        :param conc_array:      All initial concentrations at the start of the reaction step,
                                    as a Numpy array for all the chemical species, in their index order;
                                    this can be thought of (as a copy of) the "SYSTEM STATE"
        :param rxn_list:        OPTIONAL list of reactions (specified by their indices) to include in this simulation step ;
                                    EXAMPLE: [1, 3, 7]
                                    If None, do all the reactions
        :param tag_reactions:   OPTIONAL boolean indicating whether to consider a pre-determined "measure of change" in the concentrations
                                    caused by each reaction, and store the results with the reactions

        [ALL THE REMAINING ARGUMENTS ARE SPECIFIC TO INVOCATIONS USING THE VARIABLE TIME RESOLUTION]
        :param time_subdivision:        Integer with the number of subdivisions currently being used for the "full" time step
                                            of the caller function;
                                            if > 1, it means we're using an adaptive variable time resolution.
                                            (Note: the "full" time step of the calling function will be delta_time * time_subdivision)
        :param fast_threshold_fraction: The minimum relative size of the concentration baseline over its change, AS A FRACTION,
                                            for a reaction to be regarded as "Slow".  IMPORTANT: this refers to the FULL step size
        :param substep_number:          Zero-based counting of the time substeps

        :return:            The increment vector for the concentrations of ALL the chemical species
                                (whether involved in the reactions or not),
                                as a Numpy array for all the chemical species, in their index order
                            EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):   array([7. , -21.])
        """
        # The increment vector is cumulative for ALL the requested reactions
        increment_vector = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)       # One element per chemical species

        # Compute the forward and reverse "conversions" of all the applicable reactions
        delta_dict = self.compute_all_reaction_deltas(conc_array=conc_array, delta_time=delta_time, rxn_list=rxn_list)
        if 3 in self.verbose_list:
            print(f"      delta_list: {delta_dict}")


        if rxn_list is None:    # Meaning ALL reactions
            rxn_list = range(self.reaction_data.number_of_reactions())  # This will be a list of all the reaction index numbers


        # For each applicable reaction, find the needed adjustments ("deltas")
        #   to the concentrations of the reactants and products,
        #   based on the forward and reverse rates of the reaction
        for rxn_index in rxn_list:      # Consider each reaction in turn
            # TODO: maybe switch to a call to the experimental _reaction_elemental_step_SINGLE_REACTION()

            if 2 in self.verbose_list:
                print(f"      Determining the conc.'s changes as a result of rxn # {rxn_index}")

            # One element per chemical species; notice that this array is being RESET for EACH REACTION
            # TODO: instead of using a potentially very large array (mostly of zeros) for each rxn, consider a dict instead
            #       then combine then at end
            increment_vector_single_rxn = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)

            # TODO: turn into a more efficient single step, as as:
            #(reactants, products) = cls.all_reactions.unpack_terms(rxn_index)
            reactants = self.reaction_data.get_reactants(rxn_index)
            products = self.reaction_data.get_products(rxn_index)


            # If requested to evaluate the reaction "speeds"
            if tag_reactions:
                self.set_rxn_speed(rxn_index, "S")  # TENTATIVE assignment, that will be changed
                                                    #   if ANY chemical experiences significant concentration changes
                if 1 in self.verbose_list:
                    print(f"      (tentatively tagging rxn #{rxn_index} as 'S')")


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
                self.validate_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                increment_vector_single_rxn[species_index] += delta_conc


            # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
            for p in products:
                stoichiometry, species_index, order = p             # Unpack
                delta_conc = stoichiometry * delta_dict[rxn_index]  # Increment to this reaction product from the reaction being considered
                # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
                # Note: not enough to detect conc going negative from combined changes from multiple reactions!  Further testing done upstream
                self.validate_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                increment_vector_single_rxn[species_index] += delta_conc


            increment_vector += increment_vector_single_rxn     # Accumulate the increment vector across all reactions
                                                                # TODO: consider using a small dict in lieu of increment_vector_single_rxn


            # If requested to evaluate the reaction "speeds"
            if tag_reactions:
                # Mark the reaction "fast", if appropriate
                self.examine_increment_array(rxn_index=rxn_index,
                                             delta_conc_array=increment_vector_single_rxn, baseline_conc_array=conc_array,
                                             time_subdivision=time_subdivision,
                                             fast_threshold_fraction=fast_threshold_fraction)

            if self.diagnostics:
                self.save_diagnostic_data(delta_time, increment_vector_single_rxn, rxn_index, substep_number, time_subdivision)
        # END for (over rxn_list)

        return increment_vector




    def _reaction_elemental_step_SINGLE_REACTION(self, delta_time: float, conc_array: np.array, increment_vector,
                                                 rxn_index: int, tag_reactions :bool,
                                                 delta_dict):     # TODO: EXPERIMENTAL; not yet in use

        if 2 in self.verbose_list:
            print(f"      Determining the conc.'s changes as a result of rxn # {rxn_index}")

        # NOTE: instead of using a potentially very large array (mostly of zeros) for each rxn,
        #       we use a dict instead; then combine all at end
        #increment_vector_single_rxn = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)
        increment_dict_single_rxn = {}      # The key will be the elements in the reaction

        # TODO: turn into a more efficient single step, as as:
        #(reactants, products) = cls.all_reactions.unpack_terms(rxn_index)
        reactants = self.reaction_data.get_reactants(rxn_index)
        products = self.reaction_data.get_products(rxn_index)


        # If requested to evaluate the reaction "speeds"    TODO: maybe move to calling function
        if tag_reactions:
            self.set_rxn_speed(rxn_index, "S")  # TENTATIVE assignment, that will be changed
            #   if ANY chemical experiences significant concentration changes
            if 1 in self.verbose_list:
                print(f"      (tentatively tagging rxn #{rxn_index} as 'S')")


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
            # TODO: it might be more effecient to check this in bulk than with many function calls
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
        if tag_reactions:
            # Mark the reaction "fast", if appropriate
            self.examine_increment_array(rxn_index=rxn_index,
                                         delta_conc_array=increment_vector_single_rxn, baseline_conc_array=conc_array,
                                         time_subdivision=time_subdivision,
                                         fast_threshold_fraction=fast_threshold_fraction)

        if self.diagnostics:
            self.save_diagnostic_data(delta_time, increment_vector_single_rxn, rxn_index, substep_number, time_subdivision)
        '''




    def save_diagnostic_data(self, delta_time, increment_vector_single_rxn, rxn_index, substep_number, time_subdivision) -> None:
        """
        Add detailed data from current step-in-progress of the reaction simulation to the self.diagnostic_data attribute

        :param delta_time:
        :param increment_vector_single_rxn:
        :param rxn_index:
        :param substep_number:
        :param time_subdivision:
        :return:                            None
        """
        if self.diagnostic_data == {}:    # INITIALIZE the dictionary self.diagnostic_data, if needed
            for i in range(self.reaction_data.number_of_reactions()):
                self.diagnostic_data[i] = MovieTabular(parameter_name="TIME")       # One per reaction

        data_snapshot = {}
        # Add more entries to the above dictionary, starting with the Delta conc. for all the chemicals
        for index, conc in enumerate(increment_vector_single_rxn):
            data_snapshot["Delta " + self.reaction_data.get_name(index)] = conc

        data_snapshot["reaction"] = rxn_index           # TODO: now redundant because factored out into separate data frames
        data_snapshot["substep"] = substep_number
        data_snapshot["time_subdivision"] = time_subdivision
        data_snapshot["delta_time"] = delta_time
        local_system_time = self.system_time + substep_number * delta_time
        self.diagnostic_data[rxn_index].store(par=local_system_time,
                                              data_snapshot=data_snapshot)



    def criterion_fast_reaction(self, delta_conc, time_subdivision, fast_threshold_fraction,
                                baseline_conc=None) -> bool:
        """
        Apply a criterion to determine, from the given data,
        whether the originating reaction (the source of the data) needs to be classified as "Fast".
        All the passed data is for the concentration changes in 1 chemical from 1 reaction

        :param delta_conc:
        :param time_subdivision:
        :param fast_threshold_fraction:
        :param baseline_conc:           # TODO: probably phase out

        :return:                        True if the concentration change is so large (based on some criteria)
                                            that the reaction that caused it, ought to be regarded as "fast"
        """
        if self.fast_criterion_use_baseline:    # TODO: this criterion gives poor results, and will probably be eliminated
            # if abs(delta_conc) / baseline_conc > fast_threshold_fraction / time_subdivision
            # Perhaps more intuitive as shown above;
            # but, to avoid time-consuming divisions (and potential divisions by zero), re-formulated as below:
            return abs(delta_conc) * time_subdivision > fast_threshold_fraction * baseline_conc

        else:
            # Perhaps more intuitive written as:  if abs(delta_conc) > fast_threshold_fraction / time_subdivision
            return abs(delta_conc) * time_subdivision > fast_threshold_fraction



    def validate_increment(self,  delta_conc, baseline_conc: float,
                           rxn_index: int, species_index: int, delta_time) -> None:
        """
        Examine the requested concentration change given by delta_conc
        (typically, as computed by an ODE solver),
        relative to the baseline (pre-reaction) value baseline_conc,
        for the given SINGLE chemical species and SINGLE reaction.

        If the concentration change would render the concentration negative,
        raise an Exception (of custom type "ExcessiveTimeStep")

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
            print(f"******** CAUTION: negative concentration in chemical `{self.reaction_data.get_name(species_index)}`"
                  f" (resulting from reaction {self.reaction_data.single_reaction_describe(rxn_index=rxn_index, concise=True)})\n",
                  f"        upon advancing from system time t={self.system_time:,.5g} [Baseline value: {baseline_conc:.5g} ; delta conc: {delta_conc:.5g}]\n"
                  f"         It will be AUTOMATICALLY CORRECTED with a reduction in the time step size")
            raise ExcessiveTimeStep(f"The chosen time interval ({delta_time}) "
                                    f"leads to a NEGATIVE concentration of `{self.reaction_data.get_name(species_index)}` (index {species_index}) "
                                    f"from reaction {self.reaction_data.single_reaction_describe(rxn_index=rxn_index, concise=True)} (rxn # {rxn_index}): "
                                    f"MUST MAKE THE INTERVAL SMALLER!\n"
                                    f"[System Time: {self.system_time:.5g} : Baseline value: {baseline_conc:.5g} ; delta conc: {delta_conc:.5g}]")



    def examine_increment_array(self, rxn_index: int,
                                delta_conc_array: np.array,
                                time_subdivision: int,
                                fast_threshold_fraction,
                                baseline_conc_array=None
                                ) -> None:
        """
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
                                            before the last simulated reaction step or substep  TODO: phase out
        :param time_subdivision:        Integer with the number of subdivisions currently used for delta_time;
                                            used for adaptive variable time resolution
        :param fast_threshold_fraction: The minimum relative size of the concentration baseline over its change, AS A FRACTION,
                                            for a reaction to be regarded as "Slow".
                                            IMPORTANT: this refers to the FULL step size, and will be scaled for time_subdivision
                                            (e.g. a time_subdivision of 2 signifies 1/2 a time step,
                                            and will mean that just 1/2 of the change will constitute a threshold)

        :return:                        None (the equation is marked as "Fast", if appropriate, in its data structure)
        """

        if self.get_rxn_speed(rxn_index) == "F":
            # If the reaction was already marked as "Fast" then no action is taken
            if 1 in self.verbose_list:
                print(f"        Rxn # {rxn_index}: no action taken because already tagged as 'F'")

            return


        # If we get this far, the reaction came tagged as "SLOW"

        # If the reaction was (tentatively) labeled as "Slow", decide whether to flip it to "Fast"
        #   Note: the status will be used for the current simulation sub-cycle

        if 1 in self.verbose_list:
            print(f"        For Rxn # {rxn_index}, checking concentrations of chems: {self.reaction_data.get_chemicals_in_reaction(rxn_index)}")


        # Loop over just the chemicals in the given reaction (not over ALL the chemicals!)
        for i in self.reaction_data.get_chemicals_in_reaction(rxn_index):

            delta_conc = delta_conc_array[i]
            if baseline_conc_array is None:
                baseline_conc = None
            else:
                baseline_conc = baseline_conc_array[i]

            # Determine whether the reaction needs to be classified as "Fast"
            if self.criterion_fast_reaction(delta_conc=delta_conc, baseline_conc=baseline_conc,
                                            time_subdivision=time_subdivision, fast_threshold_fraction=fast_threshold_fraction):
                self.set_rxn_speed(rxn_index, "F")
                if 1 in self.verbose_list:
                    print(f"        Rxn # {rxn_index} FLIPPED TAG to 'F', "
                          f"based on a change of {delta_conc:.5g} (rel. to baseline of {baseline_conc:.6g}) "
                          f"for the conc. of `{self.reaction_data.get_name(species_index=i)}`: "
                          f"abs({100 * delta_conc / (baseline_conc+1e-09):.3g}%) > ({100 * fast_threshold_fraction:.3g}% /{time_subdivision})"
                          )   # Note: we're adding a tiny amount to baseline_conc to avoid potential divisions by zero
                return

        # END for
        # If we get thus far, the reaction is still regarded as "Slow", and its status is left unchanged

        if 1 in self.verbose_list:
            chem_list = self.reaction_data.get_chemicals_in_reaction(rxn_index)
            delta_conc_array_relevant = np.array([delta_conc_array[n] for n in chem_list])
            baseline_conc_array_relevant = np.array([baseline_conc_array[n] for n in chem_list])
            print(f"        Rxn # {rxn_index} left tagged as 'S', "
                  f"based on a change of {delta_conc_array_relevant} (rel. to baseline of {baseline_conc_array_relevant}) :\n"
                  f"        elements of abs({100 * delta_conc_array_relevant / (baseline_conc_array_relevant+1e-09)}%) are all < ({100 * fast_threshold_fraction:.3g}% /{time_subdivision})"
                  )   # Note: we're adding a tiny amount to baseline_conc to avoid potential divisions by zero



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
            rxn_list = range(self.reaction_data.number_of_reactions())

        # Process the requested reactions
        for i in rxn_list:      # Consider each reaction in turn
            delta = self.compute_reaction_delta(rxn_index=i, conc_array=conc_array, delta_time=delta_time)
            delta_dict[i] = delta

        return delta_dict



    def compute_reaction_delta(self, rxn_index: int, delta_time: float, conc_array: np.array) -> float:
        """
        For the SINGLE specified reaction, the given time interval, and the specified concentrations of chemicals,
        compute the difference of the reaction's forward and back "conversions",
        a non-standard term we're using here to refer to delta_time * (Forward_Rate − Reverse_Rate)

        For background info: https://life123.science/reactions
        What we're computing here, is referred to as:  (Δt)∗delta_forward(n)

        :param rxn_index:   An integer that indexes the reaction of interest (numbering starts at 0)
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

        # TODO: turn into a more efficient single step, as as:
        #(reactants, products, fwd_rate_constant, rev_rate_constant) = cls.all_reactions.unpack_reaction(i)
        reactants = self.reaction_data.get_reactants(rxn_index)
        products = self.reaction_data.get_products(rxn_index)
        fwd_rate_constant = self.reaction_data.get_forward_rate(rxn_index)
        rev_rate_constant = self.reaction_data.get_reverse_rate(rxn_index)

        forward_rate = fwd_rate_constant
        for r in reactants:
            stoichiometry, species_index, order = r     # Unpack the data of the reactant r
            conc = conc_array[species_index]
            #assert conc is not None, \
            #f"ReactionDynamics.compute_rate_delta(): lacking the value for the concentration of the chemical species `{self.reaction_data.get_name(species_index)}`"
            forward_rate *= conc ** order      # Raise to power

        reverse_rate = rev_rate_constant
        for p in products:
            stoichiometry, species_index, order = p     # Unpack the data of the reaction product p
            conc = conc_array[species_index]
            #assert conc is not None, \
            #f"ReactionDynamics.compute_rate_delta(): lacking the concentration value for the species `{self.reaction_data.get_name(species_index)}`"
            reverse_rate *= conc ** order     # Raise to power

        return delta_time * (forward_rate - reverse_rate)



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
                description = self.reaction_data.single_reaction_describe(rxn_index=rxn_index, concise=True)
                print(description)

            status = self.reaction_in_equilibrium(rxn_index=rxn_index, conc=conc, tolerance=tolerance, explain=explain)
            if not status:
                failures_dict = {False: [rxn_index]}

        else:
            # Check all the reactions
            status = True   # Overall status
            description_list = self.reaction_data.multiple_reactions_describe(concise=True)
            for i in range(self.reaction_data.number_of_reactions()):
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
        rxn = self.reaction_data.get_reaction(rxn_index)

        reactants = self.reaction_data.extract_reactants(rxn)     # A list of triplets
        products = self.reaction_data.extract_products(rxn)      # A list of triplets
        kF = self.reaction_data.extract_forward_rate(rxn)
        kB = self.reaction_data.extract_back_rate(rxn)

        rate_ratio = kF / kB                    # Ratio of forward/reverse reaction rates

        conc_ratio = 1.
        numerator = ""
        denominator = ""
        all_concs = []      # List of strings

        for p in products:
            # Loop over the reaction products
            species_index = self.reaction_data.extract_species_index(p)
            rxn_order = self.reaction_data.extract_rxn_order(p)

            species_name = self.reaction_data.get_name(species_index)
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
            species_index = self.reaction_data.extract_species_index(r)
            rxn_order = self.reaction_data.extract_rxn_order(r)

            species_name = self.reaction_data.get_name(species_index)
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




    #############################################################################################
    #                                                                                           #
    #                                      FOR DIAGNOSTICS                                      #
    #                                                                                           #
    #############################################################################################
    def _________FOR_DIAGNOSTICS___________(DIVIDER):  # Used to get a better structure listing in IDEs such asPycharm
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
        :param delta_arr:   Numpy array with the concentrations changes of ALL the chemicals (whether involved
                                in the reaction or not), in their index order, as a result of JUST the reaction of interest
                            TODO: maybe also accept a Panda's data frame row
        :param suppress_warning:
        :return:            True if the change in reactant/product concentrations is consistent with the
                                reaction's stoichiometry, or False otherwise
        """
        if (not suppress_warning) and (self.reaction_data.number_of_reactions() > 1):
            print(f"*** WARNING: {self.reaction_data.number_of_reactions()} reactions are present.  "
                  f"stoichiometry_checker() currently only works for 1-reaction simulations")

        self.reaction_data.assert_valid_rxn_index(rxn_index)

        reactants = self.reaction_data.get_reactants(rxn_index)
        products = self.reaction_data.get_products(rxn_index)

        # Pick (arbitrarily) the first reactant,
        # to establish a baseline change in concentration relative to its stoichiometric coefficient
        baseline_term = reactants[0]
        baseline_species = self.reaction_data.extract_species_index(baseline_term)
        baseline_stoichiometry = self.reaction_data.extract_stoichiometry(baseline_term)
        baseline_ratio =  (delta_arr[baseline_species]) / baseline_stoichiometry
        #print("baseline_ratio: ", baseline_ratio)

        for i, term in enumerate(reactants):
            if i != 0:
                species = self.reaction_data.extract_species_index(term)
                stoichiometry = self.reaction_data.extract_stoichiometry(term)
                ratio =  (delta_arr[species]) / stoichiometry
                #print(f"ratio for `{self.reaction_data.get_name(species)}`: {ratio}")
                if not np.allclose(ratio, baseline_ratio):
                    return False

        for term in products:
            species = self.reaction_data.extract_species_index(term)
            stoichiometry = self.reaction_data.extract_stoichiometry(term)
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
        if self.diagnostic_data == {}:
            print("WARNING *** In order to run stoichiometry_checker_entire_run(), "
                  "the diagnostics must be turned on, with set_diagnostics(), prior to the simulation run!")
            return False

        for rxn_index in range(self.reaction_data.number_of_reactions()):
            diagnostic_data = self.get_diagnostic_data(rxn_index=rxn_index, print_reaction=False)
            for row_index in range(len(diagnostic_data)):
                df_row = diagnostic_data.loc[row_index]     # Row in the Panda's data frame of diagnostic data
                chemical_delta_list = self._delta_names()            # EXAMPLE: ["Delta A", "Delta B", "Delta C"]
                delta = df_row[chemical_delta_list].to_numpy()      # Extract select columns from the data frame row, and turn into Numpy array
                status = self.stoichiometry_checker_from_deltas(rxn_index=rxn_index, delta_arr=delta, suppress_warning=True)
                if not status:
                    print(f"Stoichiometry NOT satisfied by reaction # {rxn_index}: "
                          f"see row # {row_index} in the diagnostic data for that reaction")
                    return False

        return True



    def _delta_names(self) -> [str]:
        """
        Return a list of strings, with names of all the chemicals,
        each prefixed by the string "Delta "
        EXAMPLE: ["Delta A", "Delta B", "Delta C"]
        :return:
        """
        chemical_list = self.reaction_data.get_all_names()      # EXAMPLE: ["A", "B", "C"]
        chemical_delta_list = ["Delta " + name
                               for name in chemical_list]
        return chemical_delta_list



    def examine_run(self, df, time_step, fast_threshold) -> None:
        """
        Analyze, primarily for debugging purposes, the data produced by a run of single_compartment_react()
        The primary focus is for diagnostic information of the adaptive variable time steps.

        CAUTION: this approach is only usable with single-reaction runs or non-overlapping multiple reactions,
                 because it looks into concentration changes of the various chemicals,
                 which in the case of multiple coupled reactions can be attributed to multiple reactions

        :param df:              A Pandas data frame, as created by single_compartment_react()
        :param time_step:       The same value (for the FULL STEP size) that was used in single_compartment_react()
        :param fast_threshold:  The same value (for the FULL STEP size) that was used in single_compartment_react()
        :return:                None
        """
        print("WARNING: The explanations below are based on an older system of using RELATIVE concentration changes, and no longer applicable to the newer approach; "
              "it'll be corrected in future versions (PLEASE DIS-REGARD FOR NOW)\n\n")

        number_datapoints = len(df)
        chemical_list = self.reaction_data.get_all_names()
        print("Examining the reaction data for the chemicals: ", chemical_list)

        if self.reaction_data.number_of_reactions() > 1:
            print(f"*** WARNING: {self.reaction_data.number_of_reactions()} reactions are present.  "
                  f"The analysis below won't be meaningful if any of the reactions are coupled!")

        for i in range(number_datapoints - 1):   # The -1 is b/c this analysis needs to look ahead to following data point
            sys_time = df.iloc[i]['SYSTEM TIME']
            print(f"\n------ Reaction step OR substep {i} (Initial SYSTEM TIME: {sys_time:,.4g}) ------")

            before = df.iloc[i][chemical_list].to_numpy()       # Concentrations at the start of this step/substep
            print("Baseline conc's at start: ", before)

            after = df.iloc[i+1][chemical_list].to_numpy()      # Concentrations at the end of this step/substep

            rel_changes = (after - before) / before * 100.      # Relative conc. changes (in %)
            largest_abs_rel_change = max(abs(rel_changes))      # Max of absolute values of % relative conc. changes

            duration = df.iloc[i+1]['SYSTEM TIME'] - df.iloc[i]['SYSTEM TIME']  # Duration of this step or substep
            divisor = round(time_step / duration)               # An integer, indicating how many of current step or substep
                                                                #   will fit in the full time interval of the simulation

            next_time = df.iloc[i+1]['SYSTEM TIME']             # Time of the next step/substep.  This is used to determine
                                                                #   if we are at the last substep of a full time interval
            time_ratio = next_time / time_step                  # If this is an integer, it means that a full time interval is coming to an end
            multiple_of_time_step = np.allclose(time_ratio, round(time_ratio))  # A Boolean that will be True iff time_ratio is an integer

            if (divisor > 1) and (not multiple_of_time_step):
                print("This is an INTERMEDIATE SUBSTEP. Therefore, no evaluations of reaction speeds")
            else:
                if divisor == 1:
                    print("This is a FULL STEP")
                else:
                    print("This is a FINAL SUBSTEP")

                print("Final conc's: ", after)
                print("Delta conc's: ", after - before)
                print("Relative conc. changes (%): ", rel_changes)
                print(f"Largest MAX rel. change in abs. value (%):  {largest_abs_rel_change:.3g}")
                print("Time subdivisions: ", divisor)

                if largest_abs_rel_change > fast_threshold / divisor:
                    print(f"Since {largest_abs_rel_change:.3g} > ({fast_threshold}% over {divisor}), the reaction is classified as *FAST*")
                else:
                    print(f"Since {largest_abs_rel_change:.3g} < ({fast_threshold}% over {divisor}), the reaction is classified as *Slow*")



    def set_diagnostics(self):
        self.diagnostics = True

    def unset_diagnostics(self):
        self.diagnostics = False



    def get_diagnostic_data(self, rxn_index: int, tail=None, print_reaction=True) -> pd.DataFrame:
        """
        Return a Pandas data frame with the diagnostic data of the requested reaction.
        Also, optionally print out a brief description of the reaction.
        Optionally, restrict the result with a start and/or end times,
        or by limiting to a specified numbers of rows at the end

        :param rxn_index:   An integer that indexes the reaction of interest (numbering starts at 0)
        :param tail:        (OPTIONAL) Number of records to consider, from the end of the diagnostic dataframe
        :param print_reaction:
        :return:            A Pandas data frame with (all or some of) the diagnostic data of the above reaction
        """
        # First, validate the reaction index
        self.reaction_data.assert_valid_rxn_index(rxn_index)

        if print_reaction:
            print("Reaction: ", self.reaction_data.single_reaction_describe(rxn_index=rxn_index, concise=True))

        return self.diagnostic_data[rxn_index].get(tail=tail)



    def diagnose_variable_time_steps(self, rel_fast_threshold) -> None:
        """
        Analyze, primarily for debugging purposes, the diagnostics data produced
        when the self.diagnostics attribute is set to True prior to running single_compartment_react().

        The primary focus is for diagnostic information of the adaptive variable time steps.

        This approach will eventually be suitable for any
        type of reaction simulations.

        TODO: CAUTION - only meant for 1 reaction.  Generalize to multiple reactions.
              For more than 1 reaction, use explain_reactions() instead

        :param rel_fast_threshold:  The same value that was used in single_compartment_react()

        :return:                None
        """
        print("WARNING: The explanations below are based on an older system of using RELATIVE concentration changes, and no longer applicable to the newer approach; "
              "it'll be corrected in future versions (PLEASE DIS-REGARD FOR NOW)\n\n")

        if not self.diagnostics:
            print("No diagnostic information is available:\nIn order to run diagnose_variable_time_steps(), "
                  "call set_diagnostics() prior to running single_compartment_react()")
            return

        assert self.reaction_data.number_of_reactions() == 1, \
            "diagnose_variable_time_steps() currently ONLY works when exactly 1 reaction is present. " \
            "For more than 1 reaction, use explain_reactions() instead"

        #diagnostic_df = self.diagnostic_data.get()  # TODO: self.diagnostic_data[rxn_index].get()

        for rxn_index in range(self.reaction_data.number_of_reactions()):
            print(f"\nDiagnostics for reaction {rxn_index}")
            diagnostic_df = self.diagnostic_data[rxn_index].get()
            number_diagnostic_points = len(diagnostic_df)

            baselines_df = self.diagnostic_data_baselines.get()

            assert len(baselines_df) == number_diagnostic_points + 1, \
                f"diagnose_variable_time_steps(): error in relative lengths of diagnostic changes Pandas frame ({number_diagnostic_points}) " \
                f"and diagnostic baselines Pandas frame ({len(baselines_df)}); the latter should be 1 longer than the former"


            chemical_list = self.reaction_data.get_all_names()
            print("Examining the reactions' diagnostic data for all the chemicals")
            chemical_delta_list = self._delta_names()
            print("The concentration increments are: ", chemical_delta_list)

            for i in range(number_diagnostic_points):
                print(f"\n---- {i} (Reaction step OR substep) ----")
                debug_time = diagnostic_df.iloc[i]['TIME']  # Note: each entry in the diagnostic_df data frame has the conc. increments
                #       for full step or substep that starts at this time
                print(f"Start time: {debug_time:.5g} (Start of Full time interval or sub-interval)")

                time_subdivision = diagnostic_df.iloc[i]['time_subdivision']
                print(f"time_subdivision: {time_subdivision}")

                baseline = baselines_df.iloc[i][chemical_list].to_numpy()
                print("Baseline conc's:", baseline)


                substep = diagnostic_df.iloc[i]['substep']
                if substep < time_subdivision - 1:
                    print("This is an INTERMEDIATE SUBSTEP. Therefore, no evaluations of reaction speeds")
                    continue

                delta = diagnostic_df.iloc[i][chemical_delta_list].to_numpy()
                print("Delta conc's:", delta)

                ratio = delta / baseline * 100.
                print("Ratios (%):", ratio)
                print(f"Max abs ratio (%): {max(abs(ratio)):.3g}")
                print(f"Comparing the above value against {rel_fast_threshold / time_subdivision}% (i.e. {rel_fast_threshold}% /{time_subdivision})")
                if max(abs(ratio)) > rel_fast_threshold/time_subdivision:
                    print("*FAST* reaction")
                else:
                    print("*Slow* reaction")



    def explain_reactions(self) -> bool:
        """
        Provide a detailed explanation of all the steps/substeps of the reactions,
        from the saved diagnostic data

        WARNING: Currently designed only for exactly 2 reactions!  TODO: generalize to any number of reactions

        TODO: allow arguments to specify the min and max reaction times during which to display the explanatory data

        :return:    True is the diagnostic data is consistent for all the steps/substeps of all the reactions,
                    or False otherwise
        """
        if not self.diagnostics:
            print("No diagnostic information is available:\nIn order to run diagnose_variable_time_steps(), "
                  "call set_diagnostics() prior to running single_compartment_react()")
            return False

        number_reactions = self.reaction_data.number_of_reactions()

        assert number_reactions == 2, \
            "explain_reactions() currently ONLY works when exactly 2 reactions are present. " \
            "Future versions will lift this restriction"

        row_baseline = 0
        row_list = [0, 0]       # TODO: generalize
        active_list = [0, 1]    # ALL the reactions.  TODO: generalize

        self._explain_reactions_helper(active_list=active_list,
                                       row_baseline=row_baseline, row_list=row_list)


        while row_baseline < len(self.diagnostic_data_baselines.get()) - 2:
            row_baseline += 1
            #if self.get_diagnostic_data(rxn_index=0).loc[row_list[0]]['time_subdivision'] == self.get_diagnostic_data(rxn_index=1).loc[row_list[1]]['time_subdivision']:
            if active_list == [0, 1]:    # ALL the reactions  TODO: generalize
                print("Normal advance (single-step across all tables)")
            else:
                print("Advance a step in the tables for the following reactions: ", active_list)

            for i in active_list:
                row_list[i] += 1

            # TODO: generalize
            time_sub0 = self.get_diagnostic_data(rxn_index=0).loc[row_list[0]]['time_subdivision']
            time_sub1 = self.get_diagnostic_data(rxn_index=1).loc[row_list[1]]['time_subdivision']
            print("time subdivisions: ", time_sub0, time_sub1)

            # TODO: generalize
            substep0 = self.get_diagnostic_data(rxn_index=0).loc[row_list[0]]['substep']
            substep1 = self.get_diagnostic_data(rxn_index=1).loc[row_list[1]]['substep']
            print("substeps: ", substep0, substep1)


            last_substep = False
            if active_list == [0]:  # TODO: generalize
                if substep0 == time_sub0 - 1:
                    print("At LAST SUBSTEP of rxn 0")
                    last_substep = True
                else:
                    print("Not at last substep of rxn 0")
            elif active_list == [1]:
                if substep1 == time_sub1 - 1:
                    print("At LAST SUBSTEP of rxn 1")
                    last_substep = True
                else:
                    print("Not at last substep of rxn 1")
            else:   # The active_list is ALL rxn's
                if time_sub0 > time_sub1:
                    active_list = [0]   # TODO: generalize
                    print("CHANGING active_list (rxns to advance in substeps) to: ", active_list)
                elif time_sub1 > time_sub0:
                    active_list = [1]
                    print("CHANGING active_list (rxns to advance in substeps) to: ", active_list)
                else:
                    print("Will be advancing all reactions in lockstep")


            status = self._explain_reactions_helper(active_list=active_list, row_baseline=row_baseline, row_list=row_list,
                                                    time_subdivision=max(time_sub0, time_sub1)) # TODO: generalize

            if last_substep:
                print("Changing active list to ALL reactions")
                active_list = [0,1]     # TODO: generalize


            if not status:
                return False    # Error termination

        # END while

        return True             # Successful termination



    def _explain_reactions_helper(self, active_list, row_baseline, row_list, time_subdivision=None) -> bool:
        """
        Helper function for explain_reactions()

        :param active_list:
        :param row_baseline:
        :param row_list:
        :param time_subdivision:
        :return:                True is the diagnostic data is consistent for this (sub)step, or False otherwise
        """
        print("-----------")
        print("ROW of baseline data: ", row_baseline)

        current_time = self.diagnostic_data_baselines.get().loc[row_baseline]['TIME']
        print(f"TIME = {current_time:.5g}")

        print("row_list: ", row_list)
        print("active_list: ", active_list)

        chemical_list = self.reaction_data.get_all_names()
        chemical_delta_list = self._delta_names()

        conc_arr_before = self.diagnostic_data_baselines.get().loc[row_baseline][chemical_list].to_numpy().astype(float)
        print("baseline concentrations: ", conc_arr_before)

        delta_cumulative = np.zeros(self.reaction_data.number_of_chemicals(),
                                    dtype=float)  # One element per chemical species

        # For each reaction
        for rxn_index in range(self.reaction_data.number_of_reactions()):
            if (rxn_index in active_list):  # If reaction is tagged as "fast"
                row = row_list[rxn_index]   # Row in the data frame for the diagnostic data on this reaction
                delta_rxn = self.get_diagnostic_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy().astype(float)
                print(f"From fast rxn {rxn_index}: delta_rxn = {delta_rxn}")
            else:                           # If reaction is tagged as "slow"
                row = row_list[rxn_index]   # Row in the data frame for the diagnostic data on this reaction
                #time_subdivision = self.get_diagnostic_data(rxn_index=rxn_index).loc[row]['time_subdivision']
                delta_rxn = self.get_diagnostic_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy().astype(float) / time_subdivision
                #delta_rxn = np.zeros(self.reaction_data.number_of_chemicals(),
                #                     dtype=float)   # This is the way it was in release beta 19
                print(f"From slow rxn {rxn_index}: delta_rxn = {delta_rxn} (using time_subdivision {time_subdivision})")


            delta_cumulative += delta_rxn
        # END for

        print("delta_cumulative: ", delta_cumulative)

        conc_after = conc_arr_before + delta_cumulative
        print("updated concentrations: ", conc_after)

        next_system_state = self.diagnostic_data_baselines.get().loc[row_baseline+1][chemical_list].to_numpy()
        print(f"concentrations from the next row ({row_baseline + 1}) of the system state: ", next_system_state)

        status = np.allclose(conc_after.astype(float), next_system_state.astype(float))
        if status:
            print("Match OK")
        else:
            print("****************************************   MISMATCH!!!  ****************************************")

        print("-----------")

        return status



    def explain_time_advance(self, return_times=False) -> Union[None, list]:
        """
        Use the saved-up diagnostic data, to print out details of the timescales of the reaction run

        If diagnostics weren't enabled ahead of calling this function, an Exception is raised

        EXAMPLE of output:
            From time 0 to 0.0168, in 42 substeps of 0.0004 (each 1/2 of full step)
            From time 0.0168 to 0.0304, in 17 FULL steps of 0.0008
            (for a grand total of 38 FULL steps)

        :param return_times:    If True, all the critical times are saved and returned as a list
        :return:                Either None, or a list of time values
        """
        assert self.diagnostics, "ReactionDynamics.explain_time_advance(): diagnostics must first be enabled; " \
                                 "use set_diagnostics() prior to the reaction run"

        df = self.diagnostic_data_baselines.get()
        assert df is not None, \
            "ReactionDynamics.explain_time_advance(): no diagnostic data found.  " \
            "Did you call set_diagnostics() prior to the reaction run?"

        n_entries = len(df)
        t = list(df["TIME"])    # List of times (simulation step points)

        if n_entries == 1:
            print(f"ReactionDynamics.explain_time_advance(): no time advance found. "
                  f"Diagnostics only contain an initial System Time of {t[0]:.3g} ;  "
                  f"did you run the reaction simulation?")
            if return_times:
                return t
            else:
                return


        # If we get thus far, we have at least 2 entries in the list "t" of times
        grad = np.diff(t)
        grad_shifted = np.insert(grad, 0, 0)  # Insert a zero at the top (shifting down the array)
        #print(grad)
        #print(grad_shifted)
        start_interval = t[0]

        critical_times = [start_interval]  # List of times to report on (start time, plus times where the steps change)

        total_number_full_steps = 0
        for i in range(1, n_entries-1):
            #print(i)
            if not np.allclose(grad[i], grad_shifted[i]):
                #print (f"   Detection at element {i} : time={t[i]}")
                #print (f"   From time {start_interval} to {t[i]}, in steps of {grad_shifted[i]}")
                primary_timestep = df.loc[i, "primary_timestep"]
                #print("ADDING FULL STEPS:", self._explain_time_advance_helper(t_start=start_interval, t_end=t[i], delta_baseline=grad_shifted[i], primary_timestep=primary_timestep))
                total_number_full_steps += self._explain_time_advance_helper(t_start=start_interval, t_end=t[i], delta_baseline=grad_shifted[i], primary_timestep=primary_timestep)
                start_interval = t[i]
                critical_times.append(t[i])

        #print (f"   From time {start_interval} to {t[-1]}, in steps of {grad_shifted[-1]}")
        primary_timestep = df.loc[n_entries-1, "primary_timestep"]
        #print("ADDING FULL STEPS:", self._explain_time_advance_helper(t_start=start_interval, t_end=t[-1], delta_baseline=grad_shifted[-1], primary_timestep=primary_timestep))
        total_number_full_steps += self._explain_time_advance_helper(t_start=start_interval, t_end=t[-1], delta_baseline=grad_shifted[-1], primary_timestep=primary_timestep)
        critical_times.append(t[-1])
        print(f"(for a grand total of the equivalent of {total_number_full_steps:,.3g} FULL steps)")

        if return_times:
            return critical_times



    def _explain_time_advance_helper(self, t_start, t_end, delta_baseline, primary_timestep)  -> int:
        """
        Using the provided data, about a group of same-size steps, create and print a description of it for the user

        :param t_start:
        :param t_end:
        :param delta_baseline:
        :param primary_timestep:
        :return:                The corresponding number of FULL steps taken (it may be fractional)
        """
        if np.allclose(t_start, t_end):
            #print(f"[Ignoring interval starting and ending at same time {t_start:.3g}]")
            return 0

        if np.allclose(delta_baseline, primary_timestep):
            n_steps = round((t_end - t_start) / delta_baseline)
            name = "step" if n_steps == 1 else "steps"  # singular vs. plural
            print(f"From time {t_start:.3g} to {t_end:.3g}, in {n_steps} FULL {name} of {delta_baseline:.3g}")
            return n_steps
        else:
            #print(f"primary_timestep/delta_baseline is:  {primary_timestep}/{delta_baseline} = {primary_timestep/delta_baseline}")
            n_steps = round((t_end - t_start) / delta_baseline)
            if n_steps == 1:    # Use wording for the singular
                print(f"From time {t_start:.3g} to {t_end:.3g}, in {n_steps} substep of {delta_baseline:.3g} (1/{round(primary_timestep/delta_baseline)} of full step)")
            else:               # Use wording for the plural
                print(f"From time {t_start:.3g} to {t_end:.3g}, in {n_steps} substeps of {delta_baseline:.3g} (each 1/{round(primary_timestep/delta_baseline)} of full step)")

            return (n_steps *  delta_baseline / primary_timestep)




    #############################################################################################
    #                                                                                           #
    #                                       GRAPHICS                                            #
    #                                                                                           #
    #############################################################################################
    def _________GRAPHICS___________(DIVIDER):  # Used to get a better structure listing in IDEs such asPycharm
        pass


    def plot_curves(self, chemicals=None, colors=None, title=None, suppress=False):
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

        :param suppress:    If True, nothing gets shown - and a plotly "Figure" object gets returned instead

        :return:            None or a plotly "Figure" object, depending on the "suppress" flag
        """
        default_colors = ['blue', 'green', 'brown', 'red', 'gray', 'orange', 'purple', 'cyan', 'darkorange', 'navy', 'darkred']

        df = self.history.get()

        if chemicals is None:
            chemicals = self.reaction_data.get_all_names()      # List of names.  EXAMPLE: ["A", "B", "H"]

        number_of_curves = len(chemicals)

        if colors is None:
            colors = default_colors[:number_of_curves]      # Pick the first default colors; TODO: rotate if needing more


        if title is None:   # If no title was specified, create one based on how many reactions are present
            number_of_rxns = self.reaction_data.number_of_reactions()
            if number_of_rxns > 2:
                title = f"Changes in concentrations for {number_of_rxns} reactions"
            elif number_of_rxns == 1:
                rxn_text = self.reaction_data.single_reaction_describe(rxn_index=0, concise=True)   # The only reaction
                title = f"Reaction `{rxn_text}` .  Changes in concentrations with time"
            else:   # Exactly 2 reactions
                rxn_text_0 = self.reaction_data.single_reaction_describe(rxn_index=0, concise=True)
                rxn_text_1 = self.reaction_data.single_reaction_describe(rxn_index=1, concise=True)
                title = f"Changes in concentration for `{rxn_text_0}` and `{rxn_text_1}`"


        fig = px.line(data_frame=df, x="SYSTEM TIME", y=chemicals,
                      title=title,
                      color_discrete_sequence = colors,
                      labels={"value":"concentration", "variable":"Chemical"})

        if not suppress:
            fig.show()
        else:
            return fig




    #############################################################################################
    #                                                                                           #
    #                                   RESULT ANALYSIS                                         #
    #                                                                                           #
    #############################################################################################
    def _________RESULT_ANALYSIS___________(DIVIDER):  # Used to get a better structure listing in IDEs such asPycharm
        pass


    def curve_intersection(self, t_start, t_end, var1, var2) -> (float, float):
        """
        Find and return the intersection of the 2 curves in the columns var1 and var2,
        in the time interval [t_start, t_end]
        If there's more than one intersection, only one - in an unpredictable choice - is returned
        TODO: the current implementation fails in cases where the 2 curves stay within some distance of each other,
              and then one curve jumps on the opposite side of the other curve, at at BIGGER distance.
              See the missed intersection at the end of experiment "reactions_single_compartment/up_regulate_1"

        :param t_start: The start of the time interval being considered
        :param t_end:   The end of the time interval being considered
        :param var1:    The name of the 1st chemical of interest
        :param var2:    The name of the 2nd chemical of interest
        :return:        The pair (time of intersection, common value)
        """
        df = self.get_history(t_start=t_start, t_end=t_end)
        row_index = abs(df[var1] - df[var2]).idxmin()   # The index of the Pandas dataframe row
                                                        #   with the smallest absolute value of difference in concentrations
        print(f"Min abs distance found at row: {row_index}")

        return num.curve_intersect_interpolate(df, row_index,
                                               x="SYSTEM TIME", var1=var1, var2=var2)
