import math
import numpy as np
import pandas as pd
import time
import plotly.express as px
import plotly.graph_objects as pgo
from typing import Union
from life123.chem_data import ChemData
from life123.diagnostics import Diagnostics
from life123.numerical import Numerical
from life123.reactions import Reactions
from life123.reaction_kinetics import ReactionKinetics, VariableTimeSteps
from life123.history import HistoryUniformConcentration, HistoryReactionRate
from life123.visualization.plotly_helper import PlotlyHelper



#############################################################################################

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


class UniformCompartment:
    """
    Used to simulate the dynamics of reactions (in a single compartment.)
    This might be thought of as a "zero-dimensional system"
    """

    def __init__(self, reactions=None, chem_data=None, names=None, preset="mid", enable_diagnostics=False):
        """
        Note: AT MOST 1 of the following 3 arguments can be passed

        :param reactions:   [OPTIONAL 1] Object of type "Reactions", with data about the reactions and the chemicals.
                                If passed, cannot pass either of the args `chem_data` nor `names` (those are both part of the "Reactions" object);
                                if not passed, the reactions can be added later, with calls to add_reaction()

        :param chem_data:   [OPTIONAL 2] Object of type "ChemData" (with data about the chemicals and their reactions).
                                If passed, cannot pass either of the args  `reactions` (an object that contains `chem_data`) nor `names`
                                (names are contained in the `chem_data` object)

        :param names:       [OPTIONAL 3] A single name, or list or tuple of names, of the chemicals;
                                providing a list can be used to make the chemicals appear in a particular desired order.
                                If passed, cannot pass either of the args `reactions` nor `chem_data` (both those object contain the chemical names)


        :param preset:      [OPTIONAL] String with code that can be adjusted make the time resolution finer or coarser;
                                it will stay in effect from now on, unless explicitly changed later

        :param enable_diagnostics:  [OPTIONAL] If True, the diagnostics mode is turned on - and will remain on unless explicitly
                                        disabled by a call to pause_diagnostics();
                                        if False (default), no action taken
        """

        self.chem_data = None       # Object of type "ChemData" (with data about the chemicals and their reactions,
                                    #                            incl. macromolecules)
        self.reactions = None       # Object ot type "Reactions" (with data about all the reactions)


        if chem_data and reactions:
            assert reactions.get_chem_data() == chem_data, \
                "BioSim1D() instantiation: the argument `reactions` is based " \
                "on a 'ChemData' object that doesn't match the one passed by the argument `chem_data`"

        if reactions is not None:
            #assert chem_data is None, \
                #"UniformCompartment instantiation: Cannot pass both `chem_data` and `reactions` as arguments (the `reactions` object already contains the `chem_data`)"
            assert names is None, \
                "UniformCompartment instantiation: Cannot pass both `names` and `reactions` as arguments (the `reactions` object already contains the `names`)"
            self.chem_data = reactions.chem_data
            self.reactions = reactions
        else:           # reactions is None
            if chem_data is None:
                self.chem_data = ChemData(names=names)      # It's ok if names is None
            else:
                assert names is None,\
                    "UniformCompartment instantiation: Cannot pass both `chem_data` and `names` as arguments (the `chem_data` object contains the `names`)"
                self.chem_data = chem_data

            self.reactions = Reactions(chem_data=self.chem_data)


        self.system_time = 0.       # Global time of the system, from initialization on

        # TODO: maybe rename "system" to "system_state", and use "system" to store a list or dict of the chemicals
        #       actually involved in this dynamic simulation
        self.system = None  # Concentration data in the single compartment we're simulating, for all the chemicals
                            # A 1-d Numpy array of the concentrations (floats), in their index order;
                            # the array size is the total number of chemical species.
                            # Each entry is the concentration of the species with that index (in the "ChemData" object)
                            # Note that this is the counterpart - with 1 less dimension - of the array by the same name
                            #       in the class BioSim1D

        self.previous_system = None # Concentration data of all the chemicals at the previous simulation step


        self.macro_system = {}      # A dict mapping macromolecule names to their counts
                                    # EXAMPLE:  {"M1": 1, "M2": 3, "M3": 1}

        self.macro_system_state = {}  # The counterpart of the system data for macro-molecules, if present
                                    # Binding fractions of the applicable transcription factors, for all the macro-molecules,
                                    # over the previous time step, indexed by macromolecule and by binding site number
                                    # EXAMPLE:   {"M1": {1: ("A", 0.2), 2: ("F", 0.93), ("P", 0.)},
                                    #             "M2": {1: ("C", 0.5), 2: ("F", 0.1)},
                                    #             "M3": {} }

                                    # For background, see: https://www.annualreviews.org/doi/10.1146/annurev-cellbio-100617-062719


        self.conc_history = HistoryUniformConcentration(active=True)    # Object used to store user-requested snapshots
                                                                        # of (some of) the chemical concentrations:
                                                                        # 'SYSTEM TIME', 'A', 'B', ..., 'comments'

        self.system_rxn_rates = {}      # Keys are the reaction indexes.  Reaction rates for the last (current) step of all reactions
                                        # EXAMPLE: {0: 0.42, 1: 4.26}

        self.rate_history = HistoryReactionRate(active=True)        # Object used to store user-requested snapshots
                                                                    # of (some of) the chemical reaction rates:
                                                                    # 'SYSTEM TIME', 'rxn0_rate', 'rxn1_rate', ...


        # FOR AUTOMATED ADAPTIVE TIME STEP SIZES
        self.adaptive_steps = VariableTimeSteps()

        if preset:
            self.adaptive_steps.use_adaptive_preset(preset)


        # The following 3 diagnostic values get reset at every run
        self.number_neg_concs = 0
        self.number_soft_aborts = 0



        # ***  FOR DIAGNOSTICS  ***

        self.verbose_list = []          # A list of integers or strings with the codes of the desired verbose checkpoints
                                        #   EXAMPLE: [1, "my_ad_hoc_tag"] to invoke sections of code marked as 1 or "my_ad_hoc_tag"
                                        #   Those sections will have entry points such as:  if "my_ad_hoc_tag" in self.verbose_list


        self.diagnostics_enabled = False  # Flag indicating whether using diagnostics

        self.diagnostics = None         # Object of class Diagnostics

        if enable_diagnostics:
            self.enable_diagnostics()       # Note: self.chem_data must be defined BEFORE this call






    #####################################################################################################

    '''                             ~   TO SET/READ SYSTEM DATA   ~                                   '''

    def ________TO_SET_AND_READ_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def get_chem_data(self) -> ChemData:
        """

        :return:    Object of type "ChemData"
        """
        return self.chem_data



    def set_conc(self, conc: Union[list, tuple, dict], snapshot=True) -> None:
        """
        Set the concentrations of some or all the chemicals

        :param conc:    EITHER
                            (1) a list or tuple of concentration values for ALL the registered chemicals,
                                in their index order
                            OR
                            (2) a dict indexed by the chemical labels, for some or all of the chemicals
                                Anything not specified will be set to zero.
                                EXAMPLE:  {"A": 12.4, "B": 0.23, "E": 2.6}

                        Note: any previous values will get over-written

        :param snapshot:[OPTIONAL] If True (default), add to the history
                            a snapshot of this state being set
        :return:        None
        """
        # TODO: more validations, incl. of individual values being wrong data type
        # TODO: also allow a Numpy array; make sure to do a copy() to it!

        if type(conc) == list or type(conc) == tuple:
            assert len(conc) == self.chem_data.number_of_chemicals(), \
                f"UniformCompartment.set_conc(): when a list or tuple is passed as the argument 'conc', " \
                f"its size must match the number of declared chemicals ({self.chem_data.number_of_chemicals()})"

            assert min(conc) >= 0, \
                f"UniformCompartment.set_conc(): values meant to be chemical concentrations cannot be negative " \
                f"(such as the passed value {min(conc)})"

            self.system = np.array(conc, dtype='d')      # float64      TODO: allow users to specify the type

        elif type(conc) == dict:
            for name, conc_value in conc.items():
                self.set_single_conc(conc=conc_value, species_name=name, snapshot=False)

        if snapshot:
            self.capture_conc_snapshot(caption="Set concentration")  # Save this operation in the history (if enabled)

        if self.diagnostics_enabled:
            # Save up the current System State, with some extra info, as "diagnostic 'concentration' data"
            system_data = self.get_conc_dict(system_data=self.system)   # The current System State, as a dict
            self.diagnostics.save_diagnostic_conc_data(system_data=system_data, system_time=self.system_time,
                                                       caption="Set concentration")



    def set_single_conc(self, conc, species_index=None, species_name=None, snapshot=True) -> None:
        """
        Set the concentrations of 1 chemical
        Note: if both species_index and species_name are provided, species_name is used     TODO: generate an error instead; and pytest this part

        :param conc:            A non-negative number with the desired concentration value
                                    of the chemical specified below.  (Any previous value will get over-written)
        :param species_index:   [OPTIONAL] An integer that indexes the chemical of interest (numbering starts at 0)
        :param species_name:    [OPTIONAL] A name for the chemical of interest.
                                    If both species_index and species_name are provided, species_name is used
                                    At least one of "species_index" and "species_name" must be specified
        :param snapshot:        [OPTIONAL] boolean: if True, add to the history
                                    a snapshot of this state being set.  Default: True
        :return:                None
        """
        # Validate the arguments
        assert conc >= 0, \
            f"UniformCompartment.set_single_conc(): chemical concentrations cannot be negative (value passed: {conc})"

        if species_name is not None:
            species_index = self.chem_data.get_index(species_name)
        elif species_index is not None:
            self.chem_data.assert_valid_chem_index(species_index)
        else:
            raise Exception("UniformCompartment.set_single_conc(): at least one "
                            "of the arguments `species_index` or `species_name` must be provided")


        if self.system is None:
            self.system = np.zeros(self.chem_data.number_of_chemicals(), dtype='d')      # float64      TODO: allow users to specify the type

        # TODO: if setting concentrations of a newly-added chemical, needs to first expand self.system
        self.system[species_index] = conc

        if snapshot:
            # Save this operation in the history (if enabled)
            self.capture_conc_snapshot(caption=f"Set concentration of `{self.chem_data.get_label(species_index)}`")



    def get_system_conc(self) -> np.array:
        """
        Retrieve the concentrations of ALL the chemicals as a Numpy array

        :return:        A Numpy array with the concentrations of ALL the chemicals,
                        in their index order
                            EXAMPLE:  array([12.3, 4.56, 0.12])
                                      The 0-th chemical has concentration 12.3, and so on...
        """
        return self.system



    def get_chem_conc(self, name: str) -> float:
        """
        Return the current system concentration of the given chemical, specified by its name.
        If no chemical by that name exists, an Exception is raised

        :param name:    The name of a chemical species
        :return:        The current system concentration of the above chemical
        """
        species_index = self.chem_data.get_index(name)
        return self.system[species_index]



    def get_conc_dict(self, species=None, system_data=None) -> dict:
        """
        Retrieve the concentrations of the requested chemicals (by default all),
        as a dictionary indexed by the chemicals' labels

        :param species:     [OPTIONAL] List or tuple of labels of the chemical species; by default, return all
        :param system_data: [OPTIONAL] A Numpy array of concentration values, in the same order as the
                                index of the chemical species; by default, use the SYSTEM DATA

        :return:            A dictionary, indexed by the chemical labels, of the concentration values;
                                EXAMPLE: {"A": 1.2, "D": 4.67}
        """
        if system_data is None:
            system_data = self.system
        else:
            assert system_data.size == self.chem_data.number_of_chemicals(), \
                f"UniformCompartment.get_conc_dict(): the argument `system_data` must be a 1-D Numpy array with as many entries " \
                f"as the declared number of chemicals ({self.chem_data.number_of_chemicals()})"


        if species is None:
            if system_data is None:
                return {}
            else:
                return {self.chem_data.get_label(index): system_data[index]
                        for index, conc in enumerate(system_data)}
        else:
            assert type(species) == list or  type(species) == tuple, \
                f"UniformCompartment.get_conc_dict(): the argument `species` must be a list or tuple" \
                f" (it was of type {type(species)})"

            conc_dict = {}
            for name in species:
                species_index = self.chem_data.get_index(name)
                conc_dict[name] = system_data[species_index]

            return conc_dict





    #####################################################################################################

    '''                                 ~  TO VISUALIZE SYSTEM  ~                                     '''

    def ________TO_VISUALIZE_SYSTEM________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def describe_state(self) -> None:
        """
        Print out various data on the current state of the system
        :return:        None
        """
        print(f"SYSTEM STATE at Time t = {self.system_time:,.8g}:")

        n_species = self.chem_data.number_of_chemicals()
        print(f"{n_species} species:")

        # Show a line of line of data for each chemical species in turn
        for species_index, name in enumerate(self.chem_data.get_all_labels()):
            if name:    # If a name was provided, show it
                name = f" ({name})"
            else:
                name = ""

            if self.system is None:
                print(f"  Species {species_index}{name}. No concentrations set yet")
            else:
                print(f"  Species {species_index}{name}. Conc: {self.system[species_index]}")

        if self.macro_system != {}:
            print("Macro-molecules, with their counts: ", self.macro_system)

        if self.macro_system_state != {}:
            print("Fractional Occupancy at the various binding sites for each macro-molecule:")
            for mm, state_dict in self.macro_system_state.items():
                state_list = [f"{a}: {b[1]} ({b[0]})" for a, b in state_dict.items()]   # EXAMPLE: ["3: 0.1 (A)", "8: 0.6 (B)"]
                state_str = " | ".join(state_list)                                      # EXAMPLE: "3: 0.1 (A) | "8: 0.6 (B)"
                print(f"     {mm} || {state_str}")

        if self.reactions.active_enzymes == set():    # If no enzymes were involved in any reaction
            print(f"Set of chemicals involved in reactions: {self.reactions.labels_of_active_chemicals()}")
        else:
            print(f"Set of chemicals involved in reactions (not counting enzymes): {self.reactions.labels_of_active_chemicals()}")
            print(f"Set of enzymes involved in reactions: {self.reactions.names_of_enzymes()}")





    #####################################################################################################

    '''                             ~  MANAGEMENT OF REACTIONS  ~                                     '''

    def ________MANAGE_REACTIONS________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def get_reactions(self):
        """
        Return all the reactions associated to this Uniform Compartment

        :return:    Object ot type "Reactions" (with data about all the reactions)
        """
        return self.reactions



    def get_single_reaction(self, i :int):
        """
        Return a single reaction, specified by its index, associated to this Uniform Compartment

        :param i:   Integer index of the desired reaction
        :return:    Object of type "ReactionGeneric"
        """
        return self.reactions.get_reaction(i)



    def clear_reactions(self) -> None:
        """
        Get rid of all reactions; start again with "an empty slate" (but still with reference
        to the same data object about the chemicals)

        # TODO: maybe offer an option to clear just one reaction, or a list of them
        # TODO: provide support for "inactivating" reactions

        :return:    None
        """
        self.reactions.clear_reactions_data()



    def add_reaction(self, **kwargs) -> int:
        """
        Register a new SINGLE chemical reaction,
        optionally including its kinetic and/or thermodynamic data.

        For details, see ChemData.add_reaction()

        :param kwargs:  Any arbitrary named arguments
        :return:        Integer index of the newly-added reaction
        """
        return self.reactions.add_reaction(**kwargs)



    def describe_reactions(self, **kwargs) -> None:
        """
        Print out a user-friendly plain-text form of ALL the reactions.

        For details, see ChemData.describe_reactions()

        :param kwargs:  Any arbitrary named arguments
        :return:        None
        """
        self.reactions.describe_reactions(**kwargs)



    def number_of_reactions(self) -> int:
        """
        Return the number of registered chemical reactions

        :return:    The number of registered chemical reactions
        """
        return self.reactions.number_of_reactions()



    def plot_reaction_network(self, graphic_component :str, unpack=False) -> None:
        """
        Send a plot of the network of reactions to the HTML log file,
        also including a brief summary of all the reactions

        EXAMPLE of usage:  plot_reaction_network("vue_cytoscape_2")

        :param graphic_component:   The name of a Vue component that accepts a "graph_data" argument,
                                        an object with the following keys
                                        'structure', 'color_mapping' and 'caption_mapping'
                                        For more details, see ChemData.prepare_graph_network()
        :param unpack:              Use True for Vue components that require their data unpacked into individual arguments;
                                        False for that accept a single data argument, named "graph_data"
        :return:                    None
        """
        self.reactions.plot_reaction_network(graphic_component=graphic_component, unpack=unpack)





    #####################################################################################################

    '''                            ~  TO PERFORM THE REACTIONS  ~                                     '''

    def ________TO_PERFORM_THE_REACTIONS________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


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
            "UniformCompartment.specify_steps(): cannot specify all 3 arguments: `total_duration`, `time_step`, `n_steps` (specify only any 2 of them)"

        assert (total_duration and time_step) or (total_duration and n_steps) or (time_step and n_steps), \
            "UniformCompartment.specify_steps(): must provide exactly 2 arguments from:  `total_duration`, `time_step`, `n_steps`"

        if not time_step:
            time_step = total_duration / n_steps

        if not n_steps:
            n_steps = math.ceil(total_duration / time_step)

        return (time_step, n_steps)     # Note: could opt to also return total_duration if there's a need for it



    def single_compartment_react(self, duration=None, target_end_time=None, stop=None,
                                 initial_step=None, n_steps=None, max_steps=None,
                                 variable_steps=True, explain_variable_steps=None,
                                 silent=False, report_interval=0.5) -> None:
        """
        Perform ALL the (previously-registered) reactions in the single compartment -
        based on the INITIAL concentrations stored in self.system

        Update the system state and the system time accordingly
        (object attributes self.system and self.system_time)

        :param duration:        The overall time advance for the reactions (it may be exceeded in case of variable steps)
        :param target_end_time: The final time at which to stop the reaction; it may be exceeded in case of variable steps
                                    If both `target_end_time` and `duration` are specified, an error will result

        :param initial_step:    The suggested size of the first step (it might be reduced automatically,
                                    in case of "hard" errors resulting from overly-large steps)

        :param stop:            Pair of the form (termination_keyword, termination_parameter), to indicate
                                    the criterion to use to stop the reaction
                                    EXAMPLES:
                                        ("conc_below", (chem_name, conc))  Stop when conc first dips below
                                        ("conc_above", (chem_name, conc))  Stop when conc first rises above

                                        TODO: add more options, such as
                                        ("before_time", t)                  Stop just before the given target time
                                        ("after_time", t)                   Stop just after the given target time
                                        ("equilibrium", tolerance)          Stop when equilibrium reached

        :param n_steps:         The desired number of steps

        :param max_steps:       [OPTIONAL] Max numbers of steps; if reached, it'll terminate regardless of any other criteria

        :param variable_steps:          [OPTIONAL] If True (default), the steps sizes will get automatically adjusted, based on thresholds
        :param explain_variable_steps:  [OPTIONAL] If not None, a brief explanation is printed about how the variable step sizes were chosen,
                                            when the System time inside that range;
                                            only applicable if variable_steps is True

        :param silent:                  [OPTIONAL] If True, less output is generated
        :param report_interval:         [OPTIONAL] How frequently, in terms of elapsed running time, in minutes,
                                            to inform the user of the current status

        :return:                        None.   The object attributes self.system and self.system_time get updated
        """

        # Default values
        capture_initial_caption = "1st reaction step"
        capture_final_caption = "last reaction step"

        # Validation
        assert self.system is not None, "UniformCompartment.single_compartment_react(): " \
                                        "the concentration values of the various chemicals must be set first"

        if variable_steps and (n_steps is not None):
            raise Exception("UniformCompartment.single_compartment_react(): if `variable_steps` is True, cannot specify `n_steps` "
                            "(because the number of steps will vary); specify `duration` or `target_end_time` instead")

        if stop is not None:
            assert type(stop) == tuple and len(stop) == 2, \
                f"UniformCompartment.single_compartment_react(): the argument `stop`, if passed, must be a pair of values"

        assert self.reactions.number_of_reactions() > 0, \
            f"UniformCompartment.single_compartment_react(): no reactions are present.  Make sure to first add them with add_reaction()"


        self.conc_history.initial_caption = capture_initial_caption     # TODO: turn into method


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
            time_step, n_steps = self.specify_steps(total_duration=duration,
                                                    time_step=initial_step,
                                                    n_steps=n_steps)
            # Note: if variable steps are requested then n_steps stops being particularly meaningful; it becomes a
            #       hypothetical value, in the (unlikely) event that the step sizes were never changed - and is only
            #       used to detect a very excessive number of actual attempted steps

            if target_end_time is None:
                if variable_steps:
                    target_end_time = self.system_time + duration
                else:
                    target_end_time = self.system_time + time_step * n_steps


        step_count = 0

        # Reset some diagnostic variables
        self.number_neg_concs = 0
        self.number_soft_aborts = 0
        self.adaptive_steps.reset_norm_usage_stats()

        # Time-related
        t_start = time.perf_counter()
        t_report = t_start
        report_interval *= 60.      # Convert to seconds


        '''
        TODO: explore a main loop of the form:
        
        try:
            [MAIN LOOP BODY]
        
        except KeyboardInterrupt:
            print("\n*** KeyboardInterrupt exception caught")
        '''


        # MAIN LOOP
        while True:
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
                    if self.get_chem_conc(chem_name) < conc_threshold:
                        break   # The concentration of the specified chemical has dropped the requested threshold
                elif termination_keyword == "conc_above":
                    chem_name, conc_threshold = termination_parameter
                    if self.get_chem_conc(chem_name) > conc_threshold:
                        break   # The concentration of the specified chemical has risen above the requested threshold

            if (not variable_steps) and (step_count == n_steps)\
                    and (target_end_time is not None) and np.allclose(self.system_time, target_end_time):
                break       # When dealing with fixed steps, catch scenarios where after performing n_steps,
                            #   the System Time is below the target_end_time because of roundoff error


            # ----------  CORE OPERATION OF MAIN LOOP  ----------
            delta_concentrations, step_actually_taken, recommended_next_step = \
                    self.reaction_step_common(delta_time=time_step,
                                              variable_steps=variable_steps, explain_variable_steps=explain_variable_steps,
                                              step_counter=step_count)

            # Update the System State
            self.previous_system = self.system.copy()
            self.system += delta_concentrations
            if min(self.system) < 0:    # Check for negative concentrations. TODO: redundant, since reaction_step_common() now does that
                print(f"***********  SYSTEM STATE ERROR: FAILED TO CATCH negative concentration "
                      f"upon advancing reactions from system time t={self.system_time:,.5g}")


            # Preserve the RATES data, as requested (part1, BEFORE updating the System Time, because reaction rates are
            # based on the *start* time of the simulation step)
            if step_count == 0:
                self.capture_rate_snapshot(force=True, step_count=0)    # Always save the initial rate
            else:
                self.capture_rate_snapshot(step_count=step_count)       # Save historical rate values (if enabled)

            # UPDATE THE SYSTEM TIME (now we're at the END of the current time step)
            self.system_time += step_actually_taken


            # Preserve the CONCENTRATION data, as requested (part2, AFTER updating the System Time, because current concentrations
            # refer to the System Time, just updated at the end of the simulation step)
            self.capture_conc_snapshot(step_count=step_count+1) # Save historical concentration values (if enabled)
                                                                # It's +1 because we save the conc. values at the END of the step


            step_count += 1

            if (n_steps is not None) and (step_count > 1000 * n_steps):  # Another approach to catch infinite loops
                raise Exception("single_compartment_react(): "
                                "the computation is taking a very large number of steps, probably from automatically trying to correct instability;"
                                " trying reducing the time_step")   # TODO: is the explanation correctly phrased?

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
            print(f"{n_steps_taken} total {step_type_str}step(s) taken in {(t_now - t_start)/60.:.2f} min")
            if variable_steps:
                if self.number_neg_concs:
                    print(f"Number of step re-do's because of negative concentrations: {self.number_neg_concs}")
                if self.number_soft_aborts:
                    print(f"Number of step re-do's because of elective soft aborts: {self.number_soft_aborts}")

                print("Norm usage:", self.adaptive_steps.norm_usage)
                print(f"System Time is now: {self.system_time:,.5g}")


        # One final snapshot, unless already taken for the last step done
        self.capture_conc_snapshot(step_count=step_count, caption=capture_final_caption, extra=True)

        # Add a caption to the very last entry in the system history
        self.conc_history.set_caption_last_snapshot(capture_final_caption)



    def reaction_step_common(self, delta_time: float, conc_array=None,
                             variable_steps=False, explain_variable_steps=None, step_counter=1) -> (np.array, float, float):
        """
        This is the common entry point for both single-compartment reactions,
        and the reaction component of reaction-diffusions in 1D, 2D and 3D.

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
        :param conc_array:      [OPTIONAL]All initial concentrations at the start of the reaction step,
                                    as a Numpy array for ALL the chemical species, in their index order.
                                    If not provided, self.system is used instead
        :param variable_steps:  If True, the step sizes will get automatically adjusted with an adaptive algorithm
        :param explain_variable_steps:  If not None, a brief explanation is printed about how the variable step sizes were chosen,
                                            when the System time inside that range;
                                            only applicable if variable_steps is True
        :param step_counter:    [OPTIONAL] Integer currently only used for diagnostics

        :return:                The triplet:
                                    1) increment vector for the concentrations of ALL the chemical species,
                                        in their index order, as a Numpy array
                                        EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):
                                            array([7. , -21.])  TODO: is this really necessary?  Maybe update self.system here?
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

        # Validate arguments
        assert self.system is not None, "UniformCompartment.reaction_step_common(): " \
                                        "the concentration values of the various chemicals must be set first"


        if explain_variable_steps:
            assert (type(explain_variable_steps) == list) and (len(explain_variable_steps) == 2), \
                "reaction_step_common(): the argument `explain_variable_steps`, if provided, must be a pair of numbers [t_start, t_end]"


        #print(f"************ At SYSTEM TIME: {self.system_time:,.4g}, calling reaction_step_common() with:")
        #print(f"             delta_time={delta_time}, system={self.system}, ")


        recommended_next_step = delta_time     # Baseline value; no reason yet to suggest a change in step size


        delta_concentrations = None

        SMALLEST_VALUE_TO_TRY = delta_time / 2000.       # Used to prevent infinite loops (only applicable in case of variable steps)

        normal_exit = False

        while delta_time > SMALLEST_VALUE_TO_TRY:       # TODO: consider moving the inside of the WHILE loop into a separate function
            try:    # We want to catch Exceptions that can arise from excessively large time steps
                    #   that lead to negative concentrations or violation of user-set thresholds ("HARD" or "SOFT" aborts)
                (delta_concentrations, recommended_next_step) = \
                        self._attempt_reaction_step(delta_time, variable_steps, explain_variable_steps, step_counter)
                normal_exit = True

                break       # IMPORTANT: this is needed because, in the absence of errors, we need to go thru the WHILE loop only once!


            # CATCH any 'ExcessiveTimeStepHard' exception raised in the loop  (i.e. a HARD ABORT)
            except ExcessiveTimeStepHard as ex:
                # Single reactions steps can fail with this error condition if the attempted time step was too large,
                # under the following scenarios:
                #       1. negative concentrations from any one reaction - caught by  validate_increment()
                #       2. negative concentration from the combined effect of multiple reactions - caught in this function
                #print("*** CAUGHT a HARD ABORT")
                self.number_neg_concs += 1
                if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                    print(ex)

                delta_time *= self.adaptive_steps.step_factors["error"]       # Reduce the excessive time step by a pre-set factor
                recommended_next_step = delta_time
                # At this point, the loop will generally try the simulation again, with a smaller step (a revised delta_time)


            # CATCH any 'ExcessiveTimeStepSoft' exception raised in the loop  (i.e. a SOFT ABORT)
            except ExcessiveTimeStepSoft as ex:
                # Single reactions steps can fail with this error condition if the attempted time step was too large,
                # under the following scenario:
                #       3. excessive norm(s) measures in the overall step - caught in this function (currently only checked in case of the variable-steps option)
                #print("*** CAUGHT a soft ABORT")
                self.number_soft_aborts += 1
                if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                    print(f"       {ex}")
                delta_time *= self.adaptive_steps.step_factors["abort"]       # Reduce the excessive time step by a pre-set factor
                recommended_next_step = delta_time
                # At this point, the loop will generally try the simulation again, with a smaller step (a revised delta_time)

        # END while


        if not normal_exit:         # i.e., if no reaction simulation took place in the WHILE loop, above
            raise Exception(f"reaction_step_common(): unable to complete the reaction step.  "
                            f"In spite of numerous automated reductions of the time step, "
                            f"it continues to lead to concentration changes that are considered excessive; "
                            f"consider reducing the original time step, and/or increasing the 'abort' thresholds with set_thresholds(). "
                            f"Current values: {self.adaptive_steps.thresholds}")

        # If we get thus far, it's the normal exit of the reaction step


        return  (delta_concentrations, delta_time, recommended_next_step)     # TODO: consider returning tentative_updated_system , since we already computed it




    def _attempt_reaction_step(self, delta_time, variable_steps, explain_variable_steps, step_counter) -> (np.array, float):
        """
        Attempt to perform a single reaction step - and then raise an Exception if it needs to be aborted,
        based on various criteria.
        If variable_steps is True, determine a new value for the "recommended next step"

        :param delta_time:              The requested time duration of the reaction step
        :param variable_steps:          If True, the step sizes will get automatically adjusted with an adaptive algorithm
        :param explain_variable_steps:  If not None, a brief explanation is printed about how the variable step sizes were chosen,
                                            when the System time inside that range;
                                            only applicable if variable_steps is True
        :param step_counter:            A pair with a time range inside which to show in the explanations about the variable step sizes;
                                            only applicable if explain_variable_steps is True

        :return:                The pair (delta_concentrations, recommended_next_step)
        """
        # *****  CORE OPERATION  *****
        delta_concentrations = self._reaction_elemental_step(delta_time=delta_time, rxn_list=None)


        if self.diagnostics_enabled:
            diagnostic_data_snapshot = {}


        recommended_next_step = delta_time       # Baseline value; no reason yet to suggest a change in step size

        if variable_steps:
            decision_data = self.adaptive_steps.adjust_timestep(n_chems=self.chem_data.number_of_chemicals(),
                                                                indexes_of_active_chemicals= self.reactions.indexes_of_active_chemicals(),
                                                                delta_conc=delta_concentrations, baseline_conc=self.system, prev_conc=self.previous_system)
            step_factor = decision_data['step_factor']
            action = decision_data['action']
            all_norms = decision_data['norms']
            applicable_norms = decision_data['applicable_norms']

            if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                if action == "abort":
                    step_status = "aborted"
                else:
                    step_status = "completed"

                print(f"\n(STEP {step_counter} {step_status}) SYSTEM TIME {self.system_time:.5g} : Examining Conc. changes "
                      f"due to tentative Δt={delta_time:.5g} ...")
                print("    Previous: ", self.previous_system)
                print("    Baseline: ", self.system)
                print("    Deltas:   ", delta_concentrations)

                if len(self.reactions.active_chemicals) < self.chem_data.number_of_chemicals():
                    print(f"    Restricting adaptive time step analysis to {len(self.reactions.active_chemicals)} "
                    f"chemicals only: {self.reactions.labels_of_active_chemicals()} , with indexes: {self.reactions.indexes_of_active_chemicals()}")

                print("    Norms:    ", all_norms)
                print("    Thresholds:    ")
                self.adaptive_steps.display_value_against_thresholds(all_norms)

                if action != "stay":    # The step is trivially 1 when the action is "stay"
                    print("    Step Factors:    ", self.adaptive_steps.step_factors)

                print(f"    => Action: '{action.upper()}'  (with step size factor of {step_factor})")


            # Abort the current step if some rate of change is deemed excessive.
            # TODO: maybe ALWAYS check this, regardless of variable-steps option
            if action == "abort":       # NOTE: this is a "strategic" abort, not a hard one from error
                msg =   f"* INFO: the tentative time step ({delta_time:.5g}) " \
                        f"leads to a value of {applicable_norms} > its ABORT threshold:\n" \
                        f"       -> will backtrack, and re-do step with a SMALLER Δt, x{step_factor:.5g} (now set to {delta_time * step_factor:.5g}) " \
                        f"[Step started at t={self.system_time:.5g}, and will rewind there]"
                #print("WARNING: ", msg)
                if self.diagnostics_enabled:
                    # Define the dict diagnostic_data_snapshot
                    diagnostic_data_snapshot['norm_A'] = all_norms.get('norm_A')
                    diagnostic_data_snapshot['norm_B'] = all_norms.get('norm_B')
                    diagnostic_data_snapshot['norm_C'] = all_norms.get('norm_C')
                    diagnostic_data_snapshot['norm_D'] = all_norms.get('norm_D')
                    diagnostic_data_snapshot['action'] = "ABORT"
                    diagnostic_data_snapshot['step_factor'] = step_factor
                    diagnostic_data_snapshot['time_step'] = delta_time
                    self.diagnostics.save_diagnostic_decisions_data(system_time=self.system_time,
                                                        data=diagnostic_data_snapshot, delta_conc_arr=delta_concentrations,
                                                        caption="excessive norm value(s)")

                    # Make a note of the abort action in all the reaction-specific diagnostics
                    self.diagnostics.annotate_abort_rxn_data("aborted: excessive norm value(s)")

                raise ExcessiveTimeStepSoft(msg)    # ABORT THE CURRENT STEP


            # Put together a recommendation to the higher-level functions, about the next best step size
            recommended_next_step = delta_time * step_factor


            if self.diagnostics_enabled:
                # Expand the dict diagnostic_data_snapshot
                diagnostic_data_snapshot['norm_A'] = all_norms.get('norm_A')    # TODO: combine all norms in 1 step
                diagnostic_data_snapshot['norm_B'] = all_norms.get('norm_B')
                diagnostic_data_snapshot['norm_C'] = all_norms.get('norm_C')
                diagnostic_data_snapshot['norm_D'] = all_norms.get('norm_D')
                diagnostic_data_snapshot['action'] = f"OK ({action})"
                diagnostic_data_snapshot['step_factor'] = step_factor
                diagnostic_data_snapshot['time_step'] = delta_time


            if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                msg = "       "

                if step_factor > 1:         # "INCREASE
                    msg +=  f"INFO: COMPLETE STEP NORMALLY and MAKE THE INTERVAL LARGER, " \
                            f"multiplied by {step_factor} (set to {recommended_next_step:.5g}) at the next round, because all norms are low"
                elif step_factor < 1:     # "DECREASE"
                    msg +=  f"INFO: COMPLETE STEP NORMALLY and MAKE THE INTERVAL SMALLER, " \
                            f"multiplied by {step_factor} (set to {recommended_next_step:.5g}) at the next round, because {applicable_norms} is high"
                else:   # "STAY THE COURSE"
                    msg +=  f"INFO: COMPLETE NORMALLY - we're inside the target range of all norms.  No change to step size."

                msg += f"\n    [The current step started at System Time: {self.system_time:.5g}, and will continue to {self.system_time + delta_time:.5g}]"
                print(msg)
        # END if variable_steps


        if self.diagnostics_enabled:
            self.diagnostics.save_diagnostic_decisions_data(system_time=self.system_time,
                                                           data=diagnostic_data_snapshot, delta_conc_arr=delta_concentrations)


        # Check whether the COMBINED delta_concentrations will make any conc negative;
        # if so, raised an "ExcessiveTimeStepHard" exception (a custom exception)
        tentative_updated_system = self.system + delta_concentrations
        if min(tentative_updated_system) < 0:
            if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                print(f"*** CAUTION: negative concentration resulting from the combined effect of all reactions, "
                      f"upon advancing reactions from system time t={self.system_time:,.5g}\n"
                      f"         It'll be AUTOMATICALLY CORRECTED with a reduction in time step size")

            # A type of HARD ABORT is detected (a negative concentration resulting from the combined effect of all reactions)
            if self.diagnostics_enabled:
                self.diagnostics.save_diagnostic_decisions_data(system_time=self.system_time,
                                                                data={"action": "ABORT",
                                                                      "step_factor": self.adaptive_steps.step_factors["error"],
                                                                      "caption": "neg. conc. from combined effect of all rxns",
                                                                      "time_step": delta_time},
                                                                delta_conc_arr=None)
                # Save up diagnostic data for ALL reactions
                self.diagnostics.save_diagnostic_aborted_rxns(system_time=self.system_time, time_step=delta_time,
                                                             caption=f"aborted: neg. conc. from combined multiple rxns")

            raise ExcessiveTimeStepHard(f"INFO: the tentative time step ({delta_time:.5g}) "
                                        f"leads to a NEGATIVE concentration of one of the chemicals: "
                                        f"\n      -> will backtrack, and re-do step with a SMALLER delta time, "
                                        f"multiplied by {self.adaptive_steps.step_factors['error']} (set to {delta_time * self.adaptive_steps.step_factors['error']:.5g}) "
                                        f"[Step started at t={self.system_time:.5g}, and will rewind there.  Baseline values: {self.system} ; delta conc's: {delta_concentrations}]")


        return  (delta_concentrations, recommended_next_step)       # Maybe also return tentative_updated_system



    def _reaction_elemental_step(self, delta_time: float, rxn_list=None) -> np.array:
        """
        Using the system concentration data of ALL the chemical species,
        do the specified SINGLE TIME STEP for ONLY the requested reactions (by default all).

        All computations are based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions (in "forward Euler" approach.)

        Return the Numpy increment vector for ALL the chemical species concentrations, in their index order
        (whether involved in these reactions or not)

        NOTES:  - the actual System Concentrations and the System Time (stored in object variables) are NOT changed
                - if any of the concentrations go negative, an Exception is raised

        :param delta_time:      The time duration of this individual reaction step - assumed to be small enough that the
                                    concentration won't vary significantly during this span.
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

        # Compute the rates ("velocities") of all the applicable reactions
        rate_dict = self.compute_all_reaction_rates(rxn_list=rxn_list)


        if rxn_list is None:    # Meaning ALL (active) reactions
            # A list of the reaction indices of all the active reactions
            rxn_list = self.reactions.active_reaction_indices()


        number_chemicals = self.chem_data.number_of_chemicals()

        # For each applicable reaction, find the needed adjustments ("deltas")
        #   to the concentrations of the reactants and products,
        #   based on the forward and reverse rates of the reaction
        for rxn_index in rxn_list:      # Consider each reaction in turn
            # TODO: maybe switch to a call to the experimental _reaction_elemental_step_SINGLE_REACTION()

            # One array entry per chemical species; notice that this array gets RESET for EACH REACTION
            # TODO: instead of using a potentially very large array (mostly of zeros) for each rxn, consider a dict instead
            #       then combine then at end
            increment_vector_single_rxn = np.zeros(number_chemicals, dtype=float)

            rxn = self.reactions.get_reaction(rxn_index)
            reactants = rxn.extract_reactants()
            products = rxn.extract_products()

            delta_rxn = rate_dict[rxn_index] * delta_time

            """
            Determine the concentration adjustments as a result of this reaction step, 
            for this individual reaction being considered
            """

            # The reactants DECREASE based on the quantity (forward reaction - reverse reaction)
            for r in reactants:
                # Unpack data from the reactant r
                species_name = rxn.extract_species_name(r)
                species_index = self.chem_data.get_index(species_name)
                if species_name == rxn.enzyme:
                    #print(f"*** SKIPPING reactant ENZYME {species_index} in reaction {rxn_index}")
                    continue    # Skip if r is an enzyme for this reaction

                stoichiometry = rxn.extract_stoichiometry(r)

                delta_conc = stoichiometry * (- delta_rxn)  # Increment to this reactant from the reaction being considered
                # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
                # Note: it's not enough to detect conc going negative from combined changes from multiple reactions!
                #       Further testing done upstream
                self.validate_increment(delta_conc=delta_conc, baseline_conc=self.system[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                increment_vector_single_rxn[species_index] += delta_conc


            # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
            for p in products:
                # Unpack data from the reactant r
                species_name = rxn.extract_species_name(p)
                species_index = self.chem_data.get_index(species_name)
                if species_name == rxn.enzyme:
                    #print(f"*** SKIPPING product ENZYME {species_index} in reaction {rxn_index}")
                    continue    # Skip if p is an enzyme for this reaction

                stoichiometry = rxn.extract_stoichiometry(p)

                delta_conc = stoichiometry * delta_rxn  # Increment to this reaction product from the reaction being considered
                # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
                # Note: it's not enough to detect conc going negative from combined changes from multiple reactions!
                #       Further testing done upstream
                self.validate_increment(delta_conc=delta_conc, baseline_conc=self.system[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                increment_vector_single_rxn[species_index] += delta_conc

            # Macro-molecule related part, if applicable    TODO: implement
            if (self.macro_system_state != {}) and (rxn.macro_enzyme is not None):
                print(f"Making adjustments for macro-molecule catalysis for reaction # {rxn_index}")
                print(f"    Macromolecule: {rxn.macro_enzyme[0]}, at site # {rxn.macro_enzyme[1]}")
                print(f"    Site occupancy at the beginning of the time step:")
                print(f"    Macromolecule count:")


            increment_vector += increment_vector_single_rxn     # Accumulate the increment vector across ALL reactions
                                                                # TODO: consider using a small dict in lieu of increment_vector_single_rxn

            if self.diagnostics_enabled:
                self.diagnostics.save_rxn_data(rxn_index=rxn_index,
                                               system_time=self.system_time, time_step=delta_time,
                                               increment_vector_single_rxn=increment_vector_single_rxn,
                                               rate=rate_dict[rxn_index])   # TODO: drop the rate?
        # END for (over rxn_list)


        self.system_rxn_rates = rate_dict

        return increment_vector



    def _reaction_elemental_step_SINGLE_REACTION(self, delta_time: float, increment_vector,
                                                 rxn_index :int,
                                                 delta_dict):     # TODO: EXPERIMENTAL; not yet in use
        """
        :param delta_time:
        :param increment_vector:
        :param rxn_index:       The integer index (0-based) to identify the reaction of interest
        :param delta_dict:
        :return:
        """
        if 2 in self.verbose_list:
            print(f"      Determining the conc.'s changes as a result of rxn # {rxn_index}")

        # NOTE: instead of using a potentially very large array (mostly of zeros) for each rxn,
        #       we use a dict instead; then combine all at end
        #increment_vector_single_rxn = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)
        increment_dict_single_rxn = {}      # The key will be the elements in the reaction

        # TODO: turn into a more efficient single step, as as:
        #(reactants, products) = cls.all_reactions.unpack_terms(rxn_index)
        reactants = self.reactions.get_reactants(rxn_index)
        products = self.reactions.get_products(rxn_index)


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
            self.validate_increment(delta_conc=delta_conc, baseline_conc=self.system[species_index],
                                    rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

            increment_dict_single_rxn[species_index] = increment_dict_single_rxn.get(species_index, 0) + delta_conc


        # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
        for p in products:
            stoichiometry, species_index, order = p             # Unpack
            delta_conc = stoichiometry * delta_dict[rxn_index]  # Increment to this reaction product from the reaction being considered
            # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
            # Note: not enough to detect conc going negative from combined changes from multiple reactions!
            #       Further testing done upstream
            # TODO: it might be more efficient to check this in bulk than with many function calls
            self.validate_increment(delta_conc=delta_conc, baseline_conc=self.system[species_index],
                                    rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

            increment_dict_single_rxn[species_index] = increment_dict_single_rxn.get(species_index, 0) + delta_conc


        for k, v in increment_dict_single_rxn.items():
            increment_vector[k] += v                # Accumulate the increment vector across all the increments from this reaction


        # TODO: need a version of examine_increment_array() and save_diagnostic_data() that accept increment_dict_single_rxn
        #       instead of increment_vector_single_rxn [also, maybe move one or both of them to calling function]
        '''  
        if self.diagnostics:
            self.save_diagnostic_data(delta_time, increment_vector_single_rxn, rxn_index)
        '''



    def validate_increment(self,  delta_conc, baseline_conc: float,
                           rxn_index :int, species_index: int, delta_time) -> None:
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
            #print(f"\n*** CAUTION: negative concentration in chemical `{self.chem_data.get_name(species_index)}` in step starting at t={self.system_time:.4g}.  "
            #      f"It will be AUTOMATICALLY CORRECTED with a reduction in time step size, as follows:")

            # A type of HARD ABORT is detected (a single reaction that, by itself, would lead to a negative concentration;
            #   while it's possible that other coupled reactions might counterbalance this - nonetheless,
            #   it's taken as a sign of excessive step size)
            if self.diagnostics_enabled:
                self.diagnostics.save_diagnostic_decisions_data(system_time=self.system_time,
                                                                data={"action": "ABORT",
                                                                      "step_factor": self.adaptive_steps.step_factors['error'],
                                                                      "caption": f"neg. conc. in {self.chem_data.get_label(species_index)} from rxn # {rxn_index}",
                                                                      "time_step": delta_time},
                                                                delta_conc_arr=None)
                self.diagnostics.save_rxn_data(rxn_index=rxn_index, system_time=self.system_time, time_step=delta_time,
                                               increment_vector_single_rxn=None,
                                               aborted=True,
                                               caption=f"aborted: neg. conc. in `{self.chem_data.get_label(species_index)}`")

            raise ExcessiveTimeStepHard(f"    INFO: the tentative time step ({delta_time:.5g}) "
                                    f"leads to a NEGATIVE concentration of `{self.chem_data.get_label(species_index)}` "
                                    f"from reaction {self.reactions.single_reaction_describe(rxn_index=rxn_index, concise=True)} (rxn # {rxn_index}): "
                                    f"\n      Baseline value: {baseline_conc:.5g} ; delta conc: {delta_conc:.5g}"
                                    f"\n      -> will backtrack, and re-do step with a SMALLER delta time, "
                                    f"multiplied by {self.adaptive_steps.step_factors['error']} (set to {delta_time * self.adaptive_steps.step_factors['error']:.5g}) "
                                    f"\n      [Step started at t={self.system_time:.5g}, and will rewind there]")



    def compute_all_reaction_rates(self, rxn_list=None) -> dict:
        """
        Compute the reaction rates, at the current chemical concentrations,
        for all the specified reaction (by default, all.)
        Return a dict with the rates indexed by the reaction number.

        For background info: https://life123.science/reactions

        :param rxn_list:    OPTIONAL list of reactions (specified by their integer index);
                                if None, do all the reactions.  EXAMPLE: [1, 3, 7]

        :return:            A dict of the reactions rates, as defined by compute_reaction_rate(),
                                at the current concentrations of chemicals in the system (as stored in self.system).
                                The dict is indexed by the reaction number, and contains as many entries as the
                                number of reactions being included in the simulation
        """
        rates_dict = {}

        if rxn_list is None:    # Meaning ALL reactions
            rxn_list = range(self.reactions.number_of_reactions())

        # Process the requested reactions
        name_mapping = self.chem_data.get_label_mapping()
        for i in rxn_list:      # Consider each desired reaction in turn
            rxn = self.reactions.get_reaction(i)
            delta = ReactionKinetics.compute_reaction_rate(rxn=rxn, conc_array=self.system, name_mapping=name_mapping)
            rates_dict[i] = delta

        return rates_dict





    #####################################################################################################

    '''                                 ~  MACROMOLECULE DYNAMICS ~                                   '''

    def ________MACROMOLECULE_DYNAMICS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def set_macromolecules(self, data=None) -> None:
        """
        Specify the macromolecules, and their counts, to be included in the system.
        The fractional occupancy is set to 0 at all binding sites of all the specified macromolecules.
        Any previous data gets over-written.

        Note: to set a single fractional occupancy value, use set_occupancy()

        :param data:    A dict mapping macromolecule names to their counts
                            EXAMPLE:  {"M1": 1, "M2": 3, "M3": 1}
                        If any of the requested macromolecules isn't registered, an Exception will be raised
                        If data=None, then the set of registered macromolecules is used,
                            and all their counts are set to 1
        :return:        None.
                        The object variables self.macro_system and self.macro_system_state get set
        """
        if data is None:
            # Use the registered macromolecules, and set all counts to 1
            data = {}
            for mm in self.chem_data.get_macromolecules():
                data[mm] = 1        # EXAMPLE, after this operation: data = {"M1": 1}

        self.macro_system = data


        self.macro_system_state = {}    # Reset
        for mm in data.keys():          # For each macromolecule in our system
            binding_sites_and_ligands = self.chem_data.get_binding_sites_and_ligands(mm)     # EXAMPLE: {1: "A", 2: "C"}
            d = {}
            for (site_number, ligand) in binding_sites_and_ligands.items():
                d[site_number] = (ligand, 0.)           # All "binding occupancy fractions" are set to 0.

            self.macro_system_state[mm] = d



    def set_occupancy(self, macromolecule, site_number: int, fractional_occupancy: float) -> None:
        """
        Set the fractional occupancy at the given binding site of the specified macromolecule,
        using the requested value.
        If the specified macromolecule hasn't yet been added to the dynamical system state,
        automatically add it with count 1

        :param macromolecule:           Name of a previously-registered macromolecule
        :param site_number:             Integer to identify a binding site on the macromolecule
        :param fractional_occupancy:    A number between 0. and 1., inclusive
        :return:                        None
        """
        assert 0. <= fractional_occupancy <= 1., \
            f"set_occupancy(): the value for the fractional occupancy must be a number between 0. and 1. " \
            f"Value given: {fractional_occupancy}"

        ligand = self.chem_data.get_ligand_name(macromolecule=macromolecule, site_number=site_number)

        # If the specified macromolecule hasn't yet been added to the dynamic system, automatically add it
        # with count 1
        if self.macro_system_state == {}:
            self.set_macromolecules({macromolecule: 1})

        self.macro_system_state[macromolecule][site_number] = (ligand, fractional_occupancy)



    def get_occupancy(self, macromolecule, site_number) -> float:
        """
        Get the fractional occupancy at the given binding site of the specified macromolecule.

        :param macromolecule:           Name of a previously-registered macromolecule
        :param site_number:             Integer to identify a binding site on the macromolecule
        :return:                        A number between 0. and 1., representing the fractional occupancy
        """
        assert self.macro_system_state != {}, \
            "get_occupancy(): The system state for macromolecules has not been set yet;  " \
            "use set_macromolecules() or set_occupancy()"

        assert macromolecule in self.macro_system_state, \
            f"get_occupancy(): No occupancy data yet set for macromolecule `{macromolecule}`"

        assert site_number in self.macro_system_state[macromolecule], \
            f"get_occupancy(): No occupancy data yet set for site number {site_number} " \
            f"of macromolecule `{macromolecule}`"

        (ligand, fractional_occupancy) = self.macro_system_state[macromolecule][site_number]
        return fractional_occupancy



    def update_occupancy(self) -> None:
        """
        Update the fractional occupancy at all binding sites,
        based on the current system concentrations of the relevant ligands

        :return:    None
        """
        for mm in self.chem_data.get_macromolecules():
            # For each macromolecule
            d = self.chem_data.get_binding_sites_and_ligands(mm)    # EXAMPLE: {1: "A", 2: "C"}
            for (site_number, ligand) in d.items():
                aff_data = self.chem_data.get_binding_site_affinity(mm, site_number)
                conc = self.get_chem_conc(ligand)
                fractional_occupancy = self.sigmoid(conc=conc, Kd=aff_data.Kd)

                self.set_occupancy(macromolecule=mm, site_number=site_number, fractional_occupancy=fractional_occupancy)



    def sigmoid(self, conc: float, Kd: float) -> float:
        """
        Return an estimate of fractional occupancy (between 0 and 1)
        on a particular binding site on a particular macromolecule,
        from the concentration of the ligand (such as a Transcription Factor)
        and its affinity to that binding site.

        A sigmoid curve is expected.

        Based on fig. 3A of the 2019 paper "Low-Affinity Binding Sites and the
        Transcription Factor Specificity Paradox in Eukaryotes"
        (https://doi.org/10.1146/annurev-cellbio-100617-062719):

            - at extremely low concentration, the occupancy is 0
            - when the concentration is 10% of Kd, the occupancy is about 0.1
            - when the concentration matches Kd, the occupancy is 1/2 by definition
            - when the concentration is 10 times Kd, the occupancy is about 0.9
            - at concentrations beyond that, the occupancy saturates to 1.0

        :param conc:    Concentration of the ligand (such as a Transcription Factor), in microMolars
        :param Kd:      Binding-side Affinity, in microMolars
        :return:        Estimated binding-site fractional occupancy : a value between
                            0. (no occupancy at all during the previous time step) and 1. (continuous occupancy)
        """
        if conc == 0:
            conc = 1e-15     # To avoid taking log of 0

        return self.logistic(x = math.log10(conc), x0 = math.log10(Kd), k = 2.1972245)



    def logistic(self, x: float, x0 = 0., k = 1.) -> float:
        """
        Compute the value of the Logistic function, in the range (0, 1), at the given point
        See: https://en.wikipedia.org/wiki/Logistic_function

        :param x:
        :param x0:
        :param k:
        :return:    The value of the Logistic function at the given point x
        """
        return 1. / ( 1 + math.exp( -k * (x-x0) ) )




    #####################################################################################################

    '''                                  ~   FOR DIAGNOSTICS   ~                                      '''

    def ________FOR_DIAGNOSTICS________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def enable_diagnostics(self):
        """
        Turn on the diagnostics mode

        :return: None
        """
        self.diagnostics_enabled = True
        if not self.diagnostics:
            self.diagnostics = Diagnostics(reactions=self.reactions)


    def pause_diagnostics(self):
        """
        Turn off the overall diagnostics mode; existing diagnostics data, if any, is left untouched

        :return:    None
        """
        self.diagnostics_enabled = False



    def get_diagnostics(self):
        """

        :return:    Object of type life123.diagnostics.Diagnostics
        """
        assert self.diagnostics is not None, \
            "get_diagnostics(): no diagnostics data is available.  Did you call enable_diagnostics() prior to running the simulation?"

        return self.diagnostics





    #####################################################################################################

    '''                                    ~   VISUALIZATION   ~                                           '''

    def ________VISUALIZATION________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def plot_history(self, chemicals=None, colors=None, title=None, title_prefix=None,
                     range_x=None, range_y=None,
                     y_label=None,
                     vertical_lines_to_add=None, show_intervals=False, show=False) -> pgo.Figure:
        """
        Using plotly, draw the plots of chemical concentration values over time, based on historical data that gets
        automatically saved when running reactions.

        Note: if this plot is to be later combined with others, use PlotlyHelper.combine_plots()
              EXAMPLE:
                    from life123 import PlotlyHelper
                    p1 = plot_history(various args, show=False)
                    p2 = plot_history(various args, show=False)
                    PlotlyHelper.combine_plots([p1, p2], other optional args)

        :param chemicals:       [OPTIONAL] Label, or list of labels, of the chemicals whose concentration changes are to be plotted;
                                    if None, then display all, in their index order
        :param colors:          [OPTIONAL] Either a single color (string with standard plotly name, such as "red"),
                                    or list of names to use, in the same order as the chemicals of the previous argument;
                                    if None, then use the registered colors (if specified),
                                    or the hardwired defaults as a last resort
        :param title:           [OPTIONAL] Title for the plot;
                                    if None, use default titles that will vary based on the # of reactions; EXAMPLES:
                                    "Changes in concentrations for 5 reactions"
                                    "Reaction `A <-> 2 B` .  Changes in concentrations with time"
                                    "Changes in concentration for `2 S <-> U` and `S <-> X`"
        :param title_prefix:    [OPTIONAL] If present, it gets prefixed (followed by ".  ") to the title,
                                    whether the title is specified by the user or automatically generated
        :param range_x:         [OPTIONAL] list of the form [t_start, t_end], to initially show only a part of the timeline.
                                    Note: it's still possible to zoom out, and see the excluded portion
        :param range_y:         [OPTIONAL] list of the form [y_min, y_max], to initially show only a part of the y values.
                                    Note: it's still possible to zoom out, and see the excluded portion
        :param y_label:          [OPTIONAL] Caption to use for the y-axis.
                                    By default, the name in `the chemicals` argument, in square brackets, if only 1 chemical,
                                    or "Concentration" if more than 1 (a legend also shown)
        :param vertical_lines_to_add: [OPTIONAL] Ignored if the argument `show_intervals` is specified.
                                    List or tuple or Numpy array or Pandas series
                                    of x-coordinates at which to draw thin vertical dotted gray lines.
                                    If the number of vertical lines is so large as to overwhelm the plot,
                                    only a sample of them is shown.
                                    Note that vertical lines, if requested, go into the plot's "layout";
                                    as a result they might not appear if this plot is later combined with another one.
        :param show_intervals:  [OPTIONAL] If True, it over-rides any value passed to the `vertical_lines` argument,
                                    and draws thin vertical dotted gray lines at all the x-coords
                                    of the data points in the saved history data;
                                    also, it adds a comment to the title.
        :param show:            If True, the plot will be shown
                                    Note: on JupyterLab, simply returning a plot object (without assigning it to a variable)
                                          leads to it being automatically shown

        :return:                A plotly "Figure" object
        """
        if chemicals is None:
            chemicals = self.chem_data.get_all_labels()      # List of the chemical labels.  EXAMPLE: ["A", "B", "H"]

        if title is None:   # If no title was specified, create a default one based on how many reactions are present
            number_of_rxns = self.reactions.number_of_reactions()
            if number_of_rxns > 2:
                title = f"Changes in concentrations for {number_of_rxns} reactions"
            elif number_of_rxns == 1:
                rxn_text = self.reactions.single_reaction_describe(rxn_index=0, concise=True)   # The only reaction
                title = f"Reaction `{rxn_text}` .  Changes in concentrations with time"
            else:   # Exactly 2 reactions
                rxn_text_0 = self.reactions.single_reaction_describe(rxn_index=0, concise=True)
                rxn_text_1 = self.reactions.single_reaction_describe(rxn_index=1, concise=True)
                title = f"Changes in concentration for `{rxn_text_0}` and `{rxn_text_1}`"

        if y_label is None:
            if type(chemicals) == str:
                y_label = f"[{chemicals}]"          # EXAMPLE:  "[A]"
            else:
                y_label = "Concentration"

        df = self.get_history()         # A Pandas dataframe that contains a column named "SYSTEM TIME"

        if colors is None:
            # TODO: switch to  colors = self.chem_data.get_registered_colors(chemicals)
            # Attempt to use the colors registered for individual chemicals, if present
            registered_colors = []
            for label in chemicals:
                stored_color = self.chem_data.get_plot_color(label)     # Will be None if no color was registered for this chemical
                registered_colors.append(stored_color)
            colors = registered_colors      # List of colors, with as many entries as the chemicals of interest;
                                            # any of the entries might be None

        return PlotlyHelper.plot_pandas(df=df, x_var="SYSTEM TIME", fields=chemicals,
                                        colors=colors, title=title, title_prefix=title_prefix,
                                        range_x=range_x, range_y=range_y,
                                        y_label=y_label, legend_header="Chemical",
                                        vertical_lines_to_add=vertical_lines_to_add, show_intervals=show_intervals, show=show)



    def plot_step_sizes(self, show_intervals=False) -> None:
        """
        Using plotly, draw the plot of the step sizes vs. time
        (only meaningful when the variable-step option was used).
        The same scale as plot_history() will be used.
        This function requires the diagnostics option to be turned on, prior to running the simulation

        :param show_intervals:  If True, will add to the plot thin vertical dotted gray lines
                                    at the time steps
        :return:                None
        """
        (transition_times, step_sizes) = self.diagnostics.explain_time_advance(return_times=True, silent=True)

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


        df = self.diagnostics.get_diagnostic_conc_data()    # Pandas dataframe with a column called "TIME"

        # Note: the step size at the final end time isn't a defined quantity - so, we'll just repeat
        #       the last value, to maintain the full x-axis size
        fig = px.line(x=new_x, y=new_y)

        if show_intervals:
            for xi in df["TIME"]:
                fig.add_vline(x=xi, line_width=1, line_dash="dot", line_color="gray")


        fig.update_layout(title='Simulation step sizes',
                          xaxis_title='SYSTEM TIME',
                          yaxis_title='Step size')

        fig.show()





    #####################################################################################################

    '''                                      ~   HISTORY   ~                                          '''

    def ________HISTORY________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def enable_history(self, kind="all", frequency=1, chem_labels=None, take_snapshot=False) -> None:
        """
        Request history capture, with the specified parameters.
        If history was already enabled, this function can be used to alter its capture parameters.

        :param kind:            [OPTIONAL] Allowed values are "all" (default),
                                    "conc" (concentration history), and "rate" (reaction rate history)
        :param frequency:       [OPTIONAL] How often to capture a snapshot of the data
                                    (default: at every simulation cycle)
        :param chem_labels:     [OPTIONAL] List of chemicals to include in the history;
                                    if None (default), include them all.
        :param take_snapshot:   If True, a snapshot of the system's current configuration is added to the history

        :return:                None
        """
        assert kind in ["all", "conc", "rate"], \
            "enable_history(): argument `kind` must be either 'all' or 'conc' or 'rate'"

        if kind in ["all", "conc"]:
            self.conc_history.enable_history(frequency=frequency, chem_labels=chem_labels)
            if take_snapshot:
                self.capture_conc_snapshot()

        if kind in ["all", "rate"]:
            self.rate_history.enable_history(frequency=frequency, chem_labels=chem_labels)
            if take_snapshot:
                self.capture_rate_snapshot()



    def capture_rate_snapshot(self, step_count=None, caption=None,
                             force=False, system_time=None) -> None:
        """

        :param step_count:
        :param caption:
        :param force:           [OPTIONAL] If True, take a snapshot regardless of step_count; default is False
        :param system_time:
        :return:                None
        """
        if not force and not self.rate_history.to_capture(step_count):
            return

        if system_time is None:
            system_time = self.system_time

        data_snapshot = {}     # rxn_rates_snapshot
        for k, v in self.system_rxn_rates.items():
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

        :param step_count:
        :param caption:
        :param extra:       [OPTIONAL] If True, it means that this is a special extra capture;
                                the capture frequency will NOT considered in the
                                decision about saving it, but an extra check will be performed
                                to make sure it's not a duplicate of the earlier capture

        :return:                    None
        """
        if not self.conc_history.to_capture(step_count, extra=extra):
            return

        data_snapshot = self.get_conc_dict(species=self.conc_history.restrict_chemicals)
        '''
           EXAMPLE of data_snapshot:
                {"A": 1.3, "B": 4.9}        
        '''
        self.conc_history.save_snapshot(step_count=step_count, system_time=self.system_time,
                                        data_snapshot=data_snapshot,
                                        caption=caption)



    def start_csv_log(self, log_file :str) -> None:
        """
        Register the specified filename for all future CSV logs.
        Any existing file by that name will get over-written.

        :param log_file:Name of a file to use for the CSV logs
        :return:        None
        """
        self.conc_history.start_csv_log(log_file)



    def get_history(self, t_start=None, t_end=None, head=None, tail=None, t=None, columns=None) -> pd.DataFrame:
        """
        Retrieve and return a Pandas dataframe with the system history that had been saved
        using capture_conc_snapshot().
        Optionally, restrict the result with a start and/or end times,
        or by limiting to a specified numbers of rows at the end

        :param t_start: [OPTIONAL] Start time in the "SYSTEM TIME" column.  Watch out for roundoff errors!
        :param t_end:   [OPTIONAL] End time.  Watch out for roundoff errors!
        :param head:    [OPTIONAL] Number of records to return,
                                   from the start of the history dataframe.
        :param tail:    [OPTIONAL] Number of records to consider, from the end of the history dataframe
        :param t:       [OPTIONAL] Individual time to pluck out from the dataframe;
                                   the row with closest time will be returned.
                                   If this parameter is specified, an extra column - called "search_value" -
                                   is inserted at the beginning of the dataframe.
                                   If either the "head" or the "tail" arguments are passed, this argument will get ignored
        :param columns: [OPTIONAL] Name, or list of names, of the column(s) to return; if not specified, all are returned.
                                   Make sure to include "SYSTEM TIME" in the list, if the time variable needs to be included

        :return:        A Pandas dataframe
        """
        #TODO: allow searches also for columns other than "SYSTEM TIME"

        # Note: the history is an object of class CollectionTabular

        df = self.conc_history.get_history().get_dataframe(head=head, tail=tail, search_val=t,
                                                           search_col="SYSTEM TIME",
                                                           val_start=t_start, val_end=t_end)

        if columns:
            assert (type(columns) == list) or (type(columns) == str), \
                "get_history(): the argument `columns`, if specified, must be a list or string"
            return df[columns]

        return df



    def get_historical_concentrations(self, row=None, t=None, df=None) -> np.array:
        """
        Typically used to retrieve a system snapshot from its history.

        Return a Numpy array with ALL the chemical concentrations (in their index order)
        from one row in the given Pandas data frame (by default, the system history);
        the row can be identified either by it row number, or by the system time.

        :param row: [OPTIONAL] Integer with the zero-based row number of the system history
                        (which is a Pandas data frame)
        :param t:   [OPTIONAL] Individual time to pluck out from the dataframe;
                        the row with closest time will be returned.
                        Exactly one of "t" and "row" must be specified
        :param df:  [OPTIONAL] A Pandas data frame with concentration information in columns that have
                        the names of the chemicals (if None, the system history is used)
        :return:    A Numpy array of floats.  EXAMPLE: array([200., 40.5])
        """
        assert row is None or t is None, "Cannot specify both arguments `row` and `t`"
        assert row is not None or t is not None, "Must specify either argument `row` and `t`"

        if df is None:
            if t:
                df = self.get_history(t=t)
            else:
                df = self.get_history()

        chem_list = self.chem_data.get_all_labels()  # List of all the chemicals' names

        if row:
            return df.loc[row][chem_list].to_numpy(dtype='float32')
        else:
            return df.iloc[0][chem_list].to_numpy(dtype='float32')



    def get_rate_history(self):
        """
        Return a Pandas dataframe with all the saved reaction rate history, for all reactions.

        Rates refer to reaction products with stoichiometric coefficient 1.
        To find the rate of change of the concentration of any product,
        multiply the reaction rate for its stoichiometric coefficient.
        For reactants, flip the signs.

        :return:    A Pandas dataframe with the following columns:
                        'SYSTEM TIME', 'rxn0_rate', 'rxn1_rate', ...
        """
        # TODO: add the option to select a subset of the reactions
        return self.rate_history.get_history().get_dataframe()




    def add_rate_to_conc_history(self, rate_name=None, new_rate_name=None):
        """
        Merge together the concentration history and a column from the reaction rate history,
        into a new dataframe.
        The original histories aren't affected.

        :param rate_name:       [OPTIONAL] Name of a single desired column from the reaction rate history
                                    EXAMPLE: "rxn4_rate"
                                    If not specified, all rates are included
        :param new_rate_name:   [OPTIONAL] New name to assign to the rate_name column above;
                                    only applicable if rate_name is specified
        :return:                A Pandas dataframe with all the concentration history,
                                    and extra columns from the reaction rate history,
                                    merged by their match of SYSTEM TIME.
                                    If either history is missing, raise an Exception
        """
        # TODO: perhaps make obsolete, by storing both histories together??

        history = self.get_history()
        rates = self.get_rate_history()

        assert (not history.empty) and (not rates.empty), \
            "add_rate_to_conc_history(): One or both the histories are missing"

        if rate_name is not None:
            rates = rates[["SYSTEM TIME", rate_name]]
            if new_rate_name:
                rates = rates.rename(columns={rate_name: new_rate_name})

        # Set the tolerance for the merges based on the size of the time values
        last_time = history["SYSTEM TIME"].iat[-1]
        tolerance = last_time * 0.000001 + 0.000001     # We use a relative and an absolute tolerance

        # Create a "left join" from the two dataframes, based on the colum 'SYSTEM TIME'
        merged_df = pd.merge_asof(history, rates, on='SYSTEM TIME', direction='nearest', tolerance=tolerance)

        # Filter to remove rows where the right DataFrame values are NaN (indicating no match)
        # This has the overall effect of producing an inner join of the two original dataframes
        # Note: the names "step_x" and "step_y" are automatically assigned by merge_asof() because of a repeated "step" column
        inner_join_df = merged_df.dropna().rename(columns={"step_x": "step_conc", "step_y": "step_rate"})

        return inner_join_df



    def add_rate_to_conc_history_OLD(self, rate_name :str, new_rate_name=None):
        """
        Merge together the concentration history and a column from the reaction rate history

        :param rate_name:       Name of the desired column from the reaction rate history
                                    EXAMPLE: "rxn1_rate"
        :param new_rate_name:   [OPTIONAL] New name for the above column
        :return:                A Pandas dataframe with all the concentration history,
                                    and an extra column from the reaction rate history
        """
        # TODO: possibly make obsolete, by storing both histories together

        history = self.get_history()
        rates = self.get_rate_history()
        assert len(history) == len(rates), \
            f"add_rate_to_conc_history(): unable to reconcile the system history data ({len(history)} rows)" \
            f"with the reaction data ({len(rates)} rows) - mismatched number of rows"

        df = history.copy()    # Duplicate the dataframe, to avoid messing up the concentration history

        if new_rate_name is None:
            new_rate_name = rate_name   # Rename column, if requested

        df[new_rate_name] = rates[rate_name]

        return df






    #####################################################################################################

    '''                                ~   RESULT ANALYSIS   ~                                        '''

    def ________RESULT_ANALYSIS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def is_in_equilibrium(self, rxn_index=None, conc=None, tolerance=1, explain=True, verbose=None) -> Union[bool, dict]:
        """
        Ascertain whether the given concentrations (by default the current System concentrations)
        are in equilibrium for the specified reactions (by default, check all reactions)

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest;
                                if None, then check all the reactions
        :param conc:        Dict with the concentrations of the species involved in the reaction(s).
                            The keys are the chemical names
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
                            If None, then use the current System concentrations instead
        :param tolerance:   Allowable relative tolerance, as a PERCENTAGE,
                                to establish satisfactory match with expected values
        :param explain:     If True, print out details about the analysis,
                                incl. the formula(s) being used to check the equilibrium
                                EXAMPLES:   "([C][D]) / ([A][B])"
                                            "[B] / [A]^2"
        :param verbose:     Alternate name for argument `explain`   TODO: DEPRECATED in favor of "explain"

        :return:            Return True if ALL the reactions are close enough to an equilibrium,
                                as allowed by the requested tolerance;
                                otherwise, return a dict of the form {False: [list of reaction indexes]}
                                for all the reactions that failed the criterion
                                EXAMPLE:  {False: [3, 6]}
        """
        # TODO: optionally display last lines in diagnostic data, if available
        # TODO: phase out "verbose" in favour of "explain"

        if verbose is not None:
            print("*** INFO: argument `verbose` is deprecated; use `explain` instead")
            explain = verbose

        if conc is None:
            conc=self.get_conc_dict()   # Use the current System concentrations, as a dict.
                                        # EXAMPLE: {'A': 23.9, 'B': 36.1}

        failures_dict = {False: []}     # 1-element dict whose value is
                                        # a list of reactions that fail to meet the criterion for equilibrium

        if rxn_index is not None:
            # Check the 1 reaction that was requested
            if explain:
                description = self.reactions.single_reaction_describe(rxn_index=rxn_index, concise=True)
                print(description)

            status = self.reaction_in_equilibrium(rxn_index=rxn_index, conc=conc, tolerance=tolerance, explain=explain)
            if not status:
                failures_dict = {False: [rxn_index]}    # Make a note that this reaction failed the equilibrium test

        else:
            # Check ALL the reactions
            status = True       # Overall status
            description_list = self.reactions.multiple_reactions_describe(concise=True)
            for rxn_index in range(self.reactions.number_of_reactions()):
                # For each reaction
                if explain:
                    print(description_list[rxn_index])

                single_status = self.reaction_in_equilibrium(rxn_index=rxn_index, conc=conc, tolerance=tolerance, explain=explain)

                if not single_status:
                    status = False
                    failures_dict[False].append(rxn_index)      # Make a note that this reaction failed the equilibrium test

        if status:
            return True
        else:
            return failures_dict



    def reaction_in_equilibrium(self, rxn_index :int, conc, tolerance, explain :bool) -> bool:
        """
        Ascertain whether the given concentrations are in equilibrium for the specified SINGLE reaction;
        return True or False, accordingly.

        Pathological case: if at least one of the reactants AND at least one of the products have zero
        concentration, then the reaction is "stuck" - and thus regarded in "equilibrium"

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest
        :param conc:        Dictionary with the concentrations of the species involved in the reaction.
                            The keys are the chemical names
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
        :param tolerance:   Allowable relative tolerance, as a PERCENTAGE,
                                to establish satisfactory match with expected values
        :param explain:     If True, print out the formula being used, as well as the concentrations (reactants, then products)
                                and some other info.
                                EXAMPLES of formulas:   "([C][D]) / ([A][B])"
                                                        "[B] / [A]^2"
        :return:            True if the given reaction is close enough to an equilibrium,
                                as allowed by the requested tolerance
        """
        rxn = self.reactions.get_reaction(rxn_index)    # Look up the requested reaction

        if np.allclose(rxn.extract_reverse_rate(), 0):
            print("reaction_in_equilibrium() currently does NOT handle irreversible reactions (with a zero reverse rate)")
            return False

        rate_ratio = rxn.extract_forward_rate() / rxn.extract_reverse_rate()  # Ratio of forward/reverse reaction rates

        result = rxn.reaction_quotient(conc=conc, explain=explain)

        # Unpack the result
        if explain:
            rxn_quotient, formula = result
        else:
            rxn_quotient = result

        if explain:
            # Prepare a concise listing from the given concentrations,
            # only including the concentrations that are applicable to this reaction
            all_applicable_concs = []
            '''
            for species_name in rxn.extract_chemicals_in_reaction(exclude_enzyme=False):
                s = f"[{species_name}] = {conc[species_name]:,.4g}"         # EXAMPLE: "[A] = 20.3"
                all_applicable_concs.append(s)
            '''
            reactants = rxn.extract_reactant_names(exclude_enzyme=False)
            for species_name in reactants:
                s = f"[{species_name}] = {conc[species_name]:,.4g}"         # EXAMPLE: "[A] = 20.3"
                all_applicable_concs.append(s)

            products = rxn.extract_product_names(exclude_enzyme=False)
            for species_name in products:
                if species_name not in reactants:           # Don't report the same concentration twice!
                    s = f"[{species_name}] = {conc[species_name]:,.4g}"         # EXAMPLE: "[B] = 0.3"
                    all_applicable_concs.append(s)

            all_applicable_concs_str = " ; ".join(all_applicable_concs)     # EXAMPLE: "[A] = 20.3 ; [B] = 0.3"


        # Handle the special case of zero concentration in some of the reactants AND in some of the products
        if np.isnan(rxn_quotient):
            if explain:
                print(f"Current concentrations: {all_applicable_concs_str}")
                print("Reaction IS in equilibrium because it can't proceed in either direction "
                      "due to zero concentrations in both some reactants and products!\n")
            return True

        if np.isinf(rxn_quotient):      # If only the denominator is 0
            if explain:
                print(f"Current concentrations: {all_applicable_concs_str}")
                print("Reaction is NOT in equilibrium because either the products contain a zero concentration\n")
            return False


        # Normal case
        status = np.allclose(rxn_quotient, rate_ratio, rtol=tolerance/100., atol=0)

        if explain:
            print(f"Current concentrations: {all_applicable_concs_str}")
            print(f"1. Ratio of reactant/product concentrations, adjusted for reaction orders: {rxn_quotient:,.6g}")
            print(f"    Formula used:  {formula}")
            print(f"2. Ratio of forward/reverse reaction rates: {rate_ratio:,.6g}")
            print(f"Discrepancy between the two values: {100 * abs(rxn_quotient - rate_ratio)/rate_ratio :,.4g} %")
            if status:
                print(f"Reaction IS in equilibrium (within {tolerance:.2g}% tolerance)\n")
            else:
                print(f"Reaction is NOT in equilibrium (not within {tolerance:.2g}% tolerance)\n")

        return status



    def find_equilibrium_conc(self, rxn_index :int) -> dict:
        """
        Determine the equilibrium concentrations that would be reached by the chemicals
        participating in the specified reversible reaction, given their current concentrations,
        IN THE ABSENCE of any other reaction.

        RESTRICTIONS:  currently limited to just aA + bB <-> cC + dD reversible reactions,
                       first-order in all chemicals (some of the terms may be missing.)
                       An Exception will be raised in all other cases!

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest
        :return:            A dictionary of the equilibrium concentrations of the
                                chemicals involved in the specified reaction
                            EXAMPLE:  {'A': 24.0, 'B': 36.0, 'C': 1.8}
        """
        #TODO: handle scenarios where kF or kR is zero

        rxn = self.reactions.get_reaction(rxn_index)    # Look up the requested reaction

        reactants = rxn.extract_reactants()
        products = rxn.extract_products()
        K = rxn.extract_equilibrium_constant()


        assert len(reactants) <= 2, \
                "find_equilibrium_conc(): Currently only implemented " \
                "for reactions with at most 2 reactants and 2 products"
        assert len(products) <= 2, \
                "find_equilibrium_conc(): Currently only implemented " \
                "for reactions with at most 2 reactants and 2 products"

        '''
        For reactions of the form aA + bB <-> cC + dD  that are first-order in all chemicals,
        the equilibrium equation is:  
                [(C0 + c*m) (D0 + d*m)] / [(A0 - a*m) (B0 - b*m)]  =  K
            
        where K is the Equilibrium constant (kF/kR),
        and the unknown m (to be solved for) is the number of "moles/liter of forward reaction"
        (i.e. the the Product concentration change)
        from the starting point to the equilibrium point
        
        For reaction terms that aren't present (for example the "D" part in the reaction A + B <-> C),
        we'll use 0 for the "stoichiometry coefficient" and 1 for the "initial concentration";
        that's a hack to make the multiplicative terms of the form X0 + x*m  (where the X's are the A,B,C,D and the x's are the a,b,c,d)
        to be identical equal to 1, and thus have no effect on the solution of the above equation
        '''
        a = 0
        b = 0
        c = 0
        d = 0
        A0 = 1
        B0 = 1
        C0 = 1
        D0 = 1


        name_map = {}   # To transform the A, B, C, D into the actual names

        # Look at the denominator of the "Reaction Quotient"
        for i, r in enumerate(reactants):
            # Loop over the reactants
            species_name =  rxn.extract_species_name(r)
            rxn_order =  rxn.extract_rxn_order(r)
            coefficient = rxn.extract_stoichiometry(r)

            assert rxn_order == 1, \
                "find_equilibrium_conc(): Currently only implemented for 1st order reactions"

            if i == 0:
                name_map["A"] = species_name
                a = coefficient
                A0 = self.get_chem_conc(species_name)
                assert A0 is not None, f"equilibrium_concentration(): unable to proceed because the " \
                                       f"concentration of `{species_name}` was not provided"
            else:
                name_map["B"] = species_name
                b = coefficient
                B0 = self.get_chem_conc(species_name)
                assert B0 is not None, f"equilibrium_concentration(): unable to proceed because the " \
                                       f"concentration of `{species_name}` was not provided"


        # Look at the numerator of the "Reaction Quotient"
        for i, p in enumerate(products):
            # Loop over the reaction products
            species_name =  rxn.extract_species_name(p)
            rxn_order =  rxn.extract_rxn_order(p)
            coefficient = rxn.extract_stoichiometry(p)

            assert rxn_order == 1, "find_equilibrium_conc(): Currently only implemented for 1st order reactions"

            if i == 0:
                name_map["C"] = species_name
                c = coefficient
                C0 = self.get_chem_conc(species_name)
                assert C0 is not None, f"equilibrium_concentration(): unable to proceed because the " \
                                       f"concentration of `{species_name}` was not provided"
            else:
                name_map["D"] = species_name
                d = coefficient
                D0 = self.get_chem_conc(species_name)
                assert D0 is not None, f"equilibrium_concentration(): unable to proceed because the " \
                                       f"concentration of `{species_name}` was not provided"


        #print("Initial values: ", A0, B0, C0, D0)
        #print("Stoichiometry coefficients: ", a, b, c, d)

        '''
        The equation we saw earlier,
                [(C0 + c*m) (D0 + d*m)] / [(A0 - a*m) (B0 - b*m)]  =  K
                
        can be expanded into a standard quadratic form for the unknown m :
        
                alpha * m**2 + beta * m + gamma = 0
                
        and then solved for m    
        '''
        alpha = (c * d - K * a * b)
        beta = d * C0 + c * D0 + K * (A0 * b + B0 * a)
        gamma = C0 * D0 - K * A0 * B0

        #print(name_map)
        #print("alpha, beta, gamma : ", alpha, beta, gamma)

        if alpha == 0:
            # The quadratic reduces to the linear equation:  beta * m + gamma = 0
            m = -gamma / beta
        else:
            sqrt_discriminant = math.sqrt(beta**2 - 4 * alpha * gamma)
            m1 = (-beta + sqrt_discriminant) / (2 * alpha)
            m2 = (-beta - sqrt_discriminant) / (2 * alpha)
            #print("m1, m2 : ", m1, m2)
            m = m1  # Let's start with one of the 2 possible solutions of the quadratic equation


        # After m "moles of forward reaction", the concentration of the reactant "A"
        # in aA + bB <-> cC + dD gets reduced by a*m . Likewise for the other terms.
        # Reaction products get increased.  Values for missing terms will be meaningless
        std_result = {"A" : A0 - a*m, "B" : B0 - b*m, "C" : C0 + c*m, "D" : D0 + d*m}

        if min(std_result.values()) < 0:    # If there's any negative value in the concentrations...
            # ...then repeat the computation using the other solution to the quadratic
            m = m2
            std_result = {"A" : A0 - a*m, "B" : B0 - b*m, "C" : C0 + c*m, "D" : D0 + d*m}


        # Let's translate our standard names A, B, C, D into the actual names,
        # and also drop any missing term
        result = {}
        for k, v in std_result.items():
            actual_name = name_map.get(k)
            if actual_name:     # Missing terms will get dropped out
                result[actual_name] = v

        return result



    def estimate_rate_constants_TODO(self, df, t, reactants, products):
        """
        TODO: not yet implemented

        EXAMPLE:  estimate_rate_constants(df, t="SYSTEM TIME", reactants=["A", "B"], products="C")
                    where df is a Pandas dataframe containing the columns "SYSTEM TIME", "A", "B" and "C"

        :param df:
        :param t:
        :param reactants:
        :param products:
        :return:
        """
        pass



    def reach_threshold(self, chem :str, threshold) -> Union[float, None]:
        """

        :param chem:        The name of the chemical of interest
        :param threshold:   A number with the concentration value being sought
        :return:            The time at which the linearly-interpolated concentration of the specified chemical
                                first reaches the given threshold;
                                or None if it never does (in which case a message is printed)
        """
        # Prepare a Pandas dataframe with 2 columns
        df = self.get_history(columns=["SYSTEM TIME", chem])

        x_intersection = Numerical.reach_threshold(df, x="SYSTEM TIME", y=chem, y_threshold=threshold)
        if x_intersection is None:
            print(f"reach_threshold(): the concentrations of `{chem}` never reaches the specified threshold of {threshold}")

        return x_intersection



    def curve_intersect(self, chem1 :str, chem2 :str, t_start, t_end, explain=False) -> (float, float):
        """
        Find and return the intersection of the concentrations of the chemicals chem1 and chem2,
        in the time interval [t_start, t_end]
        If more than one is present, only the first (smallest value of x-coord) one will be returned;
        so, the specified time interval should be narrow enough to bracket the intersection of interest

        :param chem1:   The name of the 1st chemical of interest
        :param chem2:   The name of the 2nd chemical of interest
        :param t_start: The start of the time interval being considered
        :param t_end:   The end of the time interval being considered
        :param explain: [OPTIONAL] If True, print out some details of the computation
        :return:        The pair (time of intersection, common value) ;
                            if not found, a message is printed, and None is returned
        """
        # Prepare a Pandas dataframe with 3 columns
        df = self.get_history(t_start=t_start, t_end=t_end, columns=["SYSTEM TIME", chem1, chem2])

        intersection = Numerical.curve_intersect(df, x="SYSTEM TIME", var1=chem1, var2=chem2, explain=explain)
        if intersection is None:
            print(f"curve_intersect(): No intersection detected between the concentrations of `{chem1}` and `{chem2}` "
                  f"in the time interval [{t_start} - {t_end}]")

        return intersection



    def extract_delta_concentrations(self, df, row_from :int, row_to :int, chem_list: [str]) -> np.array:
        """
        Extract the concentration changes of the specified chemical species from a Pandas dataframe
        of concentration values

        EXAMPLE:  extract_delta_concentrations(my_dataframe, 7, 8, ['A', 'B'])

        :param df:          A Pandas dataframe of concentration values (it MUST contain columns
                                with the names given in chem_list
        :param row_from:    Row number of the first row of data we're interested in
        :param row_to:      Row number of the last row of data we're interested in
        :param chem_list:   A list of names of chemicals

        :return:            A Numpy array of floats
        """
        # TODO: add validations
        from_values = df.loc[row_from][chem_list]
        to_values = df.loc[row_to][chem_list]
        return (to_values - from_values).astype(float).to_numpy(dtype='float32')
