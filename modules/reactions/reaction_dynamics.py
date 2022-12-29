import math
import numpy as np
import pandas as pd
from typing import Union
from modules.movies.movies import MovieTabular


class ReactionDynamics:
    """
    Used to simulate the dynamics of reactions (in a single compartment.)
    In the context of Life123, this may be thought of as a zero-dimensional system.
    """

    def __init__(self, reaction_data):
        """
        
        :param reaction_data:   Object of type "ReactionData" (with data about the chemicals and their reactions)
                                TODO: maybe offer an option to let the constructor instantiate that object?
        """
        self.reaction_data = reaction_data

        self.system = None  # Concentration data in the single compartment we're simulating, for all the chemicals
                            # A Numpy array of the initial concentrations of all the chemical species, in their index order
                            # 1-dimensional NumPy array of floats, of size: n_species
                            # Each entry is the concentration of the species with that index (in the "ReactionData" object)
                            # Note that this is the counterpart - with 1 less dimension - of the array by the same name
                            #       in the class BioSim1D

        self.system_time = 0.   # Global time of the system, from initialization on

        self.history = MovieTabular()   # To store user-selected snapshots of (some of) the chemical concentrations,
                                        #   whenever requested by the user.

        self.verbose_list = []          # A list of integers with the codes of the desired verbose checkpoints
                                        #   EXAMPLE: [1, 3] to invoke sections of code marked as 1 or 3
                                        #   Those sections will have entry points such as "if 1 in self.verbose_list"
                                        #   MAX code currently used: 5
        self.debug_data = MovieTabular(parameter_name="TIME")

        self.reaction_speeds = {}       # A dictionary.  EXAMPLE : { 1: "F", 4: "S", 5: "F" , 8: "S" }
                                        #   Anything missing is regarded to be "F" (Fast)




    #############################################################################################
    #                                                                                           #
    #                                ~  TO SET/READ DATA  ~                                     #
    #                                                                                           #
    #############################################################################################

    def set_conc(self, conc: Union[list, tuple], snapshot=False) -> None:
        """
        Set the concentrations of all the chemicals

        :param conc:        A list or tuple (TODO: also allow a Numpy array; make sure to do a copy() to it!)
        :param snapshot:    (OPTIONAL) boolean: if True, add to the history
                            a snapshot of this initial state being set.  Default: False
        :return:            None
        """
        # TODO: more validations, incl. of individual values being wrong data type

        assert len(conc) == self.reaction_data.number_of_chemicals(), \
            f"ReactionDynamics.set_conc(): The number of passed concentration values ({len(conc)}) " \
            f"must match the number of declared chemicals ({self.reaction_data.number_of_chemicals()})"

        assert min(conc) >= 0, \
            f"ReactionDynamics.set_conc(): Values meant to be chemical concentrations cannot be negative " \
            f"(such as the passed value {min(conc)})"


        self.system = np.array(conc)

        if snapshot:
            self.add_snapshot(caption="Initial state")



    def get_system_conc(self) -> np.array:
        """
        Retrieve the concentrations of ALL the chemicals as a Numpy array

        :return:        Either a Numpy array with the concentrations of ALL the chemicals, in their index order
                            EXAMPLE:  array([12.3, 4.56, 0.12])    The 0-th chemical has concentration 12.3, and so on...
        """
        return self.system



    def get_conc_dict(self, species=None, system_data=None) ->[]:
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
            assert system_data.size == self.reaction_data.number_of_chemicals(),\
                f"ReactionDynamics.get_conc_dict(): the argument `system_data` must be a 1-D Numpy array with as many entries " \
                f"as the declared number of chemicals ({self.reaction_data.number_of_chemicals()})"


        if species is None:
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
        Get the code about reaction speed that the given reaction is marked with

        :param rxn_index:   The index (0-based) to identify the reaction of interest
        :return:            A 1-letter string with the code "F" (for Fast) or "S" (Slow)
        """
        return self.reaction_speeds.get(rxn_index, "F")     # Any missing entry is regarded to be "F" (Fast)


    def set_rxn_speed(self, rxn_index: int, speed: str) -> None:
        """
        Set the code about reaction speed for the given reaction

        :param rxn_index:   The index (0-based) to identify the reaction of interest
        :param speed:       A 1-letter string with the code "F" (for Fast) or "S" (Slow)
        :return:            None
        """
        assert speed in ["S", "F"], "speed argument must be either 'S' or 'F'"
        self.reaction_data.assert_valid_rxn_index(rxn_index)
        self.reaction_speeds[rxn_index] = speed




    #############################################################################################
    #                                                                                           #
    #                                   TO VISUALIZE SYSTEM                                     #
    #                                                                                           #
    #############################################################################################

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

    def add_snapshot(self, species=None, caption="", time=None, system_data=None) -> None:
        """
        Preserve some or all the chemical concentrations into the history,
        linked to the passed time (by default the current System Time),
        with an optional caption.

        EXAMPLE:  add_snapshot(species=['A', 'B'])
                                  caption="Just prior to infusion")

        :param species:     (OPTIONAL) list of name of the chemical species whose concentrations we want to preserve for later use.
                                If not specified, save all
        :param caption:     (OPTIONAL) caption to attach to this preserved data
        :param time:        (OPTIONAL) time value to attach to the snapshot (default: current System Time)
        :param system_data: (OPTIONAL) a Numpy array of concentration values, in the same order as the
                                index of the chemical species; by default, use the SYSTEM DATA
                                (which is set and managed by various functions)
        :return:            None
        """

        data_snapshot = self.get_conc_dict(species=species, system_data=system_data)


        if time is None:
            time = self.system_time

        self.history.store(par=time,
                           data_snapshot=data_snapshot, caption=caption)



    def get_history(self) -> pd.DataFrame:
        """
        Retrieve and return a Pandas dataframe with the system history that had been saved

        using add_snapshot()

        :return:        a Pandas dataframe
        """
        return self.history.get()




    #############################################################################################
    #                                                                                           #
    #                                TO PERFORM THE REACTIONS                                   #
    #                                                                                           #
    #############################################################################################

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



    def single_compartment_react(self, reaction_duration=None, time_step=None, n_steps=None,
                                 snapshots=None, dynamic_step=1, fast_threshold=5) -> None:
        """
        Perform ALL the reactions in the single compartment -
        based on the INITIAL concentrations,
        which are used as the basis for all the reactions.

        Update the system state and the system time accordingly

        :param reaction_duration:  The overall time advance for the reactions (i.e. time_step * n_steps)
                                TODO: maybe also offer a "final_time" (or "stop_time") option
        :param time_step:       The size of each time step.
                                Note: if a dynamic_step > 1 is passed, then the time step will get internally
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

        :param dynamic_step:    An integer >= 1.  If > 1, individual steps may get divided by that factor,
                                    on a reaction-by-reaction basis,
                                    if that reaction has fast dynamics,
                                    or multiplied by that factor, if that reaction has slow dynamics
        :param fast_threshold:  The minimum relative size of the concentration baseline over its change, AS A PERCENTAGE,
                                for a reaction to be regarded as "Slow".  IMPORTANT: this refer to the FULL step size

        :return:                None.   The object attributes self.system and self.system_time get updated
        """
        time_step, n_steps = self.specify_steps(total_duration=reaction_duration,
                                                time_step=time_step,
                                                n_steps=n_steps)

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


        assert self.system is not None, "ReactionDynamics.single_compartment_react(): " \
                                        "the concentration values of the various chemicals must be set first"

        for i in range(n_steps):
            delta_concentrations = self.reaction_step_orchestrator(delta_time_full=time_step, conc_array=self.system,
                                                                   snapshots=snapshots,
                                                                   dynamic_step=dynamic_step, fast_threshold=fast_threshold)
            self.system += delta_concentrations
            self.system_time += time_step
            # Preserve some of the data, as requested
            if snapshots and ((i+1)%frequency == 0):
                if first_snapshot and "initial_caption" in snapshots:
                    self.add_snapshot(species=species, caption=snapshots["initial_caption"])
                    first_snapshot = False
                else:
                    self.add_snapshot(species=species)

        if snapshots and "final_caption" in snapshots:
            self.history.set_caption_last_snapshot(snapshots["final_caption"])



    def reaction_step_orchestrator(self, delta_time_full: float, conc_array,
                                   snapshots=None, dynamic_step=1, fast_threshold=5) -> np.array:
        """
        This is the common entry point for both single-compartment reactions,
        and the reaction part of reaction-diffusions in 1D, 2D and 3D.

        "Compartments" may or may not correspond to the "bins" of the higher layers;
        the calling code might have opted to merge some bins into a single "compartment"

        Using the given concentration data for all the applicable species in a single compartment,
        do a single reaction time step for ALL the reactions -
        based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        Return the increment vector for all the chemical species concentrations in the compartment

        NOTE:  the actual system concentrations are NOT changed

        :param delta_time_full: The requested time duration of the overall reaction step,
                                    which will be subdivided for the "fast" reactions, if dynamic_step is greater than 1
        :param conc_array:      All initial concentrations at the start of the reaction step,
                                    as a Numpy array for all the chemical species, in their index order;
                                    this can be thought of as the "SYSTEM STATE"
        :param snapshots:       See explanation under single_compartment_react()
        :param dynamic_step:    An integer by which to subdivide the time intervals for "fast" reactions;
                                    default 1, i.e. no subdivisions
        :param fast_threshold:  The minimum relative size of the concentration baseline over its change, AS A PERCENTAGE,
                                for a reaction to be regarded as "Slow".  IMPORTANT: this refer to the FULL step size

        :return:                The increment vector for the concentrations of ALL the chemical species,
                                    as a Numpy array for all the chemical species, in their index order
                                    EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):   array([7. , -21.])
        """
        if 5 in self.verbose_list:
            print(f"Calling reaction_step_orchestrator() with delta_time={delta_time_full}, "
                  f"conc_array={conc_array}, dynamic_step={dynamic_step}")

        assert type(dynamic_step) == int, "reaction_step_orchestrator(): the argument 'dynamic_step' must be an integer"
        assert dynamic_step >= 1, "reaction_step_orchestrator(): the argument 'dynamic_step' must be an integer greater or equal than 1"
        assert fast_threshold > 0, "reaction_step_orchestrator(): the argument 'fast_threshold' must be greater than 0"
        assert conc_array is not None, "reaction_step_orchestrator(): the argument 'conc_array' must be set to a Numpy array"

        fast_threshold_fraction = fast_threshold / 100.
        if 2 in self.verbose_list:
            print(f"************ SYSTEM TIME: {self.system_time:,.4g}")

        if dynamic_step == 1 or self.are_all_slow_rxns():
            if 2 in self.verbose_list:
                if dynamic_step == 1:
                    print("    NO adaptive variable time resolution used")
                else:
                    print("    All the reactions are SLOW")

                print(f"    Processing ALL the {self.reaction_data.number_of_reactions()} reaction(s)")

            delta_concentrations = self.single_reaction_step(delta_time=delta_time_full, time_subdivision=1, fast_threshold_fraction=fast_threshold_fraction,
                                                             conc_array=conc_array, rxn_list=None)
        else:
            # Process all the slow reactions
            slow_rxns = self.slow_rxns()
            if 2 in self.verbose_list:
                if slow_rxns == []:
                    print(f"    There are NO slow reactions")
                else:
                    print(f"    Slow reactions: {slow_rxns}")
                    print(f"    Processing SLOW reactions")

            if slow_rxns == []:
                # If there are no slow reactions
                delta_concentrations_slow = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)       # One element per chemical species
            else:
                # Process all the slow reactions
                delta_concentrations_slow = self.single_reaction_step(delta_time=delta_time_full, time_subdivision=1, fast_threshold_fraction=fast_threshold_fraction,
                                                                  conc_array=conc_array, rxn_list=slow_rxns)

            # Process all the fast reactions
            fast_rxns = self.fast_rxns()
            if 2 in self.verbose_list:
                print(f"    Fast reactions: {fast_rxns}")

            local_conc_array = conc_array.copy()   # Saved as an unchanging baseline copy
            delta_concentrations_fast = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)     # One element per chemical species
            reduced_time_step = delta_time_full / dynamic_step
            for substep in range(dynamic_step):
                if 2 in self.verbose_list:
                    print(f"    - substep: {substep} (in processing of FAST reactions)")

                incr_vector = self.single_reaction_step(delta_time=reduced_time_step, time_subdivision=dynamic_step, fast_threshold_fraction=fast_threshold_fraction,
                                                        conc_array=local_conc_array, rxn_list=fast_rxns, substep_number=substep)
                local_conc_array += incr_vector
                delta_concentrations_fast += incr_vector
                # Preserve the intermediate-state data, if requested
                if snapshots and snapshots.get("show_intermediates") and substep < dynamic_step-1:  # Skip the last one, which will be handled by the caller
                    species_to_show = snapshots.get("species")          # If not present, it will be None (meaning show all)
                    time = self.system_time + (substep+1) * reduced_time_step
                    self.add_snapshot(time=time, species=species_to_show,
                                      system_data=local_conc_array,
                                      caption=f"Interm. step, due to the fast rxns: {fast_rxns}")

            # Combine the results of all the slow reactions and all the fast reactions
            if 1 in self.verbose_list:
                print("delta_time: ", delta_time_full)
                print("    delta_concentrations_slow: ", delta_concentrations_slow)
                print("    delta_concentrations_fast: ", delta_concentrations_fast)

            delta_concentrations = delta_concentrations_slow + delta_concentrations_fast
        # END IF

        return  delta_concentrations



    def single_reaction_step(self, delta_time, conc_array, time_subdivision=1, fast_threshold_fraction=0.05,
                             rxn_list=None, substep_number=0) -> np.array:
        """
        Using the given concentration data of ALL the chemical species,
        do the specified SINGLE TIME STEP for ONLY the requested reactions (by default all).

        NOTE: this time step might be the "full" time step of the calling function,
              or a substep (if adaptive variable time resolution is being used)

        All computations are based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions (in "forward Euler" approach.)

        Return the Numpy increment vector for ALL the chemical species concentrations, in their index order
        (whether involved in these reactions or not)

        NOTE: the actual system concentrations are NOT changed

        :param delta_time:              The time duration of this individual reaction step - assumed to be small enough that the
                                            concentration won't vary significantly during this span.
                                            NOTE: this may be a "full step" or a "substep", depending on the adaptive time scale
                                                  being used by the caller function
        :param time_subdivision:        Integer with the number of subdivisions currently being used for the "full" time step
                                            of the caller function;
                                            if > 1, it means we're using an adaptive variable time resolution.
                                            (Note: the "full" time step of the calling function will be delta_time * time_subdivision)
        :param fast_threshold_fraction: The minimum relative size of the concentration baseline over its change, AS A FRACTION,
                                            for a reaction to be regarded as "Slow".  IMPORTANT: this refers to the FULL step size

        :param conc_array:              All initial concentrations at the start of the reaction step,
                                            as a Numpy array for all the chemical species, in their index order;
                                            this can be thought of as the "SYSTEM STATE"
        :param rxn_list:                OPTIONAL list of reactions (specified by index) to include in this simulation step ; EXAMPLE: [1, 3, 7]
                                            If None, do all the reactions.
        :param substep_number:          Only used for variable time resolution mode.  Zero-based counting

        :return:            The increment vector for the concentrations of ALL the chemical species,
                                as a Numpy array for all the chemical species, in their index order
                            EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):   array([7. , -21.])
        """
        increment_vector = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)       # One element per chemical species

        # Compute the forward and back "conversions" of all the applicable reactions
        delta_list = self.compute_all_reaction_deltas(conc_array=conc_array, delta_time=delta_time, rxn_list=rxn_list)
        if 3 in self.verbose_list:
            print(f"delta_list: {delta_list}")


        if rxn_list is None:    # Meaning ALL reactions
            rxn_list = range(self.reaction_data.number_of_reactions())  # This will be a list of all the reaction index numbers


        # For each applicable reaction, adjust the concentrations of the reactants and products,
        # based on the forward and back rates of the reaction
        for rxn_index in rxn_list:      # Consider each reaction in turn
            if 2 in self.verbose_list:
                print(f"      Determining the concentrations changes as a result of reaction number {rxn_index}")

            # TODO: turn into a more efficient single step, as as:
            #(reactants, products) = cls.all_reactions.unpack_terms(i)
            reactants = self.reaction_data.get_reactants(rxn_index)
            products = self.reaction_data.get_products(rxn_index)

            # Determine the concentration adjustments (and label the reaction as "Slow" or "Fast" accordingly)

            # TODO: NO NEED TO COMPUTE THE REACTION SPEEDS, EXCEPT IN MAIN STEPS and IN FINAL SUBSTEPS!!
            self.set_rxn_speed(rxn_index, "S")  # TENTATIVE assignment, that will be changed
                                                #   if ANY chemical experiences significant concentration changes

            #   The reactants decrease based on the (forward reaction - reverse reaction)
            for r in reactants:
                stoichiometry, species_index, order = r
                delta_conc = stoichiometry * (- delta_list[rxn_index])  # Increment to this reactant from the reaction being considered
                # Do a validation check to avoid negative concentrations
                self.validate_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                # If we are at the last of the substeps...  TODO: skip if caller function not using variable time resolution mode; however, this CANNOT be detected by just looking at time_subdivision == 1
                if substep_number == time_subdivision-1:
                    # Mark the reaction "fast", if appropriate
                    self.examine_increment(rxn_index=rxn_index, species_index=species_index,
                                           delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                           time_subdivision=time_subdivision,
                                           fast_threshold_fraction=fast_threshold_fraction)
                elif 2 in self.verbose_list:
                    print(f"        Skipping evaluation of reaction speed for chemical `{self.reaction_data.get_name(species_index)}` "
                          f"(substep_number = {substep_number}, time_subdivision = {time_subdivision})")

                increment_vector[species_index] += delta_conc


            #   The reaction products increase based on the (forward reaction - reverse reaction)
            for p in products:
                stoichiometry, species_index, order = p
                delta_conc = stoichiometry * delta_list[rxn_index]  # Increment to this reaction product from the reaction being considered
                # Do a validation check to avoid negative concentrations
                self.validate_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                # If we are at the last of the substeps...
                if substep_number == time_subdivision-1:
                    # Mark the reaction "fast", if appropriate
                    self.examine_increment(rxn_index=rxn_index, species_index=species_index,
                                           delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                           time_subdivision=time_subdivision,
                                           fast_threshold_fraction=fast_threshold_fraction)
                elif 2 in self.verbose_list:
                    print(f"        Skipping evaluation of reaction speed for chemical `{self.reaction_data.get_name(species_index)}` "
                          f"(substep_number = {substep_number}, time_subdivision = {time_subdivision})")

                increment_vector[species_index] += delta_conc

            if 2 in self.verbose_list:
                #TODO: split the reactions into separate data frames (one per reaction)?
                data_snapshot = self.get_conc_dict(system_data=increment_vector)    # This will be a dict of conc values
                data_snapshot["reaction"] = rxn_index
                data_snapshot["substep"] = substep_number
                data_snapshot["time_subdivision"] = time_subdivision
                self.debug_data.store(par=self.system_time, data_snapshot=data_snapshot, caption=f"delta_time: {delta_time}")
        # END for (over rxn_list)

        return increment_vector



    def validate_increment(self,  delta_conc, baseline_conc: float,
                           rxn_index: int, species_index: int, delta_time) -> None:
        """
        Examine the requested concentration change given by delta_conc
        (typically, as computed by an ode solver),
        relative to the baseline (pre-reaction) value baseline_conc,
        for the given SINGLE chemical species and SINGLE reaction.

        If the concentration change would render the concentration negative,
        raise an Exception

        :param delta_conc:              The change in concentration computed by the ode solver
                                            (for the specified chemical, in the given reaction)
        :param baseline_conc:           The initial concentration

        [The remaining arguments are ONLY USED for error printing]
        :param rxn_index:               The index (0-based) to identify the reaction of interest (ONLY USED for error printing)
        :param species_index:           The index (0-based) to identify the chemical species of interest (ONLY USED for error printing)
        :param delta_time:              The time duration of the reaction step. (ONLY USED for error printing)

        :return:                        None (an Exception is raised if a negative concentration is detected)
        """
        if (baseline_conc + delta_conc) < 0:
            raise Exception(f"The chosen time interval ({delta_time}) "
                            f"leads to a NEGATIVE concentration of the chemical species {species_index} from reaction {rxn_index}: "
                            f"must make the interval smaller!")



    def examine_increment(self, rxn_index: int, species_index: int,
                          delta_conc, baseline_conc: float,
                          time_subdivision: int,
                          fast_threshold_fraction) -> None:
        """
        Examine the requested concentration change given by delta_conc
        (typically, as computed by an ode solver),
        relative to the baseline (pre-reaction) value baseline_conc,
        for the given SINGLE chemical species and SINGLE reaction.

        If the concentration change is large, relative to its baseline value,
        then mark the given reaction as "Fast" (in its data structure);
        this will OVER-RIDE any previous annotation about the speed of that reaction.
        Note: it doesn't matter which of the chemicals in the reaction leads to this.

        :param rxn_index:               The index (0-based) to identify the reaction of interest
        :param species_index:           The index (0-based) to identify the chemical species of interest. (ONLY USED for error printing)

        :param delta_conc:              The change in concentration computed by the ode solver
                                            (for the above chemical, in the above reaction)
        :param baseline_conc:           The initial concentration

        :param time_subdivision:        Integer with the number of subdivisions currently used for delta_time;
                                            used for adaptive variable time resolution

        :param fast_threshold_fraction: The minimum relative size of the concentration baseline over its change, AS A FRACTION,
                                            for a reaction to be regarded as "Slow".
                                            IMPORTANT: this refer to the FULL step size, and will be scaled for time_subdivision
                                            (e.g. a time_subdivision of 2 signifies 1/2 a time step,
                                            and will mean that just 1/2 of the change will constitute a threshold)

        :return:                        None (the equation is marked as "Fast", if appropriate, in its data structure)
        """

        # If the reaction was (tentatively) labeled as "Slow", decide whether to flip it to "Fast"
        #   Note: the status will be used for the current simulation sub-cycle
        if self.get_rxn_speed(rxn_index) == "S":    # TODO: maybe skip this step if adaptive variable time resolution isn't being used
            #if abs(delta_conc) / baseline_conc > fast_threshold_fraction / time_subdivision # To avoid time-consuming divisions, re-formulated as below:
            if abs(delta_conc) * time_subdivision > fast_threshold_fraction * baseline_conc:
                self.set_rxn_speed(rxn_index, "F")

            if 2 in self.verbose_list:
                print(f"        Reaction # {rxn_index} determined to be '{self.get_rxn_speed(rxn_index)}', "
                      f"based on a change of {delta_conc:.6g} (relative to its baseline of {baseline_conc:.7g}) "
                      f"for the conc. of chem species # {species_index}: "
                      f"comparing abs({100 * delta_conc / baseline_conc:.3g}%) vs. ({100 * fast_threshold_fraction:.3g}% /{time_subdivision})")



    def compute_all_reaction_deltas(self, delta_time: float, conc_array=None, rxn_list=None) -> [float]:
        """
        For an explanation of the "reaction delta", see compute_reaction_delta().
        Compute the "reaction delta" for all the specified ones (by default, all.)
        Return a list with an entry for each reaction, in their index order.

        For background info: https://life123.science/reactions

        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param conc_array:  All initial concentrations at the start of the reaction step,
                                as a Numpy array for all the chemical species, in their index order;
                                this can be thought of as the "SYSTEM STATE"
        :param rxn_list:    OPTIONAL list of reactions (specified by index);
                                if None, do all the reactions.  EXAMPLE: [1, 3, 7]

        :return:            A list of the differences between forward and reverse "conversions" -
                                for explanation, see compute_reaction_delta();
                                each list has 1 entry per reaction, in the index order of the reactions
        """
        delta_list = []         # It will have 1 entry per reaction

        if rxn_list is None:    # Meaning ALL reactions
            rxn_list = range(self.reaction_data.number_of_reactions())

        # Process the requested reactions
        for i in rxn_list:      # Consider each reaction in turn
            delta = self.compute_reaction_delta(rxn_index=i, conc_array=conc_array, delta_time=delta_time)
            delta_list.append(delta)

        return delta_list



    def compute_reaction_delta(self, rxn_index: int, delta_time: float, conc_array: np.array) -> float:
        """
        For the SINGLE specified reaction, the given time interval, and the specified concentrations of chemicals,
        compute the difference of the reaction's forward and back "conversions",
        a non-standard term we're using here to refer to delta_time * (Forward_Rate − Reverse_Rate)

        For background info: https://life123.science/reactions
        What we're computing here, is referred to as:  (Δt)∗delta_forward(n)

        :param rxn_index:   An integer that indexes the reaction of interest
        :param conc_array:  All initial concentrations at the start of the reaction step,
                                as a Numpy array for all the chemical species, in their index order;
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



    def is_in_equilibrium(self, rxn_index: int, conc: dict, tolerance=0.05, explain=True) -> bool:
        """
        Ascertain whether the given concentrations are in equilibrium for the specified reaction

        :param rxn_index:   The index (0-based) to identify the reaction of interest
        :param conc:        Dict with the concentrations of the species involved in the reaction.
                            The keys are the chemical names
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
        :param tolerance:   Allowable absolute tolerance, to establish satisfactory equality
        :param explain:     If True, print out the formula being used.
                                EXAMPLES:  "([C][D]) / ([A][B])" , "[B] / [A]^2",
        :return:            True if the reaction is close enough to an equilibrium
        """
        rxn = self.reaction_data.get_reaction(rxn_index)

        reactants = self.reaction_data.extract_reactants(rxn)     # A list of triplets
        products = self.reaction_data.extract_products(rxn)       # A list of triplets
        kF = self.reaction_data.extract_forward_rate(rxn)
        kB = self.reaction_data.extract_back_rate(rxn)

        rate_ratio = kF / kB                    # Ratio of forward/reverse reaction rates

        conc_ratio = 1.
        numerator = ""
        denominator = ""
        for rxn_index, p in enumerate(products):
            species_index = self.reaction_data.extract_species_index(p)
            rxn_order = self.reaction_data.extract_rxn_order(p)

            species_name = self.reaction_data.get_name(species_index)
            species_conc = conc[species_name]
            conc_ratio *= (species_conc ** rxn_order)
            if explain:
                numerator += f"[{species_name}]"
                if rxn_order > 1:
                    numerator += f"^{rxn_order} "

        if explain and len(products) > 1:
            numerator = f"({numerator})"

        for r in reactants:
            species_index = self.reaction_data.extract_species_index(r)
            rxn_order = self.reaction_data.extract_rxn_order(r)

            species_name = self.reaction_data.get_name(species_index)
            species_conc = conc[species_name]
            conc_ratio /= (species_conc ** rxn_order)
            if explain:
                denominator += f"[{species_name}]"
                if rxn_order > 1:
                    denominator += f"^{rxn_order} "

        if explain and len(reactants) > 1:
            denominator = f"({denominator})"

        if explain:
            print(f"Ratio of forward/reverse reaction rates: {rate_ratio}")
            print(f"Ratio of reactant/product concentrations, adjusted for reaction orders: {conc_ratio:,.6g}")
            print(f"    {numerator} / {denominator}")

        return np.allclose(conc_ratio, rate_ratio, atol=tolerance)



    def stoichiometry_checker(self, rxn_index: int, conc_arr_before: np.array, conc_arr_after: np.array) -> bool:
        """
        For the indicated reaction, investigate the change in the concentration of the involved chemicals,
        to ascertain whether the change is consistent with the reaction's stoichiometry.
        See https://life123.science/reactions

        NOTE: the concentration changes in chemicals not involved in the specified reaction are ignored

        :param rxn_index:       Integer to identify the reaction of interest
        :param conc_arr_before: Numpy array with the concentrations of ALL the chemicals (whether involved
                                    in the reaction or not), in their index order, before the reaction
        :param conc_arr_after:  Same as above, but after the reaction
        :return:                True if the change in reactant/product concentrations is consistent with the
                                    reaction's stoichiometry, or False otherwise
        """
        self.reaction_data.assert_valid_rxn_index(rxn_index)

        reactants = self.reaction_data.get_reactants(rxn_index)
        products = self.reaction_data.get_products(rxn_index)

        # Pick (arbitrarily) the first reactant,
        # to establish a baseline change in concentration relative to its stoichiometric coefficient
        baseline_term = reactants[0]
        baseline_species = self.reaction_data.extract_species_index(baseline_term)
        baseline_stoichiometry = self.reaction_data.extract_stoichiometry(baseline_term)
        baseline_ratio =  (conc_arr_after[baseline_species] - conc_arr_before[baseline_species]) / baseline_stoichiometry
        #print("baseline_ratio: ", baseline_ratio)

        for i, term in enumerate(reactants):
            if i != 0:
                species = self.reaction_data.extract_species_index(term)
                stoichiometry = self.reaction_data.extract_stoichiometry(term)
                ratio =  (conc_arr_after[species] - conc_arr_before[species]) / stoichiometry
                #print(f"ratio for `{self.reaction_data.get_name(species)}`: {ratio}")
                if not np.allclose(ratio, baseline_ratio):
                    return False

        for term in products:
            species = self.reaction_data.extract_species_index(term)
            stoichiometry = self.reaction_data.extract_stoichiometry(term)
            ratio =  - (conc_arr_after[species] - conc_arr_before[species]) / stoichiometry # The minus in front is b/c we're on the other side of the eqn
            #print(f"ratio for `{self.reaction_data.get_name(species)}`: {ratio}")
            if not np.allclose(ratio, baseline_ratio):
                return False

        return True



    def debug(self, code):      # Experimental
        """
        Usage example:  if self.debug(3):
                            # DO SOMETHING
        :param code:
        :return:
        """
        if type(self.verbose_list) == list and code in self.verbose_list:
            return True

        if type(self.verbose_list) == int and code in self.verbose_list:
            return True

        return False