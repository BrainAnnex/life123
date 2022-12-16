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

        self.debug = False

        self.reaction_speeds = {}       # A dictionary.  EXAMPLE : { 1: "F", 4: "S", 5: "F" , 8: "S" }
                                        # Anything missing is assumed to be "F" (Fast)



    #############################################################################################
    #                                                                                           #
    #                                ~  TO SET/READ DATA  ~                                     #
    #                                                                                           #
    #############################################################################################

    def set_conc(self, conc: Union[list, tuple], snapshot=False) -> None:
        """
        Set the concentrations of all the chemicals

        :param conc:        A list or tuple (TODO: also allow a Numpy array)
        :param snapshot:    (OPTIONAL)
        :return:            None
        """
        # TODO: more validations, incl. of being positive

        assert len(conc) == self.reaction_data.number_of_chemicals(), \
            f"set_conc(): The number of concentrations ({len(conc)}) " \
            f"must match the number of declared chemicals ({self.reaction_data.number_of_chemicals()})"

        self.system = np.array(conc)

        if snapshot:
            self.add_snapshot(caption="Initial state")



    def get_conc(self, form="DICT"):
        """
        Retrieve the concentrations of ALL the chemicals

        :param form:    Either "ARRAY"  (EXAMPLE of returned value: array([12.3, 4.56]))
                            or "DICT"   (EXAMPLE: {"X": 12.3, "Y": 4.56}, where the keys are the names of the chemicals)
        :return:        Either a Numpy array or a dictionary
        """
        if form == "ARRAY":
            return self.system
        elif form == "DICT":
            return {self.reaction_data.get_name(index): self.system[index]
                                                                for index, conc in enumerate(self.system)}
        else:
            raise Exception(f"get_conc(): Unknown option for the `form` argument ({form}).  Allowed values are 'ARRAY' and 'DICT'")




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
        print(f"SYSTEM STATE at Time t = {self.system_time}:")

        n_species = self.reaction_data.number_of_chemicals()
        print(f"{n_species} species:")

        # Show a line of line of data for each chemical species in turn
        for species_index, name in enumerate(self.reaction_data.get_all_names()):
            if name:    # If a name was provided, show it
                name = f" ({name})"
            else:
                name = ""

            print(f"  Species {species_index}{name}.Conc: {self.system[species_index]}")




    #############################################################################################
    #                                                                                           #
    #                                       HISTORY                                             #
    #                                                                                           #
    #############################################################################################

    def add_snapshot(self, species=None, caption="") -> None:
        """
        Preserve some or all the chemical concentrations into the history, linked to the
        current System Time, with an optional caption.

        EXAMPLE:  create_snapshot(species=['A', 'B'])
                                  caption="Just prior to infusion")

        :param species: (OPTIONAL) list of name of the chemical species whose concentrations we want to preserve for later use.
                            If not specified, save all
        :param caption: (OPTIONAL) caption to attach to this preserved data
        :return:        None
        """
        if species is None:
            data_snapshot = self.get_conc(form="DICT")
        else:
            data_snapshot = {}
            for species_index, name in enumerate(species):
                data_snapshot[name] = self.system[species_index]

        self.history.store(par=self.system_time,
                           data_snapshot = data_snapshot, caption=caption)



    def get_history(self) -> pd.DataFrame:
        """
        Retrieve and return a Pandas dataframe with the system history that had been saved

        using save_snapshot()

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



    def single_compartment_react(self, total_duration=None, time_step=None, n_steps=None,
                                 dynamic_step = 1,
                                 snapshots=None) -> None:
        """
        Perform ALL the reactions in the single compartment -
        based on the INITIAL concentrations,
        which are used as the basis for all the reactions.

        NOTES: When calling this method in the context of 1D, 2D or 3D systems (as opposed to single-compartment experiments),
                    "compartments" may or may not correspond to the "bins" of the higher layers;
                    the calling code might have opted to merge some bins into a single "compartment"

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param dynamic_step:    An integer >= 1.  If > 1, individual steps may get divided by that factor,
                                    on a reaction-by-reaction basis,
                                    if that reaction has fast dynamics,
                                    or multiplied by that factor, if that reaction has slow dynamics
        :param snapshots:       OPTIONAL dict that may contain any the following keys:
                                        -"frequency" (default 1)
                                        -"species" (default all)
                                        -"initial_caption" (default blank)
                                        -"final_caption" (default blank)
                                    If provided, take a system snapshot after running a multiple
                                    of "frequency" reaction steps (default 1, i.e. at every step.)
                                    EXAMPLE: snapshots={"frequency": 2, "species": ["A", "H"]}

        :return:                None
        """
        time_step, n_steps = self.specify_steps(total_duration=total_duration,
                                                time_step=time_step,
                                                n_steps=n_steps)

        if snapshots:
            frequency = snapshots.get("frequency", 1)   # Default is 1
            species = snapshots.get("species")          # Default is None
            first_snapshot = True


        for i in range(n_steps):
            #self.reaction_step_orchestrator(delta_time=time_step, dynamic_step=dynamic_step)   # TODO: new approach, to replace the next 2 lines
            delta_concentrations = self.single_compartment_reaction_step(conc_array=self.system, delta_time=time_step, dynamic_step=dynamic_step)
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


    def mark_rxn_speed(self, rxn_index: int, speed: str) -> None:
        assert speed in ["S", "F"], "speed argument must be either 'S' or 'F'"
        self.reaction_data.assert_valid_rxn_index(rxn_index)
        self.reaction_speeds[rxn_index] = speed


    def reaction_step_orchestrator(self, delta_time: float, conc_array, dynamic_step=1) -> np.array:
        """
        TODO: This method will become the common entry point for both single-compartment reactions,
              and the reaction part of reaction-diffusions
        """

        if dynamic_step == 1 or self.are_all_slow_rxns():
            #conc_array = self.system
            delta_concentrations = self.single_reaction_step_NEW(delta_time=delta_time, conc_array=conc_array, rxn_list=None)
            #self.system += delta_concentrations
        else:
            # Process all the slow reactions
            slow_rxns = self.slow_rxns()
            #conc_array = self.system.copy()
            delta_concentrations_slow = self.single_reaction_step_NEW(delta_time=delta_time, conc_array=conc_array, rxn_list=slow_rxns)

            # Process all the fast reactions
            fast_rxns = self.fast_rxns()
            local_init_conc_array = conc_array.copy()   # Saved as an unchanging baseline copy
            #conc_array = self.system.copy()
            delta_concentrations_fast = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)       # One element per chemical species
            for i in range(dynamic_step):
                incr_vector = self.single_reaction_step_NEW(delta_time=delta_time, conc_array=local_init_conc_array, rxn_list=fast_rxns)
                local_init_conc_array += incr_vector
                delta_concentrations_fast += incr_vector

            # Combine the results of the slow and fast reactions
            delta_concentrations = delta_concentrations_slow + delta_concentrations_fast - local_init_conc_array
            # The above is a simplification of:
            # self.system = self.system + (delta_concentrations_slow - self.system) + (delta_concentrations_fast - self.system)
        # END IF

        return  delta_concentrations



    def single_reaction_step_NEW(self, delta_time: float, conc_array=None, rxn_list=None) -> np.array:
        """
        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param conc_array:  All concentrations at the start of the reaction step,
                                as a Numpy array of the initial concentrations of all the chemical species, in their index order (the "SYSTEM STATE")
        :param rxn_list:    OPTIONAL list of reactions (specified by index); EXAMPLE: [1, 3, 7]
                                If None, do all the reactions.

        :return:            The increment vector for all the chemical species concentrations
                            EXAMPLE (for a reactant and product with a 3:1 stoichiometry):   [7. , -21.]
        """
        # Compute the forward and back "conversions" of all the applicable reactions
        delta_list = self.compute_all_reaction_deltas(conc_array=conc_array, delta_time=delta_time, rxn_list=rxn_list)
        if self.debug:
            print(f"    delta_list: {delta_list}")


        increment_vector = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)       # One element per chemical species


        if rxn_list is None:    # Meaning ALL reactions
            rxn_list = range(self.reaction_data.number_of_reactions())


        # For each applicable reaction, adjust the concentrations of the reactants and products,
        # based on the forward and back rates of the reaction
        for rxn_index in rxn_list:      # Consider each reaction in turn
            self.reaction_speeds[rxn_index] = "S"   # Tentative assignment, to be changed if ANY species experiences significant concentration changes

            if self.debug:
                print(f"    adjusting the species concentrations based on reaction number {rxn_index}")

            # TODO: turn into a more efficient single step, as as:
            #(reactants, products) = cls.all_reactions.unpack_terms(i)
            reactants = self.reaction_data.get_reactants(rxn_index)
            products = self.reaction_data.get_products(rxn_index)

            # Determine the concentration adjustments

            #   The reactants decrease based on the (forward reaction - reverse reaction)
            for r in reactants:
                stoichiometry, species_index, order = r
                delta_conc = stoichiometry * (- delta_list[rxn_index])  # Increment to this reactant from the reaction being considered

                self.examine_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                       species_index=species_index, rxn_index=rxn_index, delta_time=delta_time)

                '''
                if (conc_array[species_index] + delta_conc) < 0:
                    raise Exception(f"The given time interval ({delta_time}) "
                                    f"leads to negative concentrations in reactant {species_index} in reaction {rxn_index}: make it smaller!")
                
                if delta_conc / conc_array[species_index] > 0.05:
                    self.reaction_speeds[rxn_index] = "F"
                '''
                increment_vector[species_index] += delta_conc


            #   The reaction products increase based on the (forward reaction - reverse reaction)
            for p in products:
                stoichiometry, species_index, order = p
                delta_conc = stoichiometry * delta_list[rxn_index]  # Increment to this reaction product from the reaction being considered

                self.examine_increment(delta_conc=delta_conc, baseline_conc=conc_array[species_index],
                                       species_index=species_index, rxn_index=rxn_index, delta_time=delta_time)

                '''
                if (conc_array[species_index] + delta_conc) < 0:
                    raise Exception(f"The given time interval ({delta_time}) "
                                    f"leads to negative concentrations in reaction products {species_index} in reaction {rxn_index}: make it smaller!")

                if delta_conc / conc_array[species_index] > 0.05:
                    self.reaction_speeds[rxn_index] = "F"
                '''
                increment_vector[species_index] += delta_conc
        # END for

        return increment_vector



    def examine_increment(self, delta_conc: float, baseline_conc: float, species_index: int, rxn_index: int, delta_time) -> None:
        """
        Examine the computed concentration value delta_conc, relative to the baseline (pre-reaction) value baseline_conc,
        for the given chemical species and reaction.

        :param delta_conc:
        :param baseline_conc:
        :param species_index:
        :param rxn_index:
        :param delta_time:
        :return:
        """
        THRESHOLD = 0.05
        if (baseline_conc + delta_conc) < 0:
            raise Exception(f"The given time interval ({delta_time}) "
                            f"leads to negative concentrations of the chemical species {species_index} in reaction {rxn_index}: must make it smaller!")

        if delta_conc / baseline_conc > THRESHOLD:
            self.mark_rxn_speed(rxn_index, "F")
            if self.debug:
                print(f"    Reaction number {rxn_index} marked as 'Fast'")



    def single_compartment_reaction_step(self, delta_time: float, conc_dict=None, conc_array=None, dynamic_step=1) -> np.array: # TODO: phase out, in favor of single_reaction_step_NEW()
        """
        Using the given concentration data for all the applicable species in a single compartment,
        do a single reaction time step for ALL the reactions -
        based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        Return the increment vector for all the chemical species concentrations in the compartment

        NOTES:  - the actual system concentrations are NOT changed
                - "compartments" may or may not correspond to the "bins" of the higher layers;
                  the calling code might have opted to merge some bins into a single "compartment"

        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                            concentration won't vary significantly during this span

        :param conc_dict:   Concentrations of the applicable chemicals,
                            as a dict where the key value is the chemicals index
                            EXAMPLE: {3: 16.3, 8: 0.53, 12: 1.78}
        :param conc_array:  ALTERNATE way to specify all concentrations,
                            as a Numpy array of the initial concentrations of all the chemical species, in their index order
                            NOTE: if both conc_dict and conc_array are specified, an Exception is raised
        :param dynamic_step: TODO: not yet in use

        :return:            The increment vector for all the chemical species concentrations
                            in the compartment
                            EXAMPLE (for a reactant and product with a 3:1 stoichiometry):   [7. , -21.]
        """
        if conc_array is not None:
            assert conc_dict is None,\
                "single_compartment_reaction_step(): Cannot specify both arguments conc_dict and conc_array"

            conc_dict = {}
            for index in range(len(conc_array)):
                conc_dict[index] = conc_array[index]


        # Compute the forward and back conversions of all the reactions
        delta_list = self.compute_all_reaction_deltas(conc_dict=conc_dict, delta_time=delta_time)
        if self.debug:
            print(f"    delta_list: {delta_list}")


        increment_vector = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)       # One element per chemical species

        # For each reaction, adjust the concentrations of the reactants and products,
        # based on the forward and back rates of the reaction
        for i in range(self.reaction_data.number_of_reactions()):   # For each reaction
            if self.debug:
                print(f"    adjusting the species concentrations based on reaction number {i}")

            # TODO: turn into a more efficient single step, as as:
            #(reactants, products) = cls.all_reactions.unpack_terms(i)
            reactants = self.reaction_data.get_reactants(i)
            products = self.reaction_data.get_products(i)

            # Determine the concentration adjustments

            #   The reactants decrease based on the forward reaction,
            #             and increase based on the reverse reaction
            for r in reactants:
                stoichiometry, species_index, order = r
                increment_vector[species_index] += stoichiometry * (- delta_list[i])

                if (conc_dict[species_index] + increment_vector[species_index]) < 0:    # TODO: maybe do at the very end
                    raise Exception(f"The given time interval ({delta_time}) "
                                    f"leads to negative concentrations in reactants: make it smaller!")


            #   The products increase based on the forward reaction,
            #             and decrease based on the reverse reaction
            for p in products:
                stoichiometry, species_index, order = p
                increment_vector[species_index] += stoichiometry * delta_list[i]

                if (conc_dict[species_index] + increment_vector[species_index]) < 0:    # TODO: maybe do at the very end
                    raise Exception(f"The given time interval ({delta_time}) "
                                    f"leads to negative concentrations in reaction products: make it smaller!")

        # END for

        return increment_vector



    def compute_all_reaction_deltas(self, delta_time: float, conc_dict=None, conc_array=None, rxn_list=None) -> [float]:
        """
        For an explanation of the "reaction delta", see compute_reaction_delta().
        Compute the "reaction delta" for all the reactions, or for all the specified ones.
        Return a list with an entry for each reaction

        For background info: https://life123.science/reactions

        :param conc_dict:   Concentrations of the applicable chemicals,
                            as a dict where the key value is the chemical's index
                            EXAMPLE: {3: 16.3, 8: 0.53, 12: 1.78}       # TODO: phase out
        :param conc_array:  ALTERNATE way to specify all concentrations,
                            as a Numpy array of the initial concentrations of all the chemical species, in their index order
                            NOTE: if both conc_dict and conc_array are specified, an Exception is raised
        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                            concentration won't vary significantly during this span
        :param rxn_list:    OPTIONAL list of reactions (specified by index);
                                if None, do all the reactions.  EXAMPLE: [1, 3, 7]

        :return:            A list of the differences between forward and reverse "conversions" -
                                for explanation, see compute_reaction_delta();
                                each list has 1 entry per reaction, in the index order of the reactions
        """
        if conc_array is not None:
            assert conc_dict is None, \
                "single_compartment_reaction_step(): Cannot specify both arguments conc_dict and conc_array"

            conc_dict = {}
            for index in range(len(conc_array)):
                conc_dict[index] = conc_array[index]


        delta_list = []         # It will have 1 entry per reaction

        if rxn_list is None:    # Meaning ALL reactions
            rxn_list = range(self.reaction_data.number_of_reactions())

        # Process the requested reactions
        for i in rxn_list:      # Consider each reaction in turn
            delta = self.compute_reaction_delta(rxn_index=i, conc_dict=conc_dict, delta_time=delta_time)
            delta_list.append(delta)

        return delta_list



    def compute_reaction_delta(self, rxn_index: int, conc_dict: dict, delta_time: float) -> float:
        """
        For the given time interval, the SINGLE specified reaction, and the specified concentrations of chemicals,
        compute the difference of the reaction's forward and back "conversions",
        a non-standard term we're using here to refer to delta_time * (Forward_Rate − Reverse_Rate)

        For background info: https://life123.science/reactions
        What we're computing here, is referred to as:  (Δt)∗delta_forward(n)

        :param rxn_index:   An integer that indexes the reaction of interest
        :param conc_dict:   Concentrations of the applicable chemicals,
                            as a dict where the key value is the chemical's index
                            EXAMPLE: {3: 16.3, 8: 0.53, 12: 1.78}
        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                            concentration won't vary significantly during this span

        :return:            The differences between the reaction's forward and reverse "conversions",
                            a non-standard term we're using here to refer to delta_time * (Forward_Rate − Reverse_Rate),
                            for the given reaction during the specified time span
                            TODO: also, make a note of large relative increments (to guide future time-step choices)
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
            conc = conc_dict.get(species_index)
            assert conc is not None,\
                f"ReactionDynamics.compute_rate_delta(): lacking the concentration value for the species `{self.reaction_data.get_name(species_index)}`"
            forward_rate *= conc ** order      # Raise to power

        reverse_rate = rev_rate_constant
        for p in products:
            stoichiometry, species_index, order = p     # Unpack the data of the reaction product p
            conc = conc_dict.get(species_index)
            assert conc is not None, \
                f"ReactionDynamics.compute_rate_delta(): lacking the concentration value for the species `{self.reaction_data.get_name(species_index)}`"
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
            print(f"Ratio of reactant/product concentrations, adjusted for reaction orders: {conc_ratio}")
            print(f"    {numerator} / {denominator}")

        return np.allclose(conc_ratio, rate_ratio, atol=tolerance)
