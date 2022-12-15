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



    def single_compartment_reaction_step(self, delta_time: float, conc_dict=None, conc_array=None, dynamic_step=1) -> np.array:
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
        delta_list = self.compute_all_rate_deltas(conc_dict=conc_dict, delta_time=delta_time)
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



    def compute_all_rate_deltas(self, conc_dict: dict, delta_time: float) -> list:
        """
        For an explanation of the "rate delta", see compute_rate_delta().
        Compute the "rate delta" for all the reactions.  Return a list with an entry for each reaction

        :param conc_dict:   Concentrations of the applicable chemicals,
                            as a dict where the key value is the chemicals index
                            EXAMPLE: {3: 16.3, 8: 0.53, 12: 1.78}
        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                            concentration won't vary significantly during this span

        :return:              A list of the differences between forward and reverse "conversions";
                                    each list has 1 entry per reaction, in the index order of the reactions
        """
        delta_list = []            # It will have 1 entry per reaction
        for i in range(self.reaction_data.number_of_reactions()):   # Consider each reaction in turn
            if self.debug:
                print(f"    evaluating the rates for reaction number {i}")

            delta = self.compute_rate_delta(rxn_index=i, conc_dict=conc_dict, delta_time=delta_time)
            delta_list.append(delta)

        return( delta_list )



    def compute_rate_delta(self, rxn_index: int, conc_dict: dict, delta_time: float) -> float:
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

        :return:            The differences between forward and reverse "conversions",
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
