from typing import Union, List
import math
import numpy as np


class ReactionDynamics:
    """
    TODO: In progress
    NOTE: This is the merger of the 2 former classes Reaction and Chemicals

    Carry out the reactions
    
    TODO (in-progress): add general methods to carry out reactions (currently in BioSim1D)

    """

    def __init__(self, reaction_data):
        """
        
        :param reaction_data:   Object of type "ReactionData" (with data about chemicals and reactions)
        """
        self.reaction_data = reaction_data
        self.debug = False




    #########################################################################
    #                                                                       #
    #                   TO CARRY OUT AND EXAMINE REACTIONS                  #
    #                                                                       #
    #########################################################################

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
            "specify_steps(): cannot specify all 3 arguments: `total_duration`, `time_step`, `n_steps` (should specify any 2 of them)"

        assert (total_duration and time_step) or (total_duration and n_steps) or (time_step and n_steps), \
            "specify_steps(): must provide exactly 2 arguments from:  `total_duration`, `time_step`, `n_steps`"

        if not time_step:
            time_step = total_duration / n_steps

        if not n_steps:
            n_steps = math.ceil(total_duration / time_step)

        return (time_step, n_steps)     # Note: could opt to also return total_duration is there's a need for it



    def single_compartment_reaction_step(self, conc_dict: dict, delta_time: float) -> np.array:
        """
        Using the given concentration data for all the applicable species in a single compartment,
        do a single reaction time step for ALL the reactions -
        based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        Return the increment vector for all the chemical species concentrations in the compartment

        NOTES:  - the actual system concentrations are NOT changed
                - "compartments" may or may not correspond to the "bins" of the higher layers;
                  the calling code might have opted to merge some bins into a single compartment

        :param conc_dict:           EXAMPLE: {3: 16.3, 8: 0.53, 12: 1.78}
        :param delta_time:
        :return:                    The increment vector for all the chemical species concentrations
                                    in the compartment
                                    EXAMPLE (for 2 species with a 3:1 stoichiometry):   [7. , -21.]
        """
        # Compute the forward and back conversions of all the reactions
        delta_list = self.compute_all_rate_deltas(conc_dict=conc_dict, delta_time=delta_time)
        if self.debug:
            print(f"    delta_list: {delta_list}")


        #increment_vector = np.zeros((self.reaction_data.number_of_chemicals(), 1), dtype=float)       # One element per chemical species
        increment_vector = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)       # One element per chemical species

        # For each reaction, adjust the concentrations of the reactants and products,
        # based on the forward and back rates of the reaction
        for i in range(self.reaction_data.number_of_reactions()):
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

                if (conc_dict[species_index] + increment_vector[species_index]) < 0:
                    raise Exception(f"The given time interval ({delta_time}) "
                                    f"leads to negative concentrations in reactions: make it smaller!")


            #   The products increase based on the forward reaction,
            #             and decrease based on the reverse reaction
            for p in products:
                stoichiometry, species_index, order = p
                increment_vector[species_index] += stoichiometry * delta_list[i]

                if (conc_dict[species_index] + increment_vector[species_index]) < 0:
                    raise Exception(f"The given time interval ({delta_time}) "
                                    f"leads to negative concentrations in reactions: make it smaller!")

        # END for

        return increment_vector



    def compute_all_rate_deltas(self, conc_dict: dict, delta_time: float) -> list:
        """
        For all the reactions.  Return a list with an entry for each reaction

        :param conc_dict:    EXAMPLE: {3: 16.3, 8: 0.53, 12: 1.78}
        :param delta_time:

        :return:              A list of the differences between forward and reverse conversions;
                                    each list has 1 entry per reaction, in the index order of the reactions
        """
        delta_list = []            # It will have 1 entry per reaction
        for i in range(self.reaction_data.number_of_reactions()):
        # Consider each reaction in turn
            if self.debug:
                print(f"    evaluating the rates for reaction number {i}")

            delta = self.compute_rate_delta(rxn_index=i, conc_dict=conc_dict, delta_time=delta_time)
            delta_list.append(delta)

        return( delta_list )




    def compute_rate_delta(self, rxn_index: int, conc_dict: dict, delta_time: float) -> float:
        """
        For the specified reaction, compute the difference of its forward and back "contributions",
        i.e. delta_time * reaction_rate * (concentration ** reaction_order)

        :param rxn_index:
        :param conc_dict:    EXAMPLE: {3: 16.3, 8: 0.53, 12: 1.78}
        :param delta_time:

        :return:            The differences between forward and reverse conversions
                            TODO: also, make a note of large relative increments (to guide future time-step choices)
        """

        # TODO: turn into a more efficient single step, as as:
        #(reactants, products, fwd_rate_coeff, back_rate_coeff) = cls.all_reactions.unpack_reaction(i)
        reactants = self.reaction_data.get_reactants(rxn_index)
        products = self.reaction_data.get_products(rxn_index)
        fwd_rate_coeff = self.reaction_data.get_forward_rate(rxn_index)
        back_rate_coeff = self.reaction_data.get_back_rate(rxn_index)

        delta_fwd = delta_time * fwd_rate_coeff         # TODO: save, to avoid re-computing at each bin
        for r in reactants:
            stoichiometry, species_index, order = r
            conc = conc_dict.get(species_index)
            assert conc is not None,\
                f"compute_rate_diff(): lacking the concentration value for the species `{self.reaction_data.get_name(species_index)}`"
            delta_fwd *= conc ** order      # Raise to power

        delta_back = delta_time * back_rate_coeff       # TODO: save, to avoid re-computing at each bin
        for p in products:
            stoichiometry, species_index, order = p
            conc = conc_dict.get(species_index)
            assert conc is not None, \
                f"compute_rate_diff(): lacking the concentration value for the species `{self.reaction_data.get_name(species_index)}`"
            delta_back *= conc ** order     # Raise to power

        #print(f"    delta_fwd: {delta_fwd} | delta_back: {delta_back}")
        return delta_fwd - delta_back



    def is_in_equilibrium(self, rxn_index: int, conc: dict, tolerance=0.05, explain=True) -> bool:
        """
        Ascertain whether the given concentrations are in equilibrium for the specified reaction

        :param rxn_index:   The index (0-based) to identify the reaction of interest
        :param conc:        Dict with the concentrations of the species involved in the reaction
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
