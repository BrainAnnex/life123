from typing import Union, List
import math
import numpy as np


class ReactionDynamics:
    """
    TODO: In progress
    NOTE: This is the merger of the 2 former classes Reaction and Chemicals

    Carry out the reactions
    """

    def __init__(self):
        pass



    #########################################################################
    #                                                                       #
    #                   TO CARRY OUT AND EXAMINE REACTIONS                  #
    #                                                                       #
    #########################################################################

    def specify_steps(self, total_duration=None, time_step=None, n_steps=None) -> (float, float):
        """
        Specify 2 out of 3 possible parameters, and determine the 3rd one:
            total_duration, time_step, n_steps
            Their desired relationship is total_duration = time_step * n_steps

        :param total_duration:  Float with the overall time advance (i.e. time_step * n_steps)
        :param time_step:       Float with the size of each time step
        :param n_steps:         The desired number of steps
        :return:                The pair (time_step, n_steps)
        """
        assert (not total_duration or not time_step or not n_steps), \
            "specify_steps(): cannot specify all 3 arguments: `total_duration`, `time_step`, `n_steps`"

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


        #increment_vector = np.zeros((self.chem_data.number_of_chemicals(), 1), dtype=float)       # One element per chemical species
        increment_vector = np.zeros(self.chem_data.number_of_chemicals(), dtype=float)       # One element per chemical species

        # For each reaction, adjust the concentrations of the reactants and products,
        # based on the forward and back rates of the reaction
        for i in range(self.number_of_reactions()):
            if self.debug:
                print(f"    adjusting the species concentrations based on reaction number {i}")

            # TODO: turn into a more efficient single step, as as:
            #(reactants, products) = cls.all_reactions.unpack_terms(i)
            reactants = self.get_reactants(i)
            products = self.get_products(i)

            # Determine the concentration adjustments

            #   The reactants decrease based on the forward reaction,
            #             and increase based on the reverse reaction
            for r in reactants:
                stoichiometry, species_index, order = r
                increment_vector[species_index] += stoichiometry * (- delta_list[i])

                if (conc_dict[species_index] + increment_vector[species_index]) < 0:
                    raise Exception(f"The given time interval ({delta_time}) leads to negative concentrations in reactions: make it smaller!")


            #   The products increase based on the forward reaction,
            #             and decrease based on the reverse reaction
            for p in products:
                stoichiometry, species_index, order = p
                increment_vector[species_index] += stoichiometry * delta_list[i]

                if (conc_dict[species_index] + increment_vector[species_index]) < 0:
                    raise Exception(f"The given time interval ({delta_time}) leads to negative concentrations in reactions: make it smaller!")

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
        for i in range(self.number_of_reactions()):
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
        reactants = self.get_reactants(rxn_index)
        products = self.get_products(rxn_index)
        fwd_rate_coeff = self.get_forward_rate(rxn_index)
        back_rate_coeff = self.get_back_rate(rxn_index)

        delta_fwd = delta_time * fwd_rate_coeff         # TODO: save, to avoid re-computing at each bin
        for r in reactants:
            stoichiometry, species_index, order = r
            conc = conc_dict.get(species_index)
            assert conc is not None,\
                f"compute_rate_diff(): lacking the concentration value for the species `{self.chem_data.get_name(species_index)}`"
            delta_fwd *= conc ** order      # Raise to power

        delta_back = delta_time * back_rate_coeff       # TODO: save, to avoid re-computing at each bin
        for p in products:
            stoichiometry, species_index, order = p
            conc = conc_dict.get(species_index)
            assert conc is not None, \
                f"compute_rate_diff(): lacking the concentration value for the species `{self.chem_data.get_name(species_index)}`"
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
        rxn = self.get_reaction(rxn_index)

        reactants = self.extract_reactants(rxn)     # A list of triplets
        products = self.extract_products(rxn)       # A list of triplets
        Rf = self.extract_forward_rate(rxn)
        Rb = self.extract_back_rate(rxn)

        rate_ratio = Rf / Rb                    # Ratio of forward/reverse reaction rates

        conc_ratio = 1.
        numerator = ""
        denominator = ""
        for rxn_index, p in enumerate(products):
            species_index = self.extract_species_index(p)
            rxn_order = self.extract_rxn_order(p)

            species_name = self.chem_data.get_name(species_index)
            species_conc = conc[species_name]
            conc_ratio *= (species_conc ** rxn_order)
            if explain:
                numerator += f"[{species_name}]"
                if rxn_order > 1:
                    numerator += f"^{rxn_order} "

        if explain and len(products) > 1:
            numerator = f"({numerator})"

        for r in reactants:
            species_index = self.extract_species_index(r)
            rxn_order = self.extract_rxn_order(r)

            species_name = self.chem_data.get_name(species_index)
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




    #########################################################################
    #                                                                       #
    #              SUPPORT FOR CREATION OF NETWORK DIAGRAMS                 #
    #                                                                       #
    #########################################################################

    def prepare_graph_network(self) -> dict:
        """

        :return:
        """
        return {
            # Data to define the nodes and edges of the network
            'graph': self.create_graph_network_data(),

            # Mapping the node label to its interior color
            'color_mapping': self.assign_color_mapping()
        }



    def create_graph_network_data(self) -> [{}]:
        """
        Encode the reaction data in a form suitable for visualization
        with the graph module "vue_cytoscape"

        TODO:   assign a new, separate node label to chemicals that are both reagents and product
                Ditch all None values

        :return:    A list of dictionaries.  Each dictionary must have an 'id' key with a unique value
        """
        graph_data = []
        species_in_graph = []

        # Note: the species index of the various chemicals is a UNIQUE number; so, it's suitable to be used as an ID for the nodes
        #       For the reaction nodes, use numbers from a range starting just above the end-range of the numbers for the chemicals
        next_available_id = self.chem_data.number_of_chemicals()

        for i, rxn in enumerate(self.reaction_list):
            print("\n", rxn, "\n")
            # Add a node representing the reaction
            rxn_id = next_available_id
            next_available_id += 1
            graph_data.append({'id': rxn_id, 'label': 'Reaction', 'name': 'RXN',
                               'Rf': self.extract_forward_rate(rxn), 'Rb': self.extract_back_rate(rxn)})

            # Process all products
            products = self.extract_products(rxn)
            for term in products:
                species_index = term[1]
                # Add each product to the graph as a node (if not already present)
                if species_index not in species_in_graph:
                    graph_data.append({'id': species_index, 'label': 'Product',
                                       'name': self.chem_data.get_name(species_index),
                                       'diff_rate': self.chem_data.get_diffusion_rate(species_index),
                                       'stoich': self.extract_stoichiometry(term),
                                       'rxn_order': self.extract_rxn_order(term)
                                       })
                # Append edge from "reaction node" to product
                graph_data.append({'id': next_available_id, 'source': rxn_id, 'target': species_index, 'name': 'produces'})
                next_available_id += 1

            # Process all reactants
            reactants = self.extract_reactants(rxn)
            for term in reactants:
                species_index = term[1]
                # Add each reactant to the graph as a node (if not already present)
                if species_index not in species_in_graph:
                    graph_data.append({'id': species_index, 'label': 'Reactant',
                                       'name': self.chem_data.get_name(species_index),
                                       'diff_rate': self.chem_data.get_diffusion_rate(species_index),
                                       'stoich': self.extract_stoichiometry(term),
                                       'rxn_order': self.extract_rxn_order(term)
                                       })
                # Append edge from reactant to "reaction node"
                graph_data.append({'id': next_available_id, 'source': species_index, 'target': rxn_id, 'name': 'reacts'})
                next_available_id += 1

        return graph_data



    def assign_color_mapping(self):
        return {
            'Reactant': 'neo4j_green',
            'Product': 'neo4j_red',
            'Reaction': 'neo4j_lightbrown'
        }



    #########################################################################
    #                                                                       #
    #                               PRIVATE                                 #
    #                                                                       #
    #########################################################################

    def _internal_reactions_data(self) -> None:
        """
        Print out the low-level view of the reactions data
        :return:    None
        """
        for i in range(self.number_of_reactions()):
            print(f"{i}: {self.get_reactants(i)} <-> {self.get_products(i)}   ; Fwd: {self.get_forward_rate(i)} / Back: {self.get_back_rate(i)}")



    def _parse_reaction_term(self, term: Union[int, str, tuple, list], name="term") -> (int, int, int):
        """
        Accept various ways to specify a reaction term, and return a standardized tuple form of it.
        In the tuples or lists, the 1st entry is the stoichiometry, the 2nd one is the chemical name or index,
        and the option 3rd one is the reaction order

        EXAMPLES (*assuming* that the chemical species with index 5 is called "F"):
            5       gets turned into:   (1, 5, 1)
            "F"                         (1, 5, 1)
            (3, 5)                      (3, 5, 1)
            (3, "F")                    (3, 5, 1)
            (3, 5, 2)                   (3, 5, 2)
            (3, "F", 2)                 (3, 5, 2)
            Same if lists were used in lieu of tuples

        :param term:    An integer or string or tuple or list
        :param name:    An optional nickname to refer to this term in error messages, if applicable
                            (for example, "reactant" or "product")
        :return:        A standardized tuple form, of the form (stoichiometry, species, reaction_order),
                            where all terms are integers
        """

        if type(term) == int:
            return  (1, term, 1)
        elif type(term) == str:
            return  (1, self.chem_data.get_index(term), 1)
        elif type(term) != tuple and type(term) != list:
            raise Exception(f"_parse_reaction_term(): {name} must be either an integer, or a string, or a pair or a triplet. Instead, it is {type(term)}")

        # If we get thus far, term is either a tuple or a list
        if len(term) != 3 and len(term) != 2:
            raise Exception(f"parse_reaction_term(): Unexpected length for {name} tuple/list: it should be 2 or 3. Instead, it is {len(term)}")

        stoichiometry = term[0]
        species = term[1]
        if type(species) == str:
            species = self.chem_data.get_index(species)
        elif type(species) != int:
            raise Exception(f"parse_reaction_term(): The species value must be an integer or a string. Instead, it is {species}")

        if len(term) == 2:
            return (stoichiometry, species, 1)
        else:   # Length is 3
            reaction_order = term[2]
            return (stoichiometry, species, reaction_order)
