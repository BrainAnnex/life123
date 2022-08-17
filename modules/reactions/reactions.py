from typing import Union, List
import numpy as np


class Reactions:
    """
    Data about all applicable reactions,
    including stoichiometry, reaction rates and reaction orders

    DATA STRUCTURE:
        List of reactions.
        Each reaction is a Python dictionary with 4 keys:
            "reactants"
            "products"
            "Rf"    (forward reaction rate)
            "Rb"    (back reaction rate)

        Both reactants and products are triplets consisting of (stoichiometry, species index, reaction order).
        The "reaction order" refers to the forward reaction for reactants, and the reverse reaction for products.
    """

    def __init__(self, chem_data):
        self.reaction_list = []     # List of dicts.  Each item represents a reaction, incl. its reverse
                                    # Reactions should be added by means of calls to add_reaction()

        self.chem_data = chem_data  # Object with info on the individual chemicals, incl. their names
        


    def number_of_reactions(self) -> int:
        # Return the number of registered reactions
        return len(self.reaction_list)



    def get_reaction(self, i: int) -> dict:
        """
        Return the data structure of the i-th reaction (numbering starts at 0)

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A dictionary with 4 keys ("reactants", "products", "Rf", "Rb"),
                    where "Rf" is the forward reaction rate, and "Rb" the back reaction rate
        """
        assert 0 <= i < self.number_of_reactions(), \
            f"get_reaction(): argument `i` must be in the range [0-{self.number_of_reactions()}], inclusive; " \
            f"The value passed was {i} "

        return self.reaction_list[i]



    def get_reactants(self, i: int) -> (int, int, int):
        """
        Return a triplet with details of the reactants of the i-th reaction

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A triplet (stoichiometry, species index, reaction order)
        """
        rxn = self.get_reaction(i)
        return rxn["reactants"]


    def get_reactants_formula(self, i):
        """

        :param i:   The index (0-based) to identify the reaction of interest
        :return:
        """
        rxn = self.get_reaction(i)
        reactants = rxn["reactants"]
        return self.standard_form_chem_eqn(reactants)



    def get_products(self, i: int) -> (int, int, int):
        """
        Return a triplet with details of the products of the i-th reaction

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A triplet (stoichiometry, species index, reaction order)
        """
        rxn = self.get_reaction(i)
        return rxn["products"]


    def get_products_formula(self, i):
        """

        :param i:   The index (0-based) to identify the reaction of interest
        :return:
        """
        rxn = self.get_reaction(i)
        products = rxn["products"]
        return self.standard_form_chem_eqn(products)


    def get_forward_rate(self, i: int) -> float:
        rxn = self.get_reaction(i)
        return rxn["Rf"]

    def get_reverse_rate(self, i: int) -> float:    # TODO: PHASE OUT.  Transition to name get_back_rate()
        rxn = self.get_reaction(i)
        return rxn["Rb"]

    def get_back_rate(self, i: int) -> float:   # Transitioning to this name for consistency
        rxn = self.get_reaction(i)
        return rxn["Rb"]


    def extract_reactants(self, rxn: dict) -> [(int, int, int)]:
        """
        Return a list of triplets with details of the reactants of the given reaction

        :param rxn: The data structure representing the reaction
        :return:    A list of triplets of the form (stoichiometry, species index, reaction order)
        """
        return rxn["reactants"]


    def extract_stoichiometry(self, term):
        """

        :param term:
        :return:
        """
        return term[0]

    def extract_species_index(self, term):
        return term[1]

    def extract_rxn_order(self, term):
        return term[2]



    def extract_products(self, rxn: dict) -> [(int, int, int)]:
        """
        Return a list of triplet with details of the products of the given reaction

        :param rxn: The data structure representing the reaction
        :return:    A list of triplets of the form (stoichiometry, species index, reaction order)
        """
        return rxn["products"]


    def extract_forward_rate(self, rxn: dict) -> int:
        return rxn["Rf"]

    def extract_back_rate(self, rxn: dict) -> int:
        return rxn["Rb"]



    def add_reaction(self, reactants: list, products: list, forward_rate: float, reverse_rate: float) -> None:
        """
        Add the parameters of a SINGLE reaction, including its reverse rate

        NOTE: in the next 2 parameters, if the stoichiometry and/or reaction order aren't specified, they're assumed to be 1
        :param reactants:       A list of triplets (stoichiometry, species name or index, reaction order)
        :param products:        A list of triplets (stoichiometry, species name or index, reaction order of REVERSE reaction)

        :param forward_rate:
        :param reverse_rate:
        :return:                None
        """
        assert type(reactants) == list, "add_reaction(): argument `reactants` must be a list"
        reactant_list = [self._parse_reaction_term(r, "reactant") for r in reactants]

        assert type(products) == list, "add_reaction(): argument `products` must be a list"
        product_list = [self._parse_reaction_term(r, "product") for r in products]


        rxn = {"reactants": reactant_list, "products": product_list, "Rf": forward_rate, "Rb": reverse_rate}
        self.reaction_list.append(rxn)



    def clear_reactions(self) -> None:
        """
        Get rid of all reactions; start again with "an empty slate" (but still with reference
        to the same "Chemicals" object)
        :return:    None
        """
        self.reaction_list = []



    def describe_reactions(self, concise=False, return_as_list=False) -> Union[None, List[str]]:
        """
        Print out and return a listing with a user-friendly plain-text form of all the reactions
        EXAMPLE:  ["(0) CH4 + 2 O2 <-> CO2 + 2 H2O  (Rf = 3.0 , Rb = 2.0)"]

        :param concise:         If True, less detail is shown
        :param return_as_list:  If True, in addition to being printed, a list of reaction listings gets returned
        :return:                Either None or a list of strings, based on the flag return_as_list
        """
        print("Number of reactions: ", self.number_of_reactions())

        out = []    # Output list being built (and printed out item-wise)
        for i, rxn in enumerate(self.reaction_list):
            reactants = self.extract_reactants(rxn)
            products = self.extract_products(rxn)

            left = self.standard_form_chem_eqn(reactants)     # Left side of the equation
            right = self.standard_form_chem_eqn(products)     # Right side of the equation

            rxn_description = f"{i}: {left} <-> {right}"    # Initial brief description of the reaction
            if not concise:     # Add more detail
                rxn_description += f"  (Rf = {rxn['Rf']} / Rb = {rxn['Rb']})"
                for r in reactants:
                    if r[2] > 1:
                        rxn_description += f" | {r[2]}-th order in reactant {self.chem_data.get_name(r[1])}"
                for p in products:
                    if p[2] > 1:
                        rxn_description += f" | {p[2]}-th order in product {self.chem_data.get_name(p[1])}"

            print(rxn_description)
            out.append(rxn_description)

        if return_as_list:
            return out



    def standard_form_chem_eqn(self, eqn_side: list) -> str:
        """
        Return a user-friendly form of the given side of a chemical equation.

        EXAMPLE:  turn [(1, 0, 1), (2, 10, 1)]  into  "Fe + 2 Cl"  (if species 0 is named "Fe" and species 10 is "Cl")

        :param eqn_side:    A list encoding either side of a chemical equation
        :return:            A string with a user-friendly form of a side of a chemical equation
        """
        formula_list = []
        for t in eqn_side:
            stoichiometry = t[0]
            species_index = t[1]
            # Note: the reaction order (stored in t[2]) is not used

            if stoichiometry == 1:
                term = f"{self.chem_data.get_name(species_index)}"
            else:
                term = f"{stoichiometry} {self.chem_data.get_name(species_index)}"

            formula_list.append(term)

        return " + ".join(formula_list)




    def is_in_equilibrium(self, rxn_index: int, conc: dict, tolerance=0.05, explain=True) -> bool:
        """
        Ascertain whether the given concentrations are in equilibrium for the given reaction

        :param rxn_index:   The index (0-based) to identify the reaction of interest
        :param conc:        EXAMPLE: {'A': 23.9931640625, 'B': 36.0068359375}
        :param tolerance:   Allowable absolute tolerance, to establish satisfactory equality
        :param explain:     If True, print out the formula being used.
                                EXAMPLES:  "([C][D]) / ([A][B])" , "[B] / [A]^2",
        :return:            True if the reaction is close enough to an equilibrium
        """
        rxn = self.get_reaction(rxn_index)
        reactants = self.extract_reactants(rxn) # A list of triplets
        products = self.extract_products(rxn)   # A list of triplets
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
            #([B]^2 * [X]) / [A]^^3

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
            print(f"{i}: {self.get_reactants(i)} <-> {self.get_products(i)}   ; Fwd: {self.get_forward_rate(i)} / Back: {self.get_reverse_rate(i)}")



    def _parse_reaction_term(self, term: Union[int, str, tuple, list], name="term") -> (int, int, int):
        """
        Accept various ways to specify a reaction term, and return a standardized tuple form.

        EXAMPLES (*assuming* that the species with index 5 is called "F"):
            5               (1, 5, 1)
            "F"             (1, 5, 1)
            (3, 5)          (3, 5, 1)
            (3, "F")        (3, 5, 1)
            (3, 5, 2)       (3, 5, 2)
            (3, "F", 2)     (3, 5, 2)
            Same if lists were used in lieu of tuples

        :param term:    An integer or string or tuple or list
        :param name:    An optional nickname to refer to this term in error messages, if applicable
        :return:        A standardized tuple form, of the form (stoichiometry, species, reaction_order),
                            where all terms are integers
        """

        if type(term) == int:
            return  (1, term, 1)
        elif type(term) == str:
            return  (1, self.chem_data.get_index(term), 1)
        elif type(term) != tuple and type(term) != list:
            raise Exception(f"_parse_reaction_term(): {name} must be either an integer string, or a pair or a triplet. Instead, it is {type(term)}")

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
