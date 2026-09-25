import math

from life123.visualization.py_graph_visual import PyGraphVisual
from life123.visualization.graphic_log import DisplayNetwork
from life123.reactions import ReactionDefinition, SimulationReaction
from life123.species_registry import SpeciesRegistry



class ReactionRegistry:
    """
    Manage a list of reactions, and the reaction-specific objects (defined in reactions.py file),
    such as ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition, ReactionGeneric, ReactionEnzyme, etc.

    Instances of this class are typically used by UniformCompartment objects,
    or any other high-level module that needs access to the reactions.

    A ReactionRegistry object may be shared by multiple UniformCompartment objects IF the latter
    all make use of ALL the registered reactions  (i.e. no "pick and choose" some of the reactions.)
    """

    def __init__(self, species_data=None):

        """
        :param species_data:    [OPTIONAL] Object of type "SpeciesRegistry";
                                    if not passed, it will get instantiated automatically -
                                    and then it can be obtained by means of the method get_species_data()
        """

        # TODO: consider adding to argument   "=species=None"
        """
        assert (chem_data is not None) or (labels is not None), \
            "ReactionRegistry() instantiation: exactly one of the arguments `chem_data` or `labels` must be provided"

        assert (chem_data is None) or (labels is None), \
            "ReactionRegistry() instantiation: cannot specify both the arguments `chem_data` and `labels`"

        if labels is not None:
            chem_data = SpeciesRegistry(labels=labels)
        """

        #assert chem_data is not None, \
            #"ReactionRegistry() instantiation: the arguments `chem_data` must be provided, and cannot be None"

        if species_data is None:
            self.species_data = SpeciesRegistry()
        else:
            self.species_data = species_data


        self.reaction_defn_list = []    # List of "ReactionDefinition" objects
        self.reaction_list = []         # List of "SimulationReaction" objects


        self.active_chemicals = set()   # Set of the id's of species - not counting pure catalysts - involved
                                        # in any of the registered reactions
                                        # CAUTION: the concept of "active chemical" might change in future versions, where only SOME of
                                        #          the reactions are simulated.  TODO: it might better belong to UniformCompartment





    #####################################################################################################

    '''                                          ~   TO VIEW  ~                                       '''

    def ________VIEW________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def number_of_reactions(self, include_inactive=False) -> int:
        """
        Return the number of registered chemical reactions
        (the number of DERIVED "SimulationReaction" objects)

        :param include_inactive:    [NO LONGER IN USE] If True, disabled reactions are also included
        :return:                    The number of registered chemical reactions
        """
        return len(self.reaction_list)



    def get_all_reactions(self):
        """
        Return the list of all the "simulation reactions" that have been created.
        Note: this number might be higher than the number of reactions specified by the user,
              because some overall reactions are modeled by means of multiple "simulation reactions"

        :return:    A list of "SimulationReaction" objects
        """
        return self.reaction_list



    def active_reaction_indices(self) -> [int]:
        """
        TODO: OBSOLETE

        Return a list of the reaction index numbers of all the reactions

        DEPRECATION: no active/inactive distinction is kept anymore.  All reactions are regarded as "active"

        :return:    A list of integers, to identify the active reactions by their indices
        """
        l = []
        for i, rxn in enumerate(self.reaction_list):
            #if rxn.active:
            l.append(i)

        return l



    def assert_valid_rxn_index(self, index :int) -> None:
        """
        Raise an Exception if the specified reaction index isn't valid

        :param index:   An integer that indexes the reaction of interest (numbering starts at 0)
        :return:        None
        """
        assert self.number_of_reactions() > 0, \
            f"assert_valid_rxn_index(): there are no reactions defined yet.  Use add_reaction() to add them first"

        assert (type(index) == int), \
            f"assert_valid_rxn_index(): the requested reaction index must be an integer; " \
            f"the provided value ({index}) is of type {type(index)}"

        assert 0 <= index < self.number_of_reactions(), \
            f"assert_valid_rxn_index(): the requested reaction index is not the expected range [0 to {self.number_of_reactions() - 1}], inclusive; " \
            f"the value passed was: {index} (there is no reaction whose index is {index})"



    def get_species_data(self) -> SpeciesRegistry:
        """
        Return the "SpeciesRegistry" object being used

        :return:    Object of type "SpeciesRegistry"
        """
        return self.species_data



    def get_reaction(self, i :int) -> SimulationReaction:
        """
        Return the data structure of the i-th reaction,
        in the order in which reactions were added (numbering starts at 0)

        :param i:   An integer that indexes the reaction of interest (numbering starts at 0)
        :return:    A "SimulationReaction object
        """
        self.assert_valid_rxn_index(i)

        return self.reaction_list[i]



    def get_reactants(self, i :int) -> [(int, str)]:
        """
        Return a list of pairs with details of the reactants of the i-th reaction.
        Each pair represents a "complex" of the reaction reactants

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A list of pairs of the form (stoichiometry, chem labels)
        """
        rxn = self.get_reaction(i)
        return rxn.stoichiometry.get_reactant_list()


    def get_reactants_formula(self, i :int) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A string with a user-friendly form of the left (reactants) side of the chemical reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_reactants_formula()



    def get_products(self, i :int) -> [(int, str)]:
        """
        Return a list of pairs with details of the products of the i-th reaction.
        Each pair represents a "complex" of the reaction products

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A list of pairs of the form (stoichiometry, species id)
        """
        rxn = self.get_reaction(i)
        return rxn.stoichiometry.get_product_list()


    def get_products_formula(self, i :int) -> str:
        """
        Return a string with a user-friendly form of the right (products) side of the reaction formula

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A string with a user-friendly form of the right (products) side of the chemical reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_products_formula()



    def get_forward_rate(self, i :int) -> float:
        """
        Return the value of the forward rate constant of the i-th reaction

        :param i:   The integer index (0-based) to identify the reaction of interest
        :return:    The value of the forward rate constant for the above reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_forward_rate_constant()


    def get_reverse_rate(self, i :int) -> float:
        """
        Return the value of the reverse (back) rate constant of the i-th reaction

        :param i:   The integer index (0-based) to identify the reaction of interest
        :return:    The value of the reverse (back) rate constant for the above reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_reverse_rate_constant()



    def get_species_in_reaction(self, rxn_index :int) -> set[str]:
        """
        Return the SET of the id's
        of all the species participating in the reaction with the specified index

        :param rxn_index:   An integer with the (zero-based) index to identify the reaction of interest
        :return:            A SET of indices of the id's of the species involved in the above reaction
                                Note: being a set, it's NOT in any particular order
        """
        return self.get_reaction(rxn_index).stoichiometry.get_all_species_ids()



    def get_reactions_participating_in(self, species_id :str, side :str) -> list[ReactionDefinition]:
        """
        Return a list of all the reactions that the given chemical species
        is involved in

        :param species_id:  To identify a particular species
        :param side:        Either "reagent" or "product"
        :return:            List of "ReactionDefinition" objects
        """
        assert side == "reagent" or side == "product", \
            "get_reactions_participating_in(): argument `side` must be either 'reagent' or 'product'"

        rxns_found_in = []
        for rxn_defn in self.reaction_defn_list:
            if side == "reagent":
                if species_id in rxn_defn.extract_reactant_ids():
                    rxns_found_in.append(rxn_defn)
            else:
                 if species_id in rxn_defn.extract_product_ids():
                    rxns_found_in.append(rxn_defn)


        return rxns_found_in



    def number_of_active_chemicals(self) -> int:
        """
        Return the number of all the chemicals
        involved in ANY of the registered reactions,
        but NOT counting chemicals that always appear
        in a catalytic role in all the reactions they participate in
        (if a chemical participates in a non-catalytic role in ANY reaction, it'll appear here)
        """
        return len(self.active_chemicals)



    def labels_of_active_chemicals(self, sort=False) -> list[str]:
        """
        Return a list of the id's of all the species
        involved in ANY of the registered reactions,
        but NOT counting species that always appear
        in a catalytic role in all the reactions they participate in
        (if a species participates in a non-catalytic role in ANY reaction, it'll appear here)

        The list is not in any particular order, unless sort is True.

        :param sort:    If True, the list is sorted
        :return:        A list of species id's
        """
        if not sort:
            return list(self.active_chemicals)

        return sorted(self.active_chemicals)
        #return sorted(self.active_chemicals, key=self.species_data.get_species_index)





    #####################################################################################################

    '''                             ~   TO MODIFY THE REACTIONS  ~                                    '''

    def ________MODIFY_REACTIONS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def add_elementary_reaction(self, reactants :str|list, products :str|list,
                                thermodynamic_parameters=None, temp=None,
                                kinetic_parameters=None) -> int:
        """
        Create and register a new SINGLE elementary chemical reaction,
        optionally including its kinetic and/or thermodynamic data.
        All the involved chemicals can be either previously registered, or not.

        :param reactants:       A string or pair of strings; for reactions such as 2 A -> P, pass ["A", "A"]
        :param products:        A string or pair of strings; for reactions such as A -> 2 P, pass ["P", "P"]

        :param temp:            [OPTIONAL] Temperature in Kelvins

        :return:                Integer index of the newly-added reaction
                                    (in the list self.reaction_list, stored as object variable)
        """
        # Turn `reactants` and `products` into lists, if not already lists
        if type(reactants) == str:
            reactants = [reactants]
            n_reactants = 1
        else:
            n_reactants = len(reactants)
            assert n_reactants <= 2

        if type(products) == str:
            products = [products]
            n_products = 1
        else:
            n_products = len(products)
            assert n_products <= 2

        # Determine the reaction type from the number of reactants and products
        # TODO: turn into a method
        if n_reactants == 1 and n_products == 1:
            reaction_type = "ReactionUnimolecular"
        elif n_reactants == 2 and n_products == 1:
            reaction_type = "ReactionSynthesis"
        elif n_reactants == 1 and n_products == 2:
            reaction_type = "ReactionDecomposition"
        else:
            raise Exception(f"add_elementary_reaction(): {n_reactants} reactants and {n_products} products cannot correspond to an elementary reaction")

        #print(f"add_elementary_reaction(): adding reaction of type `{reaction_type}`")

        return self.add_reaction(reactants=reactants, products=products,
                                 reaction_model="mass action",
                                 thermodynamic_parameters=thermodynamic_parameters,
                                 kinetic_parameters=kinetic_parameters,
                                 temp=temp)



    def add_reaction(self, reactants :str|list, products :str|list,
                     reaction_model=None,
                     autoregister_species=True,
                     thermodynamic_parameters=None,
                     kinetic_parameters=None,
                     temp=None) -> int:
        """
        Create and register a new SINGLE chemical reaction,
        optionally including its kinetic and/or thermodynamic data.
        All the involved chemicals can be either previously registered, or not.

        NOTE: in the reactants and products, if the stoichiometry coefficients aren't specified,
              they're assumed to be 1.

              The full structure of each term in the list of reactants and of products
              is the pair:  (stoichiometry coefficient, species id)

              EXAMPLES of formats to use for each term in the lists of the reactants and of the products:
                "F"         is taken to mean (1, "F") - default stoichiometry
                (2, "F")    is taken to mean (2, "F") - stoichiometry coefficient specified

        :param reactants:       A string or list of pairs (stoichiometry, species name),
                                    or simplified terms in various formats; for details, see above
        :param products:        A string or list of pairs (stoichiometry, species name),
                                    or simplified terms in various formats; for details, see above

        :param reaction_model:[OPTIONAL] Primarily meant for the kinetics.
                                Allowed values are "mass action", "michaelis menten",
                                "single substrate mechanism", "custom" - as detailed
                                in class ReactionModelRegistry.
                                If not provided, it default to "mass action" if the reaction
                                seems to be elementary from its stoichiometry and kinetic_parameters.

        :param temp:            [OPTIONAL] Temperature in Kelvins

        :return:                Integer index of the newly-added reaction
                                    (in the list self.reaction_list, stored as object variable)
        """
        # Merge the temp value into the rest of the thermodynamics parameters
        if temp is not None:
            if thermodynamic_parameters is None:
                thermodynamic_parameters = {"temp": temp}
            else:
                if not "temp" in thermodynamic_parameters:
                    thermodynamic_parameters["temp"] = temp
                else:
                    assert math.isclose(temp, thermodynamic_parameters["temp"]), \
                        "add_reaction(): inconsistent `temp` and `thermodynamic_parameters` arguments"

        rxn = ReactionDefinition(reactants=reactants, products=products,
                                 species_registry=self.species_data,
                                 autoregister_species=autoregister_species,
                                 reaction_model=reaction_model,
                                 thermodynamic_parameters=thermodynamic_parameters,
                                 kinetic_parameters=kinetic_parameters)

        #print(f"add_reaction(): detected reaction type `{reaction_type}`")

        return self.register_reaction(rxn_defn=rxn)



    def register_reaction(self, rxn_defn :ReactionDefinition) -> int:
        """
        Register a SINGLE chemical reaction from its "ReactionDefinition" object,
        and set all its kinetic and/or thermodynamic data from the available information,
        including the value of the temperature (stored in object variable.)

        Update the object variable `self.active_chemicals`

        :param rxn_defn:    Object of type "ReactionDefinition"

        :return:            Integer index of the newly-added reaction
                                (in the list self.reaction_list, stored as object variable)
        """
        reaction_id = len(self.reaction_list)

        rxn_defn.id = reaction_id

        self.reaction_defn_list.append(rxn_defn)

        sim_rxn_tuple = rxn_defn.sim_reactions
        for sim_rxn in sim_rxn_tuple:
            self.reaction_list.append(sim_rxn)
            involved_chemicals = sim_rxn.stoichiometry.get_all_species_ids()    # Set of species ID's
            # Update the set of "active chemicals"
            self.active_chemicals |= involved_chemicals     # Union of sets

        return reaction_id



    def clear_reactions_data(self) -> None:
        """
        Get rid of all the reactions; start again with "an empty slate" (but still with reference
        to the same data object about the chemicals and their properties)

        :return:    None
        """
        self.reaction_list = []
        self.reaction_defn_list = []
        self.active_chemicals = set()



    def _parse_reaction_term(self, term :str|tuple|list, name="term") -> (int, str):
        """
        Accept various ways to specify a reaction term, and return a standardized triplet form for it.

        NOTE: if the stoichiometry coefficient isn't specified, it defaults to 1

        In the passed tuples or lists:
            - required 1st entry is the stoichiometry
            - required 2nd entry is the chemical name

        If just a string is being passed, it is taken to be the chemical name,
        with stoichiometry coefficient of 1

        EXAMPLES:
            "F"          gets turned into:  (1, "F")   - defaults used for stoichiometry
            (2, "F")                        (2, "F")
            It's equally acceptable to use LISTS in lieu of tuples

        :param term:    A string (a chemical name)
                            OR  a pair (stoichiometry coeff, name)
        :param name:    An optional nickname, handy to refer to this term in error messages if needed
                            (for example, "reactant" or "product")
        :return:        A standardized pair of the form (stoichiometry, species_label),
                            where stoichiometry is an integer, while species_label is a string
        """
        if type(term) == str:
            return  (1, term)    # Accept simply the chemical name as a shortcut,
                                    # for when the stoichiometry coefficient and reaction order are both 1

        if type(term) != tuple and type(term) != list:
            raise Exception(f"_parse_reaction_term(): {name} must be either a string (a chemical name), "
                            f"or a pair (stoichiometry coeff, name). "
                            f"Instead, it is `{term}` (of type {type(term)})")

        # If we get thus far, term is either a tuple or a list
        assert len(term) == 2,  \
            f"_parse_reaction_term(): Unexpected length for {name} tuple/list: it should be 2. " \
            f"Instead, it is {len(term)}"

        stoichiometry = term[0]
        assert type(stoichiometry) == int, \
            f"_parse_reaction_term(): The stoichiometry coefficient must be an integer. Instead, it is {stoichiometry}"

        chem_label = term[1]
        assert type(chem_label) == str, \
                            f"_parse_reaction_term(): The chemical name must be a string. " \
                            f"Instead, it is `{chem_label}` (of type {type(chem_label)})"


        return (stoichiometry, chem_label)



    def _standardize_reaction_side(self, terms :str|list, arg_name :str):
        """

        :param terms:
        :param arg_name:
        :return:
        """
        assert terms is not None, \
            f"standardize_reaction_side(): the argument `{arg_name}` is a required one; it can't be None"

        if type(terms) == str:
            terms = [terms]
        else:
            assert type(terms) == list, \
                f"standardize_reaction_side(): the argument `{arg_name}` must be a string or a list; the passed value was {type(terms)}"

        return [self._parse_reaction_term(r, "reactant") for r in terms]   # A list of pairs




    #####################################################################################################

    '''                             ~   TO DESCRIBE THE REACTIONS  ~                                  '''

    def ________DESCRIBE_REACTIONS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def describe_reactions(self, concise=False) -> None:
        """
        Print out a user-friendly plain-text form of ALL the reactions.
        If wanting to describe just 1 reaction, use single_reaction_describe()

        EXAMPLE (not concise):
            Number of reactions: 2
            (0) CH4 + 2 O2 <-> CO2 + 2 H2O  (kF = 3.0 / kR = 2.0 / Delta_G = -1,005.13 / K = 1.5) | 1st order in all reactants & products
            (1) A + B <-> C  (kF = 5.0 / kR = 1.0 / Delta_G =  / K = 5.0) | 1st order in all reactants & products
            Chemicals involved in the above reactions: {'CH4', 'O2', 'H2O', 'A', 'B', 'C'}

        :param concise:     If True, less detail is shown
        :return:            None
        """
        #TODO: implement a Pandas version
        print(f"Number of reactions: {self.number_of_reactions()}")

        # Print a concise description of EACH REACTION IN TURN
        for description in self.multiple_reactions_describe(concise=concise):
            print(description)

        chem_labels = self.labels_of_active_chemicals(sort=True)   # Set of chem labels, sorted by species id

        # If plot colors were registered, show them alongside the chem labels
        chem_labels_with_colors = []
        for label in chem_labels:
            color = self.species_data.get_value(species_id=label, field="plot_color")
            if color:
                chem_labels_with_colors.append(f'"{label}" ({color})')
            else:
                chem_labels_with_colors.append(f'"{label}"')

        chem_labels = "{" + ", ".join(chem_labels_with_colors) + "}"


        #if self.active_enzymes == set():    # If no enzymes were involved in any reaction
        print(f"Chemicals involved in the above reactions: {chem_labels}")
        '''
        else:
            print(f"Chemicals involved in the above reactions (not counting enzymes): {chem_labels}")

            print(f"Enzymes involved in the above reactions: ")
            for enz in self.names_of_enzymes():
                pass
        '''



    def multiple_reactions_describe(self, rxn_list=None, concise=False) -> [str]:
        """
        The counterpart of single_reaction_describe() for many reactions.
        Return a list of strings, each string being a (concise or not) user-friendly plain-text form of
        each of the reactions

        :param rxn_list:    Either a list of integers, to identify the reactions of interest,
                                or None, meaning ALL reactions
        :param concise:     If True, less detail is shown
        :return:            A list of strings; each string is the description of one of the reactions
        """
        if rxn_list is None:
            rxn_list = range(self.number_of_reactions())    # Show ALL reactions, by default

        out = []    # Output list being built

        for i in rxn_list:
            description = self.single_reaction_describe(rxn_index=i, concise=concise)
            description = f"{i}: {description}"
            out.append(description)

        return out



    def single_reaction_describe(self, rxn_index: int, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the given reaction
        EXAMPLE (concise):      "CH4 + 2 O2 <-> CO2 + 2 H2O"
        EXAMPLE (not concise):  "CH4 + 2 O2 <-> CO2 + 2 H2O  (kF = 3.0 / kR = 2.0 / Delta_G = -1,005.13 / K = 1.5) | 1st order in all reactants & products"

        :param rxn_index:   Integer to identify the reaction of interest
        :param concise:     If True, less detail is shown
        :return:            A string with a description of the specified reaction
        """
        rxn_defn = self.reaction_defn_list[rxn_index]

        return rxn_defn.describe(concise)    # Invoke the individual "ReactionDefinition" object






    #####################################################################################################

    '''                          ~   FOR CREATION OF NETWORK DIAGRAMS  ~                              '''

    def ________NETWORK_DIAGRAMS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def prepare_graph_network(self) -> dict:
        """
        Prepare and return a data structure with chemical-reaction data in a network format,
        ready to be passed to the front end, for network-diagram visualization with the Cytoscape.js library
        (in the graph module "vue_cytoscape")

        A "bipartite" graph ("Petri net") representation is used for reaction, where the reaction itself - as well as the reactants
        and products - are all turned into graph vertices.

        4 parts are generated, and assembled together as a dictionary with 4 keys: 'nodes', 'edges', 'color_mapping', 'caption_mapping'

        EXAMPLE of the 'nodes' structure part of the returned object for an  A <-> B reaction:
           [{'id': 'C-0', 'labels': ['Chemical'], 'name': 'A', 'diff_rate': None},
            {'id': 'C-1', 'labels': ['Chemical'], 'name': 'B', 'diff_rate': None},

            {'id': 'R-0', 'labels': ['Reaction'], 'name': 'RXN', 'kF': 3.0, 'kR': 2.0, 'K': 1.5, 'Delta_G': -1005.13}
           ]

        EXAMPLE of the 'edges' structure part of the returned object for that same  A <-> B reaction:
           [
            {'id': 'edge-1', 'name': 'produces', 'source': 'R-0', 'target': 'C-1', 'stoich': 1, 'rxn_order': 1},
            {'id': 'edge-2', 'name': 'reacts',   'source': 'C-0', 'target': 'R-0', 'stoich': 1, 'rxn_order': 1}
           ]

        EXAMPLE of `color_mapping`:     {'Chemical': '#8DCC92', 'Reaction': '#D9C8AD'},
        EXAMPLE of `caption_mapping`:   {'Chemical': 'name', 'Reaction': 'name'}}

        :return:    A dictionary with 4 keys: 'nodes', 'edges', 'color_mapping', 'caption_mapping'
        """
        # TODO: manage shapes; e.g., use rectangles for reactions

        graph = PyGraphVisual()     # Object to facilitate data preparation for graph visualization

        # Note: the graph nodes representing Chemicals will be given an id such as "C-123" and a label "Chemical";
        #       the graph nodes representing Reactions will be given an id such as "R-456" and a label "Reaction"

        for i, rxn in enumerate(self.reaction_list):    # Consider each "SimulationReaction" object in turn
            # Add a node representing the reaction
            rxn_id = f"RXN-{i}"               # Example: "RXN-456"
            node_data = {'name': 'RXN', 'formula': rxn.standard_chemical_formula()}

            # Show the parameter of the reaction
            rxn_properties = rxn.source_object.extract_rxn_properties()
            #print(rxn_properties)

            for k,v in rxn_properties.items():
                if type(v) is float:
                    node_data[k] = f"{v:,.6g}"
                else:
                     node_data[k] = f"{v}"

            graph.add_node(node_id=rxn_id, labels='Reaction', properties=node_data)


            # Process all the PRODUCTS of this reaction
            products = rxn.stoichiometry.get_product_list()
            for stoich, species_id in products:
                #chemical_id = f"C-{self.species_data.get_species_index(species_name)}"      # Example: "C-12"
                # Use the species' "sort_order" attribute as a short unique ID
                chemical_id = f"C-{self.species_data.get_species(species_id).sort_order}"      # Example: "C-12"

                # Add each product to the graph as a node (if not already present)
                properties={'name': species_id}
                if diff := self.species_data.get_value(species_id=species_id, field="diffusion_rate"):
                    properties['diff_rate'] = diff
                graph.add_node( node_id=chemical_id, labels="Chemical",
                                properties=properties)

                # Append edge from "reaction node" to "product node"
                graph.add_edge(from_node=rxn_id, to_node=chemical_id, name="produces",
                               properties={'stoich': stoich})


            # Process all the REACTANTS of this reaction
            reactants = rxn.stoichiometry.get_reactant_list()
            for stoich, species_id in reactants:
                #chemical_id = f"C-{self.species_data.get_species_index(species_id)}"      # Example: "C-34"
                # Use the species' "sort_order" attribute as a short unique ID
                chemical_id = f"C-{self.species_data.get_species(species_id).sort_order}"   # Example: "C-34"

                # Add each reactant to the graph as a node (if not already present)
                properties={'name': species_id}
                if diff := self.species_data.get_value(species_id=species_id, field="diffusion_rate"):
                    properties['diff_rate'] = diff
                graph.add_node(node_id=chemical_id, labels="Chemical",
                               properties=properties)

                # Append edge from "reactant node" to "reaction node"
                graph.add_edge(from_node=chemical_id, to_node=rxn_id, name="reacts",
                               properties={'stoich': stoich})


        graph.assign_color_mapping(label='Chemical', color='graph_green')
        graph.assign_color_mapping(label='Reaction', color='graph_lightbrown')

        graph.assign_caption(label='Chemical', caption='name')
        graph.assign_caption(label='Reaction', caption='id')

        #print(graph)

        return graph.get_graph_data()



    def plot_reaction_network(self, log_file :str, graphic_component="vue_cytoscape_5") -> None:
        """
        Send a plot of the network of reactions to the HTML log file,
        also including a brief summary of all the reactions.

        :param log_file:            The name of the file into which to place the HTML code
                                        to create the interactive network plot.
                                        The suffix ".htm" will be added if it doesn't end with ".htm" or ".html"
                                        If the file already exists, it will get overwritten.
                                        (Note: this file will automatically include an internal reference to the JavaScript
                                        file specified in `graphic_component`)
        :param graphic_component:   The name of a Vue component that accepts a "graph_data" argument,
                                        an object with the following keys
                                        'nodes', 'edges', 'color_mapping' and 'caption_mapping
                                        For more details, see prepare_graph_network()

        :return:                    None
        """
        assert graphic_component == "vue_cytoscape_5", \
            "plot_reaction_network(): the only value supported for argument `graphic_component` is 'vue_cytoscape_5'"

        # Send a brief summary of all the reactions to the HTML log file
        header = "<h3>List of reactions:</h3>\n"

        header += "<pre>"    #"<p style='font-family: monospace; padding-left:20px'>"
        rxn_descriptions = self.multiple_reactions_describe(concise=True)
        for desc in rxn_descriptions:
            header += f"    {desc}\n"
        header += "</pre>\n"

        graph_data = self.prepare_graph_network()
        # A dictionary with 4 keys: ''nodes', 'edges', 'color_mapping' and 'caption_mapping'

        # Send a plot of the network of reactions to the HTML log file
        #GraphicLog.export_plot(graph_data, graphic_component, unpack=unpack)
        DisplayNetwork.export_plot(graph_data = graph_data,
                                   graphic_component = graphic_component,
                                   filename = log_file, caption=header)
