from life123.visualization.py_graph_visual import PyGraphVisual, PyGraphVisual_OLD
from life123.visualization.graphic_log import GraphicLog, DisplayNetwork
from life123.html_log import HtmlLog as log
from life123.reactions import ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition, ReactionEnzyme, ReactionGeneric
from life123.chem_data import ChemData



class ReactionRegistry:
    """
    Manage a list of reactions, and the reaction-specific objects (defined in reactions.py file),
    such as ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition, ReactionGeneric, ReactionEnzyme, etc.

    Instances of this class are typically used by UniformCompartment objects.
    A ReactionRegistry object may be shared by multiple UniformCompartment objects IF the latter
    all make use of ALL the registered reactions  (i.e. no "pick and choose" some of the reactions.)
    """

    def __init__(self, chem_data=None):

        """
        :param chem_data:   [OPTIONAL] Object of type "ChemData"; if not passed, it will get instantiated automatically -
                                and then it can be obtained by means of the method get_chem_data()
        """

        # TODO: consider adding to argument   "=labels=None"
        """
        assert (chem_data is not None) or (labels is not None), \
            "ReactionRegistry() instantiation: exactly one of the arguments `chem_data` or `labels` must be provided"

        assert (chem_data is None) or (labels is None), \
            "ReactionRegistry() instantiation: cannot specify both the arguments `chem_data` and `labels`"

        if labels is not None:
            chem_data = ChemData(labels=labels)
        """

        #assert chem_data is not None, \
            #"ReactionRegistry() instantiation: the arguments `chem_data` must be provided, and cannot be None"

        if chem_data is None:
            self.chem_data = ChemData()
        else:
            self.chem_data = chem_data

        self.reaction_list = []     # List of objects of the various individual reaction classes,
                                    # such as "ReactionGeneric" and "ReactionEnzyme"


        self.active_chemicals = set()   # Set of the labels of the chemicals - not counting pure catalysts - involved
                                        # in any of the registered reactions
                                        # CAUTION: the concept of "active chemical" might change in future versions, where only SOME of
                                        #          the reactions are simulated.  TODO: it might better belong to UniformCompartment

        #self.active_enzymes = set()     # Set of the labels of the enzymes (catalysts) involved
                                        # in any of the registered reactions
                                        # CAUTION: the concept of "active enzyme" might change in future versions, where only SOME of
                                        #          the reactions are simulated.  TODO: it might better belong to UniformCompartment




    def number_of_reactions(self, include_inactive=False) -> int:
        """
        Return the number of registered chemical reactions

        :param include_inactive:    If True, disabled reactions are also included
        :return:                    The number of registered chemical reactions
        """
        if include_inactive:
            return len(self.reaction_list)

        count = 0
        for rxn in self.reaction_list:
            if rxn.active:
                count += 1

        return count



    def get_all_reactions(self):
        """
        Return the list of all the reactions that have been registered

        :return:    A list of various types of reaction objects,
                        such as ReactionUnimolecular, ReactionSynthesis, etc.
        """
        return self.reaction_list



    def active_reaction_indices(self) -> [int]:
        """
        Return a list of the reaction index numbers of all the active reactions

        :return:    A list of integers, to identify the active reactions by their indices
        """
        l = []
        for i, rxn in enumerate(self.reaction_list):
            if rxn.active:
                l.append(i)

        return l



    def assert_valid_rxn_index(self, index :int) -> None:
        """
        Raise an Exception if the specified reaction index isn't valid

        :param index:   An integer that indexes the reaction of interest (numbering starts at 0)
        :return:        None
        """
        assert self.number_of_reactions() > 0, \
            f"ChemData.assert_valid_rxn_index(): there are no reactions defined yet.  Use add_reaction() to add them first"

        assert (type(index) == int), \
            f"ChemData.assert_valid_rxn_index(): the requested reaction index must be an integer; " \
            f"the provided value ({index}) is of type {type(index)}"

        assert 0 <= index < self.number_of_reactions(), \
            f"ChemData.assert_valid_rxn_index(): the requested reaction index is not the expected range [0 to {self.number_of_reactions() - 1}], inclusive; " \
            f"the value passed was: {index} (there is no reaction whose index is {index})"



    def get_chem_data(self):
        """
        :return:    Object of type "ChemData"
        """
        return self.chem_data



    def get_reaction(self, i :int):
        """
        Return the data structure of the i-th reaction,
        in the order in which reactions were added (numbering starts at 0)

        :param i:   An integer that indexes the reaction of interest (numbering starts at 0)
        :return:    An object of one of the various individual reaction classes,
                    such as "ReactionGeneric" and "ReactionEnz"
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
        return rxn.extract_reactants()


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
        :return:    A list of pairs of the form (stoichiometry, species label)
        """
        rxn = self.get_reaction(i)
        return rxn.extract_products()


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
        return rxn.extract_forward_rate()


    def get_reverse_rate(self, i :int) -> float:
        """
        Return the value of the reverse (back) rate constant of the i-th reaction

        :param i:   The integer index (0-based) to identify the reaction of interest
        :return:    The value of the reverse (back) rate constant for the above reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_reverse_rate()



    def get_chemicals_in_reaction(self, rxn_index :int) -> {int}:
        """
        Return a SET of indices (being a set, they're NOT in any particular order)
        of all the chemicals participating in the i-th reaction

        :param rxn_index:   An integer with the (zero-based) index to identify the reaction of interest
        :return:            A SET of indices of the chemicals involved in the above reaction
                                Note: being a set, it's NOT in any particular order
        """
        rxn = self.get_reaction(rxn_index)

        name_set = rxn.extract_chemicals_in_reaction()

        index_set = {self.chem_data.get_index(name) for name in name_set}

        return index_set


    def get_chemicals_indexes_in_reaction(self, rxn_index :int) -> [int]:
        """
        Return a sorted list of the indexes
        of all the chemicals participating in the i-th reaction

        :param rxn_index:   An integer with the (zero-based) index to identify the reaction of interest
        :return:            A sorted list of indices of the chemicals involved in the above reaction
        """
        index_set = self.get_chemicals_in_reaction(rxn_index)
        index_list = list(index_set)

        return sorted(index_list)



    def get_reactions_participating_in(self, chem_label :str, side :str) -> list:
        """
        Return a list of all the reactions that the given chemical species
        is involved in

        :param chem_label:  To identify a particular chemical
        :param side:        Either "reagent" or "product"
        :return:            List of various types of "Reaction" objects
        """
        assert side == "reagent" or side == "product"
        rxns_found_in = []
        for rxn in self.reaction_list:
            if side == "reagent":
                if chem_label in rxn.extract_reactant_labels():
                    rxns_found_in.append(rxn)
            else:
                 if chem_label in rxn.extract_product_labels():
                    rxns_found_in.append(rxn)


        return rxns_found_in



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



    def add_elementary_reaction(self, reactants :str|list, products :str|list,
                                reaction_type=None, temp=None,
                                **kwargs) -> int:
        """
        Create and register a new SINGLE elementary chemical reaction,
        optionally including its kinetic and/or thermodynamic data.
        All the involved chemicals can be either previously registered, or not.

        :param reactants:       A string or pair of strings; for reactions such as 2 A -> P, pass ["A", "A"]
        :param products:        A string or pair of strings; for reactions such as A -> 2 P, pass ["P", "P"]
        :param reaction_type:   A string with one of the following values:
                                    "ReactionUnimolecular", "ReactionSynthesis", "ReactionDecomposition"

        :param temp:            [OPTIONAL] Temperature in Kelvins
        :param kwargs:          [OPTIONAL] Other named arguments to pass to instantiate the various reaction
                                objects, such as `ReactionSynthesis`.
                                For list, see documentation of the classes in reactions.py

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


        if reaction_type is None:
            # Determine the reaction type from the number of reactants and products
            if n_reactants == 1 and n_products == 1:
                reaction_type = "ReactionUnimolecular"
            elif n_reactants == 2 and n_products == 1:
                reaction_type = "ReactionSynthesis"
            elif n_reactants == 1 and n_products == 2:
                reaction_type = "ReactionDecomposition"
            else:
                raise Exception(f"add_elementary_reaction(): {n_reactants} reactants and {n_products} products cannot correspond to any elementary reaction")


        match reaction_type:
            case "ReactionUnimolecular":
                rxn = ReactionUnimolecular(reactant=reactants[0], product=products[0], temp=temp, **kwargs)
            case "ReactionSynthesis":
                rxn = ReactionSynthesis(reactants=reactants, product=products[0], temp=temp, **kwargs)
            case "ReactionDecomposition":
                rxn = ReactionDecomposition(reactant=reactants[0], products=products, temp=temp, **kwargs)
            case _:
                raise Exception(f"add_elementary_reaction(): unknown reaction type ({reaction_type})")


        #print(f"add_elementary_reaction(): adding reaction of type `{reaction_type}`")

        return self.register_reaction(rxn=rxn, temp=temp)



    def add_reaction(self, reactants :str|list, products :str|list,
                     kF=None, kR=None,
                     enzyme=None, k1_F=None, k1_R=None, k2_F=None,
                     temp=None,
                     **kwargs) -> int:
        """
        Create and register a new SINGLE chemical reaction,
        optionally including its kinetic and/or thermodynamic data.
        All the involved chemicals can be either previously registered, or not.

        NOTE: in the reactants and products, if the stoichiometry coefficients aren't specified,
              they're assumed to be 1.

              The full structure of each term in the list of reactants and of products
              is the pair:  (stoichiometry coefficient, chemical label)

              EXAMPLES of formats to use for each term in the lists of the reactants and of the products:
                "F"         is taken to mean (1, "F") - default stoichiometry
                (2, "F")    is taken to mean (2, "F") - stoichiometry coefficient specified

              It's equally acceptable to use LISTS in lieu of tuples for the pairs

        :param reactants:       A string or list of pairs (stoichiometry, species name),
                                    or simplified terms in various formats; for details, see above
        :param products:        A string or list of pairs (stoichiometry, species name),
                                    or simplified terms in various formats; for details, see above

        :param temp:            [OPTIONAL] Temperature in Kelvins
        :param kwargs:          [OPTIONAL] Other named arguments to pass to instantiate the various reaction
                                objects, such as `ReactionSynthesis`.
                                For list, see documentation of the classes in reactions.py

        :return:                Integer index of the newly-added reaction
                                    (in the list self.reaction_list, stored as object variable)
        """
        #TODO: add an optional `reaction_type` argument;
        #      the dispatching might best be done elsewhere

        # Determine and assign the specific type of reaction, and instantiate an object of that class
        if enzyme is not None:
            reaction_type = "ReactionEnzyme"
            assert type(reactants) == str
            assert type(products) == str
            rxn = ReactionEnzyme(enzyme=enzyme, substrate=reactants, product=products,
                                 k1_F=k1_F, k1_R=k1_R, k2_F=k2_F, temp=temp, **kwargs)
        else:
            reactant_list = self._standardize_reaction_side(reactants, arg_name="reactants")
            product_list = self._standardize_reaction_side(products, arg_name="products")

            single_reactant = None
            if len(reactant_list) == 1 and reactant_list[0][0] == 1:    # A single reactant, with stoichiometry 1
                single_reactant = reactant_list[0][1]

            single_product = None
            if len(product_list) == 1 and product_list[0][0] == 1:      # A single product, with stoichiometry 1
                single_product = product_list[0][1]

            reaction_type = "ReactionGeneric"       # Default value, possibly changed below

            if single_reactant:    # A single reactant, with stoichiometry 1
                if single_product:      # A single product, with stoichiometry 1
                    reaction_type = "ReactionUnimolecular"
                    rxn = ReactionUnimolecular(reactant=single_reactant, product=single_product, kF=kF, kR=kR, temp=temp, **kwargs)
                elif len(product_list) == 2 and product_list[0][0] == 1 and product_list[1][0] == 1:      # Two products, both with stoichiometry 1
                    reaction_type = "ReactionDecomposition"
                    rxn = ReactionDecomposition(reactant=single_reactant, products=[product_list[0][1], product_list[1][1]],
                                                kF=kF, kR=kR, temp=temp, **kwargs)
                elif len(product_list) == 1 and product_list[0][0] == 2:      # A product with stoichiometry 2  (EXAMPLE : A <-> 2 B)
                    reaction_type = "ReactionDecomposition"
                    rxn = ReactionDecomposition(reactant=single_reactant, products=[product_list[0][1], product_list[0][1]],
                                                kF=kF, kR=kR, temp=temp, **kwargs)
            elif single_product:
                if len(reactant_list) == 2 and reactant_list[0][0] == 1 and reactant_list[1][0] == 1:      # Two reactants, both with stoichiometry 1
                    reaction_type = "ReactionSynthesis"
                    rxn = ReactionSynthesis(reactants=[reactant_list[0][1], reactant_list[1][1]],
                                            product=single_product, kF=kF, kR=kR, temp=temp, **kwargs)
                elif len(reactant_list) == 1 and reactant_list[0][0] == 2:  # A reactant with stoichiometry 2  (EXAMPLE : 2A <-> P)
                    reaction_type = "ReactionSynthesis"
                    rxn = ReactionSynthesis(reactants=[reactant_list[0][1], reactant_list[0][1]],
                                            product=single_product, kF=kF, kR=kR, temp=temp, **kwargs)

            if reaction_type == "ReactionGeneric":
                 rxn = ReactionGeneric(reactants=reactant_list, products=product_list, kF=kF, kR=kR, temp=temp, **kwargs)


        print(f"add_reaction(): detected reaction type `{reaction_type}`")

        return self.register_reaction(rxn=rxn, temp=temp)



    def register_reaction(self, rxn, temp=None) -> int:
        """
        Register a SINGLE chemical reaction from its reaction-specific object,
        and set all its kinetic and/or thermodynamic data from the available information,
        including the value of the temperature (stored in object variable.)

        All the involved chemicals can be either previously registered, or not;
        if not, they will get automatically registered.

        :param rxn: Object of one of the specific Reaction classes, such as
                        ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition,
                        ReactionEnz, ReactionGeneric
        :param temp:Temperature in degree Kelvin

        :return:    Integer index of the newly-added reaction
                        (in the list self.reaction_list, stored as object variable)
        """
        self.reaction_list.append(rxn)

        # Register any newly-encountered reactant not already registered
        # for aesthetic reasons, we'll do 1) catalyst, 2) reactants, 3) products

        #catalyst = rxn.extract_catalyst()

        #if catalyst:
            # Register the catalyst, if not already registered
            #self.chem_data.add_chemical(name=catalyst, skip_duplicates=True)

        # Register any newly-encountered reactant not already registered
        rxn_reactants = rxn.extract_reactant_labels()
        for label in rxn_reactants:
            self.chem_data.add_chemical(name=label, skip_duplicates=True)

        # Register any newly-encountered reaction product not already registered
        rxn_products = rxn.extract_product_labels()
        for label in rxn_products:
            self.chem_data.add_chemical(name=label, skip_duplicates=True)

        # Register any newly-encountered reaction intermediates not already registered
        rxn_intermediate = rxn.extract_intermediate()
        if rxn_intermediate:
            new_index = self.chem_data.add_chemical(name=rxn_intermediate, skip_duplicates=True)
            if new_index is not None:
                print(f"register_reaction() INFO: a reaction intermediates (`{rxn_intermediate}`), not explicitly registered, was automatically added to the chemical registry")

            if self.chem_data.get_diffusion_rate(chem_label=rxn_intermediate) is None:
                # Attempt to estimate the diffusion rate constant of the reaction intermediate
                D_enzyme = self.chem_data.get_diffusion_rate(chem_label=rxn_intermediate)
                D_substrate = self.chem_data.get_diffusion_rate(chem_label=rxn.substrate)
                # TODO: this might best belong elsewhere
                if (D_enzyme is not None) and (D_substrate is not None):
                    D_ES_rough_estimate = min(D_enzyme, D_substrate) * 0.9
                    print(f"register_reaction() INFO: diffusion rate for the reaction intermediates (`{rxn_intermediate}`), not yet specified, roughly estimated as {D_ES_rough_estimate}")
                    self.chem_data.set_diffusion_rate(chem_label=rxn_intermediate, diff_rate=D_ES_rough_estimate)


        involved_chemicals = set(rxn_reactants) | set(rxn_products)  # Union of sets
        if  rxn_intermediate:
             involved_chemicals = involved_chemicals | {rxn_intermediate}     # Union of sets

        #if catalyst is not None:
            #involved_chemicals = involved_chemicals - {catalyst}    # Difference between sets
            #self.active_enzymes.add(rxn.catalyst)                   # Add the new entry to a set

        # Update the set of "active chemicals"
        self.active_chemicals = self.active_chemicals | involved_chemicals  # Union of sets

        if temp is not None:
            rxn.set_thermodynamic_data(temp=temp)   # TODO: unclear if this is the best place to do this

        return len(self.reaction_list) - 1



    def clear_reactions_data(self) -> None:
        """
        Get rid of all the reactions; start again with "an empty slate" (but still with reference
        to the same data object about the chemicals and their properties)

        :return:    None
        """
        self.reaction_list = []
        self.active_chemicals = set()
        #self.active_enzymes = set()



    def inactivate_reaction(self, i :int) -> None:
        """
        Mark the i-th reaction as "inactive/disabled" (essentially, "deleted", but holding its positional
        index, to avoid a change in index in other reactions)

        TODO: Not yet supported by the dynamical modules; DON'T USE YET in simulations!

        :param i:   Zero-based index of the reaction to disable
        :return:    None
        """
        rxn = self.get_reaction(i)
        rxn.active = False

        # Re-construct self.active_chemicals and self.active_enzymes
        self.active_chemicals = set()
        #self.active_enzymes = set()

        for rxn in self.reaction_list:
            involved_chemicals = rxn.extract_chemicals_in_reaction()
            involved_chemicals = involved_chemicals - {rxn.catalyst}        # Set difference
            self.active_chemicals = self.active_chemicals.union(involved_chemicals)     # Union of sets
            #if rxn.catalyst is not None:
                #self.active_enzymes.add(rxn.catalyst)       # Add the new entry to a set





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

        chem_labels = self.labels_of_active_chemicals(sort_by_index=True)   # Set of chem labels, sorted by chemical index

        if self.chem_data.color_dict != {}:   # If plot colors were registered, show them alongside the chem labels
            chem_labels_with_colors = []
            for label in chem_labels:
                color = self.chem_data.get_plot_color(label)
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
                enzyme_color = self.chem_data.get_plot_color(enz)
                if enzyme_color:
                    print(f'  "{enz}" ({enzyme_color})')
                else:
                    print(f'  "{enz}"')
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
        rxn = self.get_reaction(rxn_index)

        return rxn.describe(concise)    # Invoke the individual reaction object



    def labels_of_active_chemicals(self, sort_by_index=False) -> [str]:
        """
        Return a list of the labels of all the chemicals
        involved in ANY of the registered reactions,
        but NOT counting chemicals that always appear
        in a catalytic role in all the reactions they participate in
        (if a chemical participates in a non-catalytic role in ANY reaction, it'll appear here)

        The list is not in any particular order, unless sort_by_index is True

        :param sort_by_index:   If True, the list is sorted by the index (order of registration)
                                    of the chemicals in it
        :return:                A set of chemical labels
        """
        if not sort_by_index:
            return list(self.active_chemicals)

        return sorted(self.active_chemicals, key=self.chem_data.get_index)



    def number_of_active_chemicals(self) -> int:
        """
        Return the number of all the chemicals
        involved in ANY of the registered reactions,
        but NOT counting chemicals that always appear
        in a catalytic role in all the reactions they participate in
        (if a chemical participates in a non-catalytic role in ANY reaction, it'll appear here)
        """
        return len(self.active_chemicals)



    def indexes_of_active_chemicals(self) -> [int]:
        """
        Return the ordered list (numerically sorted) of the INDEX numbers of all the chemicals
        involved in ANY of the registered reactions,
        but NOT counting chemicals that always appear in a catalytic role in all the reactions they
        participate in
        (if a chemical participates in a non-catalytic role in ANY reaction, it'll appear here.)

        EXAMPLE: [2, 7, 8]  if only those 3 chemicals (with indexes of, respectively, 2, 7 and 8)
                            are actively involved in ANY of the registered reactions

        CAUTION: the concept of "active chemical" might change in future versions, where only SOME of
                 the reactions are simulated
        """
        index_list = list(map(self.chem_data.get_index, self.active_chemicals))
        return sorted(index_list)





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
        graph = PyGraphVisual()     # Object to facilitate data preparation for graph visualization

        # Note: the graph nodes representing Chemicals will be given an id such as "C-123" and a label "Chemical";
        #       the graph nodes representing Reactions will be given an id such as "R-456" and a label "Reaction"

        for i, rxn in enumerate(self.reaction_list):    # Consider each REACTION in turn
            # Add a node representing the reaction
            rxn_id = f"RXN-{i}"               # Example: "RXN-456"
            node_data = {'name': 'RXN', 'formula': rxn.describe(concise=True)}

            # Show the parameter of the reaction
            rxn_properties = rxn.extract_rxn_properties()
            for k,v in rxn_properties.items():
                node_data[k] = f"{v:,.6g}"

            graph.add_node(node_id=rxn_id, labels='Reaction', properties=node_data)


            # Process all the PRODUCTS of this reaction
            products = rxn.extract_products()
            for term in products:
                species_name = rxn.extract_chem_label(term)
                chemical_id = f"C-{self.chem_data.get_index(species_name)}"      # Example: "C-12"

                # Add each product to the graph as a node (if not already present)
                properties={'name': species_name}
                if diff := self.chem_data.get_diffusion_rate(chem_label=species_name):
                    properties['diff_rate'] = diff
                graph.add_node( node_id=chemical_id, labels="Chemical",
                                properties=properties)

                # Append edge from "reaction node" to "product node"
                graph.add_edge(from_node=rxn_id, to_node=chemical_id, name="produces",
                               properties={'stoich': rxn.extract_stoichiometry(term)})


            # Process all the REACTANTS of this reaction
            reactants = rxn.extract_reactants()
            for term in reactants:
                species_name = rxn.extract_chem_label(term)
                chemical_id = f"C-{self.chem_data.get_index(species_name)}"      # Example: "C-34"

                # Add each reactant to the graph as a node (if not already present)
                properties={'name': species_name}
                if diff := self.chem_data.get_diffusion_rate(chem_label=species_name):
                    properties['diff_rate'] = diff
                graph.add_node(node_id=chemical_id, labels="Chemical",
                               properties=properties)

                # Append edge from "reactant node" to "reaction node"
                graph.add_edge(from_node=chemical_id, to_node=rxn_id, name="reacts",
                               properties={'stoich': rxn.extract_stoichiometry(term)})


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
