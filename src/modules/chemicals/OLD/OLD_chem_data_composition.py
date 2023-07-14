from typing import Union, List, NamedTuple
from src.modules.reactions.reaction import Reaction




class ChemData():
    def __init__(self, names=None, diffusion_rates=None):
        """

        :param names:           [OPTIONAL] A list with the names of the chemicals
        :param diffusion_rates: [OPTIONAL] A list or tuple with the diffusion rates of the chemicals
                                           If diffusion rates are provided, but no names given, the names will be
                                           auto-assigned as "Chemical 1", "Chemical 2", ...
        """
        self.core = ChemCore()
        self.diffusion = Diffusion(self.core)
        self.macromolecules = Macromolecules(self.core)

        self.reaction_list = []     # List of dicts.  Each item is an object of class "Reaction"

        self.temp = 298.15          # Temperature in Kelvins.  (By default, 25 C)
        # For now, assumed constant everywhere, and unvarying (or very slowly varying)


        if (names is not None) or (diffusion_rates is not None):
            self.init_chemical_data(names, diffusion_rates)


    #####################################################################################################

    '''                      ~   TO READ DATA STRUCTURES of the REACTIONS  ~                          '''

    def ________READ_RXN_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def get_index(self, name: str) -> int:
        return self.core.get_index(name)            # TODO: phase out

    def get_name(self, species_index: int) -> Union[str, None]:
        return self.core.get_name(species_index)    # TODO: phase out

    def number_of_chemicals(self) -> int:
        return self.core.number_of_chemicals()    # TODO: phase out

    def get_all_names(self) -> [Union[str, None]]:
        return self.core.get_all_names()         # TODO: phase out



    def number_of_reactions(self) -> int:
        # Return the number of registered chemical reactions
        return len(self.reaction_list)



    def assert_valid_rxn_index(self, index) -> None:
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



    def get_reaction(self, i: int) -> Reaction:
        """
        Return the data structure of the i-th reaction,
        in the order in which reactions were added (numbering starts at 0)

        :param i:   An integer that indexes the reaction of interest (numbering starts at 0)
        :return:    A dictionary with 4 keys ("reactants", "products", "kF", "kR"),
                    where "kF" is the forward reaction rate constant, and "kR" the back reaction rate constant
        """
        self.assert_valid_rxn_index(i)

        return self.reaction_list[i]



    def get_reactants(self, i: int) -> [(int, int, int)]:
        """
        Return a list of triplets with details of the reactants of the i-th reaction

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A list of triplets of the form (stoichiometry, species index, reaction order)
        """
        rxn = self.get_reaction(i)
        return rxn.extract_reactants()



    def get_reactants_formula(self, i) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula

        :param i:   The index (0-based) to identify the reaction of interest
        :return:
        """
        rxn = self.get_reaction(i)
        return rxn.extract_reactants_formula()



    def get_products(self, i: int) -> [(int, int, int)]:
        """
        Return a list of triplets with details of the products of the i-th reaction

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A list of triplets of the form (stoichiometry, species index, reaction order)
        """
        rxn = self.get_reaction(i)
        return rxn.extract_products()


    def get_products_formula(self, i) -> str:
        """
        Return a string with a user-friendly form of the right (products) side of the reaction formula

        :param i:   The index (0-based) to identify the reaction of interest
        :return:
        """
        rxn = self.get_reaction(i)
        return rxn.extract_products_formula()



    def get_forward_rate(self, i: int) -> float:
        """

        :param i:   The integer index (0-based) to identify the reaction of interest
        :return:    The value of the forward rate constant for the above reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_forward_rate()


    def get_reverse_rate(self, i: int) -> float:
        """

        :param i:   The integer index (0-based) to identify the reaction of interest
        :return:    The value of the reverse (back) rate constant for the above reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_reverse_rate()



    def get_chemicals_in_reaction(self, rxn_index: int) -> {int}:
        """
        Return a SET of indices (being a set, it's NOT in any particular order)
        of all the chemicals in the specified reaction

        :param rxn_index:   An integer with the (zero-based) index to identify the reaction of interest
        :return:            A SET of indices of the chemicals involved in the above reaction
                                Note: being a set, it's NOT in any particular order
        """
        rxn = self.get_reaction(rxn_index)

        return rxn.extract_chemicals_in_reaction()



    def get_reactions_participating_in(self, species_index: int) -> [Reaction]:
        """
        Return a list of all the reactions that the given chemical species
        is involved in

        :param species_index:
        :return:                List of "Reaction" objects
        """
        pass        # TODO: write




    #####################################################################################################

    '''                                  ~   TO SET DATA  ~                                           '''

    def ________SET_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def init_chemical_data(self, names=None, diffusion_rates=None)  -> None:
        """
        Initialize the names (if provided) and diffusion rates (if provided)
        of all the chemical species, in the given order.
        If no names are provided, the strings "Chemical 1", "Chemical 2", ..., are used

        IMPORTANT: this function can be invoked only once, before any chemical data is set.
                   To add new chemicals later, use add_chemical()

        :param names:           [OPTIONAL] List or tuple of the names of the chemical species
        :param diffusion_rates: [OPTIONAL] A list or tuple with the diffusion rates of the chemicals,
                                           in the same order as the names
        :return:                None
        """
        # Validate
        assert self.core.number_of_chemicals() == 0, \
            f"ChemData.init_chemical_data(): function can be invoked only once, before any chemical data is set"

        n_species = self.core.number_of_chemicals()

        if names:
            arg_type = type(names)
            assert arg_type == list or arg_type == tuple, \
                f"ChemData.init_chemical_data(): the `names` argument must be a list or tuple.  What was passed was of type {arg_type}"
            if n_species != 0:
                assert len(names) == n_species, \
                    f"ChemData.init_chemical_data(): the passed number of names ({len(names)}) " \
                    f"doesn't match the previously-set number of chemical species (n_species)"

        if diffusion_rates:
            arg_type = type(diffusion_rates)
            assert arg_type == list or arg_type == tuple, \
                f"ChemData.init_chemical_data(): the diffusion_rates must be a list or tuple.  What was passed was of type {arg_type}"
            if n_species != 0:
                assert len(diffusion_rates) == n_species, \
                    f"ChemData.init_chemical_data(): the passed number of diffusion rates ({len(diffusion_rates)}) " \
                    f"doesn't match the previously-set number of chemical species ({n_species})"

        if diffusion_rates and names:
            assert len(names) == len(diffusion_rates), \
                f"ChemData.init_chemical_data(): the supplied names and diffusion_rates " \
                f"don't match in number ({len(names)} vs. {len(diffusion_rates)})"


        # Populate the data structure
        if names is None:
            if diffusion_rates is None:   # The strings "Chemical 1", "Chemical 2", ..., will be used as names
                for i in range(n_species):
                    assigned_name = f"Chemical {i+1}"
                    self.core.add_chemical_species(assigned_name)
                    #self.chemical_data.append({"name": assigned_name})
                    #self.name_dict[assigned_name] = i
            else:
                for i, diff in enumerate(diffusion_rates):
                    assigned_name = f"Chemical {i+1}"
                    self.diffusion.assert_valid_diffusion(diff)
                    self.core.add_chemical_species(assigned_name)
                    self.diffusion.set_diffusion_rate(name=assigned_name, diff_rate=diff)
                    #self.chemical_data.append({"name": assigned_name, "diff": diff})
                    #self.name_dict[assigned_name] = i

        else:   # names is not None
            for i, chem_name in enumerate(names):
                assert type(chem_name) == str, \
                    f"ChemData.init_chemical_data(): all the names must be strings.  The passed value ({chem_name}) is of type {type(chem_name)}"
                if diffusion_rates is None:
                    self.core.add_chemical_species(chem_name)
                    #self.chemical_data.append({"name": chem_name})
                else:
                    diff = diffusion_rates[i]
                    self.diffusion.assert_valid_diffusion(diff)
                    self.core.add_chemical_species(chem_name)
                    self.diffusion.set_diffusion_rate(name=chem_name, diff_rate=diff)
                    #self.chemical_data.append({"name": chem_name, "diff": diff})

                #self.name_dict[chem_name] = i

        #self.n_species = len(self.chemical_data)



    def add_chemical(self, name: str, diffusion_rate=None, note=None) -> None:
        """
        Register a new chemical species, with a name and (optionally) a diffusion rate.

        :param name:            Name of the chemical species to add
        :param diffusion_rate:  [OPTIONAL] Floating-point number with the diffusion rate (in water) of this chemical
        :param note:            [OPTIONAL] Note to attach to the chemical
        :return:                None
        """
        assert type(name) == str, \
            f"ChemData.add_chemical(): a chemical's name must be provided, as a string value.  " \
            f"What was passed was of type {type(name)}"

        if diffusion_rate:
            self.diffusion.assert_valid_diffusion(diffusion_rate)
            self.diffusion.set_diffusion_rate(name=name, diff_rate=diffusion_rate)

        self.core.add_chemical_species(name=name, note=note)



    def set_temp(self, temp, units="K") -> None:
        """
        Specify the temperature of the environment
        (for now assumed uniform everywhere)

        :param temp:    Temperature, in Kelvins, or None
        :param units:   Not yet implemented
        :return:        None
        """
        self.temp = temp



    def add_reaction(self, reactants: Union[int, str, list], products: Union[int, str, list],
                     forward_rate=None, reverse_rate=None,
                     delta_H=None, delta_S=None, delta_G=None) -> Reaction:
        """
        Add the parameters of a SINGLE reaction, optionally including kinetic and/or thermodynamic data.
        The involved chemicals must be already registered - use add_chemical() if needed.

        NOTE: in the reactants and products, if the stoichiometry and/or reaction order aren't specified,
              they're assumed to be 1.
              Their full structure is the triplet (stoichiometry coefficient, name, reaction order)

        EXAMPLES of formats for each item of the reactants and products
        (*assuming* that the chemical species with index 5 is called "F"):
                    "F"         gets turned into:   (1, 5, 1)
                    (3, "F")                        (3, 5, 1)
                    (3, "F", 2)                     (3, 5, 2)
                    It's equally acceptable to use LISTS in lieu of tuples

        :param reactants:       A list of triplets (stoichiometry, species name or index, reaction order),
                                    or simplified terms in various formats; for details, see above.
                                    If not a list, it will get turned into one
        :param products:        A list of triplets (stoichiometry, species name or index, reaction order of REVERSE reaction),
                                    or simplified terms in various formats; for details, see above.
                                    If not a list, it will get turned into one
        :param forward_rate:    [OPTIONAL] Forward reaction rate constant
        :param reverse_rate:    [OPTIONAL] Reverse reaction rate constant
        :param delta_H:         [OPTIONAL] Change in Enthalpy (from reactants to products)
        :param delta_S:         [OPTIONAL] Change in Entropy (from reactants to products)
        :param delta_G:         [OPTIONAL] Change in Free Energy (from reactants to products)
        :return:                Object of type "Reaction"
                                (note: the object variable self.reaction_list gets appended to)
        """

        rxn = Reaction(self, reactants, products, forward_rate, reverse_rate,
                       delta_H, delta_S, delta_G)

        self.reaction_list.append(rxn)

        return rxn



    def clear_reactions_data(self) -> None:
        """
        Get rid of all reactions; start again with "an empty slate" (but still with reference
        to the same data object about the chemicals)

        :return:    None
        """
        self.reaction_list = []





    #####################################################################################################

    '''                                ~   TO DESCRIBE THE DATA  ~                                    '''

    def ________DESCRIBE_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def describe_reactions(self, concise=False) -> None:
        """
        Print out a user-friendly plain-text form of all the reactions.
        If wanting to describe just 1 reaction, use single_reaction_describe()

        EXAMPLE (not concise):
            Number of reactions: 2 (at temp. 25 C)
            (0) CH4 + 2 O2 <-> CO2 + 2 H2O  (kF = 3.0 / kR = 2.0 / Delta_G = -1,005.13 / K = 1.5) | 1st order in all reactants & products"
            (1) A + B <-> C  (kF = 5.0 / kR = 1.0 / Delta_G =  / K = 5.0) | 1st order in all reactants & products"

        :param concise:     If True, less detail is shown
        :return:            None
        """
        print(f"Number of reactions: {self.number_of_reactions()} (at temp. {self.temp - 273.15:,.4g} C)")
        for description in self.multiple_reactions_describe(concise=concise):
            print(description)



    def multiple_reactions_describe(self, rxn_list=None, concise=False) -> [str]:
        """
        The counterpart of single_reaction_describe() for many reactions

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

        return rxn.describe(concise)




    #####################################################################################################

    '''                          ~   FOR CREATION OF NETWORK DIAGRAMS  ~                              '''

    def ________NETWORK_DIAGRAMS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def prepare_graph_network(self) -> dict:
        """

        :return:    A dictionary with 2 keys: 'graph' and 'color_mapping'
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

        :return:    A list of dictionaries.  Each dictionary must have an 'id' key with a unique value.
                    EXAMPLE, for an  A <-> B reaction:
                       [{'id': 0, 'label': 'Chemical', 'name': 'A', 'diff_rate': None, 'stoich': 1, 'rxn_order': 1},
                        {'id': 1, 'label': 'Chemical', 'name': 'B', 'diff_rate': None, 'stoich': 1, 'rxn_order': 1},

                        {'id': 2, 'label': 'Reaction', 'name': 'RXN', 'kF': 3.0, 'kR': 2.0, 'K': 1.5, 'Delta_G': -1005.13},

                        {'id': 3, 'name': 'produces', 'source': 2, 'target': 1},
                        {'id': 4, 'name': 'reacts', 'source': 0, 'target': 2}
                       ]
        """
        graph_data = []
        species_in_graph = []

        # Note: the species index of the various chemicals is a UNIQUE number; so, it's suitable to be used as an ID for the nodes
        #       For the reaction nodes, use numbers from a range starting just above the end-range of the numbers for the chemicals
        next_available_id = self.core.number_of_chemicals()

        for i, rxn in enumerate(self.reaction_list):    # Consider each reaction in turn
            # Add a node representing the reaction
            rxn_id = next_available_id
            next_available_id += 1
            node_data = {'id': rxn_id, 'label': 'Reaction', 'name': 'RXN'}

            rxn_properties = rxn.extract_rxn_properties()
            for k,v in rxn_properties.items():
                node_data[k] = f"{v:,.6g}"
                #'kF': self.extract_forward_rate(rxn), 'kR': self.extract_back_rate(rxn)})

            graph_data.append(node_data)

            # Process all products
            products = rxn.extract_products()
            for term in products:
                species_index = term[1]
                # Add each product to the graph as a node (if not already present)
                if species_index not in species_in_graph:
                    graph_data.append({'id': species_index, 'label': 'Chemical',
                                       'name': self.core.get_name(species_index),
                                       'diff_rate': self.diffusion.get_diffusion_rate(species_index),
                                       'stoich': rxn.extract_stoichiometry(term),
                                       'rxn_order': rxn.extract_rxn_order(term)
                                       })
                # Append edge from "reaction node" to product
                graph_data.append({'id': next_available_id, 'source': rxn_id, 'target': species_index, 'name': 'produces'})
                next_available_id += 1

            # Process all reactants
            reactants = rxn.extract_reactants()
            for term in reactants:
                species_index = term[1]
                # Add each reactant to the graph as a node (if not already present)
                if species_index not in species_in_graph:
                    graph_data.append({'id': species_index, 'label': 'Chemical',
                                       'name': self.core.get_name(species_index),
                                       'diff_rate': self.diffusion.get_diffusion_rate(species_index),
                                       'stoich': rxn.extract_stoichiometry(term),
                                       'rxn_order': rxn.extract_rxn_order(term)
                                       })
                # Append edge from reactant to "reaction node"
                graph_data.append({'id': next_available_id, 'source': species_index, 'target': rxn_id, 'name': 'reacts'})
                next_available_id += 1

        return graph_data



    def assign_color_mapping(self):
        return {
            'Chemical': 'neo4j_green',
            #'Product': 'neo4j_red',
            'Reaction': 'neo4j_lightbrown'
        }





########################################################################################################

class ChemCore():
    """
    Data about all the chemicals and (if applicable) reactions,
    including:
        - names
        - diffusion rates
        - macro-molecules Binding Site Affinities (for Transcription Factors)
        - reaction data (see class "Reaction", in "reaction_data.py")


    Note: for now, the temperature is assumed constant everywhere, and unvarying (or very slowly varying)
    """

    def __init__(self):
        """
        If chemical names and their diffusion rates are both provided, they must have the same count,
        and appear in the same order.
        It's ok not to pass any data, and later add it.
        Reactions, if applicable, need to be added later by means of calls to add_reaction()
        Macro-molecules, if applicable, need to be added later

        """

        self.chemical_data = []     # Basic data for all chemicals, *except* water and macro-molecules
        # EXAMPLE: [{"name": "A"} ,
        #           {"name": "B", "note": "some note"}]
        # The position in the list is referred to as the "index" of that chemical
        # TODO: maybe use a Pandas dataframe

        self.name_dict = {}         # To map assigned names to their positional index (in the ordered list of chemicals);
        # this is automatically set and maintained
        # EXAMPLE: {"A": 0, "B": 1, "C": 2}




    #####################################################################################################

    '''                       ~   TO READ DATA STRUCTURES of the CHEMICALS  ~                         '''

    def ________READ_CHEM_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def number_of_chemicals(self) -> int:
        """
        Return the number of registered chemicals - exclusive of water and of macro-molecules

        :return:
        """
        return len(self.chemical_data)



    def assert_valid_species_index(self, species_index: int) -> None:
        """
        Raise an Exception if the specified species_index (meant to identify a chemical) isn't valid

        :param species_index:   An integer that indexes the chemical of interest (numbering starts at 0)
        :return:                None
        """
        n_species = self.number_of_chemicals()
        assert (type(species_index) == int) and (0 <= species_index < n_species), \
            f"The requested species index ({species_index}) is not the expected integer the range [0 - {n_species - 1}], inclusive"



    def get_index(self, name: str) -> int:
        """
        Return the index of the chemical species with the given name.
        If not found, an Exception is raised

        :param name:    Name of the chemical species of interest
        :return:        The index of the species with the given name
        """
        index =  self.name_dict.get(name)
        assert index is not None, \
            f"ChemData.get_index(): No chemical species named `{name}` was found"

        return index



    def get_name(self, species_index: int) -> Union[str, None]:
        """
        Return the name of the species with the given index.
        If no name was assigned, return None.

        :param species_index:   An integer (starting with zero) corresponding to the
                                    original order with which the chemical species were first added
        :return:                The name of the species with the given index if present,
                                    or None if not
        """
        self.assert_valid_species_index(species_index)

        return self.chemical_data[species_index].get("name")    # If "name" is not present, None will be returned



    def get_all_names(self) -> [Union[str, None]]:
        """
        Return a list with the names of all the chemical species, in their index order.
        If any is missing, None is used     TODO: raise an Exception instead

        :return:    A list of strings
        """
        return [c.get("name", None)
                for c in self.chemical_data]



    def add_chemical_species(self, name, note=None) -> int:     # TODO: test
        """

        :param name:
        :param note:
        :return:
        """
        index = len(self.chemical_data)     # The next available index number (list position)
        self.name_dict[name] = index

        if note is None:
            self.chemical_data.append({"name": name})
        else:
            self.chemical_data.append({"name": name, "note": note})

        return index




########################################################################################################

class Diffusion():

    def __init__(self, core):
        self.diffusion_rates = {}   # EXAMPLE: {"A": 6.4, "B": 12.0}
        self.core = core



    def assert_valid_diffusion(self, diff) -> None:
        """
        Raise an Exception if the specified diffusion value isn't valid

        :param diff:    Diffusion rate
        :return:        None
        """
        assert type(diff) == float or type(diff) == int, \
            f"Diffusion.assert_valid_diffusion(): The value for the diffusion rate ({diff}) is not a number (float or int)"

        assert diff >= 0., \
            f"Diffusion.assert_valid_diffusion(): the diffusion rate ({diff}) cannot be negative"



    def get_diffusion_rate(self, species_index: int) -> Union[str, None]:
        """
        Return the diffusion rate of the chemical species with the given index.
        If no value was assigned, return None.
        TODO: also accept lookup by name

        :param species_index:   An integer (starting with zero) corresponding to the
                                    original order with which the chemical species were first added
        :return:                The value of the diffusion rate for the species with the given index if present,
                                    or None if not
        """
        self.core.assert_valid_species_index(species_index)

        name = self.core.get_name(species_index)

        return self.diffusion_rates.get(name)      # If not present, None is returned



    def get_all_diffusion_rates(self) -> list:
        """
        Return a list of the diffusion rates of all the chemicals,
        in the order of their indexes.

        If any value is missing, None is used for it

        :return:    A list of numbers with the diffusion rates
        """
        return [self.diffusion_rates.get(name) for name in self.core.get_all_names()]
        # If any value is not present, None is used for it



    def missing_diffusion_rate(self) -> bool:
        """
        Determine whether any of the registered chemicals has a missing diffusion rates

        :return:    True if any of the diffusion rates (for all the registered chemicals) is missing;
                        False otherwise
        """
        if len(self.diffusion_rates) < self.core.number_of_chemicals():
            return True

        for name, diff in self.diffusion_rates.items():
            if diff is None:
                return True     # TODO: this should never occur

        return False



    def set_diffusion_rate(self, name: str, diff_rate) -> None:
        """
        Set the diffusion rate of the given chemical species (identified by its name)

        :param name:        Name of a chemical species
        :param diff_rate:   Diffusion rate (in water) for the above chemical
        :return:            None
        """
        self.assert_valid_diffusion(diff_rate)
        self.diffusion_rates[name] = diff_rate





########################################################################################################


class ChemicalAffinity(NamedTuple):
    # Used for binding to macromolecules (e.g. Transcription Factors to DNA)
    chemical: str
    affinity: float



class Macromolecules():
    """
    Modeling of large molecules (such as DNA) with multiple binding sites (for example,
    for Transcription Factors)
    """

    def __init__(self, core):
        """

        :param core:
        """
        self.core = core

        self.macro_molecules = []   # List of names.  EXAMPLE: ["M1", "M2"]
        # The position in the list is referred to as the "index" of that macro-molecule
        # Names will be enforced to be unique

        self.binding_sites = {}     # A dict whose keys are macromolecule names.  The values are in turn dicts, indexed by site number.
        # EXAMPLE:
        #       {"M1": {1: ChemicalAffinity("A", 2.4), 2: ChemicalAffinity("C", 5.1)},
        #        "M2": {1: ChemicalAffinity("C", 9.1), 2: ChemicalAffinity("B", 0.3), 3: ChemicalAffinity("A", 1.8), 4: ChemicalAffinity("C", 2.3)}
        #        }
        #       where "M1", "M2" are macro-molecules, and "A", "B", "C" are bulk chemicals (such as transcription factors),
        #           all previously-declared;
        #           the various ChemicalAffinity's are NamedTuples (objects) storing a bulk-chemical name and its affinity at that site.

        # Info on Binding Site Affinities : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6787930/

        #       OLD: {"M1": {"A": 2.4, "B": 853.} }       # Alt: [("A", 2.4), ("B", 853.)]


    def add_macromolecules(self, names: Union[str, List[str]]) -> None:
        """
        Register one or more macromolecule species, specified by their name(s)
        Note: this is a register of names, NOT dynamical information about counts of macromolecules in the system
              (that's the domain of the class ReactionDynamics)

        :param names:   A string, or list of strings, with the name(s) of the macromolecule(s)
        :return:        None.  The object variable self.macro_molecules will get modified.
        """
        if type(names) == str:  # If a single name was given
            names = [names]

        for m in names:
            if m in self.macro_molecules:
                print(f"WARNING: Macromolecule `{m}` was already registered.  Skipped...")    # Redundant attempt to re-register an existing name
            else:
                self.macro_molecules.append(m)  # Grow the list of registered macromolecules



    def get_macromolecules(self) -> [str]:
        """
        Return a list of the names of all the registered macromolecules

        :return:    A (possibly empty) list of names of all the registered macromolecules
        """
        return self.macro_molecules



    def set_binding_site_affinity(self, macromolecule: str, site_number: int, chemical: str, affinity) -> None:
        """
        Set the values of the binding affinity of the given macromolecule, at the indicated site on it,
        for the specified chemical species.

        Any previously-set value (for that macromolecule/site_number/chemical) will get over-written.

        Only 1 chemical can be associated to a particular site of a given macromolecule; attempting
        to associate another one will result in error.

        NOTE: at present, no allowance is made for the fact that if the macromolecule is already bound
              to chemical "X" then its affinity of "Y" might be different than in the absence of "X"

        :param macromolecule:   Name of a macromolecule; if not previously-declared,
                                    it will get added to the list of registered macromolecules
        :param site_number:     Integer to identify a binding site on the macromolecule
        :param chemical:        Name of a previously-declared (bulk) chemical;
                                    if not found, an Exception will be raised       TODO: inconsistent with "macromolecule" arg
        :param affinity:        A number, in units of concentration
        :return:                None
        """
        assert chemical in self.core.get_all_names(), \
            f"set_binding_site_affinity(): no chemical named `{chemical}` found; use add_chemical() first"

        assert type(site_number) == int, \
            f"set_binding_site_affinity(): the argument `site_number` must be an integer"

        if macromolecule not in self.macro_molecules:
            self.add_macromolecules([macromolecule])

        if self.binding_sites.get(macromolecule) is None:
            self.binding_sites[macromolecule] = {}          # Initialize the dict of binding sites for this macromolecule

        binding_data = self.binding_sites[macromolecule]    # This will be a dict whose key is the site_number

        if site_number in binding_data:
            existing_affinity_data = binding_data[site_number]
            if existing_affinity_data.chemical != chemical:
                raise Exception(f"set_binding_site_affinity(): "
                                f"site number {site_number} of macromolecule `{macromolecule}` was previously associated to chemical `{existing_affinity_data.chemical}` "
                                f"(attempting to set an affinity value for chemical `{chemical}`)")

        binding_data[site_number] = ChemicalAffinity(chemical=chemical, affinity=affinity)



    def get_binding_site_affinity(self, macromolecule: str, site_number: int) -> ChemicalAffinity:
        """
        Return the value of the binding affinity of the given macromolecule.
        If no value was previously set, an Exception is raised

        :param macromolecule:   Name of a macromolecule; if not found, an Exception will get raised
        :param site_number:     Integer to identify a binding site on the macromolecule
        :return:                The NamedTuple (chemical, affinity)
                                    if no value was previously set, an Exception is raised
        """
        assert macromolecule in self.macro_molecules, \
            f"get_binding_site_affinity(): no macromolecule named `{macromolecule}` found"

        assert type(site_number) == int, \
            f"get_binding_site_affinity(): the argument `site_number` must be an integer"

        binding_data = self.binding_sites.get(macromolecule)    # This will be a dict whose key is the site_number.  (None if not found)
        assert binding_data is not None, "not found"

        chem_affinity = binding_data.get(site_number)           # Will be None if not found
        assert chem_affinity is not None, "not found"

        return chem_affinity



    def get_binding_sites(self, macromolecule) -> [int]:
        """
        Get a list of all the binding-site numbers of the given macromolecule.
        If the requested macromolecule isn't registered, an Exception will be raised

        :param macromolecule:   The name of a macromolecule
        :return:                A (possibly empty) list of integers, representing the binding-site numbers.
                                    EXAMPLE: [1, 2]
        """
        binding_site_info = self.binding_sites.get(macromolecule)   # EXAMPLE: {1: ChemicalAffinity("A", 2.4), 2: ChemicalAffinity("C", 5.1)}
        #   will be None if no binding sites are found
        if binding_site_info is None:
            assert macromolecule in self.get_macromolecules(), \
                f"get_binding_sites(): the requested macromolecule ({macromolecule}) isn't registered; use add_macromolecules()"
            return []
        return list(binding_site_info)                              # EXAMPLE: [1, 2]



    def get_binding_sites_and_ligands(self, macromolecule) -> dict:
        """
        Return a mapping (python dict) from binding-site number to ligand species, for the given macromolecule
        If the requested macromolecule isn't registered, an Exception will be raised

        :param macromolecule:   The name of a macromolecule
        :return:                A dict whose keys are binding-site numbers and values are their respetive ligands
                                    EXAMPLE: {1: "A", 2: "C"}
        """
        binding_site_info = self.binding_sites.get(macromolecule)   # EXAMPLE: {1: ChemicalAffinity("A", 2.4), 2: ChemicalAffinity("C", 5.1)}
        #   will be None if no binding sites are found
        if binding_site_info is None:
            assert macromolecule in self.get_macromolecules(), \
                f"get_binding_sites(): the requested macromolecule ({macromolecule}) isn't registered; use add_macromolecules()"
            return {}

        d = {}
        for (site_number, affinity_obj) in binding_site_info.items():
            ligand = affinity_obj.chemical
            d[site_number] = ligand

        return d        # EXAMPLE: {1: "A", 2: "C"}



    def reset_macromolecule(self, macromolecule):
        """

        :param macromolecule:
        :return:
        """
        if self.binding_sites.get(macromolecule) is not None:
            del self.binding_sites[macromolecule]



    def clear_macromolecules(self):
        """
        Reset all macromolecules to their original state

        :return:    None
        """
        self.macro_molecules = []
        self.binding_sites = {}
