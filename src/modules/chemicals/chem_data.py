from typing import Union, List, NamedTuple, Set
from src.modules.reactions.reaction import Reaction
from src.modules.visualization.py_graph_visual.py_graph_visual import PyGraphVisual
from src.modules.html_log.html_log import HtmlLog as log
from src.modules.visualization.graphic_log import GraphicLog
import pandas as pd



class ChemCore:
    """
    Core data about the chemical species, such as their names and indexes (position in their listing)

    Note: end users will generally utilize the class ChemData, which extends this one
    """

    def __init__(self):

        self.chemical_data = [] # Basic data for all chemicals, *except* water and macro-molecules
                                # Each list entry represents 1 chemical,
                                # and is a dict required to contain the key "name", plus any other keys
                                # EXAMPLE: [{"name": "A"} ,
                                #           {"name": "B", "note": "some note"}]
                                # The position in this list is referred to as the "INDEX" of that chemical
                                # TODO: maybe use a Pandas dataframe

        self.name_dict = {}     # To map assigned names to their positional index (in the ordered list of chemicals, self.chemical_data);
                                # this is automatically set and maintained
                                # EXAMPLE: {"A": 0, "B": 1, "C": 2}

        self.active_chemicals = set()   # Set of the INDEXES of the chemicals - not counting catalysts - involved
                                        # in any of the registered reactions

        self.active_enzymes = set()     # Set of the INDEXES of the enzymes (catalysts) involved in any of the registered reactions



    def number_of_chemicals(self) -> int:
        """
        Return the number of registered chemicals - exclusive of water and of macro-molecules

        :return:    The number of registered chemicals - exclusive of water and of macro-molecules
        """
        return len(self.chemical_data)



    def assert_valid_species_index(self, species_index: int) -> None:
        """
        Raise an Exception if the specified species_index (meant to identify a chemical) isn't valid

        :param species_index:   An integer that indexes the chemical of interest (numbering starts at 0)
        :return:                None
        """
        n_species = self.number_of_chemicals()
        assert type(species_index) == int,  f"The specified species index ({species_index}) must be an integer: " \
                                            f"instead, it has type {type(species_index)}"

        assert 0 <= species_index < n_species, \
            f"The specified species index ({species_index}) is not in the expected range [0 - {n_species - 1}], inclusive.  " \
            f"I.e., there is no chemical species assigned to this index"



    def get_index(self, name: str) -> int:
        """
        Return the index of the chemical species with the given name.
        If not found, an Exception is raised

        :param name:    Name of the chemical species of interest
        :return:        The index of the species with the given name
        """
        index = self.name_dict.get(name)
        assert index is not None, \
            f"ChemData.get_index(): No chemical species named `{name}` was found"

        return index



    def get_name(self, species_index: int) -> str:
        """
        Return the name of the species with the given index.

        :param species_index:   An integer (starting with zero) corresponding to the
                                    original order with which the chemical species were first added
        :return:                The name of the species with the given index.
                                    If missing or blank, an Exception is raised
        """
        self.assert_valid_species_index(species_index)

        name = self.chemical_data[species_index].get("name")    # If "name" is not present, None will be returned

        assert name, \
            f"get_name(): A chemical species with the requested index ({species_index}) is present, but it lacks a name"

        return name



    def get_all_names(self) -> [str]:
        """
        Return a list with the names of all the chemical species, in their index order.
        If any is missing or blank, an Exception instead

        :return:    A list of strings with the chemical names,
                        in their registered index order
        """
        all_names = []
        for i, c in enumerate(self.chemical_data):
            name = c.get("name", None)
            assert name, f"get_all_names(): missing or blank chemical name in index position {i}"
            all_names.append(name)

        return all_names



    def all_chemicals(self) -> pd.DataFrame:
        """
        Returns a Pandas dataframe with all the known information
        about all the registered chemicals (not counting macro-molecules),
        in their index order

        :return:    A Pandas dataframe
        """
        return pd.DataFrame(self.chemical_data)



    def add_chemical(self, name :str, note=None, **kwargs) -> int:
        """
        Register a new chemical species, with a name
        and (optionally) :
            - a note
            - any other named argument(s) that the user wishes to store (i.e. arbitrary named arguments)

        EXAMPLE:  add_chemical("P1", note = "my note about P1", full_name = "protein P1")

        Note: if also wanting to set the diffusion rate in a single function call,
              use ChemData.add_chemical_with_diffusion() instead

        :param name:    Name of the new chemical species to register
        :param note:    (OPTIONAL) string to attach to the given chemical
        :param kwargs:  (OPTIONAL) dictionary of named arguments
        :return:        The integer index assigned to the newly-added chemical
        """
        assert type(name) == str, \
            f"add_chemical(): a chemical's name must be provided, as a string value.  " \
            f"What was passed ({name}) was of type {type(name)}"

        assert name != "", \
            f"add_chemical(): the chemical's name cannot be a blank string"

        index = len(self.chemical_data)     # The next available index number (list position)
        self.name_dict[name] = index

        if note is None:
            d = {"name": name}
        else:
            d = {"name": name, "note": note}

        d.update(kwargs)        # Merge dictionary kwargs into d

        self.chemical_data.append(d)

        return index




###############################################################################################################
###############################################################################################################


class Diffusion(ChemCore):
    """
    Extends its parent class, to manage diffusion-related data

    End users will generally utilize the class ChemData, which extends this one
    """

    def __init__(self):

        super().__init__()          # Invoke the constructor of its parent class

        self.diffusion_rates = {}   # Values for the diffusion rates, indexed by chemical name
                                    # EXAMPLE: {"A": 6.4, "B": 12.0}



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
        self.assert_valid_species_index(species_index)

        name = self.get_name(species_index)

        return self.diffusion_rates.get(name)      # If not present, None is returned



    def get_all_diffusion_rates(self) -> list:
        """
        Return a list of the diffusion rates of all the chemicals,
        in the order of their indexes.

        If any value is missing, None is used for it

        :return:    A list of numbers with the diffusion rates
        """
        return [self.diffusion_rates.get(name) for name in self.get_all_names()]
        # If any value is not present, None is used for it



    def missing_diffusion_rate(self) -> bool:
        """
        Determine whether any of the registered chemicals has a missing diffusion rates

        :return:    True if any of the diffusion rates (for all the registered chemicals) is missing;
                        False otherwise
        """
        if len(self.diffusion_rates) < self.number_of_chemicals():
            return True

        for name, diff in self.diffusion_rates.items():
            if diff is None:
                return True     # TODO: this should never occur

        return False



    def set_diffusion_rate(self, name :str, diff_rate) -> None:
        """
        Set the diffusion rate of the given chemical species (identified by its name)

        :param name:        Name of a chemical species
        :param diff_rate:   Diffusion rate (in water) for the above chemical
        :return:            None
        """
        self.assert_valid_diffusion(diff_rate)
        self.diffusion_rates[name] = diff_rate




###############################################################################################################
###############################################################################################################


class AllReactions(Diffusion):
    """
    Extends its parent class, to manage reaction-related data

    End users will generally utilize the class ChemData, which extends this one
    """

    def __init__(self):

        super().__init__()          # Invoke the constructor of its parent class

        self.reaction_list = []     # List of dicts.  Each item is an object of class "Reaction"

        self.temp = 298.15          # Temperature in Kelvins.  (By default, 25 C)
        # For now, assumed constant everywhere, and unvarying (or very slowly varying)



    def number_of_reactions(self) -> int:
        """
        Return the number of registered chemical reactions

        :return:
        """
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



    def add_chemical_with_diffusion(self, name :str, diffusion_rate, note=None, **kwargs) -> int:
        """
        Register a new chemical species, with a name, a diffusion rate (in water),
        and (optionally) :
            - a note
            - any other named argument(s) that the user wishes to store (i.e. arbitrary named arguments)

        EXAMPLE:  add_chemical("P1", diffusion_rate = 0.1, note = "my note about P1", full_name = "protein P1")

        Note: if no diffusion is to be set, use the method add_chemical()

        :param name:            Name of the chemical species to add
        :param diffusion_rate:  Floating-point number with the diffusion rate (in water) of this chemical
        :param note:            [OPTIONAL] Note to attach to the chemical
        :param kwargs:          [OPTIONAL] Any arbitrary named arguments
        :return:                The integer index assigned to the newly-added chemical
        """
        # Validate the diffusion_rate, if present; other arguments get validated by add_chemical()
        self.assert_valid_diffusion(diffusion_rate)

        index = self.add_chemical(name=name, note=note, **kwargs)

        self.set_diffusion_rate(name=name, diff_rate=diffusion_rate)

        return index



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
                     delta_H=None, delta_S=None, delta_G=None) -> int:
        """
        Register a new SINGLE chemical reaction,
        optionally including its kinetic and/or thermodynamic data.
        All the involved chemicals must be already registered - use add_chemical() if needed.

        NOTE: in the reactants and products, if the stoichiometry coefficients aren't specified,
              they're assumed to be 1.
              The reaction orders, if not specified, are assumed to be equal to their corresponding
              stoichiometry coefficients.

              The full structure of each term in the list of reactants and of products
              is the triplet:  (stoichiometry coefficient, name, reaction order)

              EXAMPLES of formats to use for each term in the lists of the reactants and of the products:
                "F"         is taken to mean (1, "F", 1) - default stoichiometry and reaction order
                (2, "F")    is taken to mean (2, "F", 2) - stoichiometry coefficient used as default for reaction order
                (2, "F", 1) means stoichiometry coefficient 2 and reaction order 1 - no defaults invoked
              It's equally acceptable to use LISTS in lieu of tuples for the pair or triplets

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

        :return:                Integer index of the newly-added reaction
                                 (in the object variable self.reaction_list)
        """
        rxn = Reaction(self, reactants, products, forward_rate, reverse_rate,
                       delta_H, delta_S, delta_G)
        self.reaction_list.append(rxn)

        involved_chemicals = rxn.extract_chemicals_in_reaction(exclude_enzyme=True)

        self.active_chemicals = self.active_chemicals.union(involved_chemicals)     # Union of sets
        if rxn.enzyme is not None:
            self.active_enzymes.add(rxn.enzyme)       # Add the new entry to a set

        return len(self.reaction_list) - 1



    def clear_reactions_data(self) -> None:
        """
        Get rid of all reactions; start again with "an empty slate" (but still with reference
        to the same data object about the chemicals)

        :return:    None
        """
        self.reaction_list = []



    #####################################################################################################

    '''                             ~   TO DESCRIBE THE REACTIONS  ~                                  '''

    def ________DESCRIBE_RXNS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def describe_reactions(self, concise=False) -> None:
        """
        Print out a user-friendly plain-text form of ALL the reactions.
        If wanting to describe just 1 reaction, use single_reaction_describe()

        EXAMPLE (not concise):
            Number of reactions: 2 (at temp. 25 C)
            (0) CH4 + 2 O2 <-> CO2 + 2 H2O  (kF = 3.0 / kR = 2.0 / Delta_G = -1,005.13 / K = 1.5) | 1st order in all reactants & products
            (1) A + B <-> C  (kF = 5.0 / kR = 1.0 / Delta_G =  / K = 5.0) | 1st order in all reactants & products
            Set of chemicals involved in the above reactions: {'CH4', 'O2', 'H2O', 'A', 'B', 'C'}

        :param concise:     If True, less detail is shown
        :return:            None
        """
        print(f"Number of reactions: {self.number_of_reactions()} (at temp. {self.temp - 273.15:,.4g} C)")
        for description in self.multiple_reactions_describe(concise=concise):
            print(description)

        if self.active_enzymes == set():    # If no enzymes were involved in any reaction
            print(f"Set of chemicals involved in the above reactions: "
                  f"{self.names_of_active_chemicals()}")
        else:
            print(f"Set of chemicals involved in the above reactions (not counting enzymes): "
                  f"{self.names_of_active_chemicals()}")
            print(f"Set of enzymes involved in the above reactions: "
                  f"{self.names_of_enzymes()}")



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

        return rxn.describe(concise)



    def names_of_active_chemicals(self) -> Set[str]:
        """
        Return a set of the names of all the chemicals involved in ANY of the registered reactions,
        NOT counting catalysts
        """
        name_set = set()
        for ac_index in self.active_chemicals:
            name_set.add(self.get_name(ac_index))

        return name_set



    def names_of_enzymes(self) -> Set[str]:
        """
        Return a set of the names of all the enzymes involved in ANY of reactions
        """
        name_set = set()
        for e_index in self.active_enzymes:
            name_set.add(self.get_name(e_index))

        return name_set



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

        EXAMPLE of the graph structure part of the returned object for an  A <-> B reaction:
           [{'id': 'C-0', 'labels': ['Chemical'], 'name': 'A', 'diff_rate': None},
            {'id': 'C-1', 'labels': ['Chemical'], 'name': 'B', 'diff_rate': None},

            {'id': 'R-0', 'labels': ['Reaction'], 'name': 'RXN', 'kF': 3.0, 'kR': 2.0, 'K': 1.5, 'Delta_G': -1005.13},

            {'id': 'edge-1', 'name': 'produces', 'source': 'R-0', 'target': 'C-1', 'stoich': 1, 'rxn_order': 1},
            {'id': 'edge-2', 'name': 'reacts',   'source': 'C-0', 'target': 'R-0', 'stoich': 1, 'rxn_order': 1}
           ]

        :return:    A dictionary with 3 keys: 'structure', 'color_mapping', 'caption_mapping'
        """
        graph = PyGraphVisual()

        # Note: the graph nodes representing Chemicals will be given an id such as "C-123" and a label "Chemical";
        #       the graph nodes representing Reactions will be given an id such as "R-456" and a label "Reaction"

        for i, rxn in enumerate(self.reaction_list):    # Consider each REACTION in turn
            # Add a node representing the reaction
            rxn_id = f"R-{i}"               # Example: "R-456"
            node_data = {'name': 'RXN'}

            rxn_properties = rxn.extract_rxn_properties()
            for k,v in rxn_properties.items():
                node_data[k] = f"{v:,.6g}"

            graph.add_node(node_id=rxn_id, labels='Reaction', data=node_data)


            # Process all the PRODUCTS of this reaction
            products = rxn.extract_products()
            for term in products:
                species_index = term[1]
                chemical_id = f"C-{species_index}"      # Example: "C-12"
                # Add each product to the graph as a node (if not already present)
                graph.add_node( node_id=chemical_id, labels="Chemical",
                                data={'name': self.get_name(species_index),
                                      'diff_rate': self.get_diffusion_rate(species_index)
                                      })

                # Append edge from "reaction node" to "product node"
                graph.add_edge(from_node=rxn_id, to_node=chemical_id, name="produces",
                               data={'stoich': rxn.extract_stoichiometry(term),
                                     'rxn_order': rxn.extract_rxn_order(term)
                                     })


            # Process all the REACTANTS of this reaction
            reactants = rxn.extract_reactants()
            for term in reactants:
                species_index = term[1]
                chemical_id = f"C-{species_index}"      # Example: "C-34"
                # Add each reactant to the graph as a node (if not already present)
                graph.add_node(node_id=chemical_id, labels="Chemical",
                               data={'name': self.get_name(species_index),
                                     'diff_rate': self.get_diffusion_rate(species_index)
                                     })

                # Append edge from "reactant node" to "reaction node"
                graph.add_edge(from_node=chemical_id, to_node=rxn_id, name="reacts",
                               data={'stoich': rxn.extract_stoichiometry(term),
                                     'rxn_order': rxn.extract_rxn_order(term)
                                     })


        graph.assign_color_mapping(label='Chemical', color='graph_green')
        graph.assign_color_mapping(label='Reaction', color='graph_lightbrown')

        graph.assign_caption(label='Chemical', caption='name')
        graph.assign_caption(label='Reaction', caption='name')

        #print(graph)

        return graph.serialize()



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
        # Send a brief summary of all the reactions to the HTML log file
        log.write("List of reactions:", style=log.h3, newline=False, also_print=False)

        rxn_descriptions = self.multiple_reactions_describe(concise=True)
        for desc in rxn_descriptions:
            log.write(desc, indent=4, also_print=False)

        #log.blank_line()

        graph_data = self.prepare_graph_network()
        # A dictionary with 3 keys: 'structure', 'color_mapping' and 'caption_mapping'

        # Send a plot of the network of reactions to the HTML log file
        GraphicLog.export_plot(graph_data, graphic_component, unpack=unpack)





###############################################################################################################
###############################################################################################################


class ChemicalAffinity(NamedTuple):
    """
    Used for binding of ligands to macromolecules (e.g. Transcription Factors to DNA)
    """
    chemical: str   # Name of ligand
    Kd: float       # Dissociation constant; inversely related to binding affinity
                    # Note: dissociation constants are for now assumed to be constant,
                    #       regardless of what other (nearby) sites are occupied by ligands



class Macromolecules(AllReactions):
    """
    Extends its parent class to manage modeling of large molecules (such as DNA)
    with multiple binding sites (for example, for Transcription Factors)

    End users will generally utilize the class ChemData, which extends this one
    """

    def __init__(self):

        super().__init__()          # Invoke the constructor of its parent class


        self.macromolecules = []    # List of names.  EXAMPLE: ["M1", "M2"]
                                    # The position in the list is referred to as the
                                    #   "index" of that macro-molecule
                                    # Names are enforced to be unique

        self.binding_sites = {}     # A dict whose keys are macromolecule names.
                                    # The values are in turn dicts, indexed by binding-site number.
        # EXAMPLE:
        #       {"M1": {1: ChemicalAffinity("A", 2.4), 2: ChemicalAffinity("C", 5.1)},
        #        "M2": {1: ChemicalAffinity("C", 9.1), 2: ChemicalAffinity("B", 0.3),
        #               3: ChemicalAffinity("A", 1.8), 4: ChemicalAffinity("C", 2.3)}
        #        }
        #       where "M1", "M2" are macro-molecules,
        #           and "A", "B", "C" are bulk chemicals (such as transcription factors),
        #           all previously-declared;
        #           the various ChemicalAffinity's are NamedTuples (objects)
        #           storing a ligand name and its dissociation constant at that site.

        # TODO: maybe make a new class for a SINGLE macromolecule (akin to what done for reactions)

        # Info on Binding Site Affinities : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6787930/




    def add_macromolecules(self, names: Union[str, List[str]]) -> None:
        """
        Register one or more macromolecule species, specified by their name(s)
        Note: this is a register of names, NOT of dynamical information
              about counts of macromolecules in the system (which is the domain of the class ReactionDynamics)

        :param names:   A string, or list of strings, with the name(s) of the macromolecule(s)
        :return:        None.  The object attribute self.macro_molecules will get modified
        """
        if type(names) == str:  # If a single name was given, turn it into a list
            names = [names]

        for m in names:
            if m in self.macromolecules:
                # Warn of redundant attempt to re-register an existing macromolecule name
                print(f"WARNING: Macromolecule `{m}` was already registered.  Skipped...")
            else:
                self.macromolecules.append(m)  # Grow the list of registered macromolecules



    def get_macromolecules(self) -> [str]:
        """
        Return a list of the names of all the registered macromolecules

        :return:    A (possibly empty) list of the names of all the registered macromolecules
        """
        return self.macromolecules



    def set_binding_site_affinity(self, macromolecule: str, site_number: int, ligand: str, Kd) -> None:
        """
        Set the values of the binding affinity of the given macromolecule, at the indicated site on it,
        for the specified chemical species.

        Any previously-set value (for that macromolecule/site_number/chemical) will get over-written.

        IMPORTANT: Only 1 chemical (ligand type) can be associated to a particular site of a given macromolecule; attempting
        to associate another one will result in error.  In case multiple ligands can bind to the same physical site on
        the macromolecule, simply assign multiple site numbers for each of them

        NOTE: at present, no allowance is made for the fact that if the macromolecule is already bound
              to chemical "X" then its affinity of "Y" might be different than in the absence of "X"

        :param macromolecule:   Name of a macromolecule; if not previously-declared,
                                    it will get added to the list of registered macromolecules
        :param site_number:     Unique integer to identify a binding site on the macromolecule
        :param ligand:          Name of a previously-declared (bulk) chemical;
                                    if not found, an Exception will be raised       TODO: inconsistent with "macromolecule" arg
        :param Kd:              Dissociation constant, in units of concentration (typically microMolar).
                                    Note that the dissociation constant is inversely proportional to the binding affinity
        :return:                None
        """
        assert ligand in self.get_all_names(), \
            f"set_binding_site_affinity(): no chemical named `{ligand}` found; use add_chemical() first"

        assert type(site_number) == int, \
            f"set_binding_site_affinity(): the argument `site_number` must be an integer"

        if macromolecule not in self.macromolecules:
            self.add_macromolecules([macromolecule])

        if self.binding_sites.get(macromolecule) is None:
            self.binding_sites[macromolecule] = {}          # Initialize the dict of binding sites for this macromolecule

        binding_data = self.binding_sites[macromolecule]    # This will be a dict whose key is the site_number

        if site_number in binding_data:
            existing_affinity_data = binding_data[site_number]
            if existing_affinity_data.chemical != ligand:
                raise Exception(f"set_binding_site_affinity(): "
                                f"site number {site_number} of macromolecule `{macromolecule}` was previously associated to chemical `{existing_affinity_data.chemical}` "
                                f"(attempting to set an affinity value for chemical `{ligand}`)")

        binding_data[site_number] = ChemicalAffinity(chemical=ligand, Kd=Kd)



    def get_binding_site_affinity(self, macromolecule: str, site_number: int) -> ChemicalAffinity:
        """
        Return the value of the binding affinity of the given macromolecule.
        If no value was previously set, an Exception is raised

        :param macromolecule:   Name of a macromolecule; if not found, an Exception will get raised
        :param site_number:     Integer to identify a binding site on the macromolecule
        :return:                The NamedTuple (ligand name, dissociation constant)
                                    if no value was previously set, an Exception is raised
        """
        assert macromolecule in self.macromolecules, \
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
                                                                    #   It will be None if no binding sites are found
        if binding_site_info is None:
            assert macromolecule in self.get_macromolecules(), \
                f"get_binding_sites(): the requested macromolecule ({macromolecule}) isn't registered; use add_macromolecules()"
            return {}

        d = {}
        for (site_number, affinity_obj) in binding_site_info.items():
            ligand = affinity_obj.chemical
            d[site_number] = ligand

        return d        # EXAMPLE: {1: "A", 2: "C"}



    def get_ligand_name(self, macromolecule: str, site_number: int) -> str:
        """
        Return the name of the ligand associated to the specified site
        on the given macromolecule.
        If not found, an Exception is raised

        :param macromolecule:   The name of a macromolecule
        :param site_number:     Integer to identify a binding site on the macromolecule
        :return:                The name of the ligand (chemical species)
        """
        binding_site_info = self.binding_sites.get(macromolecule)   # EXAMPLE: {1: ChemicalAffinity("A", 2.4), 2: ChemicalAffinity("C", 5.1)}
                                                                    #   It will be None if no binding sites are found
        if binding_site_info is None:
            assert macromolecule in self.get_macromolecules(), \
                f"get_ligand_name(): the requested macromolecule (`{macromolecule}`) isn't registered; use add_macromolecules()"

            raise Exception(f"get_ligand_name(): no binding sites are defined on macromolecule {macromolecule}")

        ligand_data = binding_site_info.get(site_number)        # EXAMPLE:  ChemicalAffinity("A", 2.4)
                                                                # It will be None if this binding site isn't present
        assert ligand_data is not None, \
            f"get_ligand_name(): binding site {site_number} is not found on macromolecule {macromolecule}"

        return  ligand_data.chemical



    def show_binding_affinities(self) -> None:
        """
        Print out the Dissociation Constant for each Binding Site in each Macromolecule

        :return:    None
        """
        for mm in self.get_macromolecules():
            #binding_sites_and_ligands = self.get_binding_sites_and_ligands(mm)
            print(mm, " :")
            for site_number in self.get_binding_sites(mm):
                aff = self.get_binding_site_affinity(macromolecule=mm, site_number=site_number)
                print(f"   Site {site_number} - Kd (dissociation const) for {aff.chemical} : {aff.Kd}")



    def reset_macromolecule(self, macromolecule) -> None:
        """
        Erase all data for the specified macromolecule

        :param macromolecule:   The name of a macromolecule
        :return:                None
        """
        if self.binding_sites.get(macromolecule) is not None:
            del self.binding_sites[macromolecule]



    def clear_macromolecules(self) -> None:
        """
        Reset all macromolecules to their original state

        :return:    None
        """
        self.macromolecules = []
        self.binding_sites = {}




###############################################################################################################
###############################################################################################################


class ChemData(Macromolecules):
    """
    Data about all the chemicals and (if applicable) reactions,
    including:
        - names
        - diffusion rates
        - macro-molecules Binding Site Affinities (for Transcription Factors)
        - reaction data (see also class "Reaction", in "reaction_data.py")


    Notes:  - for now, the temperature is assumed constant everywhere, and unvarying (or very slowly varying)

            - we're using a "daisy chain" of classes extending the previous one, starting from ChemCore
              and ending in this user-facing class:
                    ChemCore <- Diffusion <- AllReactions <- Macromolecules <- ChemData
    """
    def __init__(self, names=None, diffusion_rates=None):
        """
        If chemical names and their diffusion rates are both provided, they must have the same count,
        and appear in the same order.
        It's ok to avoid passing any data at instantiation, and later add it.
        Reactions, if applicable, need to be added later by means of calls to add_reaction()
        Macro-molecules, if applicable, need to be added later

        :param names:           [OPTIONAL] A single name, or list or tuple of names, of the chemicals
        :param diffusion_rates: [OPTIONAL] A list or tuple with the diffusion rates of the chemicals
                                           If diffusion rates are provided, but no names given, the names will be
                                           auto-assigned as "Chemical 1", "Chemical 2", ...
        TODO: allow a way to optionally pass macromolecules as well
        """
        super().__init__()       # Invoke the constructor of its parent class

        # Initialize with the passed data, if provided
        if (names is not None) or (diffusion_rates is not None):
            self.init_chemical_data(names, diffusion_rates)




    def init_chemical_data(self, names=None, diffusion_rates=None)  -> None:
        """
        Initialize the names (if provided) and diffusion rates (if provided)
        of all the chemical species, in the given order.
        If no names are provided, the strings "Chemical 1", "Chemical 2", ..., are used

        IMPORTANT: this function can be invoked only once, before any chemical data is set.
                   To add new chemicals later, use add_chemical()

        :param names:           [OPTIONAL] A single name, or list or tuple of names, of the chemicals
        :param diffusion_rates: [OPTIONAL] A list or tuple with the diffusion rates of the chemicals,
                                           in the same order as the names
        :return:                None
        """
        n_species = self.number_of_chemicals()

        # Validate
        assert n_species == 0, \
            f"ChemData.init_chemical_data(): this function can be invoked only once, before any chemical data is set"


        if names:
            arg_type = type(names)
            if arg_type == str:
                names = [names]     # Turn a single name into a list
            else:
                assert arg_type == list or arg_type == tuple, \
                    f"ChemData.init_chemical_data(): the `names` argument must be a list or tuple, or a string for a single name.  " \
                    f"What was passed was of type {arg_type}"

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
                    self.add_chemical(assigned_name)
            else:
                for i, diff in enumerate(diffusion_rates):
                    assigned_name = f"Chemical {i+1}"
                    self.assert_valid_diffusion(diff)
                    self.add_chemical(assigned_name)
                    self.set_diffusion_rate(name=assigned_name, diff_rate=diff)

        else:   # names is not None
            for i, chem_name in enumerate(names):
                assert type(chem_name) == str, \
                    f"ChemData.init_chemical_data(): all the names must be strings.  The passed value ({chem_name}) is of type {type(chem_name)}"
                if diffusion_rates is None:
                    self.add_chemical(chem_name)
                else:
                    diff = diffusion_rates[i]
                    self.assert_valid_diffusion(diff)
                    self.add_chemical(chem_name)
                    self.set_diffusion_rate(name=chem_name, diff_rate=diff)
