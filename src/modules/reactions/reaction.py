from typing import Union
#from reaction_data import ReactionData
#from src.modules.reactions.reaction_data import ReactionData
import math
import numpy as np


class Reaction:
    """
    Data about a SINGLE reaction,
    including:
        - stoichiometry
        - kinetic data (reaction rates, reaction orders)
        - thermodynamic data (temperature, changes in enthalpy/entropy/Gibbs Free Energy)
        - list of involved enzymes


    (Note: this will eventually be stored in a Neo4j graph database)

    Each reaction contains:
            "reactants"
            "products"
            "kF"    (forward reaction rate constant)
            "kR"    (reverse reaction rate constant)
            "K"     (equilibrium constant - from either kinetic or thermodynamic data; if both present, they must match up!)
            "Delta_H" (change in Enthalpy: Enthalpy of Products - Enthalpy of Reactants)
            "Delta_S" (change in Entropy)
            "Delta_G" (change in Gibbs Free Energy)
                        Note - at constant temperature T :  Delta_G = Delta_H - T * Delta_S
                        Equilibrium constant = exp(-Delta_G / RT)
            "enzymes" (list of the indices of the chemical species that appear as catalysts in the reaction)

        Each Reactant and each Product is a triplet of the form: (stoichiometry, species index, reaction order).
        The "reaction order" refers to the forward reaction for reactants, and the reverse reaction for products.
    """

    def __init__(self, rxn_data):
        """
        """
        self.rxn_data = rxn_data
        
        self.reactants = None
        self.products = None
        self.kF = None
        self.kR = None
        self.Delta_H = None
        self.Delta_S = None
        self.Delta_G = None
        self.K = None
        self.enzymes = None

        self.R = 8.314462           # Ideal gas constant, in units of Joules / (Kelvin * Mole)





#####################################################################################################

    '''                       ~   TO READ DATA STRUCTURE of the REACTION  ~                           '''

    def ________READ_RXN_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def extract_reactants(self) -> [(int, int, int)]:
        """
        Return a list of triplets with details of the reactants of the given reaction

        :return:    A list of triplets of the form (stoichiometry, species index, reaction order)
        """
        return self.reactants



    def extract_reactants_formula(self) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula

        :return:
        """
        reactants = self.extract_reactants()
        return self._standard_form_chem_eqn(reactants)



    def extract_products(self) -> [(int, int, int)]:
        """
        Return a list of triplet with details of the products of the given reaction

        :return:    A list of triplets of the form (stoichiometry, species index, reaction order)
        """
        return self.products



    def extract_products_formula(self) -> str:
        """
        Return a string with a user-friendly form of the right (products) side of the reaction formula

        :return:
        """
        products = self.extract_products()
        return self._standard_form_chem_eqn(products)



    def extract_forward_rate(self) -> float:
        """

        :return:    The value of the forward rate constant for the above reaction
        """
        return self.kF


    def extract_reverse_rate(self) -> float:
        """

        :return:    The value of the reverse (back) rate constant for the above reaction
        """
        return self.kR



    def extract_stoichiometry(self, term: (int, int, int)) -> int:
        """
        Return the stoichiometry coefficient, from a reaction term

        :param term:
        :return:
        """
        return term[0]

    def extract_species_index(self, term: (int, int, int)) -> int:
        """
        Return the index of the chemical species, from a reaction term

        :param term:
        :return:
        """
        return term[1]

    def extract_rxn_order(self, term: (int, int, int)) -> int:
        """
        Return the reaction order, from a reaction term

        :param term:
        :return:
        """
        return term[2]



    def extract_rxn_properties(self) -> {}:
        """
        Create a dictionary with the numerical properties of the given reaction
        (skipping any lists or None values)
        For example, reaction rates, Delta G, equilibrium constant

        :return:    EXAMPLE: {'kF': 3.0, 'kR': 2.0, 'Delta_G': -1005.1305052750387, 'K': 1.5}
        """
        properties = {}

        if self.kF is not None:
            properties['kF'] = self.kF

        if self.kR is not None:
            properties['kR'] = self.kR

        if self.Delta_G is not None:
            properties['Delta_G'] = self.Delta_G

        if self.Delta_G is not None:
            properties['Delta_G'] = self.Delta_G

        if self.K is not None:
            properties['K'] = self.K

        return properties



    def extract_chemicals_in_reaction(self) -> {int}:
        """
        Return a SET of indices (being a set, it's NOT in any particular order)
        of all the chemicals in the specified reaction

        :return:            A SET of indices of the chemicals involved in the above reaction
                                Note: being a set, it's NOT in any particular order
        """
        chem_set = set()    # Running set being built

        reactants = self.extract_reactants()
        products = self.extract_products()

        for r in reactants:
            species_index = self.extract_species_index(r)
            chem_set.add(species_index)

        for p in products:
            species_index = self.extract_species_index(p)
            chem_set.add(species_index)

        return chem_set




    #####################################################################################################

    '''                                  ~   TO SET DATA  ~                                           '''

    def ________SET_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def define_reaction(self, reactants: Union[int, str, tuple, list], products: Union[int, str, tuple, list],
                        forward_rate=None, reverse_rate=None,
                        Delta_H=None, Delta_S=None, Delta_G=None) -> None:
        """
        TODO: how to handle if called multiple times    (maybe always do a reset?  Or let this be part of the instantiation?)

        Add the parameters of a SINGLE reaction, optionally including kinetic and/or thermodynamic data.
        The involved chemicals must be already registered - use add_chemical() if needed.

        NOTE: in the reactants and products, if the stoichiometry and/or reaction order aren't specified,
              they're assumed to be 1.
              Their full structure is the triplet (stoichiometry coefficient, name, reaction order)

        EXAMPLES of formats for reactants and products (*assuming* that the chemical species with index 5 is called "F"):
                    "F"         gets turned into:   (1, 5, 1)
                    (3, "F")                        (3, 5, 1)
                    (3, "F", 2)                     (3, 5, 2)
                    It's equally acceptable to use LISTS in lieu of tuples

        :param reactants:       A list of triplets (stoichiometry, species name or index, reaction order),
                                    or simplified terms in various formats; for details, see above
        :param products:        A list of triplets (stoichiometry, species name or index, reaction order of REVERSE reaction),
                                    or simplified terms in various formats; for details, see above
        :param forward_rate:    [OPTIONAL] Forward reaction rate constant
        :param reverse_rate:    [OPTIONAL] Reverse reaction rate constant
        :param Delta_H:         [OPTIONAL] Change in Enthalpy (from reactants to products)
        :param Delta_S:         [OPTIONAL] Change in Entropy (from reactants to products)
        :param Delta_G:         [OPTIONAL] Change in Free Energy (from reactants to products)
        :return:                None
        """
        reactant_list = [self._parse_reaction_term(r, "reactant") for r in reactants]   # A list of triples of integers
        product_list = [self._parse_reaction_term(r, "product") for r in products]      # A list of triples of integers

        # TODO: use a more sophisticated approach to catch identical reaction sides even if
        #       terms are reshuffled
        assert reactant_list != product_list, \
            f"Reaction.define_reaction(): the reactants and the products can't all be identical! " \
            f"Internal structure: {reactant_list}"

        # Locate any enzymes (catalysts) that might be present
        enzyme_list = []
        for reactant in reactant_list:
            if reactant in product_list:
                enzyme_list.append(self.extract_species_index(reactant))

        number_enzymes = len(enzyme_list)
        if number_enzymes == len(reactant_list) or number_enzymes == len(product_list):
            raise Exception(f"Reaction.define_reaction(): all the terms in the reaction appear to be enzymes!  "
                            f"Enzymes: {[self.rxn_data.get_name(e) for e in enzyme_list]}")


        self.reactants = reactant_list
        self.products = product_list
        self.enzymes = enzyme_list


        # Process kinetic data, if available
        if forward_rate:
            self.kF = forward_rate

        if reverse_rate:
            self.kR = reverse_rate
            if forward_rate:
                # If all the kinetic data is available...
                equil_const = forward_rate / reverse_rate    # ...compute the equilibrium constant (from kinetic data)
                self.K = equil_const
                if self.rxn_data.temp:
                    self.Delta_G = - self.R * self.rxn_data.temp * math.log(equil_const)   # the change in Gibbs Free Energy


        # Process thermodynamic data, if available
        if Delta_H is not None:
            self.Delta_H = Delta_H

        if Delta_S is not None:
            self.Delta_S = Delta_S
            if (self.rxn_data.temp is not None) and (Delta_H is not None):
                # If all the thermodynamic data is available...
                Delta_G = Delta_H - self.rxn_data.temp * Delta_S                  # ...compute the change in Gibbs Free Energy

                if self.Delta_G is not None:      # If already set from kinetic data, make sure that the two match!
                    assert np.allclose(Delta_G, self.Delta_G), \
                        f"define_reaction(): Kinetic data (leading to Delta_G {self.Delta_G}) " \
                        f"is inconsistent with thermodynamic data (leading to Delta_G {Delta_G})"
                else:
                    # The kinetic data was incomplete; fill in the missing parts from the thermodynamic data
                    self.Delta_G = Delta_G
                    equil_const = math.exp(- Delta_G / (self.R * self.rxn_data.temp))    # Compute the equilibrium constant (from the thermodynamic data)
                    self.K = equil_const
                    # If only one of the Forward or Reverse rates was provided, compute the other one
                    if (self.kF is None) and (self.kR is not None):
                        self.kF = equil_const * self.kR
                    if (self.kR is None) and (self.kF is not None):
                        self.kR = self.kF / equil_const

        # TODO: add more combinations of arguments supplied
        elif (Delta_G is not None) and (self.rxn_data.temp is not None):
            # The kinetic data was incomplete; fill in the missing parts from the thermodynamic data
            self.Delta_G = Delta_G
            equil_const = math.exp(- Delta_G / (self.R * self.rxn_data.temp))    # Compute the equilibrium constant (from the thermodynamic data)
            self.K = equil_const
            # If only one of the Forward or Reverse rates was provided, compute the other one
            if (self.kF is None) and (self.kR is not None):
                self.kF = equil_const * self.kR
            if (self.kR is None) and (self.kF is not None):
                self.kR = self.kF / equil_const





    #####################################################################################################

    '''                                ~   TO DESCRIBE THE DATA  ~                                    '''

    def ________DESCRIBE_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def describe(self, rxn_index=0, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction
        EXAMPLE (concise):      "CH4 + 2 O2 <-> CO2 + 2 H2O"
        EXAMPLE (not concise):  "(0) CH4 + 2 O2 <-> CO2 + 2 H2O  (kF = 3.0 / kR = 2.0 / Delta_G = -1,005.13 / K = 1.5) | 1st order in all reactants & products"

        :param rxn_index:   Integer associated to this reaction.  Only used if concise is False
        :param concise:     If True, less detail is shown
        :return:            A string with a description of the specified reaction
        """
        reactants = self.extract_reactants()
        products = self.extract_products()

        left = self._standard_form_chem_eqn(reactants)       # Left side of the equation, as a user-friendly string
        right = self._standard_form_chem_eqn(products)       # Right side of the equation

        if concise:
            return f"{left} <-> {right}"        # Concise start point for a description of the reaction

        # If we get this far, we're looking for a more detailed description

        rxn_description = f"{rxn_index}: {left} <-> {right}"        # Start point for a description of the reaction
        details = []
        rxn_properties = self.extract_rxn_properties()
        for k,v in rxn_properties.items():
            details.append(f"{k} = {v:,.6g}")          # EXAMPLE: "kF = 3"

        rxn_description += "  (" + ' / '.join(details) + ")"    # EXAMPLE: "  (kF = 3 / kR = 2 / Delta_G = -1,005.13)"

        high_order = False      # Is there any term whose order is greater than 1 ?
        for r in reactants:
            if r[2] > 1:
                rxn_description += f" | {r[2]}-th order in reactant {self.rxn_data.get_name(r[1])}"
                high_order = True
        for p in products:
            if p[2] > 1:
                rxn_description += f" | {p[2]}-th order in product {self.rxn_data.get_name(p[1])}"
                high_order = True

        if not high_order:
            rxn_description += " | 1st order in all reactants & products"

        return rxn_description



    #####################################################################################################

    '''                                    ~   PRIVATE  ~                                             '''

    def ________PRIVATE________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def _standard_form_chem_eqn(self, eqn_side: list) -> str:
        """
        Return a user-friendly form of a side of a chemical equation.

        EXAMPLE:  turn [(1, 0, 1), (2, 8, 1)]  into  "Fe + 2 Cl"  (if species 0 is named "Fe" and species 8 is "Cl")

        Note: the reaction order is NOT used

        :param eqn_side:    A list encoding either side of a chemical equation
        :return:            A string with a user-friendly form of a side of a chemical equation
        """
        assert type(eqn_side) == list, \
            f"Reaction._standard_form_chem_eqn(): the argument must be a list (it was of type {type(eqn_side)})"

        formula_list = []
        for t in eqn_side:
            stoichiometry = self.extract_stoichiometry(t)
            species_index = self.extract_species_index(t)

            if stoichiometry == 1:
                term = f"{self.rxn_data.get_name(species_index)}"
            else:
                term = f"{stoichiometry} {self.rxn_data.get_name(species_index)}"

            formula_list.append(term)

        return " + ".join(formula_list)


    def _parse_reaction_term(self, term: Union[int, str, tuple, list], name="term") -> (int, int, int):
        """
        Accept various ways to specify a reaction term, and return a standardized tuple form of it.
        In the tuples or lists:
            - optional 1st entry is the stoichiometry
            - required entry is the chemical name
            - optional 3rd one is the reaction order

        EXAMPLES (*assuming* that the chemical species with index 5 is called "F"):
            "F"          gets turned into:  (1, 5, 1)
            (3, "F")                        (3, 5, 1)
            (3, "F", 2)                     (3, 5, 2)
            It's equally acceptable to use LISTS in lieu of tuples

        :param term:    A string (a chemical name)
                            OR  a pair (stoichiometry coeff, name)
                            OR  a triplet (stoichiometry coeff, name, reaction order)
        :param name:    An optional nickname to refer to this term in error messages, if applicable
                            (for example, "reactant" or "product")
        :return:        A standardized tuple form, of the form (stoichiometry, species index, reaction_order),
                            where all terms are integers
        """
        if type(term) == str:
            return  (1, self.rxn_data.get_index(term), 1)    # Accept simply the chemical name as a shortcut
                                                    # for when the stoichiometry coefficient and reaction order are both 1
        elif type(term) != tuple and type(term) != list:
            raise Exception(f"_parse_reaction_term(): {name} item must be either a string (a chemical name), "
                            f"or a pair (stoichiometry coeff, name) or a triplet (stoichiometry coeff, name, reaction order). "
                            f"Instead, it is `{term}` (of type {type(term)})")

        # If we get thus far, term is either a tuple or a list
        if len(term) != 3 and len(term) != 2:
            raise Exception(f"_parse_reaction_term(): Unexpected length for {name} tuple/list: it should be 2 or 3. Instead, it is {len(term)}")

        stoichiometry = term[0]
        species = term[1]
        if type(species) == str:
            species = self.rxn_data.get_index(species)
        else:
            raise Exception(f"_parse_reaction_term(): The chemical name must be a string. Instead, it is {species} (of type {type(species)})")

        if len(term) == 2:
            return (stoichiometry, species, 1)
        else:   # Length is 3
            reaction_order = term[2]
            return (stoichiometry, species, reaction_order)