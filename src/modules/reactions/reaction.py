from typing import Union, Set, Tuple
import numpy as np
from src.modules.reactions.thermodynamics import ThermoDynamics


class Reaction:
    """
    Data about a SINGLE reaction,
    including:
        - stoichiometry
        - kinetic data (reaction rates, reaction orders)
        - thermodynamic data (temperature, changes in enthalpy/entropy/Gibbs Free Energy)
        - list of involved enzymes


    (Note: this data will eventually be stored in a Neo4j graph database)

    Each reaction contains:
            "reactants"
            "products"
            "kF"    (forward reaction rate constant)
            "kR"    (reverse reaction rate constant)
            "enzyme" (the index of a chemical that catalyzes this reaction)
            "macro_enzyme" (The pair (macromolecule name, binding site number))
            plus thermodynamic data: delta_H, delta_S, delta_G, K (equilibrium constant) -
                                     for details see class "ThermoDynamics"

        Each Reactant and each Product is a triplet of the form: (stoichiometry, species index, reaction order).
        The "reaction order" in that triplet refers to the forward reaction for reactants, and the reverse reaction for products.
        Note that reactants and products might be catalysts
    """


    def __init__(self, chem_data, reactants: Union[int, str, list], products: Union[int, str, list],
                 forward_rate=None, reverse_rate=None,
                 delta_H=None, delta_S=None, delta_G=None):
        """
        Add the parameters of a SINGLE reaction, optionally including kinetic and/or thermodynamic data.
        The involved chemicals must be already registered - use add_chemical() if needed.

        NOTE: in the reactants and products, if the stoichiometry and/or reaction order aren't specified,
              they're assumed to be 1.
              Their full structure is the triplet (stoichiometry coefficient, name, reaction order)

        EXAMPLES of formats for each item in the lists of the reactants and products
        (*assuming* that the chemical species with index 5 is called "F"):
                    "F"         gets turned into:   (1, 5, 1)
                    (3, "F")                        (3, 5, 1)
                    (3, "F", 2)                     (3, 5, 2)
                    It's equally acceptable to use LISTS in lieu of tuples

        :param chem_data:       Object of type "ReactionData"
        :param reactants:       A list of triplets (stoichiometry, species name, reaction order),
                                    or simplified terms in various formats; for details, see above.
                                    If not a list, it will get turned into one
        :param products:        A list of triplets (stoichiometry, species name, reaction order of REVERSE reaction),
                                    or simplified terms in various formats; for details, see above.
                                    If not a list, it will get turned into one
        :param forward_rate:    [OPTIONAL] Forward reaction rate constant
        :param reverse_rate:    [OPTIONAL] Reverse reaction rate constant
        :param delta_H:         [OPTIONAL] Change in Enthalpy (from reactants to products)
        :param delta_S:         [OPTIONAL] Change in Entropy (from reactants to products)
        :param delta_G:         [OPTIONAL] Change in Free Energy (from reactants to products)

        """
        self.chem_data = chem_data
        
        self.reactants = None
        self.products = None
        self.kF = forward_rate
        self.kR = reverse_rate
        self.delta_H = delta_H
        self.delta_S = delta_S
        self.delta_G = delta_G
        self.K = None
        self.enzyme = None          # The INDEX of a chemical that catalyzes this reaction, if applicable
                                    #   Note: enzymes are automatically extracted from the reaction formula
        self.macro_enzyme = None    # The pair (macromolecule name, binding site number)
                                    #   EXAMPLE: ("M2", 5)          TODO: maybe turn into a data object

        assert reactants is not None, "Reaction(): the argument `reactants` is a required one; it can't be None"
        if type(reactants) != list:
            reactants = [reactants]


        assert products is not None, "Reaction(): the argument `products` is a required one; it can't be None"
        if type(products) != list:
            products = [products]


        reactant_list = [self._parse_reaction_term(r, "reactant") for r in reactants]   # A list of triples of integers
        product_list = [self._parse_reaction_term(r, "product") for r in products]      # A list of triples of integers

        # Catch identical reaction sides, even if terms are reshuffled
        assert set(reactant_list) != set(product_list), \
            f"Reaction(): the two sides of the reaction can't be identical! " \
            f"Same reactants and products: {self._standard_form_chem_eqn(reactant_list)}"


        # Locate any enzymes (catalysts) that might be present - though for now a warning is issued if more than 1
        enzyme_list = []
        for reactant in reactant_list:
            if reactant in product_list:
                enzyme_list.append(self.extract_species_index(reactant))

        number_enzymes = len(enzyme_list)

        if number_enzymes == len(reactant_list) or number_enzymes == len(product_list):
            raise Exception(f"Reaction(): all the terms in the reaction appear to be enzymes!  "
                            f"Enzymes: {[self.chem_data.get_name(e) for e in enzyme_list]}")


        self.reactants = reactant_list
        self.products = product_list

        if number_enzymes >= 1:
            self.enzyme = enzyme_list[0]    # In the irregular scenarios that there appear to be multiple enzymes, only one
                                            #   is considered, and a warning is printed out (the other apparent enzyme
                                            #   will be treated as any other reagent/product)
        if number_enzymes > 1:
            print(f"Reaction(): WARNING - the reaction appears to have multiple enzymes:"
                  f" {[self.chem_data.get_name(e) for e in enzyme_list]}")


        # Process the kinetic and thermodynamic data, and update various object attributes accordingly
        self._set_kinetic_and_thermodynamic(forward_rate=forward_rate, reverse_rate=reverse_rate,
                                            delta_H=delta_H, delta_S=delta_S, delta_G=delta_G, temp=self.chem_data.temp)



    def set_macro_enzyme(self, macromolecule: str, site_number: int) -> None:
        """
        Register that the given macromolecule, at the given site on it,
        catalyzes this reaction

        :param macromolecule:   Name of macromolecule acting as a catalyst
        :param site_number:     Integer to identify a binding site on the above macromolecule
        :return:                None
        """
        self.macro_enzyme = (macromolecule, site_number)




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
        EXAMPLE: "CH4 + 2 O2"

        :return:    A string with one side of a chemical reaction
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

        :return:    A string with one side of a chemical reaction
        """
        products = self.extract_products()
        return self._standard_form_chem_eqn(products)



    def extract_forward_rate(self) -> float:
        """

        :return:    The value of the forward rate constant for this reaction
        """
        return self.kF


    def extract_reverse_rate(self) -> float:
        """

        :return:    The value of the reverse (back) rate constant for this reaction
        """
        return self.kR



    def unpack_for_dynamics(self) -> tuple:
        """
        A convenient unpacking meant for dynamics simulations
        that need the reactants, the products, and the forward and reverse rate constants

        :return:    A 4-element tuple, containing:
                        (reactants , products , forward rate constant , reverse rate constant)
                        Note: both reactants products are lists of triplets
        """
        return (self.reactants, self.products, self.kF, self.kR)



    def extract_stoichiometry(self, term: (int, int, int)) -> int:
        """
        Return the stoichiometry coefficient, from a reaction term

        :param term:    A triplet of integers representing a reaction term
        :return:        An integer with the stoichiometry coefficient
        """
        return term[0]

    def extract_species_index(self, term: (int, int, int)) -> int:
        """
        Return the index of the chemical species, from a reaction term

        :param term:    A triplet of integers representing a reaction term
        :return:        An integer with the index of the chemical species in the term
        """
        return term[1]

    def extract_rxn_order(self, term: (int, int, int)) -> int:
        """
        Return the reaction order, from a reaction term

        :param term:    A triplet of integers representing a reaction term
        :return:        An integer with the reaction order for this term
        """
        return term[2]



    def extract_rxn_properties(self) -> {}:
        """
        Create a dictionary with the numerical properties of the given reaction
        (skipping any lists or None values)
        Possible values include:
            forward and reverse reaction rates, ΔH, ΔS, ΔG, K (equilibrium constant)

        :return:    EXAMPLE: {'kF': 3.0, 'kR': 2.0, 'delta_G': -1005.1305052750387, 'K': 1.5}
        """
        properties = {}

        if self.kF is not None:
            properties['kF'] = self.kF

        if self.kR is not None:
            properties['kR'] = self.kR

        if self.delta_H is not None:
            properties['delta_H'] = self.delta_H

        if self.delta_S is not None:
            properties['delta_S'] = self.delta_S

        if self.delta_G is not None:
            properties['delta_G'] = self.delta_G

        if self.K is not None:
            properties['K'] = self.K

        return properties



    def extract_chemicals_in_reaction(self, exclude_enzyme=False) -> Set[int]:
        """
        Return a SET of indices (being a set, it's NOT in any particular order)
        identifying all the chemicals appearing in this reaction.
        Optionally, exclude any that participate in a catalytic role
        (appearing identically on both sides of the reaction)

        :param exclude_enzyme:  If True, any enzyme, if present, won't be included
        :return:                A SET of indices of the chemicals involved in this reaction
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

        if exclude_enzyme:
            chem_set = chem_set - {self.enzyme}     # Difference between sets

        return chem_set




    #####################################################################################################

    '''                               ~   TO DESCRIBE THE DATA  ~                                     '''

    def ________DESCRIBE_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def describe(self, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction
        EXAMPLE (concise):      "CH4 + 2 O2 <-> CO2 + 2 H2O"
        EXAMPLE (not concise):  "CH4 + 2 O2 <-> CO2 + 2 H2O  (kF = 3.0 / kR = 2.0 / Delta_G = -1,005.13 / K = 1.5) | 1st order in all reactants & products"

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

        rxn_description = f"{left} <-> {right}"        # Start point for a description of the reaction
        details = []
        rxn_properties = self.extract_rxn_properties()
        for k,v in rxn_properties.items():
            details.append(f"{k} = {v:,.5g}")          # EXAMPLE: "kF = 3"

        rxn_description += "  (" + ' / '.join(details) + ")"    # EXAMPLE: "  (kF = 3 / kR = 2 / Delta_G = -1,005.13)"


        # If an ENZYME is involved, show it
        if self.enzyme is not None:
            rxn_description += f" | Enzyme: {self.chem_data.get_name(self.enzyme)}"

        if self.macro_enzyme is not None:
            rxn_description += f" | Macromolecule Enzyme: {self.macro_enzyme[0]}, at site # {self.macro_enzyme[1]}"


        # Show details about the ORDER of the reaction
        high_order = False      # Is there any term whose order is greater than 1 ?
        for r in reactants:
            if r[2] > 1:
                rxn_description += f" | {r[2]}-th order in reactant {self.chem_data.get_name(r[1])}"
                high_order = True
        for p in products:
            if p[2] > 1:
                rxn_description += f" | {p[2]}-th order in product {self.chem_data.get_name(p[1])}"
                high_order = True

        if not high_order:
            rxn_description += " | 1st order in all reactants & products"

        return rxn_description




    #####################################################################################################

    '''                                     ~   ANALYSIS  ~                                           '''

    def ________ANALYSIS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def reaction_quotient(self, conc, explain=False) -> Union[np.double, Tuple[np.double, str]]:
        """
        Compute the "Reaction Quotient" (aka "Mass–action Ratio"),
        given the concentrations of chemicals involved in this reaction

        :param conc:        Dictionary with the concentrations of the species involved in the reaction.
                            The keys are the chemical names
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
        :param explain:     If True, it also returns the math formula being used for the computation
                                EXAMPLES:   "([C][D]) / ([A][B])"
                                            "[B] / [A]^2"

        :return:            If explain is False, return value for the "Reaction Quotient" (aka "Mass–action Ratio");
                                if True, return a pair with that quotient and a string with the math formula that was used.
                                Note that the reaction quotient is a Numpy scalar that might be np.inf or np.nan
        """
        numerator = np.double(1)    # The product of all the concentrations of the reaction products (adjusted for reaction order)
        denominator = np.double(1)  # The product of all the concentrations of the reactants (also adjusted for reaction order)

        numerator_text = ""      # First part of the the textual explanation
        denominator_text = ""    # Second part of the the textual explanation


        # Computer the numerator of the "Reaction Quotient"
        for p in self.products:
            # Loop over the reaction products
            species_index =  self.extract_species_index(p)
            rxn_order =  self.extract_rxn_order(p)

            species_name =  self.chem_data.get_name(species_index)
            species_conc = conc.get(species_name)
            assert species_conc is not None, f"reaction_quotient(): unable to proceed because the " \
                                             f"concentration of `{species_name}` was not provided"

            numerator *= (species_conc ** rxn_order)
            if explain:
                numerator_text += f"[{species_name}]"
                if rxn_order > 1:
                    numerator_text += f"^{rxn_order} "

        if explain and len(self.products) > 1:
            numerator_text = f"({numerator_text})"  # In case of multiple terms, enclose them in parenthesis


        # Computer the denominator of the "Reaction Quotient"
        for r in self.reactants:
            # Loop over the reactants
            species_index =  self.extract_species_index(r)
            rxn_order =  self.extract_rxn_order(r)

            species_name = self.chem_data.get_name(species_index)
            species_conc = conc.get(species_name)
            assert species_conc is not None, f"reaction_quotient(): unable to proceed because the " \
                                             f"concentration of `{species_name}` was not provided"

            denominator *= (species_conc ** rxn_order)
            if explain:
                denominator_text += f"[{species_name}]"
                if rxn_order > 1:
                    denominator_text += f"^{rxn_order} "

        if explain and len(self.reactants) > 1:
            denominator_text = f"({denominator_text})"  # In case of multiple terms, enclose them in parenthesis


        quotient = numerator / denominator      # It might be np.inf or np.nan

        if explain:
            formula = f"{numerator_text} / {denominator_text}"
            return (quotient, formula)

        return quotient




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
                term = f"{self.chem_data.get_name(species_index)}"
            else:
                term = f"{stoichiometry} {self.chem_data.get_name(species_index)}"

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
            return  (1, self.chem_data.get_index(term), 1)  # Accept simply the chemical name as a shortcut
                                                            # for when the stoichiometry coefficient and reaction order are both 1
        elif type(term) != tuple and type(term) != list:
            raise Exception(f"_parse_reaction_term(): {name} item must be either a string (a chemical name), "
                            f"or a pair (stoichiometry coeff, name) or a triplet (stoichiometry coeff, name, reaction order). "
                            f"Instead, it is `{term}` (of type {type(term)})")

        # If we get thus far, term is either a tuple or a list
        if len(term) != 3 and len(term) != 2:
            raise Exception(f"_parse_reaction_term(): Unexpected length for {name} tuple/list: it should be 2 or 3. Instead, it is {len(term)}")

        stoichiometry = term[0]
        assert type(stoichiometry) == int, \
            f"_parse_reaction_term(): The stoichiometry coefficient, if provided, must be an integer. Instead, it is {stoichiometry}"

        species = term[1]
        if type(species) == str:
            species = self.chem_data.get_index(species)
        else:
            raise Exception(f"_parse_reaction_term(): The chemical name must be a string. Instead, it is {species} (of type {type(species)})")

        if len(term) == 2:
            return (stoichiometry, species, 1)
        else:   # Length is 3
            reaction_order = term[2]
            return (stoichiometry, species, reaction_order)



    def _set_kinetic_and_thermodynamic(self, forward_rate, reverse_rate,
                                       delta_H, delta_S, delta_G, temp) -> None:
        """
        Set all the kinetic and thermodynamic data derivable - directly or indirectly - from the passed arguments,
        storing it in object attributes.
        Raise an Exception if any inconsistency is detected.

        :return:    None
        """
        self.kF = forward_rate
        self.kR = reverse_rate
        self.delta_H = delta_H
        self.delta_S = delta_S
        self.delta_G = delta_G


        # Process kinetic data, if available
        #       (extracting thermodynamic data when feasible)
        if (self.kF is not None) and (self.kR is not None):
            # If all the kinetic data is available...
            self.K = self.kF / self.kR    # ...compute the equilibrium constant (from kinetic data)

            if temp:
                # If the temperature is set, compute the change in Gibbs Free Energy
                delta_G_kinetic = ThermoDynamics.delta_G_from_K(K = self.K, temp = temp)
                if self.delta_G is None:
                    self.delta_G = delta_G_kinetic
                else:   # If already present (passed as argument), make sure that the two match!
                    assert np.allclose(delta_G_kinetic, self.delta_G), \
                        f"_set_kinetic_and_thermodynamic(): Kinetic data (leading to Delta_G={delta_G_kinetic}) " \
                        f"is inconsistent with the passed value of Delta_G={self.delta_G})"


        if (self.delta_H is not None) and (self.delta_S is not None) and (temp is not None):
            # If all the thermodynamic data (possibly except delta_G) is available...

            # Compute the change in Gibbs Free Energy from delta_H and delta_S, at the current temperature
            delta_G_thermo = ThermoDynamics.delta_G_from_enthalpy(delta_H = self.delta_H, delta_S = self.delta_S, temp = temp)

            if self.delta_G is None:
                self.delta_G = delta_G_thermo
            else:  # If already present (passed as argument or was set from kinetic data), make sure that the two match!
                if not np.allclose(delta_G_thermo, self.delta_G):
                    if delta_G is not None:
                        raise Exception(f"_set_kinetic_and_thermodynamic(): thermodynamic data (leading to Delta_G={delta_G_thermo}) "
                                        f"is inconsistent with the passed value of delta_G={delta_G})")
                    else:
                        raise Exception(f"_set_kinetic_and_thermodynamic(): thermodynamic data (leading to Delta_G={delta_G_thermo}) "
                                        f"is inconsistent with kinetic data (leading to Delta_G={self.delta_G})")


        if self.delta_G is not None:
            if (self.K is None) and (temp is not None):
                # If the temperature is known, compute the equilibrium constant (from the thermodynamic data)
                # Note: no need to do it if self.K is present, because we ALREADY handled that case
                self.K = ThermoDynamics.K_from_delta_G(delta_G = self.delta_G, temp = temp)

                # If only one of the Forward or Reverse rates was provided, compute the other one
                if (self.kF is None) and (self.kR is not None):
                    self.kF = self.K * self.kR
                if (self.kR is None) and (self.kF is not None):
                    self.kR = self.kF / self.K

            if temp is not None:
                # If either Enthalpy or Entropy is missing, but the other one is known, compute the missing one
                if (self.delta_H is None) and (self.delta_S is not None):
                    self.delta_H = ThermoDynamics.delta_H_from_gibbs(delta_G=self.delta_G, delta_S=self.delta_S, temp=temp)
                elif (self.delta_H is not None) and (self.delta_S is None):
                    self.delta_S = ThermoDynamics.delta_S_from_gibbs(delta_G=self.delta_G, delta_H=self.delta_H, temp=temp)
