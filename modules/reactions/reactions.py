

class Reactions:
    """
    Data about all applicable reactions,
    including stoichiometry, reaction rates and reaction orders
    """

    def __init__(self, chem_data):
        self.reaction_list = []     # List of dicts.  Each item represents a reaction, incl. its reverse
                                    # Reactions should be added by means of calls to add_reaction()

        self.chem_data = chem_data  # Object with info on the individual chemicals, incl. their names
        

    def number_of_reactions(self) -> int:
        # Return the number of registered reactions
        return len(self.reaction_list)


    def get_reaction(self, i) -> dict:
        return self.reaction_list[i]


    def get_reactants(self, i) -> (int, int, int):
        """

        :param i:
        :return:        A triplet (stoichiometry, species index, reaction order)
        """
        rxn = self.get_reaction(i)
        return rxn["reactants"]

    def get_reactants_formula(self, i):
        rxn = self.get_reaction(i)
        reactants = rxn["reactants"]
        return self.standard_form(reactants)


    def get_products(self, i) -> (int, int, int):
        """

        :param i:
        :return:        A triplet (stoichiometry, species index, reaction order)
        """
        rxn = self.get_reaction(i)
        return rxn["products"]

    def get_products_formula(self, i):
        rxn = self.get_reaction(i)
        products = rxn["products"]
        return self.standard_form(products)


    def get_forward_rate(self, i) -> float:
        rxn = self.get_reaction(i)
        return rxn["Rf"]

    def get_reverse_rate(self, i) -> float:
        rxn = self.get_reaction(i)
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
        assert type(reactants) == list, "ERROR: argument `reactants` must be a list"
        reactant_list = [self._parse_reaction_term(r, "reactants") for r in reactants]

        assert type(products) == list, "ERROR: argument `products` must be a list"
        product_list = [self._parse_reaction_term(r, "products") for r in products]

        rxn = {"reactants": reactant_list, "products": product_list, "Rf": forward_rate, "Rb": reverse_rate}
        self.reaction_list.append(rxn)



    def clear_reactions(self) -> None:
        """
        Get rid of all reactions; start again with "an empty slate" (but still with reference
        to the same "Chemicals" object)
        :return:
        """
        self.reaction_list = []



    def describe_reactions(self, concise=False) -> list:
        """
        Print out and return a listing with a friendly form of all the reactions
        EXAMPLE:  ["(0) A <=> B  (Rf = 3.0 , Rb = 2.0)"]

        :param concise:
        :return:
        """
        print("Number of reactions: ", self.number_of_reactions())
        out = []
        for i, rxn in enumerate(self.reaction_list):
            left = self.standard_form(rxn["reactants"])
            right = self.standard_form(rxn["products"])
            rxn_description = f"{i}: {left} <-> {right}"
            if not concise:
                rxn_description += f"  (Rf = {rxn['Rf']} / Rb = {rxn['Rb']})"

            print(rxn_description)
            out.append(rxn_description)

        return out



    def standard_form(self, eqn_side: list) -> str:
        """
        Return a friendly form of the given side of a chemical equation.

        EXAMPLE:  turn [(1, 0, 1), (2, 1, 1)]  into  "Fe + 2 Cl"  (if species 0 is named "Fe" and species 1 "Cl")

        :param eqn_side:    A list encoding either side of a chemical equation
        :return:
        """
        formula_list = []
        for t in eqn_side:
            stoichiometry = t[0]
            species_index = t[1]
            # Note: the reaction order (t[2]) is not used

            if stoichiometry == 1:
                term = f"{self.chem_data.get_name(species_index)}"
            else:
                term = f"{stoichiometry} {self.chem_data.get_name(species_index)}"

            formula_list.append(term)

        return " + ".join(formula_list)




    ###################       PRIVATE             ###################


    def _parse_reaction_term(self, term, name="term"):
        """
        EXAMPLES (assuming that the species with index 6 is called "E"):
            5               (1, 5, 1)
            "F"             (1, 5, 1)
            (3, 5)          (3, 5, 1)
            (3, "F")        (3, 5, 1)
            (3, 5, 2)       (3, 5, 2)
            (3, "F", 2)     (3, 5, 2)
            Same if lists were used in lieu of tuples

        :param term:
        :param name:
        :return:
        """

        if type(term) == int:
            return  (1, term, 1)
        elif type(term) == str:
            return  (1, self.chem_data.get_index(term), 1)
        elif type(term) != tuple and type(term) != list:
            raise Exception(f"{name} must be either an integer string, or a pair or a triplet. Instead, it is {type(term)}")

        # If we get thus far, term is either a tuple or a list
        if len(term) != 3 and len(term) != 2:
            raise Exception(f"Unexpected length for {name} tuple/list: it should be 2 or 3. Instead, it is {len(term)}")

        stoichiometry = term[0]
        species = term[1]
        if type(species) == str:
            species = self.chem_data.get_index(species)
        elif type(species) != int:
            raise Exception(f"The species value must be an integer or a string. Instead, it is {species}")

        if len(term) == 2:
            return (stoichiometry, species, 1)
        else:   # Length is 3
            reaction_order = term[2]
            return (stoichiometry, species, reaction_order)
