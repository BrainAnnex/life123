from life_1D.bio_sim_1d import BioSim1D as bio

class Reactions:
    """
    Data about all applicable reactions
    """

    def __init__(self):
        self.reaction_list = []      # List of dicts.  Each item represents a reaction, incl. its reverse


    def number_of_reactions(self):
        return len(self.reaction_list)


    def get_reaction(self, i) -> dict:
        return self.reaction_list[i]

    def get_reactants(self, i):
        rxn = self.get_reaction(i)
        return rxn["reactants"]

    def get_reactant_names(self, i):
        rxn = self.get_reaction(i)
        reactants = rxn["reactants"]
        return self.standard_form(reactants)


    def get_products(self, i):
        rxn = self.get_reaction(i)
        return rxn["products"]

    def get_product_names(self, i):
        rxn = self.get_reaction(i)
        products = rxn["products"]
        return self.standard_form(products)


    def get_forward_rate(self, i):
        rxn = self.get_reaction(i)
        return rxn["Rf"]

    def get_back_rate(self, i):
        rxn = self.get_reaction(i)
        return rxn["Rb"]


    def add_reaction(self, reactants: list, products: list, forward_rate: float, back_rate: float) -> None:
        """
        Add the parameters of a SINGLE reaction, including its reverse rate

        :param reactants:
        :param products:
        :param forward_rate:
        :param back_rate:
        :return:                None
        """
        rxn = {"reactants": reactants, "products": products, "Rf": forward_rate, "Rb": back_rate}
        self.reaction_list.append(rxn)



    def standard_form(self, eqn_side: list) -> str:
        """
        Return a friendly form of the given side of a chemical equation.

        EXAMPLE:  turn [(0, 1), (1, 5)]  into  "Fe + 5 Cl"  (if species 0 is named "Fe" and species 1 "Cl")

        :param eqn_side:    A list encoding either side of a chemical equation
        :return:
        """
        formula_list = []
        for t in eqn_side:
            species_index, stoichiometry = t
            if stoichiometry == 1:
                term = f"{bio.get_name(species_index)}"
            else:
                term = f"{stoichiometry} {bio.get_name(species_index)}"

            formula_list.append(term)

        return " + ".join(formula_list)
