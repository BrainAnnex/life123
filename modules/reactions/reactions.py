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

    def get_products(self, i):
        rxn = self.get_reaction(i)
        return rxn["products"]

    def get_forward_rate(self, i):
        rxn = self.get_reaction(i)
        return rxn["Rf"]

    def get_back_rate(self, i):
        rxn = self.get_reaction(i)
        return rxn["Rb"]


    def add_reaction(self, reactants: list, products: list, forward_rate: float, back_rate: float):
        rxn = {"reactants": reactants, "products": products, "Rf": forward_rate, "Rb": back_rate}
        self.reaction_list.append(rxn)
