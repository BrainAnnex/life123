from __future__ import annotations      # To facilitate type annotations
import numpy as np
from typing import Set, Mapping
from dataclasses import dataclass, field



@dataclass(frozen=True)
class Stoichiometry:
    """
    Managing the SIGNED stoichiometric coefficients,
    for all species in the reaction.

    Reaction reagents have negative signs, and products have positive signs.

    Catalysts (which would have a net coefficient of zero) are stored separately
    """

    # EXAMPLES below are for the reaction A + E -> B + E

    vector: Mapping[str, int]
    # Signed stoichiometric reaction vector.
    # Maps species id -> signed stoichiometric coefficient.
    # Reactants are negative; products are positive.
    # Catalysts are not represented here, since their net coefficient is zero.
    # The signed coefficients ν_i, given a set of species X_i,
    #       allow the reaction to be expressed as : ∑i ν_i X_i = 0

    # EXAMPLE: {"A": -1, "B": 1}

    catalysts: list[str] = field(default_factory=list)
    # EXAMPLE: ["E"]    (list of species that have a net coefficient of zero in the reaction vector)


    def __post_init__(self) -> None:
        """
        Automatically invoked by the constructor, just before it terminates.

        Enforce non zero values in argument `coefficients`

        :return: None
        """
        for sp, coeff in self.vector.items():
            assert coeff != 0, f"Stoichiometry instantiation: coefficients cannot be zero (species \"{sp}\") - " \
            f"if you want to specify an enzyme/catalyst, pass it as a list to the argument `catalysts`"



    def to_dict(self) -> dict:
        """
        Return a dictionary form of the stoichiometry,
        mapping all species id's to their SIGNED stoichiometric coefficients.
        IMPORTANT: catalysts, if any, carry a net signed stoichiometric coefficient of zero

        :return:    A dictionary with the mapping of the species id's
                    to their signed stoichiometric coefficients
        """
        d = dict(self.vector).copy()  # Clone the dictionary
        for cat in self.catalysts:
            d[cat] = 0

        return d



    def get_reaction_vector(self) -> dict:
        """
        Following Martin Feinberg's "Foundations of Chemical Reaction Network Theory",
        we define the "reaction vector" of a reaction y -> y' (where y any y' are vectors)
        as:  y' - y

        IMPORTANT: Note that catalysts, if any, are NOT included.

        The component of (y′ − y) corresponding to species s is  y′_s − y_s,
        i.e the difference between the stoichiometric coefficient of s in the product complex y′ (the right-hand side of the reaction)
        and its stoichiometric coefficient in the reactant complex y (the left-hand side of the equation).
        This difference is the net number of molecules of s
        produced with each occurrence of the reaction y → y′.

        Other terms for reaction vectors: "reaction difference vector", or "reaction increment"
        ~~~
        EXAMPLE: for reaction  A + E -> 2P + Q + E
                 the reactant complex y is:     A + E
                 while product complex y′ is:   2P + Q + E
                 and the corresponding reaction vector y′ - y is:  2P + Q - A
                 The non-zero components of the reaction vector,
                 written as a mapping, are: {"A": -1, "P": 2, "Q": 1}
        ~~~

        :return:    The non-zero components of the reaction vector,
                        written as a dict mapping of species id to its component value
        """
        # Note that catalysts, if any, are NOT included
        return dict(self.vector)



    def get_reactant_list(self) -> list:
        """
        Return all reactants as a list of pairs
        of the form (stoichiometry coefficient, species id).
        Catalysts, if any, are also included.
        ~~~
        EXAMPLE: {"A": -2, "B":- 2, "P": 3, "E": 0}
                 would return to [(2, "A"), (2, "B"), (1, "E")]
        :return:
        """
        # TODO: consider returning a set instead of a list
        return [  (-v, k) for k,v in self.vector.items() if v < 0]  \
                + [(1, c) for c in self.catalysts]



    def get_reactant_ids(self, exclude_catalysts=False) -> Set[str]:
        """
        Return the set of all reactant species id's

        :param exclude_catalysts:
        :return:
        """
        s = { k for k,v in self.vector.items()  if v < 0 }    #  Set construction

        if not exclude_catalysts:
            s |= { c for c in self.catalysts}   # Set union

        return s


    def get_product_list(self) -> list:
        """
        Return all products as a list of pairs
        of the form (stoichiometry coefficient, species id).
        Catalysts, if any, are also included.
        ~~~
        EXAMPLE: {"A": -2, "B":- 2, "P": 3, "E": 0}
                 would return to [(3, "P"), (1, "E")]
        :return:
        """
        # TODO: consider returning a set instead of a list
        return [   (v, k) for k,v in self.vector.items() if v > 0] \
                + [(1, c) for c in self.catalysts]



    def get_product_ids(self, exclude_catalysts=False) -> Set[str]:
        """
        Return the set of all product species id's

        :param exclude_catalysts:
        :return:
        """
        s = { k for k,v in self.vector.items()  if v > 0 }    #  Set construction

        if not exclude_catalysts:
            s |= { c for c in self.catalysts}   # Set union

        return s



    def get_all_species_ids(self) -> Set[str]:
        """
        Return the SET of the id's of ALL the species appearing in this reaction

        :return:    A SET of the id's of the species involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        # Use set construction and set union
        return    { k for k,_ in self.vector.items() }  \
                | { c for c in self.catalysts}



    def get_reaction_complexes(self) -> tuple:
        """
        It separates the reactants and products in a pair;
        basically (the left-hand side, the right-hand side).
        Following Martin Feinberg's "Foundations of Chemical Reaction Network Theory",
        a reaction is:  reactant complex -> product complex
        ~~~
        EXAMPLE: for reaction  A + E -> 2P + Q + E
                 the reactant complex is:     A + E
                 while product complex is:   2P + Q + E
                 So, this function will return  ( {"A": 1, "E": 1} ,  {"P": 2, "Q": 1, "E": 1} )
        ~~
        :return:    The pair (reactants complex , products complex)
        """
        d = self.vector
        reactants = {k:-v  for k,v in d.items() if v < 0}
        for cat in self.catalysts:
            reactants[cat] = 1

        products = {k:v    for k,v in d.items() if v > 0}
        for cat in self.catalysts:
            products[cat] = 1

        return (reactants, products)



    def _standard_form_complex(self, complex :dict[str]) -> str:
        """
        Return a user-friendly form of a "complex" (a side of a chemical equation)

        EXAMPLE:  turn {"Fe": 1,  "Cl": 2]  into  "Fe + 2 Cl"

        :param complex:     A dictionary encoding either side of a chemical equation
        :return:            A string with a user-friendly form of a side of a chemical equation
        """
        formula_list = []
        for species_name, stoichiometry in complex.items():

            if stoichiometry == 1:
                term = species_name
            else:
                term = f"{stoichiometry} {species_name}"

            formula_list.append(term)

        return " + ".join(formula_list)



    def standard_chemical_formula(self, reversible=False) -> str:
        """
        Return,  as a string, a user-friendly plain-text form of the reaction

        :return:
        """
        reactants, products = self.get_reaction_complexes()

        arrow = " <-> " if reversible else " -> "

        return self._standard_form_complex(reactants) + arrow + self._standard_form_complex(products)



    def reaction_pattern(self) -> tuple[int, int, int]:
        """
        Return the triplet (n_reactants, n_products, n_catalysts)
        Catalysts, if present, are counted separated, and NOT included under either "reactants" or "products"
        ~~~
        EXAMPLES:
            The reaction A + B -> C      gives (2, 1, 0)
            The reaction A + E -> B + E  gives (1, 1, 1)
        ~~~

        :return:    A triplet of integers (n_reactants, n_products, n_catalysts)
        """
        values = self.vector.values()     # EXAMPLE: dict_values([-1, 1])

        negative_count = sum(v < 0 for v in values)
        positive_count = sum(v > 0 for v in values)
        zero_count     = len(self.catalysts)

        return (negative_count, positive_count, zero_count)



    def is_elementary_like(self) -> bool:
        """
        Return True if the reaction is of the form A->B, A->B+C, A->2B, A+B->C, or 2A->B
        IMPORTANT:  this stoichiometry check is neither sufficient nor strictly necessary for
                    a reaction to be elementary.
                    This function is meant as a coarse check - for example, for the purpose of
                    providing plausible defaults for common reaction types

        :return:    True if the stoichiometry is consistent with typical elementary reactions,
                        or False otherwise
        """
        if len(self.vector) > 3 or self.catalysts != []:
            return False

        # If we get thus far, we have at most 3 terms in the overall reaction - and no catalysts

        coeffs = sorted(self.vector.values())
        print(coeffs)
        patterns = [
                        [-1, 1],        # unimolecular rearrangement/isomerization
                        [-1, 1, 1],     # decomposition
                        [-1, 2],        # decomposition
                        [-1, -1, 1],    # synthesis
                        [-2, 1]         # synthesis
                ]

        return coeffs in patterns



    def consistency_checker(self, conc_before :dict, conc_after :dict) -> None:
        """
        Investigate the change in the concentration of the species involved in the reaction,
        to ascertain whether the change is consistent with the reaction's stoichiometry.

        In other words, a "stoichiometric-ratio consistency check"

        More formally stated, for a single reaction,
        with a reaction vector ν = (ν_A, ν_B, ...)  , where the coefficients are signed numbers,
        and the observed concentration change Delta c = conc_after - conc_before ,
        the question is whether there exists a single scalar "reaction extent" ξ such that
        Δc = ξ ν

        :param conc_before: A dict that maps species `id` to its initial concentration
        :param conc_after:  A dict that maps species `id` to its final concentration
        :return:            None.  Raise an Exception if the change in reactant/product concentrations
                                is not consistent with the reaction's stoichiometry
        """
        assert len(conc_before) == len(self.vector), \
            f"consistency_checker(): the argument `conc_before` must contain exactly the same keys " \
            f"as the species in this reaction: {list(self.vector.keys())}"

        assert len(conc_after) == len(self.vector), \
            f"consistency_checker(): the argument `conc_after` must contain exactly the same keys " \
            f"as the species in this reaction: {list(self.vector.keys())}"

        delta_conc = {}
        for sp,conc in conc_before.items():
            assert sp in self.vector, \
                f"consistency_checker(): the species \"{sp}\" appearing in the argument `conc_before` " \
                f"is not part of the reaction: {self.vector}"

            assert sp in conc_after, \
                f"consistency_checker(): the species \"{sp}\" appearing in the argument `conc_before` " \
                f"is not found among the species provided to the argument `conc_after`: {conc_after}"

            delta_conc[sp] = conc_after[sp] - conc_before[sp]

        #print("delta_conc: ", delta_conc)

        # Example: if the reaction's signed stoichiometric coefficient are
        #          {"A": -1, "B": 1, "E": 0}
        #          then a delta_conc of {"A": -10, "B": 10, "E": 0}
        #          is consistent because of multiple of the reaction vector
        ratio = None
        # Loop over all the values of delta_conc,
        # and make sure they are all the same multiple of the corresponding reaction coefficients
        for sp,conc in delta_conc.items():
            if self.vector[sp] == 0:
                continue
            if ratio is None:
                ratio = delta_conc[sp] / self.vector[sp]
                #print("ratio: ", ratio)
            else:
                assert np.allclose(delta_conc[sp], ratio * self.vector[sp]), \
                    f"consistency_checker(): the delta concentration {delta_conc} " \
                    f"is incompatible with the reaction's stoichiometry of {self.vector}"



    @classmethod
    def from_reactants_products(cls, reactants :list|str, products :list|str) -> Stoichiometry:
        """
        This serves as an "alternate constructor / factory" for the Stoichiometry dataclass.
        Notice that this is a CLASS method.
        ~~~
        EXAMPLE:
            stoich_obj = Stoichiometry.from_reactants_products(reactants="A", products=["B", (2,"C")])

        :param reactants:   A string,
                                or a list whose elements are either strings or pairs of the form (stoichiometry, species id);
                                whenever strings are used, the stoichiometry coefficient is taken to be 1
        :param products:    A string,
                                or a list whose elements are either strings or pairs of the form (stoichiometry, species id);
                                whenever strings are used, the stoichiometry coefficient is taken to be 1
        :return:            A dataclass object of type "Stoichiometry"
        """
        assert reactants is not None, \
            "from_reactants_products(): the argument `reactants` is a required one"
        if type(reactants) == str:
            reactants = [reactants]
        else:
            assert type(reactants) is list, \
                "from_reactants_products(): the argument `reactants` must be a list or a string"

        assert products is not None, \
            "from_reactants_products(): the argument `products` is a required one"
        if type(products) == str:
            products = [products]
        else:
            assert type(products) is list, \
                "from_reactants_products(): the argument `products` must be a list or a string"


        coeffs = {}     # Signed stoichiometric coefficients
        r = {}          # Standardized reactants: used for error checking
        p = {}          # Standardized products:  used for error checking

        for term in reactants:
            if type(term) is str:       # EXAMPLE : "R"
                c, species = 1, term
            else:                       # EXAMPLE : (2, "R")
                assert (type(term) is tuple) and (len(term) == 2), \
                    "from_reactants_products(): List elements in argument `reactants` must be strings or pairs of the form (int, string)"
                c, species = term

            assert type(c) is int, \
                "from_reactants_products(): List elements in argument `reactants` must be strings or pairs of the form (int, string)"
            assert type(species) is str, \
                "from_reactants_products(): List elements in argument `reactants` must be strings or pairs of the form (int, string)"

            coeffs[species] = coeffs.get(species, 0) - c    # Accumulate the sum of the SIGNED stoichiometric coefficients for this species
            r[species] = r.get(species, 0) + c              # Accumulate the UN-signed coefficients for the reactants


        for term in products:
            if type(term) is str:       # EXAMPLE : "P"
                c, species = 1, term
            else:                       # EXAMPLE : (2, "P")
                assert (type(term) is tuple) and (len(term) == 2), \
                    "List elements in argument `products` must be strings or pairs of the form (int, string)"
                c, species = term

            assert type(c) is int, \
                "from_reactants_products(): List elements in argument `reactants` must be strings or pairs of the form (int, string)"
            assert type(species) is str, \
                "from_reactants_products(): List elements in argument `reactants` must be strings or pairs of the form (int, string)"

            coeffs[species] = coeffs.get(species, 0) + c    # Accumulate the sum of the SIGNED stoichiometric coefficients for this species
            p[species] = p.get(species, 0) + c              # Accumulate the UN-signed coefficients for the products


        vector = {k: v for k,v in coeffs.items() if v != 0}
        catalysts = [k for k,v in coeffs.items() if v == 0]


        # Catch identical reaction sides, even if terms are reshuffled
        if r == p:      # Alternatively:    if vector == {} and catalysts != []
            lhs = [f"{v} {k}" for k,v in r.items()]             # EXAMPLE : ["2 A" , "3 B"]
            rhs = [f"{v} {k}" for k,v in p.items()]             # EXAMPLE : ["2 A" , "3 B"]
            msg = f"{' + '.join(lhs)} -> {' + '.join(rhs)}"     # EXAMPLE: "2 A + 3 B  -> 2 A + 3 B"
            raise Exception(f"from_reactants_products(): the two sides of the reaction can't be identical! \"{msg}\"")


        # Instantiate and return a Stoichiometry dataclass
        return cls(
            vector=vector,
            catalysts=catalysts,
        )
