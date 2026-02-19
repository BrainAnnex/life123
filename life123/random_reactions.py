from life123.chem_data import ChemData
from life123.reaction_registry import ReactionRegistry
import numpy as np



class RandomReactionNetwork:
    """
    Generate and manage random networks of reactions
    """

    def __init__(self, n_chems :int, n_rxns :int, relative_rxn_prob=None, seed=None):
        """
        :param n_chems:             Integer >=3, with the desired number of chemicals
        :param n_rxns:              Number of desired reactions (1 or higher)
        :param relative_rxn_prob:   [OPTIONAL] A 3-element list with the relative probability
                                        of each of the 3 reaction choices that are currently available:
                                        [ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition]
                                        Probabilities must add up to 1;
                                        if not specified, a uniform distribution is used.
                                        EXAMPLE:  [0.2, 0.5, 0.3]
        """
        RXN_CHOICES = ["ReactionSynthesis", "ReactionDecomposition"]  # Available reaction classes
        #RXN_CHOICES = ["ReactionUnimolecular", "ReactionSynthesis", "ReactionDecomposition"]

        assert type(n_chems) == int, \
            f"RandomReactionNetwork Instantiation: the argument `n_chem` must be an integer, not {type(n_chems)}"

        assert n_chems >= 3, \
            f"RandomReactionNetwork Instantiation: the number of chemicals must be at least 3 (value passed was {n_chems})"

        assert type(n_rxns) == int, \
            f"RandomReactionNetwork Instantiation: the argument `n_rxns` must be an integer, not {type(n_rxns)}"

        #assert n_rxns >= 1, \
            #f"RandomReactionNetwork Instantiation: the value passed to argument `n_rxns` ({n_rxns}) must be at least 1"


        self.chem_data = ChemData(n_chems=n_chems)        # Object of type "ChemData", with info on the individual chemicals

        all_chems = self.chem_data.get_all_labels()
        # Auto-generated strings "A", "B", ..., "Z", "Z2", "Z3"

        self.registry = ReactionRegistry(chem_data=self.chem_data)      # Object of type ReactionRegistry,
                                                                        # to manage a list of reactions, and the reaction-specific objects

        # A Generator object with Numpy’s default BitGenerator
        self.rng = np.random.default_rng(seed=seed)

        # Select each of the reactions types (using the relative probabilities, if available)
        rxn_type = self.rng.choice(a=RXN_CHOICES, size=n_rxns, replace=True, p=relative_rxn_prob)  # List of class names
        # EXAMPLE: [ReactionUnimolecular, ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition]

        self.standard_species_enthalpy = {}


        MAX_ATTEMPTS = 6        # To avoid infinite loops
        for r in rxn_type:
            # For each reaction to create
            #print(f"Attempting to add reaction of type: {r}")
            reactants, products = self._assign_chems_to_rxn(r, all_chems)
            #print(f"    Tentatively considering reactants: {reactants} | products: {products}")

            n_attempts = 1
            # Prevent re-adding reactions that identically occurred before
            while n_attempts <= MAX_ATTEMPTS and self.already_used(reactants, products):
                reactants, products = self._assign_chems_to_rxn(r, all_chems)
                #print(f"    the previous combination had already been used.  Tentatively considering new "
                #      f"reactants: {reactants} | products: {products}")
                n_attempts += 1
            assert n_attempts <= MAX_ATTEMPTS, \
                f"Giving up after {n_attempts} attempts.  It seems that you are requesting such a large number of reactions ({n_rxns}), " \
                f"relative to the {n_chems} available chemicals, that random reactions would be repeating"


            delta_enthalpy = self.random_reaction_enthalpy(reactants, products)
            i = self.registry.add_elementary_reaction(reactants=reactants, products=products, delta_H=delta_enthalpy)
            print(f"-- Added reaction {i}, with reactants: {reactants} | products: {products} | delta_enthalpy: {delta_enthalpy}")
            print(self.registry.get_reaction(i).describe(concise=False))



    def _assign_chems_to_rxn(self, reaction_type :str, all_chem_labels :[str]) -> (list, list):
        """
        From the pool of available chemicals, randomly assign reactants and products to
        the specified type of reaction.

        Note: validation checks are assumed to be done in the calling function

        :param reaction_type:    One of the following reaction types:
                                    "ReactionUnimolecular", "ReactionSynthesis", "ReactionDecomposition"
        :param all_chem_labels: All the chemical labels for the random network we're building
        :return:                The pair (reactants, products),
                                    where each element is a list of chemical labels
                                    EXAMPLE (for a synthesis reaction):    (['C', 'C'], ['A'])
        """
        n_reactants = 1
        n_products = 1

        if reaction_type == "ReactionDecomposition":
            n_products = 2      # The only scenario leading to 2 products
        elif reaction_type == "ReactionSynthesis":
            n_reactants = 2     # The only scenario leading to 2 reactants

        # Random selection of reactants are with replacement, to allow reactions such as 2 A -> B
        reactants = self.rng.choice(a=all_chem_labels, size=n_reactants, replace=True)

        # Likewise, random selection of products are with replacement, to allow reactions such as A -> 2 B
        # However, do not pick the products from the pool of reactants that were used in this reaction!
        available_pool = list(set(all_chem_labels) - set(reactants))    # EXAMPLE: ['X', 'Y']
        #print(available_pool)          # Note: the order the elements in available_pool isn't guaranteed,
                                        #       and might vary, even with runs with the same random seed
        products = self.rng.choice(a=available_pool, size=n_products, replace=True)

        reactants = reactants.tolist()  # Converts both the container (numpy.ndarray) and the element types (numpy.str_),
                                        # so as to get a python list of strings
        products = products.tolist()

        return reactants, products



    def already_used(self, reactants :[str], products :[str]) -> bool:
        """
        Check whether the given set of reactants, and the given set of products,
        have already been used in a previously-generated random reaction,
        or in its reverse reaction.

        The stoichiometry is not taken into account because, for example, if we already
        utilized D -> C , we do not want to generate D -> 2 C ,
        and thus regard it as "already used".

        ~~~ EXAMPLES ~~~
        If we previous generated  A -> B  then:
            already_used(self, reactants="A", products="B") will be True
            already_used(self, reactants="B", products="A") will be True
            already_used(self, reactants=["A", "A'], products="B") will be True

        :param reactants:   A list of the labels of the reactants in the reaction to look up
        :param products:    A list of the labels of the products in the reaction to look up
        :return:            True if found, or False if not
        """
        reactant_set = set(reactants)
        product_set = set(products)

        for rxn in self.registry.reaction_list:
            if set(rxn.extract_reactant_labels()) == reactant_set \
                    and set(rxn.extract_product_labels()) == product_set:
                return True

            # Check the reverse reaction as well
            if set(rxn.extract_reactant_labels()) == product_set \
                    and set(rxn.extract_product_labels()) == reactant_set:
                return True

        return False



    def random_species_enthalpy(self, sigma = 40.4145, n=None):
        """
        A random standard enthalpy value to assign to each hypothetical chemical.
        Essentially, a essentially a generalized, randomized version of enthalpy of formation
        for n chemical species.

        The values are generated from a normal distribution all centered at zero,
        and with the passed value for the standard deviation sigma.

        The default value for sigma has the effect of making the derived SD
        of the sum of 3 random values to be about
        sqrt(3) * sigma = 70 kJ/mol
        For example, the sum of 3 random values comes up in determining
        the reaction enthalpy changes in reactions such as  A + B <-> C ,
        as computed by random_reaction_enthalpy()
        [But in the case of reactions such as 2 A <-> C, the reaction enthalpy change
        would be sqrt(5) * sigma = 90.4 kJ/mol]

        :param sigma:   [OPTIONAL] The standard deviation of the (zero-centered)
                            normal distribution used to generate all values.
                            Default: 40.4145 kJ/mol
        :param n:       [OPTIONAL] If not provided, return a single number;
                            otherwise, return a list of n values.
                            All values are kJ/mol
        :return:
        """
        return self.rng.normal(loc=0, scale=sigma, size=n)



    def get_species_enthalpy(self, label :str):
        """
        Look up the standard enthalpy value to assign to the given chemical;
        if not found, assign and store a new one

        :param label:
        :return:
        """
        if value := self.standard_species_enthalpy.get(label):
            return value
        value = self.random_species_enthalpy()   # Assign a new value
        self.standard_species_enthalpy[label] = value
        print(f"Standard species enthalpy of {value} assigned to `{label}`")
        return value



    def random_reaction_enthalpy(self, reactants :[str], products :[str]):
        """
        Attempt to provide a plausible regime for non-enzymatic small-molecule aqueous reactions,
        and thermodynamically consistent with all the previous random reactions added so far.

        It aims to randomly generate reaction ΔH values that are normally-distributed, with a mean of zero,
        and a standard deviation of roughly 70-90 kJ/mol

        With σ ≈ 70 kJ/mol for ΔH:
            ~68% of reactions fall within ±70 kJ/mol    (one standard deviation from the mean)
            ~95% within ±140 kJ/mol                     (2 standard deviations from the mean)

        That is the plausible range for non-enzymatic small-molecule chemistry.  [TODO: FURTHER VERIFY]

        So statistically:
            ✔ You’ll get many moderate reactions
            ✔ Some strongly exothermic/endothermic ones

        For a given elementary reaction, specified by its reactants and products,
        compute the Sum_i[v_i * H0_i] ,
        where i ranges over all the species in the reaction,
        H0_i is the standard enthalpy of that species,
        and v_i is the stoichiometric coefficient (negative for reactants,
        and positive for products)

        :param reactants:
        :param products:
        :return:            A random value for the change in enthalpy of the given reaction,
                                normally-distributed with mean 0 and σ ≈ 70 kJ/mol,
                                thermodynamically consistent with all the previous random reactions added so far
        """
        reaction_enthalpy = 0

        if len(reactants) == 1:
            h = self.get_species_enthalpy(reactants[0])
            reaction_enthalpy -= h
        elif reactants[0] == reactants[1]:
            h = self.get_species_enthalpy(reactants[0])
            reaction_enthalpy -= 2 * h
        else:
            h_1 = self.get_species_enthalpy(reactants[0])
            h_2 = self.get_species_enthalpy(reactants[1])
            reaction_enthalpy = reaction_enthalpy - h_1 - h_2


        if len(products) == 1:
            h = self.get_species_enthalpy(products[0])
            reaction_enthalpy += h
        elif products[0] == products[1]:
            h = self.get_species_enthalpy(products[0])
            reaction_enthalpy += 2 * h
        else:
            h_1 = self.get_species_enthalpy(products[0])
            h_2 = self.get_species_enthalpy(products[1])
            reaction_enthalpy = reaction_enthalpy + h_1 + h_2


        return reaction_enthalpy
