from life123.species_registry import Species, SpeciesRegistry
from life123.reaction_registry import ReactionRegistry
from life123.reactions import SimulationReaction, ReactionDefinition
from life123.stoichiometry import Stoichiometry
from life123.thermodynamics import ThermoDynamics
from typing import NamedTuple
import math
import numpy as np



class ReactionTerm(NamedTuple):
    # EXAMPLE: (2, "R", 20)     used to indicate a reaction term 2R , with a total atom count of 20 (each species R contains 10 atoms)
    coeff: int
    species_id: str
    atom_count: int


class RandomReactionNetwork:
    """
    Generate and manage random networks of reactions that are thermodynamically possible,
    and biologically plausible
    """

    def __init__(self, n_species :int, n_rxns :int, seed=None,
                 temp=298.15, thermodynamic_ruggedness = 2.0, verbose=False):
        """
        :param n_species:           Integer >=3, with the desired number of chemicals species
        :param n_rxns:              Number of desired reactions (ok to use 0 for testing)
        :param seed:                [OPTIONAL] A (large) integer, to guarantee consistent random runs
        :param temp:                [OPTIONAL] Temperature in Kelvins.  By default, 298.15 K
        :param thermodynamic_ruggedness: [OPTIONAL] It affects the magnitudes of the random reactions delta Gibbs energy values,
                                            by means of picking the standard deviation of the normal distribution
                                            used to generate the random values for the SPECIES enthalpy.
                                            Very rough guidelines:
                                                ruggedness	σ(H)	    behavior
                                                    0.5	    10	    metabolism-like
                                                    1	    20	    moderate
                                                    1.5	    30	    broad chemistry
                                                    2	    40	    rugged chemistry
                                            For example, a value of 40.4145 for sigma has the effect of
                                            making the derived SD for the REACTION delta_H value to be sqrt(3) * sigma = 70 kJ/mol
                                            For more details, see random_species_enthalpy()

        :param verbose:             [OPTIONAL] Some extra printout if True
        """
        self.species_sigma_entropy = thermodynamic_ruggedness * 40.4145 / 2

        self.kB_over_h = 2.08366191233e10 # K-1 s-1
                                # Boltzmann constant / Plank constant
                                #       Boltzmann constant = 1.380649 × 10-23 J K-1
                                #       Boltzmann constant = 1.380649 × 10-23 J K-1
                                #       Plank constant = 6.62607015 × 10-34 J s

        self.RXN_CHOICES = ["ReactionSynthesis", "ReactionDecomposition"]  # NOT IN CURRENT USE
        # Available reaction classes


        assert type(n_species) == int, \
            f"RandomReactionNetwork Instantiation: the argument `n_species` must be an integer, not {type(n_species)}"

        assert n_species >= 3, \
            f"RandomReactionNetwork Instantiation: the number of species must be at least 3 (value passed was {n_species})"

        assert type(n_rxns) == int, \
            f"RandomReactionNetwork Instantiation: the number of reactions specified in `n_rxns` must be an integer, not {type(n_rxns)}"


        self.n_species = n_species
        self.n_rxns = n_rxns
        self.species_registry = SpeciesRegistry(n_species=n_species)    # Object of type "SpeciesRegistry", with info on the individual species
                                                                        # They will get named "A", "B", "C", ...


        self.reaction_registry = ReactionRegistry(species_registry=self.species_registry)   # Object of type ReactionRegistry,
                                                                                    # to manage the randomly-generated reactions


        self.rng = np.random.default_rng(seed=seed)      # A Generator object with Numpy’s default BitGenerator

        # Select each of the reactions types for the reactions (using the relative probabilities, if available)
        #rxn_types = self.rng.choice(a=self.RXN_CHOICES, size=n_rxns, replace=True, p=relative_rxn_prob)  # List of class names
        # EXAMPLE: [ReactionSynthesis, ReactionDecomposition, ReactionDecomposition, ReactionSynthesis]

        self.standard_species_enthalpy = {}             # A dictionary of random species enthalpies assigned to each species.
                                                        # EXAMPLE: { 'A': -10.474684757627816, 'B': 3.4399997584256923}

        self.temp = temp                                # System temperature in Kelvins

        #self.setup_random_reactions(rxn_types=rxn_types, verbose=verbose)


    """
    def _generate_reaction_candidate(self) -> str:
        
        Randomly pick a reaction type from the list of choices in `self.RXN_CHOICES`
        The object variable `self.relative_rxn_prob`, if set, determines the relative probabilities

        :return:    A string that is one of the entries of self.RXN_CHOICES ;
                        EXAMPLE: "ReactionSynthesis"
        
        rxn_type = self.rng.choice(a=self.RXN_CHOICES, size=1, p=self.relative_rxn_prob)
        return rxn_type
    """

    def foo_OLD(self):
        # Reaction R1 + R2 <-> P
        r, p = self._assign_species_to_rxn(self.species_registry.get_all_species_ids())

        r1, r2 = r

        n_r1 = self.get_number_atoms(r1)
        n_r2 = self.get_number_atoms(r2)
        n_p = self.get_number_atoms(p)

        # Get the total number of None values appearing in those variables
        none_count = sum(x is None for x in (n_r1, n_r2, n_p))

        if none_count == 0:
            # Atom counts are present in all variables; only thing we can do is to verify if compatible
            return self.validate_atom_count(reactant_counts=(n_r1, n_r2), product_counts=(n_p,))
        elif none_count == 1:
            self.assign_single_atom_count(reactants=r, reactant_counts=(n_r1, n_r2),
                                    products=(p,), product_counts=(n_p,))
        elif none_count == 2:
            self.assign_double_atom_count(reactants=r, reactant_counts=(n_r1, n_r2),
                                    products=(p,), product_counts=(n_p,))
        else:
            self.assign_triple_atom_count(reactants=r,
                                    products=(p,))


    def assign_triple_atom_count(self, reactants :tuple,
                                products :tuple,) -> bool:

        r1, r2 = reactants
        p, _ = products

        # All values are unset: we have complete freedom to set them as we please
        self.set_number_atoms(r1, atom_count=16)
        self.set_number_atoms(r2, atom_count=16)
        self.set_number_atoms(p, atom_count=32)
        return True


    def assign_single_atom_count(self, reactants :tuple, reactant_counts :tuple,
                                products :tuple, product_counts :tuple) -> bool:
        """
        Exactly one count is missing

        :param reactants:
        :param reactant_counts:
        :param products:
        :param product_counts:
        :return:                True if successful, or False if not
        """
        r1, r2 = reactants
        n_r1, n_r2 = reactant_counts
        p, _ = products
        n_p, _ = product_counts

        if r1 == r2:
            # Reaction of type 2 R <-> P
            if n_p is None:
                # Only n_p is None
                self.set_number_atoms(p, atom_count = 2*n_r1)
                return True
            else:
                # Only n_r1 (and the identical n_r2) is None
                if n_p % 2 != 0:    # If odd
                    return False    # Impossible scenario
                else:
                    half_product_count = n_p / 2
                    self.set_number_atoms(r1, atom_count=half_product_count)
                    self.set_number_atoms(r2, atom_count=half_product_count)
                    return True

        else:
            # Reaction of type R1 + R2 <-> P , with R1 != R2
            if n_p is None:
                self.set_number_atoms(p, atom_count=n_r1+n_r2)
                return True
            elif n_r1 is None:
                if n_r2 >= n_p:
                    return False
                self.set_number_atoms(r1, atom_count=n_p-n_r2)
                return True
            else:
                # n_r2 is None
                if n_r1 >= n_p:
                    return False
                self.set_number_atoms(r1, atom_count=n_p-n_r1)
                return True



    def assign_double_atom_count(self, reactants :tuple, reactant_counts :tuple,
                          products :tuple, product_counts :tuple):
        """
        Exactly 2 values are missing

        :param reactants:
        :param reactant_counts:
        :param products:
        :param product_counts:
        :return:                True if successful, or False if not
        """
        r1, r2 = reactants
        n_r1, n_r2 = reactant_counts
        p, _ = products
        n_p, _ = product_counts

        if r1 == r2:
            # Reaction of type 2 R <-> P
            if n_r1 is None and n_p is None:
                # All values are unset: we have complete freedom to set them as we please
                self.set_number_atoms(r1, atom_count=16)
                self.set_number_atoms(r2, atom_count=16)
                self.set_number_atoms(p, atom_count=32)
                return True

        else:
            # Reaction of type R1 + R2 <-> P , with R1 != R2
            if n_p is None:
                if n_r1 is None:
                    n_r1 = self.set_number_atoms(r1, atom_count=16)
                if n_r2 is None:
                    n_r2 = self.set_number_atoms(r2, atom_count=16)
                self.set_number_atoms(p, atom_count=n_r1+n_r2)
                return True
            elif n_r1 is None:
                if n_r2 is None:
                    n_r2 = self.set_number_atoms(r2, atom_count=16)
                if n_p is None:
                    n_p = self.set_number_atoms(r2, atom_count=16)
                self.set_number_atoms(r1, atom_count=n_p-n_r2)
                return True
            else:
                # n_r2 is None
                if n_r1 is None:
                    n_r1 = self.set_number_atoms(r1, atom_count=16)
                if n_p is None:
                    n_p = self.set_number_atoms(p, atom_count=32)
                self.set_number_atoms(r2, atom_count=n_p-n_r1)
                return True



    def validate_atom_count(self, reactant_counts :tuple, product_counts :tuple):
        """

        :param reactant_counts:
        :param product_counts:
        :return:
        """
        return reactant_counts[0] + reactant_counts[1] == product_counts[0]   # TODO: generalize


    def set_number_atoms(self, species_id: str, atom_count :int) -> int:
        """
        Annotate the Species dataclass associated to the given species id,
        to record its specified atom count.

        The value is turned into a string, and stored as the "annotation" property.

        :param species_id:
        :param atom_count:
        :return:            The atom count passed to it
        """
        species = self.species_registry.get_species(species_id)
        species.annotation = str(atom_count)
        return atom_count



    def get_number_atoms(self, species_id: str) -> int:
        """
        Look up the number of atoms set for the given species id.
        If not set, return zero

        :param species_id:
        :return:
        """
        field_value = self.species_registry.get_species(species_id).annotation
        if field_value == "":   # Not set
            return 0

        return int(field_value)


    def big(self):
        # Reaction R1 + R2 <-> P
        rxn = self._assign_species_to_rxn()     # Object of type SimulationReaction

        self._fulfill_reaction(rxn.stoichiometry)



    def solve_lhs_LONG(self, unset_reactants, total_reactant_atoms, total_product_atoms) -> bool:
        """
        In general, this would be an equation of the form SUM[a_i * x_i] + total_reactant_atoms = total_product_atoms
        i.e.   SUM[a_i * x_i]  = (total_product_atoms - total_reactant_atoms)
        but the only scenarios we're covering for the lhs are :  x0 + x1  or  2 x0

        :param unset_reactants:
        :param total_reactant_atoms:
        :param total_product_atoms:
        :return:
        """
        number_unset_reactants = len(unset_reactants)

        diff = total_product_atoms - total_reactant_atoms
        if diff <= 0:
            return False        # No solution exists

        if number_unset_reactants == 1:     # We already handled the case where it's zero
            # Scenario:  a0 * x0 = diff
            b0 = unset_reactants[0].coeff
            if diff % b0 != 0:
                return False    # No integer solution exists
            x0 = diff // b0
            self.set_number_atoms(unset_reactants[0][1], x0)
            return True
        elif number_unset_reactants == 2:
            # Scenario:  a0 * x0 + a1 * x1 = diff, but we only support x0 + x1 = diff
            assert unset_reactants[0].coeff == 1, "Unsupported reaction stoichiometry"
            assert unset_reactants[1].coeff == 1, "Unsupported reaction stoichiometry"
            # We have x0 + x1 = diff
            if diff < 2:
                return False    # No solution exists
            x0 = diff // 2      # Floor division (discard the decimal and rounds down)
            x1 = diff - x0
            self.set_number_atoms(unset_reactants[0].species_id, x0)
            self.set_number_atoms(unset_reactants[1].species_id, x1)
            return True
        else:
            raise Exception("Unsupported reaction stoichiometry")



    def solve_lhs(self, a_list, alpha, beta) -> None|tuple:
        """
        In general, this would be an equation of the form SUM[a_i * x_i] + alpha = beta
        i.e.   SUM[a_i * x_i]  = (beta - alpha)
        but the only scenarios we're covering for the summation are :  x0 + x1  or  a0*x0

        :param a_list: Only acceptable values are [1, 1] and [2]
        :param alpha:
        :param beta:
        :return:
        """
        a_length = len(a_list)

        assert a_length > 0, "this scenario (empty `a_list`) should be handled by the calling function"

        diff = beta - alpha
        if diff <= 0:
            return None        # No solution exists

        if a_length == 1:
            # Scenario:  a0 * x0 = diff
            b0 = a_list[0]
            if diff % b0 != 0:
                return None     # No integer solution exists
            x0 = diff // b0
            return (x0, )

        elif a_length == 2:
            # Scenario:  a0 * x0 + a1 * x1 = diff, but we only support x0 + x1 = diff
            assert a_list == [1, 1], f"Unsupported reactant stoichiometry: {a_list}"
            # We have x0 + x1 = diff
            if diff < 2:
                return None     # No solution exists
            x0 = diff // 2      # Floor division (discard the decimal and rounds down)
            x1 = diff - x0
            return (x0, x1)
        else:
            raise Exception(f"Unsupported reactant stoichiometry: {a_list}")



    def solve_both_sides(self, a_list, alpha, b_list, beta) -> None|tuple:
        """
        In general, this would be an equation of the form:
                    SUM[a_i * x_i] + alpha = SUM[b_i * y_i] + beta
        i.e.:       SUM[a_i * x_i] - SUM[b_i * y_i] = diff  , where diff = (beta - alpha)

        We only support:
            x0 + x1 + alpha = y0 + beta
            2 x0    + alpha = y0 + beta

        :param a_list: Only acceptable values are [1, 1] and [2]
        :param alpha:
        :param b_list: Only acceptable value is [1]
        :param beta:
        :return:
        """
        assert len(a_list) > 0, "this scenario (empty `a_list`) should be handled by the calling function"
        assert len(b_list) > 0, "this scenario (empty `b_list`) should be handled by the calling function"

        diff = beta - alpha

        # We only support:
        #             x0 + x1 + alpha = y0 + beta
        #             2 x0    + alpha = y0 + beta
        assert a_list == [1,1] or a_list == [2], f"Unsupported reactant stoichiometry: {a_list}"
        assert b_list == [1], f"Unsupported product stoichiometry: {b_list}"


        if len(a_list) == 1:
            # Scenario: 2 x0  + alpha = y0 + beta
            # ie.:      y0 = 2 x0 - diff
            # if we set x0 = |diff| + 1 , it will be positive and > 0;
            # then: y0 = 2 x0 - diff = 2 (|diff| + 1) - diff = (2 |diff| - diff) + 2 will be positive and > 0
            x0 = abs(diff) + 1
            y0 = 2 * x0 - diff
            #self.set_number_atoms(unset_reactants[0].species_id, x0)
            #self.set_number_atoms(unset_products[0].species_id, y0)
            return  ( (x0, ), (y0, ) )

        else:   # If we get here, `a_list` is of size 2
            # Scenario: x0 + x1 + alpha = y0 + beta
            # i.e.:     x0 + x1 - y0 = diff
            if diff == 0:
                # Equation becomes  x0 + x1 - y0 = 0 ; one solution is x0=16, x1=16, y0=32
                x0 = 16
                x1 = 16
                y0 = 32
            elif diff > 0:
                # Equation is x0 + x1 - y0 = diff, where diff > 0 ; can simply set all vars to diff
                x0 = diff
                x1 = diff
                y0 = diff
            else:   # Case diff < 0
                # Equation is x0 + x1 - y0 = -eps, with eps = -diff
                # Can be re-written as:  y0 = (x0 + x1) + eps
                # Since eps > 0, can set x0=eps, x1=eps, y0=3 eps
                eps = -diff
                x0 = eps
                x1 = eps
                y0 = 3*eps

            #self.set_number_atoms(unset_reactants[0].species_id, x0)
            #self.set_number_atoms(unset_reactants[1].species_id, x1)
            #self.set_number_atoms(unset_products[0].species_id, y0)
            return ( (x0,x1), (y0, ) )



    def solve_rhs_UNNECESSARY(self, alpha, b_list, beta) -> None | tuple:
        """
        In general, this would be an equation of the form:   alpha = SUM[b_i * y_i] + beta
        i.e.   SUM[b_i * y_i]  = (alpha - beta)
        but the only scenario we're covering for the rhs are :  y0

        :param alpha:
        :param b_list:
        :param beta:
        :return:
        """

        diff =  alpha - beta
        if diff <= 0:
            return None        # No solution exists

        if number_unset_products == 1:     # We already handled the case where it's zero
            # Scenario:  b0 * y0 = diff
            b0 = unset_products[0].coeff
            assert b0 == 1, "Unsupported reaction stoichiometry"
            # We have y0 = diff
            y0 = diff
            #self.set_number_atoms(unset_products[0].species_id, diff)
            return (y0, )
        else:
            raise Exception("Unsupported reaction stoichiometry")



    def _fulfill_reaction(self, stoic :Stoichiometry) -> bool:
        """
        Only supports reactions of type  R1 + R2 <-> P  or  2 R <-> P

        :param stoic:   Object of type "Stoichiometry"
        :return:        True if possible, or False otherwise
        """
        reactants = stoic.get_reactant_list(exclude_catalysts=True) # EXAMPLE: [(1, "R1"), (1, "R2")]
                                                                    # EXAMPLE: [(2, "R")]

        products = stoic.get_product_list(exclude_catalysts=True)   # EXAMPLE: [(1, "P")]


        reactants_with_counts = [ReactionTerm( coeff=r[0], species_id=r[1], atom_count=r[0] * self.get_number_atoms(r[1]) )
                                    for r in reactants]

        products_with_counts = [ReactionTerm( coeff=p[0], species_id=p[1], atom_count=p[0] * self.get_number_atoms(p[1]) )
                                    for p in products]
        print(f"   reactants_with_counts: {reactants_with_counts} \n   products_with_counts: {products_with_counts}")

        total_reactant_atoms = sum(r.atom_count for r in reactants_with_counts)
        total_product_atoms  = sum(p.atom_count for p in products_with_counts)
        print(f"   total_reactant_atoms: {total_reactant_atoms} | total_product_atoms: {total_product_atoms}")

        unset_reactants = [ r  for r in reactants_with_counts if r.atom_count == 0]     # EXAMPLE: [(2, "R", 0)]
        unset_products =  [ p  for p in products_with_counts if p.atom_count == 0]     # EXAMPLE: [(1, "P", 0)]
        print(f"   unset_reactants: {unset_reactants} \n   unset_products: {unset_products}")

        unset_reactant_coefficients = [ r.coeff  for r in reactants_with_counts if r.atom_count == 0]       # EXAMPLE: [1, 1]
        unset_product_coefficients =  [ p.coeff  for p in products_with_counts if p.atom_count == 0]        # EXAMPLE: [2, 1]
        print(f"   unset_reactant_coefficients: {unset_reactant_coefficients} \n   unset_product_coefficients: {unset_product_coefficients}")

        number_unset_reactants = len(unset_reactants)
        number_unset_products =  len(unset_products)



        if number_unset_reactants == 0 and number_unset_products == 0:
            # No variables to adjust; the only thing we can do is to verify a match
            return total_reactant_atoms == total_product_atoms


        if number_unset_reactants > 0 and number_unset_products == 0:
            # In general, this would be an equation of the form SUM[a_i * x_i] + total_reactant_atoms = total_product_atoms
            #   i.e.   SUM[a_i * x_i]  = (total_product_atoms - total_reactant_atoms)
            #   but the only scenarios we're covering for the lhs are :  x0 + x1  or  2 x0
            result = self.solve_lhs(a_list=unset_reactant_coefficients, alpha=total_reactant_atoms, beta=total_product_atoms)
            if result is None:
                return False
            for i, x in enumerate(result):
                self.set_number_atoms(unset_reactants[i].species_id, x)
            return True


        if number_unset_reactants == 0 and number_unset_products > 0:
            # In general, this would be an equation of the form total_reactant_atoms = SUM[b_i * y_i] + total_product_atoms
            #   i.e.   SUM[b_i * y_i]  = (total_reactant_atoms - total_product_atoms)
            #   but the only scenario we're covering for the rhs are :  y0
            result = self.solve_lhs(a_list=unset_product_coefficients, alpha=total_product_atoms, beta=total_reactant_atoms)
            if result is None:
                return False
            for i, x in enumerate(result):
                self.set_number_atoms(unset_products[i].species_id, x)
            return True


        if number_unset_reactants > 0 and number_unset_products > 0:
            # In general, this would be an equation of the form:
            #           SUM[a_i * x_i] + total_reactant_atoms = SUM[b_i * y_i] + total_product_atoms
            result = self.solve_both_sides(a_list=unset_reactant_coefficients, alpha=total_reactant_atoms,
                                           b_list=unset_product_coefficients, beta=total_product_atoms)
            if result is None:
                return False

            r, p = result
            for i, x in enumerate(r):
                self.set_number_atoms(unset_reactants[i].species_id, x)
            for i, x in enumerate(p):
                self.set_number_atoms(unset_products[i].species_id, x)
            return True


        raise Exception("Uncaught scenario")    # TODO: delete



    def _assign_species_to_rxn(self) -> SimulationReaction:
        """
        From the pool of available species, randomly assign reactants and products to
        the standard reaction  R + S <-> P  (note that R and S at times get assigned to the same species,
        and so the reaction form becomes 2 R <-> P)

        #Along the way, assign compatible number of atoms to each reaction, and modify the species id
        #to be of the form 'X (n)' , where n is the total number of atoms

        Note: validation checks are assumed to be done in the calling function

        :return:    The corresponding "SimulationReaction" object
                       The pair (reactants, products),
                                    where each element is a list of species id's
                                    EXAMPLE (for a synthesis reaction):    (['C', 'C'], ['A'])
        """
        n_reactants = 2     # Note: they might end up being the same species
        n_products = 1

        # All the species id's available for the random network we're building
        all_species_ids = self.species_registry.get_all_species_ids()

        # Random selection of reactants are WITH replacement, to allow reactions such as 2 A -> B
        reactants = self.rng.choice(a=all_species_ids, size=n_reactants, replace=True)

        # Likewise, random selection of products are with replacement, to allow reactions such as A -> 2 B
        # However, do NOT pick the products from the pool of reactants that were used in this reaction!
        available_pool = [c for c in all_species_ids
                          if c not in reactants]         # EXAMPLE: ['X', 'Y']
        #print(available_pool)

        products = self.rng.choice(a=available_pool, size=n_products, replace=True)

        reactants = reactants.tolist()  # Converts both the container (numpy.ndarray) and the element types (numpy.str_),
                                        # so as to get a python list of strings
        products = products.tolist()

        rxn = ReactionDefinition(reactants=reactants, products=products, reaction_model="mass action",
                                 species_registry=self.species_registry)
        return rxn.get_simulation_reactions(0)



    def _assign_species_to_rxn_OLD(self, all_species_ids :[str]) -> (list, list):
        """
        From the pool of available species, randomly assign reactants and products to
        the standard reaction  R + S <-> P  (note that R and S at times get assigned to the same species,
        and so the reaction form becomes 2 R <-> P)

        Along the way, assign compatible number of atoms to each reaction, and modify the species id
        to be of the form 'X (n)' , where n is the total number of atoms

        Note: validation checks are assumed to be done in the calling function

        :param all_species_ids: All the species id's for the random network we're building
        :return:                The pair (reactants, products),
                                    where each element is a list of species id's
                                    EXAMPLE (for a synthesis reaction):    (['C', 'C'], ['A'])
        """
        n_reactants = 2     # Note: they might end up being the same species
        n_products = 1

        # Random selection of reactants are WITH replacement, to allow reactions such as 2 A -> B
        reactants = self.rng.choice(a=all_species_ids, size=n_reactants, replace=True)

        # Likewise, random selection of products are with replacement, to allow reactions such as A -> 2 B
        # However, do NOT pick the products from the pool of reactants that were used in this reaction!
        available_pool = [c for c in all_species_ids
                          if c not in reactants]         # EXAMPLE: ['X', 'Y']
        #print(available_pool)

        products = self.rng.choice(a=available_pool, size=n_products, replace=True)

        reactants = reactants.tolist()  # Converts both the container (numpy.ndarray) and the element types (numpy.str_),
                                        # so as to get a python list of strings
        products = products.tolist()

        return reactants, products



    def _atomic_accounting(self, reactants: list[str], products: list[str]):
        """

        :param reactants:   EXAMPLE: ['C', 'C']
        :param products:    EXAMPLE: ['A']
        :return:
        """
        pass


    def setup_random_reactions(self, rxn_types :list[str], verbose=False) -> None:
        """
        Create a network of random reactions, and populate with them the registry `self.reaction_data`

        :param rxn_types:   List of names of reaction classes.
                                EXAMPLE: ["ReactionSynthesis", "ReactionDecomposition", "ReactionDecomposition"]
        :param verbose:
        :return:            None
        """
        MAX_ATTEMPTS = 6        # To avoid infinite loops

        all_chems = self.species_registry.get_all_species_ids()             # Get all the id's of the species we're using
        # Auto-generated strings "A", "B", ..., "Z", "Z2", "Z3", ...

        k_diffusion_limit = 1e9  # Smoluchowski diffusion limit for bimolecular reactions

        for r in rxn_types:
            # For each reaction to create.  EXAMPLE: "ReactionSynthesis"
            #print(f"Attempting to add reaction of type: {r}")
            reactants, products = self._assign_species_to_rxn(r, all_chems)
            #print(f"    Tentatively considering reactants: {reactants} | products: {products}")

            n_attempts = 1
            # Prevent re-adding reactions that identically occurred before
            while n_attempts <= MAX_ATTEMPTS and self._already_used(reactants, products):
                reactants, products = self._assign_species_to_rxn(r, all_chems)
                #print(f"    the previous combination had already been used.  Tentatively considering new "
                #      f"reactants: {reactants} | products: {products}")
                n_attempts += 1
            assert n_attempts <= MAX_ATTEMPTS, \
                f"Giving up after {n_attempts} attempts.  It seems that you are requesting such a large number of reactions ({self.n_rxns}), " \
                f"relative to the {self.n_species} available chemicals, that random reactions would be repeating"


            delta_enthalpy = self.random_reaction_enthalpy(reactants, products)
            delta_entropy = self.random_reaction_entropy(reaction_type=r)

            delta_gibbs_energy = ThermoDynamics.gibbs_energy_from_enthalpy_entropy(delta_H=delta_enthalpy,
                                                                                   delta_S=delta_entropy,
                                                                                   temp=self.temp)

            forward_activation_gibbs_energy = self.forward_activation_gibbs_energy_BEP(delta_G=delta_gibbs_energy)
            kF = self.rate_constant_from_activation_gibbs_energy(activation_delta_G=forward_activation_gibbs_energy,
                                                                temp=self.temp)
            if r == "ReactionSynthesis":
                kF = min(kF, k_diffusion_limit)     # Apply an upper bound to reactions of the type A+B -> C
            elif r == "ReactionDecomposition":
                pass        # TODO: apply a diffusion limit to kR

            i = self.reaction_registry.add_elementary_reaction(reactants=reactants, products=products,
                                                               thermodynamic_parameters={"delta_H": delta_enthalpy, "delta_S": delta_entropy},
                                                               #delta_H=delta_enthalpy, delta_S=delta_entropy,
                                                               kinetic_parameters={"kF": kF},
                                                               temp=self.temp)

            if verbose:
                print(f"-- Added reaction {i}, with reactants: {reactants} | products: {products} | delta_enthalpy: {delta_enthalpy}")
                print(self.reaction_registry.get_reaction(i).describe(detail="full"))



    def get_reaction_data(self) -> ReactionRegistry:
        """
        Extract and return an object containing all the reactions

        :return:    Object of type "ReactionRegistry"
        """
        return self.reaction_registry



    def get_all_random_reactions(self) -> list:
        """

        :return:     List of objects of type "SimulationReaction"
        """
        return self.reaction_registry.get_all_reactions()



    def _already_used(self, reactants :[str], products :[str]) -> bool:
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

        :param reactants:   A list of the reactant species id's in the reaction to look up
        :param products:    A list of the products species id's in the reaction to look up
        :return:            True if found, or False if not
        """
        reactant_set = set(reactants)
        product_set = set(products)

        for rxn in self.reaction_registry.reaction_list:
            if rxn.stoichiometry.get_reactant_ids() == reactant_set \
                    and set(rxn.stoichiometry.get_product_ids()) == product_set:
                return True

            # Check the reverse reaction as well
            if rxn.stoichiometry.get_reactant_ids() == product_set \
                    and set(rxn.stoichiometry.get_product_ids()) == reactant_set:
                return True

        return False



    def random_species_enthalpy(self, sigma=None, n=None) -> np.ndarray | float:
        """
        Generate a random "standard enthalpy" value to assign to each of our chemical species.
        Essentially, a species standard enthalpy is a generalized, randomized version
        of enthalpy of formation, for the n chemical species.

        The values are generated from a normal distribution all centered at zero,
        and with the given value for the standard deviation sigma.

        EXAMPLE: a value of 40.4145 for sigma has the effect of making the derived SD
        of the sum of 3 random values to be about
        sqrt(3) * sigma = 70 kJ/mol
        For example, the sum of 3 random values comes up in determining
        the reaction enthalpy changes in reactions such as  A + B <-> C ,
        as computed by random_reaction_enthalpy()
        [But in the case of reactions such as 2 A <-> C, the reaction enthalpy change
        would be sqrt(2**2 + 1) * sigma = sqrt(5) * sigma = 90.4 kJ/mol]

        :param sigma:   [OPTIONAL] The standard deviation of the (zero-centered)
                            normal distribution used to generate all the values.
                            Default: 40.4145 kJ/mol
        :param n:       [OPTIONAL] The number of desired values.
                            If None (default) a single value is generated
        :return:        If n isn't None, return a Numpy Array with n numbers;
                            otherwise, return a single float.
                            All values are kJ/mol
        """
        if sigma is None:
            sigma = self.species_sigma_entropy

        return self.rng.normal(loc=0, scale=sigma, size=n)



    def get_species_enthalpy(self, label :str):
        """
        Look up the standard enthalpy value to assign to the given chemical;
        if not found, assign and store a new random one

        :param label:   To identify a chemical species
        :return:        A standard enthalpy value already assigned, or just assigned,
                            to the given chemical
        """
        if value := self.standard_species_enthalpy.get(label):
            return value
        value = self.random_species_enthalpy()   # Assign a new value
        self.standard_species_enthalpy[label] = value
        #print(f"Standard species enthalpy of {value} assigned to `{label}`")
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

        So, statistically:
            - You'll get many moderate reactions
            - Some strongly exothermic/endothermic ones

        For a given elementary reaction, specified by its reactants and products,
        compute the Sum_i[v_i * H0_i] ,
        where i ranges over all the species in the reaction,
        H0_i is the standard enthalpy of that species,
        and v_i is the stoichiometric coefficient (negative for reactants,
        and positive for products)

        :param reactants:   A list of the labels of the reactants
        :param products:    A list of the labels of the reaction products
        :return:            A random value for the change in enthalpy of the given reaction, in kJ/mol,
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



    def random_reaction_entropy(self, reaction_type :str) -> float:
        """
        Attempt to provide plausible random entropy values
        for non-enzymatic small-molecule aqueous reactions.

        Current starting point (to further verify):

        For decomposition reactions:
            ΔS ~ Normal with mean=40 and σ=20 ,  J/(mol·K)
            Note that 95% of the values will fall between 0 and 80 (i.e within 2σ of the mean)
            Further clipped to within the range [0, 100]

        For synthesis reactions:
            ΔS ~ Normal with mean=-40 and σ=20 ,  J/(mol·K)
            Note that 95% of the values will fall between -80 and 0 (i.e within 2σ of the mean)
            Further clipped to within the range [-100, 0]

        At biological temperature T ≈ 300K, our distribution centers of 40 J/mol·K will correspond to entropy values of:
            300 * 40 J/mol = 12 kJ/mol

        :param reaction_type:One of the following reaction types:
                                "ReactionSynthesis", "ReactionDecomposition"
        :return:            A random value for the change in entropy of the given reaction,
                                in J/(mol·K)
        """
        if reaction_type == "ReactionDecomposition":
            entropy = self.rng.normal(loc=40, scale=20)    # Mean and SD
            if entropy < 0:
                return 0.
            if entropy > 100:
                return 100
        elif reaction_type == "ReactionSynthesis":
            entropy = self.rng.normal(loc=-40, scale=20)    # Mean and SD
            if entropy < -100:
                return -100.
            if entropy > 0:
                return 0
        else:
            raise Exception("random_reaction_entropy(): only `ReactionDecomposition` and `ReactionSynthesis` are currently allowed")

        return entropy



    def forward_activation_gibbs_energy_normal(self, delta_G :float, mean=70., sigma=15.) -> float:
        """
        Draw values for the forward activation free energy ΔGf‡ from a normal distribution.

        We also enforce :  ΔGf‡ >= 0   and  ΔGf‡ >= ΔG

        Tentatively (to be further verified and validated),
        we use the following default parameter values:

            mean=70. kJ/mol
            sigma=15. kJ/mol

        This is a coarse way to generate random forward activation free energies,
        because no consideration is given to the reaction free energy ΔG, other than
        enforcing ΔGf‡ >= ΔG.

        For a more sophisticated approach, see forward_activation_free_energy_BEP()

        :param delta_G: The reaction free energy ΔG, in kJ/mol
        :param mean:    The standard deviation, in kJ/mol, of the Gaussian values
        :param sigma:   The standard deviation, in kJ/mol, of the Gaussian values

        :return:        An estimate, based on a Normal distribution,
                            of a plausible forward activation free energy, in kJ/mol,
                            for a random reaction
        """
        forward_activation_delta_G = self.rng.normal(loc=mean, scale=sigma)    # Mean and SD

        return max(forward_activation_delta_G, max(0, delta_G))



    def forward_activation_gibbs_energy_BEP(self, delta_G :float, alpha=0.5, beta=70., sigma=12.) -> float:
        """
        Using the Brønsted-Evans-Polanyi relation (BEP),
        the forward activation Gibbs free energy ΔGf‡ is modeled by:
            ΔGf‡ = α ΔG + β + ϵ
        where
            ΔG reaction free energy
            α and β are two parameters
            ϵ is noise (because real reactions are not perfectly correlated)

        We also enforce :  ΔGf‡ >= 0   and  ΔGf‡ >= ΔG

        Tentatively (to be further verified and validated),
        we use the following default parameter values:

            α=0.5
            β=70 kJ/mol
            ϵ modeled with a Normal distribution N(0, 12^2) , i.e. σ=12 kJ/mol

        See: "Molecular Driving Forces - Statistical Thermodynamics", by Dill & Bromberg (2nd edn, 2011)

        :param delta_G: The reaction free energy ΔG, in kJ/mol
        :param alpha:   The multiplicative parameter of the Brønsted-Evans-Polanyi relation (BEP) model
        :param beta:    The multiplicative parameter of the BEP model, in 70 kJ/mol
        :param sigma:   The standard deviation, in kJ/mol, of the Gaussian noise to add to the estimates

        :return:        An estimate, based on the Brønsted-Evans-Polanyi relation (BEP),
                            of a plausible forward activation free energy, in kJ/mol,
                            for the given reaction free energy ΔG
        """
        eps = self.rng.normal(loc=0, scale=sigma)    # Mean and SD
        forward_activation_delta_G = alpha*delta_G + beta + eps

        return max(forward_activation_delta_G, max(0, delta_G))     # Disregard the PyCharm complaint



    def rate_constant_from_activation_gibbs_energy(self, activation_delta_G :float, temp :float) -> float:
        """
        Using the Eyring equation from transition state theory,
        estimate the kinetic rate constant from the activation free energy.

        We are working under the assumption that the transmission coefficient κ = 1
        (the fraction of the flux through the transition state that proceeds to the product
        without recrossing the transition state)

        For more info: https://en.wikipedia.org/wiki/Eyring_equation

        Note that, at T=300 K, every additional 5.7434 kJ/mol (R T ln10) of activation free energy
        leads to a x10 decrease in the rate constant

        :param activation_delta_G:  Activation free energy (forward or reverse), in kJ/mol
        :param temp:                System's temperature, in degree Kelvins
        :return:                    The kinetic rate constant, in 1/sec, for either the forward or reverse reaction,
                                        depending on the type of activation_delta_G
        """
        exponential = math.exp(- activation_delta_G / (ThermoDynamics.R * temp))

        # Note: self.kB_over_h * temp is about 6.25e12 s-1 at 300 K,
        #       and is called "attempt frequency"
        return self.kB_over_h * temp * exponential
