from life123.species_registry import SpeciesRegistry
from life123.reaction_registry import ReactionRegistry
from life123.thermodynamics import ThermoDynamics
import math
import numpy as np



class RandomReactionNetwork:
    """
    Generate and manage random networks of reactions that are thermodynamically possible,
    and biologically plausible
    """

    def __init__(self, n_species :int, n_rxns :int, relative_rxn_prob=None, seed=None,
                temp=298.15, thermodynamic_ruggedness = 2.0, verbose=False):
        """
        :param n_species:           Integer >=3, with the desired number of chemicals species
        :param n_rxns:              Number of desired reactions (ok to use 0 for testing)
        :param relative_rxn_prob:   [OPTIONAL] A 3-element list with the relative probability
                                        of each of the 3 reaction choices that are currently available:
                                        [ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition]
                                        Probabilities must add up to 1;
                                        if not specified, a uniform distribution is used.
                                        EXAMPLE:  [0.2, 0.5, 0.3]
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
                                            A value of 40.4145 for sigma has the effect of
                                            making the derived SD for the REACTION delta_H value to be sqrt(3) * sigma = 70 kJ/mol
                                            For more details, see random_species_enthalpy()

        :param verbose:             [OPTIONAL] Some extra printout if True
        """
        self.species_sigma_entropy = thermodynamic_ruggedness * 40.4145 / 2

        self.kB_over_h = 2.08366191233e10 # K-1 s-1
                                # Boltzmann constant / Plank constant
                                # Boltzmann constant = 1.380649 × 10-23 J K-1
                                # Plank constant = 6.62607015 × 10-34 J s

        RXN_CHOICES = ["ReactionSynthesis", "ReactionDecomposition"]  # Available reaction classes
        #RXN_CHOICES = ["ReactionUnimolecular", "ReactionSynthesis", "ReactionDecomposition"]

        assert type(n_species) == int, \
            f"RandomReactionNetwork Instantiation: the argument `n_chem` must be an integer, not {type(n_species)}"

        assert n_species >= 3, \
            f"RandomReactionNetwork Instantiation: the number of chemical species must be at least 3 (value passed was {n_species})"

        assert type(n_rxns) == int, \
            f"RandomReactionNetwork Instantiation: the argument `n_rxns` must be an integer, not {type(n_rxns)}"


        self.n_species = n_species
        self.n_rxns = n_rxns
        self.species_data = SpeciesRegistry(n_species=n_species)        # Object of type "SpeciesRegistry", with info on the individual chemical species


        self.reaction_data = ReactionRegistry(species_data=self.species_data)   # Object of type ReactionRegistry,
                                                                                # to manage a list of reactions, and the reaction-specific objects


        self.rng = np.random.default_rng(seed=seed)      # A Generator object with Numpy’s default BitGenerator

        # Select each of the reactions types for the reactions (using the relative probabilities, if available)
        rxn_types = self.rng.choice(a=RXN_CHOICES, size=n_rxns, replace=True, p=relative_rxn_prob)  # List of class names
        # EXAMPLE: [ReactionSynthesis, ReactionDecomposition, ReactionDecomposition, ReactionSynthesis]

        self.standard_species_enthalpy = {}

        self.temp = temp

        self.setup_random_reactions(rxn_types=rxn_types, verbose=verbose)



    def setup_random_reactions(self, rxn_types :list[str], verbose=False) -> None:
        """
        Create a network of randon reactions

        :param rxn_types:   List of names of reaction classes.
                                EXAMPLE: ["ReactionSynthesis", "ReactionDecomposition", "ReactionDecomposition"]
        :param verbose:
        :return:            None
        """
        MAX_ATTEMPTS = 6        # To avoid infinite loops

        all_chems = self.species_data.get_all_species_ids()             # Get all the id's of the species we're using
        # Auto-generated strings "A", "B", ..., "Z", "Z2", "Z3", ...

        k_diffusion_limit = 1e9  # Smoluchowski diffusion limit for bimolecular reactions

        for r in rxn_types:
            # For each reaction to create.  EXAMPLE: "ReactionSynthesis"
            #print(f"Attempting to add reaction of type: {r}")
            reactants, products = self._assign_chems_to_rxn(r, all_chems)
            #print(f"    Tentatively considering reactants: {reactants} | products: {products}")

            n_attempts = 1
            # Prevent re-adding reactions that identically occurred before
            while n_attempts <= MAX_ATTEMPTS and self._already_used(reactants, products):
                reactants, products = self._assign_chems_to_rxn(r, all_chems)
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

            i = self.reaction_data.add_elementary_reaction(reactants=reactants, products=products,
                                                           delta_H=delta_enthalpy, delta_S=delta_entropy,
                                                           kF=kF,
                                                           temp=self.temp)

            if verbose:
                print(f"-- Added reaction {i}, with reactants: {reactants} | products: {products} | delta_enthalpy: {delta_enthalpy}")
                print(self.reaction_data.get_reaction(i).describe(concise=False))



    def get_reaction_data(self):
        """
        Extract and return all the reactions

        :return:    Object of type "ReactionRegistry"
        """
        return self.reaction_data



    def _assign_chems_to_rxn(self, reaction_type :str, all_chem_labels :[str]) -> (list, list):
        """
        From the pool of available chemicals, randomly assign reactants and products to
        the specified type of reaction.

        Note: validation checks are assumed to be done in the calling function

        :param reaction_type:   One of the following reaction types:
                                    "ReactionSynthesis", "ReactionDecomposition"
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
        else:
            raise Exception(f"Only `ReactionDecomposition` and `ReactionSynthesis` are currently allowed; "
                            f"you requested `{reaction_type}`")


        # Random selection of reactants are WITH replacement, to allow reactions such as 2 A -> B
        reactants = self.rng.choice(a=all_chem_labels, size=n_reactants, replace=True)

        # Likewise, random selection of products are with replacement, to allow reactions such as A -> 2 B
        # However, do NOT pick the products from the pool of reactants that were used in this reaction!
        available_pool = [c for c in all_chem_labels
                             if c not in reactants]         # EXAMPLE: ['X', 'Y']
        #print(available_pool)

        products = self.rng.choice(a=available_pool, size=n_products, replace=True)

        reactants = reactants.tolist()  # Converts both the container (numpy.ndarray) and the element types (numpy.str_),
                                        # so as to get a python list of strings
        products = products.tolist()

        return reactants, products



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

        :param reactants:   A list of the labels of the reactants in the reaction to look up
        :param products:    A list of the labels of the products in the reaction to look up
        :return:            True if found, or False if not
        """
        reactant_set = set(reactants)
        product_set = set(products)

        for rxn in self.reaction_data.reaction_list:
            if set(rxn.extract_reactant_labels()) == reactant_set \
                    and set(rxn.extract_product_labels()) == product_set:
                return True

            # Check the reverse reaction as well
            if set(rxn.extract_reactant_labels()) == product_set \
                    and set(rxn.extract_product_labels()) == reactant_set:
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

        So statistically:
            ✔ You’ll get many moderate reactions
            ✔ Some strongly exothermic/endothermic ones

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
        Attempt to provide plausible random entropy values for non-enzymatic small-molecule aqueous reactions.

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
        Using the Brønsted–Evans–Polanyi relation (BEP),
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
        :param alpha:   The multiplicative parameter of the Brønsted–Evans–Polanyi relation (BEP) model
        :param beta:    The multiplicative parameter of the BEP model, in 70 kJ/mol
        :param sigma:   The standard deviation, in kJ/mol, of the Gaussian noise to add to the estimates

        :return:        An estimate, based on the Brønsted–Evans–Polanyi relation (BEP),
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
