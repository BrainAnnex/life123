from __future__ import annotations      # To facilitate type annotations
import numpy as np
from typing import Set, Mapping
from dataclasses import dataclass, field, fields, asdict
from life123.thermodynamics import ThermoDynamics
from life123.reaction_kinetics import ReactionKinetics
from life123.species_registry import Species, SpeciesRegistry
from life123.units import show_standard_units, convert, K, C




@dataclass(frozen=True)
class Stoichiometry:
    """
    Mapping species `id` to its signed stoichiometric coefficient,
    for all species in the reaction.
    Reaction reagents have negative signs, and products have positive signs.
    (including species that have a net coefficient of zero, such as catalysts that appear on both sides)

    EXAMPLE:  {"A": -1, "B": 1, "E": 0}   for the reaction A + E -> B + E
    """
    coefficients: Mapping[str, int]



    def reaction_pattern(self) -> tuple[int, int, int]:
        """
        Return the triplet (n_reactants, n_products, n_catalysts)
        Catalysts, if present, are counted separated, and NOT included under either "reactants" or "products"
        ~~~
        EXAMPLES:
        {"A": -1, "B": -1, "C": 1} , i.e. the reaction A + B -> C, leads to (2, 1, 0)

        {"A": -1, "B": 1, "E": 0}, i.e. the reaction A + E -> B + E, leads to (1, 1, 1)
        ~~~

        :return:    A triplet of integers
        """
        values = self.coefficients.values()     # EXAMPLE: dict_values([-1, 1])

        negative_count = sum(v < 0 for v in values)
        positive_count = sum(v > 0 for v in values)
        zero_count     = sum(v == 0 for v in values)

        return (negative_count, positive_count, zero_count)



    def to_dict(self) -> dict:
        """
        Return a dictionary form of the stoichiometry.

        :return:    A dictionary with the mapping of the species id's
                    to their signed stoichiometric coefficients
        """
        d = asdict(self)
        return d.get("coefficients")



    def consistency_checker(self, conc_before :dict, conc_after :dict) -> None:
        """

        :param conc_before:
        :param conc_after:
        :return:            None.  Raise an Exception if the change in reactant/product concentrations is consistent with the
                                    reaction's stoichiometry
        """
        assert len(conc_before) == len(self.coefficients), \
            f"consistency_checker(): the argument `conc_before` must contain exactly the same keys " \
            f"as the species in this reaction: {list(self.coefficients.keys())}"

        assert len(conc_after) == len(self.coefficients), \
            f"consistency_checker(): the argument `conc_after` must contain exactly the same keys " \
            f"as the species in this reaction: {list(self.coefficients.keys())}"

        delta_conc = {}
        for sp,conc in conc_before.items():
            assert sp in self.coefficients, \
                f"consistency_checker(): the species \"{sp}\" appearing in the argument `conc_before` " \
                f"is not part of the reaction: {self.coefficients}"

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
            if self.coefficients[sp] == 0:
                continue
            if ratio is None:
                ratio = delta_conc[sp] / self.coefficients[sp]
                #print("ratio: ", ratio)
            else:
                assert np.allclose(delta_conc[sp], ratio * self.coefficients[sp]), \
                    f"consistency_checker(): the delta concentration {delta_conc} " \
                    f"is incompatible with the reaction's stoichiometry of {self.coefficients}"



    def get_reaction_vector(self) -> {}:
        """
        Following Martin Feinberg's "Foundations of Chemical Reaction Network Theory",
        we define the "reaction vector" of a reaction y -> y' (where y any y' are vectors)
        as:  y' - y

        The component of (y′ − y) corresponding to species s is  y′_s − y_s,
        i.e the difference between the stoichiometric coefficient of s in the product complex y′ (the right-hand side of the reaction)
        and its stoichiometric coefficient in the reactant complex y (the left-hand side of the equation).
        This difference is the net number of molecules of s
        produced with each occurrence of the reaction y → y′.

        ~~~
        EXAMPLE: for reaction  A + E -> 2P + Q + E
                 the reactant complex y is:     A + E
                 while product complex y′ is:   2P + Q + E
                 and the corresponding reaction vector y′ - y is:  2P + Q - A
                 The non-zero components of the reaction vector, written as a mapping, are: {"A": -1, "P": 2, "Q": 1}
        ~~~

        :return:    The non-zero components of the reaction vector,
                        written as a dict mapping of species id to its component value
        """
        # Form a new dict from the dict returned by get_signed_stoichiometric_coefficients(),
        # omitting all terms with a zero value
        d = {k:v
                for k,v in self.coefficients.items()
                if v != 0}

        return d



    def get_split_reaction_vector(self) -> tuple:
        """
        Similar to get_reaction_vector(), but it separately returns the reactants and products in a pair

        :return:    The pair (reactants portion , products portion)
        """
        coeffs = self.coefficients
        reactants = {k:v  for k,v in coeffs.items() if v < 0}
        products = {k:v   for k,v in coeffs.items() if v > 0}
        return (reactants, products)





#######################################################################################################################


class Kinetics:

    def __init__(self, law :str, parameters=None):
        """
        EXAMPLES of `law`:
            "mass action"
            "Michaelis-Menten"
            "custom"            [user-supplied Python function]
            "Hill"              [not yet supported]
            "enzyme inhibition" [hypothetical future extension; not available]


        :param law:
        :param parameters:
        """
        AVAILABLE_RATE_LAWS = ["mass action", "Michaelis-Menten", "custom"]

        if law is not None:
            assert law in AVAILABLE_RATE_LAWS, \
                f"Kinetics instantiation: the value passed to the `law` argument (\"{law}\") " \
                f"is not one the allowed values: {AVAILABLE_RATE_LAWS}"


        if parameters is None:
            parameters = {}


        self.law = law
        self.kinetic_rate_function = None
        self.parameters = {}

        self.set_parameters(parameters)



    def set_parameters(self, parameters :dict) -> None:
        """

        :param parameters:
        :return:            None
        """
        if parameters is None:
            parameters = {}


        if self.law == "mass action":

            # validate that at most only "kR" and/or "kF" were passed
            allowed_keys = {"kR", "kF"}
            extra_keys = set(parameters.keys()) - allowed_keys

            if extra_keys:
                raise ValueError(f"set_parameters(): Unexpected parameter keys: {', '.join(extra_keys)}")

            kF = parameters.get("kF", self.parameters.get("kF"))    # Over-write if passed
            kR = parameters.get("kR", self.parameters.get("kR"))    # Over-write if passed
            K = None

            if kF is None:
                kF = 0

            if kR is None:
                kR = 0

            if np.allclose(kR, 0):
                reversible = False
            else:
                reversible = True

            if reversible:
                # If we have a reversible reaction that follows mass-action kinetics
                if (kF is not None) and (kR is not None) and (not np.allclose(kR, 0)):
                    K = kF / kR    # Kinetic parameter ratio
            else:
                # If we have an IR-reversible reaction that follows mass-action kinetics
                assert not kR, \
                    f"Irreversible reactions with mass-action kinetics " \
                    f"cannot have a value for the reverse rate constant (kR = {kR})"
                kR = 0


            self.parameters = {"kF": kF, "kR": kR, "reversible": reversible, "K": K}



    def to_dict(self) -> dict:
        """
        Return a dictionary form of the dataclass.
        Unset fields are omitted

        :return:    A dictionary populated with the public fields of this data class
        """
        properties = {"kinetics_type": self.law}

         # Only include the fields that were set
        if self.parameters.get("kF") is not None:
            properties['kF'] = self.parameters.get("kF")

        if self.parameters.get("kR") is not None:
            properties['kR'] = self.parameters.get("kR")

        if self.parameters.get("reversible") is not None:
            properties['reversible'] = self.parameters.get("reversible")

        return properties



    def set_rate_constants_from_equilibrium_constant(self, K :float|int) -> None:
        """
        Set, as needed, a missing reaction rate constant (kF or kR)
        from the other one and the given equilibrium constant K.
        If all values already exist, and an inconsistency is detected, an Exception will be raised.

        Note: the reaction's equilibrium constant and its kinetic rate constants are
              in the relationship K = kF / kR for any reaction that follows "mass-action kinetics",
              i.e. whose reaction rates are proportional to the product of the reactants’ concentrations
              raised to their stoichiometric coefficients

        :param K:   The reaction's equilibrium constant
        :return:    None
        """
        assert K is not None, \
            "set_rate_constants_from_equilibrium_constant(): missing value for argument `K`"


        if self.law != "mass action":
            return

        kF = self.parameters.get("kF")
        kR = self.parameters.get("kR")

        if (not kR) and (kF is not None) and (not np.allclose(K, 0)):
            kR = kF / K
            self.parameters["kR"] = kR
            if not np.allclose(kR, 0):
                self.parameters["reversible"] = True
            return

        if (not kF) and (kR is not None):
            self.parameters["kF"] = K * kR
            return

        if (kF is not None) and (kR is not None) and (not np.allclose(kR, 0)):
            assert np.allclose(K, kF / kR), \
                f"set_rate_constants_from_equilibrium_constant(): values for kR ({kR}) and kR ({kR}) already exist, " \
                f"and are inconsistent with the passed value of K ({K})"



    def extract_intermediate(self) -> str|None:
        """
        Return the name of the reaction intermediate species,
        or None if there's no intermediate

        :return:
        """
        if self.law == "Michaelis-Menten":
            return "TBA"        # TODO: FIX!

        return None



    def set_rate_function(self, f) -> None:
        """
        Set the function used to estimate the reaction rate (aka "velocity"),
        at the start of the time step.

        :param f:   A function that takes the following args:
                        reactant_terms :[(int, str)]
                        product_terms :[(int, str)],
                        kF :float, kR :float,
                        conc_dict :dict
                    and return a float
                    EXAMPLE:  ReactionKinetics.compute_rate_mass_action_kinetics
                              # Generalized "standard rate law"

        :return:    None
        """
        self.kinetic_rate_function = f






###################################################################################################################


@dataclass(slots=True)      # Note: (slots=True) has the effect of prohibiting non-listed fields,
                            #       and of making the class more efficient
class ReactionThermodynamics:
    """
    Thermodynamic data belonging to a particular reaction
    """
    delta_H: float | None = None
    delta_S: float | None = None
    delta_G: float | None = None
    K: float | None = None
    temp: float | None = None



    def to_dict(self) -> dict:
        """
        Return a dictionary form of the dataclass.
        Unset fields are omitted

        :return:    A dictionary populated with the public fields of this data class
        """
        d = asdict(self)
        result = {}
        for k, v in d.items():
            if v is not None:
                result[k] = v    # Only include the fields that were set

        return result



    def set_temperature(self, temp :float) -> None:
        """
        Set all the thermodynamic data derivable from the given temperature,
        and all stored thermodynamic data.
        Raise an Exception if any inconsistency is detected.

        :param temp:    System temperature in Kelvins.  For now, assumed constant everywhere,
                            and unvarying (or very slowly varying).
                            If the temp gradually changes, periodically call this method.
        :return:        None
        """
        # Process the thermodynamic data, and update various object attributes accordingly
        thermo_data = ThermoDynamics.extract_thermodynamic_data(K=self.K,
                                                  delta_H=self.delta_H, delta_S=self.delta_S, delta_G=self.delta_G,
                                                  temp=temp)

        #print(f"thermo_data : {thermo_data}")
        self.K = thermo_data["K"]
        self.delta_H = thermo_data["delta_H"]
        self.delta_S = thermo_data["delta_S"]
        self.delta_G = thermo_data["delta_G"]






###################################################################################################################


class Reaction:
    """

    """

    def __init__(self, reactants :str|list|tuple, products :str|list|tuple,
                 species_registry :SpeciesRegistry, name=None, autoregister_species=False,
                 active=True,
                 delta_H=None, delta_S=None, temp=None,
                 kinetics_type=None, kinetic_parameters=None):
        """

        :param reactants:   A list/tuple of terms that are either species id's (with implied stoichiometry 1),
                                or pairs (stoichiometry coefficient , species id).
                                If not a list, it will first get turned into one
        :param products:    A list/tuple of terms that are either chemicals labels (with implied stoichiometry 1),
                                or pairs (stoichiometry coefficient , chemical label).
                                If not a list, it will first get turned into one
        :param species_registry:
        :pamar auto_register_species:

        :param active:

        :param delta_H:     [OPTIONAL] Change in Enthalpy (from reactants to products), in kJ/mol
        :param delta_S:     [OPTIONAL] Change in Entropy (from reactants to products), in Joules/(mol·K)
        :param temp:        [OPTIONAL]

        :param kinetics_type:[OPTIONAL]
        """
        self.name = name

        self.active = active

        self.thermodynamics = None
        self.kinetics = None

        self.reactants = None       # A list of pairs (stoichiometry, chemical label)
        self.products = None        # A list of pairs (stoichiometry, chemical label)
        self.stoichiometry = None   # A "Stoichiometry" object
                                    #   mapping species `id` to its signed stoichiometric coefficient,
                                    #   for all species in the reaction

        self.analytic_solution_family = None    # Available values: "ONE_TO_ONE", "ONE_TO_TWO", "TWO_TO_ONE"
        self.elementary = None
        self.reaction_category = None


        assert reactants is not None, "Reaction() instantiation: the argument `reactants` is a required one"
        if type(reactants) == str:
            reactants = [reactants]

        assert products is not None, "Reaction() instantiation: the argument `products` is a required one"
        if type(products) == str:
            products = [products]


        # Normalize the elements of each list to be (int, str) pairs; i.e. turn any single string "X" into the pair (1, "X")
        reactant_list = [(1, r) if type(r) == str else r
                            for r in reactants]   # A list of pairs
        product_list =  [(1, p) if type(p) == str else p
                            for p in products]   # A list of pairs

        # Catch identical reaction sides, even if terms are reshuffled
        assert set(reactant_list) != set(product_list), \
            f"Reaction(): the two sides of the reaction can't be identical! " \
            f"Same reactant and product complexes: \"{self._standard_form_chem_eqn(reactant_list)}\""


        # Check whether all the species in the reaction are registered ones
        for _, s_id in reactant_list:
            if not species_registry.species_exists(s_id):
                if autoregister_species:
                    species_registry.add_species(id=s_id)
                else:
                    raise Exception(f'No species with id "{s_id}" exists in the species registry')

        for _, s_id in product_list:
            if not species_registry.species_exists(s_id):
                if autoregister_species:
                    species_registry.add_species(id=s_id)
                else:
                    raise Exception(f'No species with id "{s_id}" exists in the species registry')


        c = self.get_signed_stoichiometric_coefficients(reactants=reactant_list, products=product_list)
        self.stoichiometry = Stoichiometry(c)


        # TODO: move to a separate function
        # EXAMPLE: {"A": -2, "B":- 2, "P": 3, "E": 0}
        #   would lead to self.reactants = [(2, "A"), (2, "B"), (1, "E")]
        #              and self.products = [(3, "P"), (1, "E")]
        self.reactants = [(-v, k) for k,v in c.items() if v < 0] + [(1, k) for k,v in c.items() if v == 0]
        self.products =  [(v, k)  for k,v in c.items() if v > 0] + [(1, k) for k,v in c.items() if v == 0]




        #########   Process the kinetic data   #########

        self.elementary = self._detect_elementary_reaction(kinetics_type)
        if self.elementary:
            kinetics_type = "mass action"


        self.kinetics = Kinetics(law=kinetics_type, parameters=kinetic_parameters)


        #self.reaction_type = self._determine_reaction_type()
        #print(f"detected reaction type `{self.reaction_type}`")

        self.reaction_category = self._determine_reaction_category()
        #print(f"detected reaction category `{self.reaction_category}`")


        self.analytic_solution_family = self._determine_analytic_solution_family()



        #########   Process the thermodynamic data   #########

        self.thermodynamics = ReactionThermodynamics(delta_H=delta_H, delta_S=delta_S,
                                                     K=self.kinetics.parameters.get("K"), temp=temp)

        if temp is not None:
            self.set_thermodynamic_data(temp)




    def _detect_elementary_reaction(self, kinetics_type) -> bool:
        """

        :return:
        """
        if kinetics_type is not None:
            if kinetics_type != "mass action":
                return False

        r, p, c = self.stoichiometry.reaction_pattern()     # number of reactants, products, catalysts

        if c > 0:      # If enzymes were involved
            return False

        if r == 1 and p == 1:
            return True

        if r == 1 and p == 2:
            return True

        if r == 2 and p == 1:
            return True

        return False



    def _determine_reaction_category(self) -> str:
        """

        :return:
        """
        r, p, c = self.stoichiometry.reaction_pattern()     # number of reactants, products, catalysts

        if c > 0:
            return "Enzymatic"

        # TODO: switch to using signed terms
        if (r == 1 and self.reactants[0][0] == 1) and (p == 1 and self.products[0][0] == 1):
            # Reaction is of the type A <-> B               {"A": -1, "B": 1}
            return "Unimolecular rearrangement/isomerization"

        if (r == 1 and self.reactants[0][0] == 1) \
                and (p == 2 and self.products[0][0] == 1  and self.products[1][0] == 1):
            # Reaction is of the type A <-> B + C           {"A": -1, "B": 1, "C": 1}
            return "Unimolecular decomposition"

        if (r == 1 and self.reactants[0][0] == 1) and (p == 1 and self.products[0][0] == 2):
            # Reaction is of the type A <-> 2 B             {"A": -1, "B": 2}
            return "Unimolecular decomposition"

        if (r == 2 and self.reactants[0][0] == 1 and self.reactants[1][0] == 1) \
            and (p == 1 and self.products[0][0] == 1):
            # Reaction is of the type A + B <-> C           {"A": -1, "B": -1, "C": 1}
            return "Bimolecular synthesis"

        if (r == 1 and self.reactants[0][0] == 2) and (p == 1 and self.products[0][0] == 1):
            # Reaction is of the type 2 A <-> C             {"A": -2, "C": 1}
            return "Bimolecular synthesis"

        return "General one-step"



    def _determine_reaction_type(self):
        """

        :return:
        """
        reactant_list = self.reactants
        product_list = self.products

        single_reactant = None
        if len(reactant_list) == 1 and reactant_list[0][0] == 1:    # A single reactant, with stoichiometry 1
            single_reactant = reactant_list[0][1]

        single_product = None
        if len(product_list) == 1 and product_list[0][0] == 1:      # A single product, with stoichiometry 1
            single_product = product_list[0][1]

        reaction_type = "ReactionGeneric"       # Default value, possibly changed below

        if single_reactant:    # A single reactant, with stoichiometry 1
            if single_product:      # A single product, with stoichiometry 1
                reaction_type = "ReactionUnimolecular"
                return reaction_type
            elif len(product_list) == 2 and product_list[0][0] == 1 and product_list[1][0] == 1:      # Two products, both with stoichiometry 1
                reaction_type = "ReactionDecomposition"
                return reaction_type
            elif len(product_list) == 1 and product_list[0][0] == 2:      # A product with stoichiometry 2  (EXAMPLE : A <-> 2 B)
                reaction_type = "ReactionDecomposition"
                return reaction_type
        elif single_product:
            if len(reactant_list) == 2 and reactant_list[0][0] == 1 and reactant_list[1][0] == 1:      # Two reactants, both with stoichiometry 1
                reaction_type = "ReactionSynthesis"
                return reaction_type
            elif len(reactant_list) == 1 and reactant_list[0][0] == 2:  # A reactant with stoichiometry 2  (EXAMPLE : 2A <-> P)
                reaction_type = "ReactionSynthesis"
                return reaction_type

        if reaction_type == "ReactionGeneric":
             return reaction_type



    def get_signed_stoichiometric_coefficients(self, reactants :list[tuple], products :list[tuple]) -> dict:
        """
        Return the sums of all the stoichiometric coefficients for each species in this reaction.
        The reactants get negative values, and the products positive ones

        EXAMPLE: for reaction  A + E -> 2P + Q + E
        it would return {"A": -1, "P": 2, "Q": 1, "E": 0}

        Those signed coefficients ν_i, given a set of species X_i,
        allow the reaction to be expressed as : ∑i ν_i X_i = 0

        :return:    A dictionary mapping the id's of the species in this reaction
                        to their SIGNED stoichiometric coefficients in this reaction
        """
        # TODO: maybe move to class Stoichiometry (and turn it from dataclass to regular class)
        coeffs = {}

        for c, species in reactants:        # Example: (2, "A")
            coeffs[species] = coeffs.get(species, 0) - c    # Accumulate the sum of the stoichiometric coefficients for this species

        for c, species in products:         # Example: (1, "P")
            coeffs[species] = coeffs.get(species, 0) + c    # Accumulate the sum of the stoichiometric coefficients for this species

        return coeffs



    def extract_rxn_properties(self) -> dict:
        """
        Create a dictionary with the numerical properties of the given reaction
        (skipping any None values)
        Possible values include:
            - forward and reverse reaction rates (kR and kR, respectively)
            - ΔH, ΔS, ΔG,
            - K (equilibrium constant)

        :return:    EXAMPLE: {'kF': 3.0, 'kR': 2.0, 'delta_G': -1005.130505, 'K': 1.5}
        """
        thermo_properties = self.thermodynamics.to_dict()
        kinetic_properties = self.kinetics.to_dict()

        return thermo_properties | kinetic_properties   # Combine the two dictionaries



    def set_thermodynamic_data(self, temp :float) -> None:
        """
        Set all the thermodynamic data derivable from the given temperature,
        and all stored kinetic and thermodynamic data.
        Raise an Exception if any inconsistency is detected.

        :param temp:    System temperature in Kelvins.  For now, assumed constant everywhere,
                            and unvarying (or very slowly varying).
                            If the temp gradually changes, periodically call this method.
        :return:        None
        """
        # Process the thermodynamic data, and update various object attributes accordingly
        if temp is not None:
            self.thermodynamics.set_temperature(temp)

        if self.thermodynamics.K is not None:
            self.kinetics.set_rate_constants_from_equilibrium_constant(K=self.thermodynamics.K)



    def extract_intermediate(self) -> str|None:
        """
        Return the name of the reaction intermediate species,
        or None if there's no intermediate

        :return:
        """
        return self.kinetics.extract_intermediate()



    def describe(self, concise=False) -> str:
        """
        Return,  as a string, a user-friendly plain-text form of the reaction

        :param concise:     If True, less detail is shown
        :return:            A string with a description of this reaction
        """
        if self.kinetics.law == "mass action":
            left = self._standard_form_chem_eqn(self.reactants)       # Left side of the equation, as a user-friendly string
            right = self._standard_form_chem_eqn(self.products)       # Right side of the equation

            if self.kinetics.parameters["reversible"]:
                rxn_description = f"{left} <-> {right}"
            else:
                rxn_description = f"{left} -> {right}"

            if concise:
                return rxn_description      # Minimalist description

            # If we get this far, we're looking for a more detailed description
            rxn_description += "  "
            if self.elementary:
                rxn_description += "Elementary "

            rxn_description += self.reaction_category + " reaction\n       "

            rxn_properties = self.extract_rxn_properties()     # A dict
            rxn_description += self.format_reaction_details(rxn_properties)

            return rxn_description


        return "TBA" # TODO: expand



    def extract_reactant_ids(self) -> [str]:
        """
        Return the list of ALL the reactant id's in this reaction
        (including any catalysts, if applicable)

        :return:    A list of unique chemical labels,
                        in the order they appeared in when this reaction was first defined
        """
        return [t[1] for t in self.reactants]


    def extract_reactants(self) -> list[(int, str)]:
        """
        Return a list of pairs with details of the reactants of the given reaction,
        incl. their stoichiometry and species id

        :return:    A list of pairs of the form (stoichiometry coefficient, species id)
        """
        return self.reactants


    def extract_reactants_formula(self) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula
        (aka the reactant "complex")

        :return:    A string with the left (reactant) side of the reaction formula
        """
        return self._standard_form_chem_eqn(self.reactants)



    def extract_product_ids(self) -> [str]:
        """
        Return the list of ALL the product id's in this reaction
        (including any catalysts, if applicable)

        :return:    A list of unique chemical labels,
                        in the order they appeared in when this reaction was first defined
        """
        return [t[1] for t in self.products]


    def extract_products(self) -> list[(int, str)]:
        """
        Return a list of pairs with details of the products of the given reaction,
        incl. their stoichiometry and species id

        :return:    A list of pairs of the form (stoichiometry coefficient, species id)
        """
        return self.products


    def extract_products_formula(self) -> str:
        """
        Return a string with a user-friendly form of the right (products) side of the reaction formula
        (aka the product "complex")

        :return:    A string with the right (product) side of the reaction formula
        """
        return self._standard_form_chem_eqn(self.products)



    def extract_species_in_reaction(self) -> Set[str]:
        """
        Return a SET of the id's of ALL the species appearing in this reaction

        :return:    A SET of the id's of the species involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return set(self.extract_reactant_ids()) | set(self.extract_product_ids())   # Union of sets



    def reaction_quotient(self, conc, explain=False) -> np.double | tuple[np.double, str]:
        """
        Compute the "Reaction Quotient" (aka "Mass–action Ratio"),
        given the concentrations of chemicals involved in this reaction.

        Note: this implementation only covers reactions that have "mass action" kinetics

        :param conc:        Dictionary with the concentrations of the species involved in the reaction.
                            The keys are the chemical labels
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
        :param explain:     If True, it also returns the math formula being used for the computation
                                EXAMPLES:   "([C][D]) / ([A][B])"
                                            "[B] / [A]^2"

        :return:            If explain is False, return value for the "Reaction Quotient" (aka "Mass–action Ratio");
                                if True, return a pair with that quotient and a string with the math formula that was used.
                                Note that the reaction quotient is a Numpy scalar that might be np.inf or np.nan
        """
        assert self.kinetics.law == "mass action", \
            "reaction_quotient(): only implemented for reactions that have \"mass action\" kinetics"

        return ReactionKinetics.compute_reaction_quotient(reactant_data=self.reactants,
                                                        product_data=self.products,
                                                        conc=conc, explain=explain)



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the generic reaction,
        determine its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step.

        :param conc_dict:   A dict mapping specie id's to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """

        if self.kinetics.law == "mass action":
            return ReactionKinetics.compute_rate_elementary(reactants = self.extract_reactant_ids(),
                                                            products = self.extract_product_ids(),
                                                            kF = self.kinetics.parameters["kF"], kR=self.kinetics.parameters["kR"],
                                                            reversible=self.kinetics.parameters["reversible"],
                                                            conc_dict=conc_dict)

        function_to_call = self.kinetics.rate_function
        assert function_to_call is not None, \
            f"determine_reaction_rate(): no kinetic rate function was provide for the reaction `{self.describe(concise=True)}` isn't set; " \
            f"make sure to first call set_rate_function()"
        #print(f"determine_reaction_rate() - function being invoked to determine the reaction's rate: `{function_to_call.__name__}()`")

        return function_to_call(reactant_terms=self.reactants, product_terms=self.products,
                                kF = self.kinetics.parameters["kF"], kR=self.kinetics.parameters["kR"],
                                conc_dict=conc_dict)                        # Carry out the function call



    def step_simulation(self, delta_time, conc_dict :dict, exact=False) -> (dict, float):
        """
        Simulate the generic reaction, over the specified time interval.
        The forward Euler method is used

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentrations won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :param exact:       Only available if this reaction type has a known analytical solution

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn     The mapping of chemical labels
                                                                    to their concentration CHANGES
                                                                    during this step
                                - rxn_rate                      The reaction rate ("velocity") for this reaction
                                EXAMPLE of increment_dict_single_rxn: {"B": -1.3, "F": 2.9, "D": -1.6}
        """
        # TODO: move out of this modules
        increment_dict_single_rxn = {}      # The keys are the species id's,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)


        if exact:
            if self.analytic_solution_family == "ONE_TO_ONE":
                r = self.reactants[0][1]           # EXAMPLE: "R"
                p = self.products[0][1]            # EXAMPLE: "P"

                R0 = conc_dict[r]
                P0 = conc_dict[p]
                # Compute the respective increments of R0 and P0
                if self.kinetics.parameters["reversible"]:
                    delta_p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=self.kinetics.parameters["kF"], kR=self.kinetics.parameters["kR"],
                                                                                     A0=R0, P0=P0, t=delta_time, incremental=True)
                else:
                    delta_p = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=self.kinetics.parameters["kF"],
                                                                                       A0=R0, P0=P0, t=delta_time, incremental=True)

                # Work out the stoichiometry for all the species
                increment_dict_single_rxn = {r: -delta_p, p: delta_p}
                return (increment_dict_single_rxn, rxn_rate)
            else:
                raise Exception("step_simulation(): no exact analytical solution is available for this reaction type")



        # If we get thus far, exact=False

        # In the "forward Euler" approximation, the following rate is taken to remain unvaried during the entire (small) time step
        delta_rxn = rxn_rate * delta_time      # forward reaction - reverse reaction


        reactants = self.reactants      # A list of pairs of the form (stoichiometry coefficient, species id))
        products = self.products        # A list of pairs of the form (stoichiometry coefficient, species id))


        """
        Determine the concentration adjustments as a result of this reaction step, 
        for this individual reaction being considered
        """

        # The reactants DECREASE based on the quantity delta_rxn
        for stoichiometry, species_id in reactants:         # Unpack data from each reactant
            delta_conc = stoichiometry * (- delta_rxn)      # Increment to this reactant from the reaction being considered

            increment_dict_single_rxn[species_id] = increment_dict_single_rxn.get(species_id,0) + delta_conc


        # The reaction products INCREASE based on the quantity delta_rxn
        for stoichiometry, species_id in products:      # Unpack data from each product
            delta_conc = stoichiometry * delta_rxn      # Increment to this reaction product from the reaction being considered

            increment_dict_single_rxn[species_id] = increment_dict_single_rxn.get(species_id,0) + delta_conc


        assert len(increment_dict_single_rxn) == len(self.extract_species_in_reaction())  # TODO: temporary check to eventually drop

        return (increment_dict_single_rxn, rxn_rate)



    def find_equilibrium_conc(self, conc_dict :dict) -> dict:
        """
        Determine the equilibrium concentrations that would be eventually reached
        by the species in this reaction,
        given their initial concentrations,
        IN THE ABSENCE of any other reaction.

        :param conc_dict:   A dict mapping species id's to their initial concentrations,
                                for all the species involved in this reaction
                                EXAMPLE:  {"X": 4.3, "Y": 1.5, "F": 31.6, "G": 3.6}

        :return:            A dict mapping the above chemical id's to their equilibrium concentrations
        """
        reactants = self.reactants
        products = self.products

        if self.kinetics.law != "mass action":
            assert self.kinetics.kinetic_rate_function == ReactionKinetics.compute_rate_first_order, \
                "find_equilibrium_conc(): for reactions that don't exhibit mass-action kinetics, " \
                "it's only implemented when the kinetic rate function is `ReactionKinetics.compute_rate_first_order` \n" \
                "if that's the case, make sure to first invoke:   set_rate_function(ReactionKinetics.compute_rate_first_order)"

            # To conform with ReactionKinetics._compute_equilibrium_conc_first_order(),
            # we'll express the reaction in the form aA + bB <-> bC + dD

            assert len(reactants) == 2, \
                f"find_equilibrium_conc(): for reactions that don't exhibit mass-action kinetics, " \
                f"this function is only implemented for case with 2 reactants (we have {len(reactants)})"

            assert len(products) == 2, \
                f"find_equilibrium_conc(): for reactions that don't exhibit mass-action kinetics, " \
                f"this function is only implemented for case with 2 products (we have {len(products)})"

            r1, r2 = reactants      # Each value is a pair (stoichiometry coefficient, species id)
            p1, p2 = products       # Each value is a pair (stoichiometry coefficient, species id)

            A0 = conc_dict.get(r1[1])
            assert A0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                   f"concentration of the reactant `{r1[1]}` was not provided"

            B0 = conc_dict.get(r2[1])
            assert B0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                   f"concentration of the reactant `{r2[1]}` was not provided"


            C0 = conc_dict.get(p1[1])
            assert C0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                   f"concentration of the product `{p1[1]}` was not provided"

            D0 = conc_dict.get(p2[1])
            assert D0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                   f"concentration of the product `{p2[1]}` was not provided"


            eq_dict = ReactionKinetics._compute_equilibrium_conc_first_order(kF=self.kinetics.parameters["kF"], kR=self.kinetics.parameters["kR"],
                                                                             a=r1[0], b=r2[0],
                                                                             p=p1[0], q=p1[0],
                                                                             A0=A0, B0=B0, P0=C0, Q0=D0)

            # eq_dict contains the keys "A", "B", "P", "Q";
            # translate the standard names A, B, P, Q into the actual names, and also drop any missing term
            return  {r1[1]: eq_dict["A"], r2[1]: eq_dict["B"],
                     p1[1]: eq_dict["P"], p2[1]: eq_dict["Q"]}


        """
        # If we get thus far, we have MASS-ACTION kinetics
        """
        assert len(reactants) <= 2, \
                f"find_equilibrium_conc(): for reactions that exhibit mass-action kinetics, " \
                f"it's only implemented when there are no more than 2 reactants (number of reactants is {len(reactants)})"
        assert len(products) <= 2, \
                f"find_equilibrium_conc(): for reactions that exhibit mass-action kinetics, " \
                f"it's only implemented when there are no more than 2 products (number of products is {len(products)})"


        # To conform with functions available in ReactionKinetics,
        # we'll express the reaction in the form   A + B <-> P + Q   or  2 A <-> P  or  A <-> 2 P   (no other solver is available)

        standard_names = ["A", "B", "P", "Q"]
        coeffs = [0, 0, 0, 0]       # For general reaction A + B <-> P + Q
        concs  = [0., 0., 0., 0.]   # For A0, B0, P0, Q0

        vec_r, vec_p = self.stoichiometry.get_split_reaction_vector()
        #print("reaction vectors (reactants) (products): ", vec_r, vec_p)

        name_map = {}   # To map standard names into actual species names


        # Process reactants (which will set the first element or first two, in coeffs and concs)
        index = 0
        for sp, c in vec_r.items():
            coeffs[index] = -c      # Negative of stoichiometric coefficient because it's a reactant
            conc = conc_dict.get(sp)
            assert conc is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                       f"concentration of the reactant `{sp}` was not provided"
            concs[index] = conc
            std_name = standard_names[index]
            name_map[std_name] = sp
            index += 1

        # Process products (which will set the next element or two, in coeffs and concs)
        index = 2
        for sp, c in vec_p.items():
            coeffs[index] = c
            conc = conc_dict.get(sp)
            assert conc is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                       f"concentration of the product `{sp}` was not provided"
            concs[index] = conc
            std_name = standard_names[index]
            name_map[std_name] = sp
            index += 1


        # Unpack
        a, b, p, q = coeffs
        A0, B0, P0, Q0 = concs

        """
        print(f"coeffs: {coeffs} | concs: {concs} | name_map: {name_map}")
        print(f"a: {a} | b: {b} | p: {p} | q: {q}")
        print(f"A0: {A0} | B0: {B0} | P0: {P0} | Q0: {Q0}")
        print(f"kF: {self.kinetics.parameters['kF']} | kR: {self.kinetics.parameters['kR']}")
        print(self.analytic_solution_family)
        """
        
        if (self.analytic_solution_family == "ONE_TO_ONE"):
            # Reaction is of the form A <-> P
            eq_dict = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=self.kinetics.parameters["kF"],
                                                                            kR=self.kinetics.parameters["kR"],
                                                                            A0=A0, P0=P0)

        elif (self.analytic_solution_family == "TWO_TO_ONE") and (a == 1):
            # Reaction is of the form A + B <-> P
            eq_dict = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=self.kinetics.parameters["kF"],
                                                                            kR=self.kinetics.parameters["kR"],
                                                                            A0=A0, B0=B0, P0=P0)

        elif (self.analytic_solution_family == "TWO_TO_ONE") and (a == 2):
            # Reaction is of the form 2 A <-> P
            eq_dict = ReactionKinetics.compute_equilibrium_conc_elementary_synthesis(kF=self.kinetics.parameters["kF"],
                                                                                     kR=self.kinetics.parameters["kR"],
                                                                                     A0=A0, P0=P0)

        elif (self.analytic_solution_family == "ONE_TO_TWO") and (p == 1):
            # Reaction is of the form A <-> P + Q
            eq_dict = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=self.kinetics.parameters["kF"],
                                                                            kR=self.kinetics.parameters["kR"],
                                                                            A0=A0, P0=P0, Q0=Q0)

        elif (self.analytic_solution_family == "ONE_TO_TWO") and (p == 2):
            # Reaction is of the form A <-> 2 P
            eq_dict = ReactionKinetics.compute_equilibrium_conc_elementary_decomposition(kF=self.kinetics.parameters["kF"],
                                                                                         kR=self.kinetics.parameters["kR"],
                                                                                         A0=A0, P0=P0)

        else:
            raise Exception(f"find_equilibrium_conc(): Not implemented for this reaction type ({self.analytic_solution_family})")

        """                                                                      
        eq_dict = ReactionKinetics._compute_equilibrium_conc_first_order(kF=self.kinetics.parameters["kF"], kR=self.kinetics.parameters["kR"],
                                                                         a=a, b=b, p=c, q=d,
                                                                         A0=A0, B0=B0, P0=C0, Q0=D0)
        """
        #print("eq_dict", eq_dict)

        # eq_dict contains the keys "A", "B", "P", "Q";
        # Translate the standard names A, B, P, Q into the actual species names, and also drop any missing term
        converted_dict = {}
        for k, v in eq_dict.items():
            converted_name = name_map[k]
            converted_dict[converted_name] = v

        return converted_dict






    #####################################################################################################

    '''                                    ~   PRIVATE  ~                                             '''

    def ________PRIVATE________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def _standard_form_chem_eqn(self, eqn_side :list[tuple]) -> str:
        """
        Return a user-friendly form of a "complex" (a side of a chemical equation)

        EXAMPLE:  turn [(1, "Fe"), (2, "Cl")]  into  "Fe + 2 Cl"

        :param eqn_side:    A list encoding either side of a chemical equation
        :return:            A string with a user-friendly form of a side of a chemical equation
        """
        assert type(eqn_side) == list, \
            f"Reaction._standard_form_chem_eqn(): the argument must be a list (it was of type {type(eqn_side)})"

        formula_list = []
        for t in eqn_side:
            stoichiometry, species_name = t
            #stoichiometry = self.extract_stoichiometry(t)
            #species_name = self.extract_species(t)

            if stoichiometry == 1:
                term = species_name
            else:
                term = f"{stoichiometry} {species_name}"

            formula_list.append(term)

        return " + ".join(formula_list)



    def _determine_analytic_solution_family(self) -> str|None:
        """

        :return:
        """
        if self.kinetics.law != "mass action":
            return None

        r, p, c = self.stoichiometry.reaction_pattern()     # number of reactants, products, catalysts

        if c > 0:      # If enzymes were involved
            return None

        # TODO: switch to using signed terms
        if (r == 1 and self.reactants[0][0] == 1) and (p == 1 and self.products[0][0] == 1):
            # Reaction is of the type A <-> B               {"A": -1, "B": 1}
            return "ONE_TO_ONE"

        if (r == 1 and self.reactants[0][0] == 1) \
                and (p == 2 and self.products[0][0] == 1  and self.products[1][0] == 1):
            # Reaction is of the type A <-> B + C           {"A": -1, "B": 1, "C": 1}
            return "ONE_TO_TWO"

        if (r == 1 and self.reactants[0][0] == 1) and (p == 1 and self.products[0][0] == 2):
            # Reaction is of the type A <-> 2 B             {"A": -1, "B": 2}
            return "ONE_TO_TWO"

        if (r == 2 and self.reactants[0][0] == 1 and self.reactants[1][0] == 1) \
            and (p == 1 and self.products[0][0] == 1):
            # Reaction is of the type A + B <-> C           {"A": -1, "B": -1, "C": 1}
            return "TWO_TO_ONE"

        if (r == 1 and self.reactants[0][0] == 2) and (p == 1 and self.products[0][0] == 1):
            # Reaction is of the type 2 A <-> C             {"A": -2, "C": 1}
            return "TWO_TO_ONE"

        return None



    def format_reaction_details(self, rxn_properties :dict) -> str:
        """
        Format and return a string with some details about the parameters of this reaction,
        contained in the passed dictionary.
        Any property named "temp" gets converted from degree K to C.

        :param rxn_properties:  A dictionary with numerical properties of interest for the reaction
                                    EXAMPLE: {'kF': 3.0, 'kR': 2.0, 'delta_G': 1.2345, 'K': 1.5}

        :return:                A string with some details about the parameters of this reaction
                                    EXAMPLE: "  (kF = 3 | kR = 2 | delta_G = 1.2345 kJ/mol | Temp = 25 C)"
        """
        print("rxn_properties: ", rxn_properties)
        details = []    # Running list of strings with each of the individual details

        for k,v in rxn_properties.items():
            if k == "temp":
                single_detail = f"Temp = {convert(v, from_unit=K, to_unit=C):,.4g} C"
                # EXAMPLE: "Temp = 25 C"
                details.append(single_detail)
                continue

            if type(v) == str:
                single_detail = f"{k} = '{v}'"
            elif type(v) == bool:
                single_detail = f"{k} = {v}"
            else:
                single_detail = f"{k} = {v:,.5g}"   # EXAMPLES: "kF = 3"
                                                    #           "delta_G = 1.2345"
            units = show_standard_units(k)
            if units is not None:
                single_detail += " " + units        # EXAMPLE: "delta_G = 1.2345 kJ/mol"

            details.append(single_detail)


        description = ""

        #if temp:
        #    details.append(f"Temp = {convert(temp, from_unit=K, to_unit=C):,.4g} C")          # EXAMPLE: "Temp = 25 C"

        if details:     # If there is any data
            description = "  (" + ' | '.join(details) + ")"   # EXAMPLE: "  (kF = 3 | kR = 2 | delta_G = 1.2345 kJ/mol)"

        return description