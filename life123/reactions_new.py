from __future__ import annotations      # To facilitate type annotations
import numpy as np
from typing import Union, Set, Tuple, Mapping
from dataclasses import dataclass, field, fields, asdict
from life123.thermodynamics import ThermoDynamics
from life123.reaction_kinetics import ReactionKinetics
from life123.species_registry import Species, SpeciesRegistry
from life123.units import show_standard_units, convert, K, C




@dataclass(frozen=True)
class Stoichiometry:
    """
    Mapping species `id` to its signed stoichiometric coefficient,
    for all species in the reaction
    (including species that have a net coefficient of zero, such as catalysts that appear on both sides)
    """
    coefficients: Mapping[str, int]



    def reaction_pattern(self) -> tuple[int, int, int]:
        """
        Return the triplet (n_reactants, n_products, n_catalysts)
        Catalysts, if present, are NOT counted under either "reactants" or "products"
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

        if law == "mass action":
            if parameters.get("reversible") is None:
                parameters["reversible"] = True      # Set  default

            kF = parameters.get("kF")
            kR = parameters.get("kR")

            parameters["kF"] = kF   # Might be None
            parameters["kR"] = kR   # Might be None

            if parameters.get("reversible"):
                # If we have a reversible reaction that follows mass-action kinetics
                if (kF is not None) and (kR is not None) and (not np.allclose(kR, 0)):
                    parameters["K"] = kF / kR    # Kinetic parameter ratio

            else:
                # If we have an IR-reversible reaction that follows mass-action kinetics
                assert not parameters.get("kR"), \
                    f"Irreversible reactions with mass-action kinetics " \
                    f"cannot have a value for the reverse rate constant (kR = {parameters.get('kR')})"
                parameters["kR"] = 0


        self.law = law
        self.parameters = parameters



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

        if (kR is None) and (kF is not None):
            self.parameters["kR"] = kF / K
            return

        if (kF is None) and (kR is not None):
            self.parameters["kF"] = K * kR
            return

        if (kF is not None) and (kR is not None):
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

    def __init__(self, reactants :str|list|tuple, products :str|list|tuple, species_registry :SpeciesRegistry, name=None,
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

        self.analytic_solution_family = None
        self.elementary = None
        self.reaction_category = None


        assert reactants is not None, "ReactionGeneric(): the argument `reactants` is a required one; it can't be None"
        if type(reactants) == str:
            reactants = [reactants]

        assert products is not None, "ReactionGeneric(): the argument `products` is a required one; it can't be None"
        if type(products) == str:
            products = [products]


        # Normalize the elements of each list to be (int, str) pairs
        reactant_list = [(1, r) if type(r) == str else r
                            for r in reactants]   # A list of pairs
        product_list =  [(1, p) if type(p) == str else p
                            for p in products]   # A list of pairs

        # Catch identical reaction sides, even if terms are reshuffled
        assert set(reactant_list) != set(product_list), \
            f"Reaction(): the two sides of the reaction can't be identical! " \
            f"Same reactant and product complexes: \"{self._standard_form_chem_eqn(reactant_list)}\""


        # Verify that all the species in the reaction are registered ones
        for _, s_id in reactant_list:
            assert species_registry.species_exists(s_id), \
                f"No species with id \"{s_id}\" exists in the species registry"

        for _, s_id in product_list:
            assert species_registry.species_exists(s_id), \
                f"No species with id \"{s_id}\" exists in the species registry"


        self.reactants = reactant_list
        self.products = product_list


        c = self.get_signed_stoichiometric_coefficients()
        self.stoichiometry = Stoichiometry(c)



        #########   Process the kinetic data   #########

        self.elementary = self._detect_elementary_reaction(kinetics_type)
        if self.elementary:
            kinetics_type = "mass action"


        self.kinetics = Kinetics(law=kinetics_type, parameters=kinetic_parameters)


        self.reaction_type = self._determine_reaction_type()
        print(f"detected reaction type `{self.reaction_type}`")

        self.reaction_category = self._determine_reaction_category()
        print(f"detected reaction category `{self._determine_reaction_category}`")


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

        if r == 1 and p == 1:
            return "Unimolecular rearrangement/isomerization"

        if r == 1 and p == 2:
            return "Unimolecular decomposition"

        if r == 2 and p == 1:
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



    def get_signed_stoichiometric_coefficients(self) -> dict:
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
        coeffs = {}

        for c, species in self.reactants:
            coeffs[species] = coeffs.get(species, 0) - c    # Accumulate the sum of the stoichiometric coefficients for this species

        for c, species in self.products:
            coeffs[species] = coeffs.get(species, 0) + c    # Accumulate the sum of the stoichiometric coefficients for this species

        return coeffs



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
                for k,v in self.get_signed_stoichiometric_coefficients().items()
                if v != 0}

        return d



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

            rxn_description += self.reaction_category + "\n       "

            rxn_properties = self.extract_rxn_properties()     # A dict
            rxn_description += self.format_reaction_details(rxn_properties)

            return rxn_description


        return "TBA" # TODO: expand




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

        if r == 1 and p == 1:
            return "ONE_TO_ONE"

        if r == 1 and p == 2:
            return "ONE_TO_TWO"

        if r == 2 and p == 1:
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