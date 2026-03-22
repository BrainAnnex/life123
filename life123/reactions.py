# classes: "ReactionCommon", "ReactionOneStep", "ReactionUnimolecular",
#          "ReactionSynthesis", "ReactionDecomposition", "ReactionEnzyme",
#          and "ReactionGeneric"

from typing import Union, Set, Tuple
import numpy as np
from life123.thermodynamics import ThermoDynamics
from life123.reaction_kinetics import ReactionKinetics



###################################################################################################################

class ReactionCommon:
    """
    Foundational base class for ALL individual reactions.

    Typically NOT instantiated by the user.

    Individual reactions classes need to provide the following methods:

        extract_stoichiometry()
        extract_chem_label()

        reaction_details()
        extract_intermediate()
        extract_rxn_properties()
        set_thermodynamic_data()
        extract_forward_rate()
        extract_reverse_rate()
        extract_equilibrium_constant()

        describe()

        extract_reactant_labels()
        extract_reactants()

        extract_product_labels()
        extract_products()

        extract_chemicals_in_reaction()

        extract_rxn_properties()     A dict

        step_simulation()
    """

    def __init__(self, active=True, temp=None):
        """
        :param active:  [NOT YET IMPLEMENTED] If False, the reaction is regarded as inactive
        :param temp:    In Kelvins
        """
        self.active = active        # TODO: not yet in use
        self.temp = temp            # In Kelvins




    def reaction_details(self, rxn_properties :dict) -> str:
        """
        Return a string with some details about the parameters of this reaction

        :param rxn_properties:  A dictionary with numerical properties of interest for the reaction
                                    EXAMPLE: {'kF': 3.0, 'kR': 2.0, 'delta_G': -1005.13, 'K': 1.5}

        :return:                A string with some details about the parameters of this reaction
                                    EXAMPLE: "  (kF = 3 / kR = 2 / Delta_G = -1,005.13 / Temp = 25 C)"
        """
        details = []
        #rxn_properties = self.extract_rxn_properties()
        for k,v in rxn_properties.items():
            details.append(f"{k} = {v:,.5g}")          # EXAMPLE: "kF = 3"

        description = ""

        if self.temp:
            details.append(f"Temp = {self.temp - 273.15:,.4g} C")          # EXAMPLE: "Temp = 25 C"

        if details:
            description = "  (" + ' / '.join(details) + ")"   # EXAMPLE: "  (kF = 3 / kR = 2 / Delta_G = -1,005.13)"

        return description



    def extract_stoichiometry(self, term :(int, str)) -> int:
        """
        Return the stoichiometry coefficient, from a reaction "TERM"

        :param term:    A pair (int, str) representing a reaction term
        :return:        An integer with the stoichiometry coefficient
        """
        try:
            return term[0]
        except Exception:
            raise Exception(f"extract_stoichiometry(): the argument must be a pair consisting of an integer and a string; "
                            f"the passed value was {term}")


    def extract_chem_label(self, term :(int, str)) -> str:
        """
        Return the label of the chemical species, from a reaction "TERM"

        :param term:    A pair (int, str) representing a reaction term
        :return:        The name of the chemical species in the term
        """
        try:
            return term[1]
        except Exception:
            raise Exception(f"extract_chem_label(): the argument must be a pair consisting of an integer and a string; "
                            f"the passed value was {term}")





###################################################################################################################


class ReactionElementary(ReactionCommon):
    """
    Base class for all reactions that can be modeled kinetically as happening in 1 step
    (i.e. with no intermediaries).

    Typically NOT instantiated by the user.
    """
    def __init__(self, reversible=True, kF=None, kR=None,
                 delta_H=None, delta_S=None, delta_G=None, **kwargs):
        """
        :param reversible:
        :param kF:
        :param kR:
        :param delta_H:
        :param delta_S:
        :param delta_G:
        :param kwargs:
        """

        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        self.reversible = reversible

        if not self.reversible:
            assert not kR, \
                f"ReactionElementary instantiation: irreversible reactions " \
                f"cannot have a value for the reverse rate constant (kR = {kR})"
            kR = 0

        self.kF = kF                # Forward rate constant
        self.kR = kR                # Reverse rate constant
        self.K = None               # Forward–reverse rate constant ratio (Kinetic parameter ratio)

        # Process the kinetic  data
        self.equilibrium_constant_from_kinetic_data(kF=kF, kR=kR) # This will set self.K if possible

        self.delta_H = delta_H
        self.delta_S = delta_S
        self.delta_G = delta_G

        if self.temp is not None:
            # Process the thermodynamic data, and update various object attributes accordingly
            thermo_data = ThermoDynamics.extract_thermodynamic_data(K=self.K,
                                                                delta_H=delta_H, delta_S=delta_S, delta_G=delta_G,
                                                                temp=self.temp)
            #print(f"thermo_data : {thermo_data}")
            self.K = thermo_data["K"]
            self.delta_H = thermo_data["delta_H"]
            self.delta_S = thermo_data["delta_S"]
            self.delta_G = thermo_data["delta_G"]

            self.set_rate_constants_from_equilibrium_constant(K=self.K)




    def equilibrium_constant_from_kinetic_data(self, K=None, kF=None, kR=None) -> None:
        """
        True for Elementary reactions
        (and, more generally, for any reaction that follows mass-action kinetics)

        :param K:   The reaction's equilibrium constant
        :param kF:  The forward reaction's reaction rate constant
        :param kR:  The reverse reaction's reaction rate constant
        :return:    None
        """
        #print(f"In equilibrium_constant_from_kinetic_data() : K={K}, kF={kF}, kR={kR}")
        if (self.K is None) and (kF is not None) and (kR is not None) and (not np.allclose(self.kR, 0)):
            self.K = kF / kR
            return

        if (self.kR is None) and (kF is not None) and (K is not None):
            self.kR = kF / K
            return

        if (self.kF is None) and (K is not None) and (kR is not None):
            self.kF = K * kR



    def set_rate_constants_from_equilibrium_constant(self, K) -> None:
        """
        Set, as needed, a missing reaction rate constant (kF or kR)
        from the other one and the given equilibrium constant K.
        If all values already exist, and an inconsistency is detected, an Exception will be raised.

        Note: the reaction's equilibrium constant and its kinetic rate constants are
              in the relationship K = kF / kR for any
              (and, more generally, for any reaction that follows mass-action kinetics)

        :param K:   The reaction's equilibrium constant
        :return:    None
        """
        if self.K is None:
            return

        if (self.kR is None) and (self.kF is not None):
            self.kR = self.kF / K
            return

        if (self.kF is None) and (self.kR is not None):
            self.kF = K * self.kR
            return

        if (self.kF is not None) and (self.kR is not None):
            assert np.allclose(K, self.kF / self.kR), \
                f"set_rate_constants_from_equilibrium_constant(): values for kR and kR already exist, " \
                f"and are inconsistent with the passed value of K ({K})"



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



    def set_thermodynamic_data(self, temp :float) -> None:
        """
        Set all the thermodynamic data derivable from the given temperature,
        and all previously passed kinetic and thermodynamic data.
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

        self.set_rate_constants_from_equilibrium_constant(K=self.K)



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


    def extract_equilibrium_constant(self) -> float:
        """

        :return:    The value of the equilibrium constant for this reaction
        """
        return self.K


    def extract_intermediate(self) -> str | None:
        """
        Return the name of the reaction intermediate species,
        or None if there's no intermediate
        :return:
        """
        return None   # There are no intermediates






###################################################################################################################

class ReactionUnimolecular(ReactionElementary):
    """
    Elementary reaction of type A <-> P, of first order in `A` and `P`
    """

    def __init__(self, reactant :str, product :str, **kwargs):
        """
        :param reactant:    Label of the reaction's single reactant
        :param product:     Label of the reaction's single product (cannot be the same as the reactant)
        :param kwargs:      Other named arguments to pass thru to the parent class
        """
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        assert type(reactant) == str, "ReactionUnimolecular instantiation: argument `reactant` must be a string"
        assert type(product) == str, "ReactionUnimolecular instantiation: argument `product` must be a string"
        assert reactant != product, \
            f"ReactionUnimolecular instantiation: the `reactant` and the `product` cannot be identical (`{reactant}`)"

        self.reactant = reactant
        self.product = product




    def describe(self, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction

        :param concise:     If True, less detail is shown
        :return:            A string with a description of this reaction
        """
        description = f"{self.reactant} <-> {self.product}"

        if not concise:
            if self.reversible:
                description += "  Elementary Unimolecular reaction"
            else:
                description += "  Elementary Unimolecular Irreversible reaction"

            description += self.reaction_details(self.extract_rxn_properties())

        return description



    def extract_reactant_labels(self) -> [str]:
        """
        Return the list of ALL the reactant labels in this reaction

        :return:    A list of ALL the reactant labels in this reaction
        """
        return [self.reactant]


    def extract_reactants(self) -> [(int, str)]:
        """
        Return a list of pairs with details of the reactants of the given reaction,
        incl. their stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
        """
        return [(1, self.reactant)]


    def extract_reactants_formula(self) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        return self.reactant



    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels ALL the products of this reaction

        :return:
        """
        return [self.product]


    def extract_products(self) -> [(int, str)]:
        """
        Return a list of pairs with details of the products of the given reaction,
        incl. their stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
        """
        return [(1, self.product)]


    def extract_products_formula(self) -> str:
        """
        Return a string with a user-friendly form of the right (product) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        return self.product



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the labels of ALL the chemicals appearing in this reaction

        :return:    A SET of the labels of the chemicals involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return {self.reactant, self.product}



    def reaction_quotient(self, conc, explain=False) -> np.double | Tuple[np.double, str]:
        """
        Compute the "Reaction Quotient" (aka "Mass–action Ratio"),
        given the concentrations of chemicals involved in this reaction

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
        return ThermoDynamics.compute_reaction_quotient(reactant_data=[(1, self.reactant)], product_data=[(1, self.product)],
                                                        conc=conc, explain=explain)



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the unimolecular reaction,
        determine its current reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate".

        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """
        return ReactionKinetics.compute_rate_elementary(reactants = [self.reactant], products=[self.product],
                                                        kF = self.kF, kR=self.kR, reversible=self.reversible,
                                                        conc_dict=conc_dict)



    def step_simulation(self, delta_time, conc_dict :dict, exact=False) -> (dict, float):
        """
        Simulate the unimolecular reaction R <-> P, over the specified time interval,
        using either the exact analytical solution, or the "Forward Euler" approximation method

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentrations won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"R": 1.5, "P": 31.6}
        :param exact:       [OPTIONAL] If True, use the exact analytical solution;
                                if False (default), use the "Forward Euler" approximation method

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical label to their concentration CHANGES
                                                            during this step
                                                            EXAMPLE: {"R": -0.2, "P": 0.2}
                                                                     (meaning [R] decreases by 0.2 and [P] increases by the same amount)

                                - rxn_rate                  is the reaction rate ("velocity") for this reaction
                                                            (rate of change of the product), at the START of the simulation step
                                                            EXAMPLE: 3.5
        """
        increment_dict_single_rxn = {}      # The keys are the chemical labels,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        r = self.reactant           # EXAMPLE: "R"
        p = self.product            # EXAMPLE: "P"

        if exact:
            R0 = conc_dict[r]
            P0 = conc_dict[p]
            # Compute the respective increments of R0 and P0
            if self.reversible:
                delta_p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=self.kF, kR=self.kR,
                                                                                 A0=R0, P0=P0, t=delta_time, incremental=True)
            else:
                delta_p = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=self.kF,
                                                                                   A0=R0, P0=P0, t=delta_time, incremental=True)

            increment_dict_single_rxn = {r: -delta_p, p: delta_p}
            return (increment_dict_single_rxn, rxn_rate)


        # If we get thus far, exact=False

        # In the "forward Euler" approximation, the following rate is taken to remain unvaried during the entire (small) time step
        delta_rxn = rxn_rate * delta_time   # forward reaction - reverse reaction

        # Determine the concentration adjustments as a result of this reaction step

        # The reactant DECREASES based on the quantity delta_rxn
        delta_conc = - delta_rxn    # Increment to this reactant from the reaction step
        increment_dict_single_rxn[r] = delta_conc


        # The reaction product INCREASES based on the quantity delta_rxn
        delta_conc = delta_rxn      # Increment to this reaction product from the reaction step
        increment_dict_single_rxn[p] = delta_conc

        return (increment_dict_single_rxn, rxn_rate)



    def find_equilibrium_conc(self, conc_dict :dict) -> dict:
        """
        Determine the equilibrium concentrations that would be eventually reached
        by the chemicals of this reaction,
        given the specified concentrations,
        IN THE ABSENCE of any other reaction.

        :param conc_dict:   A dict mapping chemical labels to their initial concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6}

        :return:            A dict mapping the above chemical labels to their equilibrium concentrations
        """
        # To conform with ReactionKinetics.compute_equilibrium_conc,
        # we'll express the reaction in the form A <-> C
        A0 = conc_dict.get(self.reactant)
        assert A0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the reactant `{self.reactant}` was not provided"

        C0 = conc_dict.get(self.product)
        assert C0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the produce `{self.product}` was not provided"

        eq_dict = ReactionKinetics.compute_equilibrium_conc_first_order(kF=self.kF, kR=self.kR,
                                                                        a=1, p=1,
                                                                        A0=A0, P0=C0)

        # eq_dict contains the keys "A", "B", "P", "Q";
        # translate the standard names A, B, P, Q into the actual names, and also drop any missing term
        return  {self.reactant: eq_dict["A"], self.product: eq_dict["P"]}





#######################################################################################################################

class ReactionSynthesis(ReactionElementary):
    """
    Bimolecular reactions of type A + B <-> P,
    of first order for each participating chemical.
    """

    def __init__(self, reactants :(str, str), product :str, **kwargs):
        """
        :param reactants:   Pair with the labels of the 2 reactants (possibly the same)
        :param product:     Label of the reaction's single product (cannot be the same as any of the reactants)
        :param kwargs:      Other named arguments to pass thru to the parent class
        """
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        assert type(reactants) == list or type(reactants) == tuple, \
            "ReactionSynthesis instantiation: argument `reactants` must be a list or tuple"
        assert len(reactants) == 2, \
            "ReactionSynthesis instantiation: argument `reactants` must be a pair"
        assert type(product) == str, "ReactionSynthesis instantiation: argument `product` must be a string"

        (r1, r2) = reactants
        assert (r1 != product) and (r2 != product), \
            "ReactionSynthesis instantiation: the `product` cannot be identical to any of the reactants"


        self.reactant_1 = r1
        self.reactant_2 = r2
        self.product = product




    def describe(self, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction

        :param concise:     If True, less detail is shown
        :return:            A string with a description of this reaction
                            EXAMPLES of concise version:
                                A + B <-> P
                                2 C <-> P
        """
        description = f"{self.extract_reactants_formula()}  <-> {self.product}"

        if not concise:
            if self.reversible:
                description += "  Elementary Synthesis reaction"
            else:
                description += "  Elementary Synthesis Irreversible reaction"

            description += self.reaction_details(self.extract_rxn_properties())

        return description



    def extract_reactant_labels(self) -> [str]:
        """
        Return the list of ALL the reactant labels in this reaction

        :return:    The list of ALL the reactant labels in this reaction
        """
        return [self.reactant_1, self.reactant_2]


    def extract_reactants(self) -> [(int, str)]:
        """
        Return a list of pairs with details of the reactants of the given reaction,
        incl. their stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
        """
        if self.reactant_1 != self.reactant_2:
            return [(1, self.reactant_1) , (1, self.reactant_2)]

        return [(2, self.reactant_1)]


    def extract_reactants_formula(self) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        if self.reactant_1 != self.reactant_2:
            return f"{self.reactant_1} + {self.reactant_2}"

        return f"2 {self.reactant_1}"



    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels ALL the products of this reaction

        :return:
        """
        return [self.product]


    def extract_products(self) -> [(int, str)]:
        """
        Return a list of pairs with details of the products of the given reaction,
        incl. their stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
        """
        return [(1, self.product)]


    def extract_products_formula(self) -> str:
        """
        Return a string with a user-friendly form of the right (product) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        return self.product



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the labels of ALL the chemicals appearing in this reaction

        :return:    A SET of the labels of the chemicals involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return {self.reactant_1, self.reactant_2, self.product}



    def reaction_quotient(self, conc, explain=False) -> np.double | Tuple[np.double, str]:
        """
        Compute the "Reaction Quotient" (aka "Mass–action Ratio"),
        given the concentrations of chemicals involved in this reaction

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
        return ThermoDynamics.compute_reaction_quotient(reactant_data=[(1, self.reactant_1) , (1, self.reactant_2)],
                                                        product_data= [(1, self.product)],
                                                        conc=conc, explain=explain)



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the synthesis reaction,
        determine its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step.

        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """
        return ReactionKinetics.compute_rate_elementary(reactants = [self.reactant_1, self.reactant_2],
                                                        products=[self.product],
                                                        kF = self.kF, kR=self.kR, reversible=self.reversible,
                                                        conc_dict=conc_dict)



    def step_simulation(self, delta_time, conc_dict :dict, exact=False) -> (dict, float):
        """
        Simulate the synthesis reaction A + B <-> C, over the specified time interval,
        using the "Forward Euler" method

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentrations won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"A": 1.5, "B": 31.6, "C": 19.9}
        :param exact:       [OPTIONAL] If True, use the exact analytical solution;
                                if False (default), use the "Forward Euler" approximation method

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical labels to their concentration CHANGES
                                                            during this step
                                - rxn_rate                  is the reaction rate ("velocity") for this reaction
                                                            (rate of change of the product)
                                EXAMPLE of increment_dict_single_rxn: "A": -1.3, "B": 2.9, "C": -1.6
        """
        increment_dict_single_rxn = {}      # The keys are the chemical labels,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        if exact and self.reversible:   # TODO: take care of the ir-reversible case
            A0 = conc_dict[self.reactant_1]    # TODO: Look into whether this approach will work for 2A -> C
            B0 = conc_dict[self.reactant_2]
            C0 = conc_dict[self.product]

            increment_triplet = ReactionKinetics.exact_advance_synthesis_reversible(kF=self.kF, kR=self.kR,
                                                                                    A0=A0, B0=B0, P0=C0, t=delta_time, incremental=True)

            increment_dict_single_rxn = {self.reactant_1: increment_triplet[0], self.reactant_2: increment_triplet[1],
                                         self.product: increment_triplet[2]}
            return (increment_dict_single_rxn, rxn_rate)


        # If we get thus far, exact=False

        # In the "forward Euler" approximation, the following rate is taken to remain unvaried during the entire (small) time step
        delta_rxn = rxn_rate * delta_time   # forward reaction - reverse reaction


        # Determine the concentration adjustments as a result of this reaction step:

        # The reactants DECREASE based on the quantity delta_rxn
        if self.reactant_1 != self.reactant_2:
            for r in [self.reactant_1, self.reactant_2]:
                # EXAMPLE of r: "A"
                # stoichiometry = 1
                delta_conc = - delta_rxn    # Increment to this reactant from the reaction step
                increment_dict_single_rxn[r] = delta_conc
        else:
            # The left side of the reaction is of the form 2 R
            # stoichiometry = 2
            delta_conc = - 2 * delta_rxn    # Increment to this reactant from the reaction step
            increment_dict_single_rxn[self.reactant_1] = delta_conc


        # The reaction product INCREASES based on the quantity delta_rxn
        p = self.product            # EXAMPLE: "F"
        # stoichiometry = 1
        delta_conc = delta_rxn      # Increment to this reaction product from the reaction step
        increment_dict_single_rxn[p] = delta_conc

        return (increment_dict_single_rxn, rxn_rate)



    def find_equilibrium_conc(self, conc_dict :dict) -> dict:
        """
        Determine the equilibrium concentrations that would be eventually reached
        by the chemicals of this reaction,
        given the specified concentrations,
        IN THE ABSENCE of any other reaction.

        :param conc_dict:   A dict mapping chemical labels to their initial concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6}

        :return:            A dict mapping the above chemical labels to their equilibrium concentrations
        """
        # To conform with ReactionKinetics.compute_equilibrium_conc,
        # we'll express the reaction in the form A + B <-> C
        A0 = conc_dict.get(self.reactant_1)
        assert A0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                       f"concentration of the reactant `{self.reactant_1}` was not provided"

        B0 = conc_dict.get(self.reactant_2)
        assert B0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the reactant `{self.reactant_2}` was not provided"

        C0 = conc_dict.get(self.product)
        assert C0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the produce `{self.product}` was not provided"

        if self.reactant_1 == self.reactant_2:
            eq_dict = ReactionKinetics.compute_equilibrium_conc_first_order(kF=self.kF, kR=self.kR,
                                                                            a=2, b=2, p=1,
                                                                            A0=A0, B0=B0, P0=C0)
        else:
            eq_dict = ReactionKinetics.compute_equilibrium_conc_first_order(kF=self.kF, kR=self.kR,
                                                                            a=1, b=1, p=1,
                                                                            A0=A0, B0=B0, P0=C0)

        # eq_dict contains the keys "A", "B", "P", "Q";
        # translate the standard names A, B, P, Q into the actual names, and also drop any missing term
        return  {self.reactant_1: eq_dict["A"], self.reactant_2: eq_dict["B"], self.product: eq_dict["P"]}






#######################################################################################################################

class ReactionDecomposition(ReactionElementary):
    """
    Bimolecular reactions of type A <-> B + C
    of first order for each participating chemical.
    """

    def __init__(self, reactant :str, products :(str, str), **kwargs):
        """
        :param reactant:    Label of the reaction's single reactant (cannot be the same as any of the products)
        :param products:    Pair with the labels of the 2 reaction products (possibly the same)
        :param kwargs:      Other named arguments to pass thru to the parent class
        """
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        assert type(products) == list or type(products) == tuple, \
            "ReactionDecomposition instantiation: argument `reactants` must be a list or tuple"
        assert len(products) == 2, \
            "ReactionDecomposition instantiation: argument `reactants` must be a pair"
        assert type(reactant) == str, \
            "ReactionDecomposition instantiation: argument `product` must be a string"

        (p1, p2) = products
        assert (p1 != reactant) and (p2 != reactant), \
            "ReactionDecomposition instantiation: the `reactant` cannot be identical to any of the reaction products"


        self.reactant = reactant
        self.product_1 = p1
        self.product_2 = p2




    def describe(self, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction

        :param concise:     If True, less detail is shown
        :return:            A string with a description of this reaction
                            EXAMPLES of concise version:
                                A <-> P + Q
                                D <-> 2 P
        """
        description = f"{self.reactant} <-> {self.extract_products_formula()}"

        if not concise:
            if self.reversible:
                description += "  Elementary Decomposition reaction"
            else:
                description += "  Elementary Decomposition Irreversible reaction"

            description += self.reaction_details(self.extract_rxn_properties())

        return description



    def extract_reactant_labels(self) -> [str]:
        """
        Return the list of ALL the reactant labels in this reaction

        :return:
        """
        return [self.reactant]


    def extract_reactants(self) -> [(int, str)]:
        """
        Return a list of pairs with details of the reactants of the given reaction,
        incl. their stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
        """
        return [(1, self.reactant)]


    def extract_reactants_formula(self) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        return self.reactant



    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels ALL the products of this reaction

        :return:
        """
        return [self.product_1, self.product_2]


    def extract_products(self) -> [(int, str)]:
        """
        Return a list of pairs with details of the products of the given reaction,
        incl. their stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
        """
        if self.product_1 != self.product_2:
            return [(1, self.product_1), (1, self.product_2)]

        return [(2, self.product_1)]


    def extract_products_formula(self) -> str:
        """
        Return a string with a user-friendly form of the right (product) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        if self.product_1 != self.product_2:
            return f"{self.product_1} + {self.product_2}"

        return f"2 {self.product_1}"



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the labels of ALL the chemicals appearing in this reaction

        :return:    A SET of the labels of the chemicals involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return {self.reactant, self.product_1, self.product_2}



    def reaction_quotient(self, conc, explain=False) -> np.double | Tuple[np.double, str]:
        """
        Compute the "Reaction Quotient" (aka "Mass–action Ratio"),
        given the concentrations of chemicals involved in this reaction

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
        return ThermoDynamics.compute_reaction_quotient(reactant_data=[(1, self.reactant)],
                                                        product_data=[(1, self.product_1) , (1, self.product_2)],
                                                        conc=conc, explain=explain)



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the decomposition reaction,
        determine its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step.

        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """
        return ReactionKinetics.compute_rate_elementary(reactants = [self.reactant],
                                                        products=[self.product_1, self.product_2],
                                                        kF = self.kF, kR=self.kR, reversible=self.reversible,
                                                        conc_dict=conc_dict)



    def step_simulation(self, delta_time, conc_dict :dict, exact=False) -> (dict, float):
        """
        Simulate the decomposition reaction A <-> B + C, over the specified time interval,
        using the "Forward Euler" method

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentrations won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"A": 1.5, "B": 31.6, "C": 19.9}

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical labels to their concentration CHANGES
                                                            during this step
                                - rxn_rate                  is the reaction rate ("velocity") for this reaction
                                                            (rate of change of either of the products)
                                EXAMPLE of increment_dict_single_rxn: "A": -1.3, "B": 2.9, "C": -1.6
        """
        increment_dict_single_rxn = {}      # The keys are the chemical labels,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        delta_rxn = rxn_rate * delta_time      # forward reaction - reverse reaction


        # Determine the concentration adjustments as a result of this reaction step:

        # The reactant DECREASES based on the quantity delta_rxn
        r = self.reactant            # EXAMPLE: "F"
        # stoichiometry = 1
        delta_conc = - delta_rxn    # Increment to this reactant from the reaction step
        increment_dict_single_rxn[r] = delta_conc


        # The reaction products INCREASE based on the quantity delta_rxn
        if self.product_1 != self.product_2:
            for p in [self.product_1, self.product_2]:
                # EXAMPLE of p: "B"
                # stoichiometry = 1
                delta_conc = delta_rxn      # Increment to this reaction product from the reaction step
                increment_dict_single_rxn[p] = delta_conc
        else:
            # The right side of the reaction is of the form 2 P
            # stoichiometry = 2
            delta_conc = 2 * delta_rxn      # Increment to this reaction product from the reaction step
            increment_dict_single_rxn[self.product_1] = delta_conc


        return (increment_dict_single_rxn, rxn_rate)



    def find_equilibrium_conc(self, conc_dict :dict) -> dict:
        """
        Determine the equilibrium concentrations that would be eventually reached
        by the chemicals of this reaction,
        given the specified concentrations,
        IN THE ABSENCE of any other reaction.

        :param conc_dict:   A dict mapping chemical labels to their initial concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6}

        :return:            A dict mapping the above chemical labels to their equilibrium concentrations
        """
        # To conform with ReactionKinetics.compute_equilibrium_conc,
        # we'll express the reaction in the form A <-> C + D
        A0 = conc_dict.get(self.reactant)
        assert A0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                       f"concentration of the reactant `{self.reactant}` was not provided"

        C0 = conc_dict.get(self.product_1)
        assert C0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the reactant `{self.product_1}` was not provided"

        D0 = conc_dict.get(self.product_2)
        assert D0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the produce `{self.product_2}` was not provided"

        if self.product_1 == self.product_2:
            eq_dict = ReactionKinetics.compute_equilibrium_conc_first_order(kF=self.kF, kR=self.kR,
                                                                            a=1, p=2, q=2,
                                                                            A0=A0, P0=C0, Q0=D0)
        else:
            eq_dict = ReactionKinetics.compute_equilibrium_conc_first_order(kF=self.kF, kR=self.kR,
                                                                            a=1, p=1, q=1,
                                                                            A0=A0, P0=C0, Q0=D0)

        # eq_dict contains the keys "A", "B", "P", "Q";
        # translate the standard names A, B, P, Q into the actual names, and also drop any missing term
        return  {self.reactant: eq_dict["A"], self.product_1: eq_dict["P"], self.product_2: eq_dict["Q"]}





###################################################################################################################

class ReactionEnzyme(ReactionCommon):
    """
    Data about a SINGLE enzyme-catalyzed reaction that can be modeled kinetically as:

    E + S <-> ES -> E + P

        E : Enzyme
        S : Substrate
        ES: Intermediate Enzyme-Substrate complex
        P : Product
    """

    def __init__(self, enzyme :str, substrate :str, product :str, intermediate=None,
                 k1_F=None, k1_R=None, k2_F=None, k2_R=None,
                 kM=None, kcat=None, **kwargs):
        """
        :param enzyme:      The label for the chemical acting as enzyme
        :param substrate:   The reactant
        :param product:     The final reaction product
        :param intermediate:[OPTIONAL] The label for the reaction intermediate;
                                if not specified, the labels of the enzyme and substrate get concatenated.
                                EXAMPLE: "ES"
        :param k1_F:        [OPTIONAL]The forward reaction rate of the 1st part of the reaction
        :param k1_R:        [OPTIONAL]The reverse reaction rate of the 1st part of the reaction
        :param k2_F:        [OPTIONAL]The forward reaction rate of the 2nd part of the reaction
        :param k2_R:        [NOT IN USE]
        :param kM:          [OPTIONAL]"Michaelis constant"
        :param kcat:        [OPTIONAL]"Catalytic rate constant" aka "Turnover number" aka "Collective rate constant"
                            (equal to k2_F)
        """
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        assert type(enzyme) == str, "ReactionEnzyme instantiation: argument `enzyme` must be a string"
        assert type(substrate) == str, "ReactionEnzyme instantiation: argument `substrate` must be a string"
        assert type(product) == str, "ReactionEnzyme instantiation: argument `product` must be a string"
        # TODO: assert that the strings are valid labels (not just empty space, etc); create function in ChemData

        assert (substrate != product), \
            "ReactionEnzyme instantiation: the `substrate` cannot be the same as the `product`"

        self.enzyme = enzyme
        self.substrate = substrate


        self.intermediate = intermediate
        if intermediate is None:
            # Concatenate the labels of the enzyme and the substrate.  EXAMPLE: "ES"
            self.intermediate = enzyme + substrate

        self.product = product

        self.k1_F = k1_F
        self.k1_R = k1_R
        self.k2_F = k2_F
        self.k2_R = k2_R

        self.kM = kM
        self.kcat = kcat

        self.vmax = None

        if all(v is not None for v in [k1_F, k1_R, k2_F]):
            self.kM = (k2_F + k1_R) / k1_F
            if kM is not None:
                assert np.allclose(self.kM, kM), \
                    f"ReactionEnzyme instantiation: inconsistent arguments.  " \
                    f"The passed `kM` value ({kM}) doesn't match the value ({self.kM}) inferred from the given reaction rate constants"

        if k2_F is not None:
            self.kcat = k2_F
            if kcat is not None:
                assert np.allclose(self.kcat, kcat), \
                    f"ReactionEnzyme instantiation: inconsistent arguments.  " \
                    f"The passed `kcat` value ({kcat}) doesn't match the value ({self.kcat}) of the given `k2_F` reaction rate constants"




    def describe(self, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction

        :param concise:     If True, less detail is shown
        :return:            A string with a description of this reaction
        """
        description = f"{self.enzyme} + {self.substrate} <-> {self.intermediate} -> {self.enzyme} + {self.product}"

        if not concise:
            description += "  (Enzymatic reaction)"
            description += f"  (k1_F = {self.k1_F:,.5g} / k1_R = {self.k1_R:,.5g} / k2_F = {self.k2_F:,.5g} / kM = {self.kM:,.5g} / kcat = {self.kcat:,.5g}"
            if self.temp:
                 description += f" / Temp = {self.temp - 273.15:,.4g} C"

            description += ")"
            # TODO: add more thermodynamic data, if available

        return description



    def extract_reactant_labels(self) -> [str]:
        """
        Return the list of ALL the reactant labels in this reaction,
        (including the enzyme)

        :return:
        """
        return [self.substrate, self.enzyme]


    def extract_reactants(self) -> [(int, str)]:
        """
        Return a list of pairs with details of ALL the reactants of the given reaction
        (including the enzyme)
        Each pair contains stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
        """
        return [(1, self.substrate), (1, self.enzyme)]


    def extract_reactants_formula(self) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        return f"{self.substrate} + {self.enzyme}"



    def extract_intermediate(self) -> str | None:
        """
        Return the name of the reaction intermediate species,
        or None if there's no intermediate
        :return:
        """
        return self.intermediate



    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels of ALL the products of this reaction,
        (including the enzyme)

        :return:
        """
        return [self.product, self.enzyme]


    def extract_products(self) -> [(int, str, int)]:
        """
        Return a list of pairs with details of the products of the given reaction,
        incl. their stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
        """
        return [(1, self.product), (1, self.enzyme)]


    def extract_products_formula(self) -> str:
        """
        Return a string with a user-friendly form of the right (product) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        return f"{self.product} + {self.enzyme}"



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the labels of ALL the chemicals appearing in this reaction,
        including enzyme and intermediate

        :return:    A SET of the labels of the chemicals involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return {self.enzyme, self.substrate, self.intermediate, self.product}



    def extract_rxn_properties(self) -> {}:
        """
        Create a dictionary with the numerical properties of the given reaction
        """
        return {'k1_F': self.k1_F, 'k1_R': self.k1_R, 'k2_F': self.k2_F}    # TODO: possibly add more



    def compute_vmax(self, E_tot :float) -> float:
        """

        :param E_tot:   Total Enzyme concentration (bound and unbound enzyme);
                            at times referred to as E0
        :return:        The maximal reaction rate
                            (the asymptote of the rate, as the Substrate concentration grows large relative to Enzyme concentration)
        """
        assert self.kcat is not None, \
            "compute_vmax(): missing value for kcat"

        self.vmax = self.kcat * E_tot
        return self.vmax



    def compute_rate(self, S_conc :float) -> float:
        """
        Based on the traditional Michaelis-Menten model.
        The argument may also be Numpy array.

        :param S_conc:  The concentration [S] of the (free) Substrate,
                            or a Numpy array of concentrations
        :return:        The corresponding reaction rate, in terms of production of the product P
        """
        # TODO: possibly move to ReactionKinetics
        return self.vmax * S_conc / (self.kM + S_conc)



    def compute_rate_morrison(self, S_tot :float, E_tot :float) -> float:
        """
        Based on the Morrison model.
        Especially useful in scenarios with high concentrations of enzyme.
        The arguments may also be Numpy arrays.

        Reference: eqn 7.32 on page 124 of "Analysis of Enzyme Reaction Kinetics, Vol. 1",
                   by F. Xavier Malcata, Wiley, 2023

        :param S_tot:   The total concentration of free Substrate and Substrate bound to Enzyme
                            (i.e. [S] + [ES])
        :param E_tot:   Total Enzyme concentration (bound and unbound enzyme);
                            at times referred to as E0
        :return:        The corresponding reaction rate, in terms of production of the product P
        """
        # TODO: possibly move to ReactionKinetics
        S_over_E = S_tot / E_tot

        kM_over_E = self.kM / E_tot

        radicand = (1 + S_over_E + kM_over_E)**2 - 4 * S_over_E

        term = 1 + S_over_E + kM_over_E - np.sqrt(radicand)

        return 0.5 * self.vmax * term



    def compute_k1_forward(self, kM, kcat, k1_reverse, verbose=False):
        """
        Compute and return the value for k1_forward, given kM, kcat and k1_reverse.
        Note that this is a linear affine transformation : k1_forward = k1_reverse * (1 / kM) + (kcat / kM)

        :param kM:
        :param kcat:
        :param k1_reverse:
        :param verbose:
        :return:
        """
        k1_forward = (k1_reverse + kcat) / kM
        if verbose:
                K = k1_forward / k1_reverse
                print(f"k1_forward: {k1_forward} , K (k1_f / k1_r) = {K}")

        return k1_forward



    def compute_k1_reverse(self, kM, kcat, k1_forward: Union[float, np.ndarray], verbose=False):
        """
        Compute and return the value for k1_reverse, given kM, kcat and k1_forward
        Note that this is a linear affine transformation : k1_reverse = k1_forward * kM  - kcat

        :param kM:
        :param kcat:
        :param k1_forward:
        :param verbose:
        :return:
        """
        # Verify that the combination of given parameter is physically possible
        min_value_k1_f = self.min_k1_forward(kM, kcat)
        if type(k1_forward) == np.ndarray:
            assert (k1_forward >= min_value_k1_f).all(), \
                f"compute_k1_reverse(): given the specified kM ({kM}) and kcat ({kcat}), some of the k1_forward values " \
                f"are not physically meaningful, as they would lead to a negative value for k1_reverse!  " \
                f"The minimum valid value for k1_forward is {min_value_k1_f}"
        else:
            assert k1_forward >= min_value_k1_f, \
                f"compute_k1_reverse(): the given values for kM ({kM}), kcat ({kcat}) and k1_forward ({k1_forward}) " \
                f"are not physically meaningful, as they would lead to a negative value for k1_reverse!  " \
                f"The minimum valid value for k1_forward is {min_value_k1_f}"

        k1_reverse = k1_forward * kM - kcat
        if verbose:
            if np.allclose(k1_reverse, 0):
                print (f"k1_reverse: {k1_reverse} , K (k1_f / k1_r) = INFINITE")
            else:
                K = k1_forward / k1_reverse
                print(f"k1_reverse: {k1_reverse} , K (k1_f / k1_r) = {K}")

        return k1_reverse



    def min_k1_forward(self, kM :float, kcat :float) -> float:
        """
        Return the minimum physically-possible value for k1_forward,
        for the given kinetic parameters kM and kcat

        :param kM:
        :param kcat:
        :return:
        """
        return kcat / kM



    def set_thermodynamic_data(self, temp :float) -> None:
        """
        Set all the thermodynamic data derivable from the given temperature,
        and all previously passed kinetic and thermodynamic data.
        Raise an Exception if any inconsistency is detected.

        :param temp:    System temperature in Kelvins.
        """
        pass    # TODO



    def step_simulation(self, delta_time, conc_dict :dict, exact=False) -> (dict, float):
        """
        Simulate the enzymatic reaction E + S <-> ES -> E + P, over the specified time interval,
        using the "Forward Euler" method

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentrations won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"E": 1.5, "S": 31.6, "ES": 0.4, "P": 19.9}

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical labels to their concentration CHANGES
                                                            during this step
                                - rxn_rate pair             one value for each elementary reaction
                                EXAMPLE of increment_dict_single_rxn: {"E": 0, "S": -2.9, "ES": 0.1, "P": 2.8}
        """

        increment_dict_single_rxn = {}      # The keys are the chemical labels,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations,
        # for the 2 parts of this reaction:
        # 1) E + S <-> ES   (synthesis reaction of the "intermediate" species ES)
        # 2) ES -> E + P    (irreversible decomposition reaction)

        rxn_rate_1 = ReactionKinetics.compute_rate_elementary(reactants = [self.enzyme, self.substrate],
                                                              products=[self.intermediate],
                                                              kF = self.k1_F, kR=self.k1_R, reversible=True,
                                                              conc_dict=conc_dict)

        rxn_rate_2 = ReactionKinetics.compute_rate_elementary(reactants = [self.intermediate],
                                                              products=[self.enzyme, self.product],
                                                              kF = self.k2_F, kR=0, reversible=False,
                                                              conc_dict=conc_dict)


        # PART 1 - Determine the concentration adjustments as a result of the 1st reaction:  E + S <-> ES

        delta_rxn = rxn_rate_1 * delta_time         # forward reaction - reverse reaction

        # The reactants (enzyme and substrate) DECREASE based on the quantity delta_rxn  (stoichiometry = 1)
        increment_dict_single_rxn[self.enzyme] = - delta_rxn
        increment_dict_single_rxn[self.substrate] = - delta_rxn

        # The reaction product (the intermediate species) INCREASES based on the quantity delta_rxn (stoichiometry = 1)
        increment_dict_single_rxn[self.intermediate] = delta_rxn


        # PART 2 - Determine the concentration adjustments as a result of the 2nd reaction:  ES -> E + P

        delta_rxn = rxn_rate_2 * delta_time         # forward reaction only (since irreversible)

        # The reactant (the intermediate species) DECREASES based on the quantity delta_rxn (stoichiometry = 1)
        increment_dict_single_rxn[self.intermediate] -= delta_rxn  # Notice the "-=", because the intermediate also occurred in previous reaction

        # The reaction products (enzyme and product) INCREASE based on the quantity delta_rxn (stoichiometry = 1)
        increment_dict_single_rxn[self.enzyme] += delta_rxn        # Notice the "+=", because the enzyme also occurred in previous reaction
        increment_dict_single_rxn[self.product] = delta_rxn


        return (increment_dict_single_rxn,
                (rxn_rate_1, rxn_rate_2))





###################################################################################################################

class ReactionGeneric(ReactionCommon):
    """
    Data about a generic SINGLE reaction of the most general type,
    with arbitrary number of reactants and products,
    arbitrary stoichiometry,
    and arbitrary kinetic reaction orders for each participating chemical.

    IMPORTANT:  Use only for reactions whose dynamics can be satisfactorily modeled with kinetics involving
                the chemical concentrations raised to given powers.
                Do *NOT* use ReactionGeneric to model typical biochemistry enzymatic reactions
                of the type E + S -> E + P   (use the "ReactionEnzyme" class instead!)

    NOTE:   Simpler, more specific classes are available for Unimolecular Reactions (A <-> B),
            Synthesis Reactions (A + B <-> C), and Decomposition Reactions (A <-> B + C)

    Included:
        - stoichiometry
        - kinetic data (reaction rates, reaction orders)
        - thermodynamic data (temperature, changes in enthalpy/entropy/Gibbs Free Energy)
        - list of involved catalysts

    Each reaction contains:
            "reactants"
            "products"
            "kF"    (forward reaction rate constant)
            "kR"    (reverse reaction rate constant)
            "enzyme" (the index of a chemical that catalyzes this reaction)
            "macro_enzyme" (The pair (macromolecule name, binding site number))
            plus thermodynamic data: delta_H, delta_S, delta_G, K (equilibrium constant) -
                                     for details see class "ThermoDynamics"

    Internally, each Reactant and each Product is a pair of the form: (stoichiometry coefficient, species label).
    Note that any reactant and products might act as catalyst.
    """

    def __init__(self, reactants :str|list, products :str|list,
                 reversible=True, kF=None, kR=None,
                 delta_H=None, delta_S=None, delta_G=None,
                 kinetic_rate_function=ReactionKinetics.compute_rate_mass_action_kinetics,
                 **kwargs):
        """
        Create the structure for a new SINGLE chemical reaction,
        optionally including its kinetic and/or thermodynamic data.
        All the involved chemicals must be already registered - use add_chemical() if needed.

        NOTE: in the reactants and products, if the stoichiometry coefficients aren't specified,
              they're assumed to be 1.
              The reaction orders, if not specified, are assumed to be equal to their corresponding
              stoichiometry coefficients.

              The full structure of each term in the list of reactants and of products
              is the pair:  (stoichiometry coefficient, chemical label

              EXAMPLES of formats to use for each term in the lists of the reactants and of the products:
                "F"         is taken to mean (1, "F) - default stoichiometry coefficient of 1

              It's equally acceptable to use LISTS in lieu of tuples for the pairs

        :param reactants:   A list of terms that are either chemicals labels (with implied stoichiometry 1),
                                or pairs (stoichiometry coefficient , chemical label).
                                If not a list, it will first get turned into one
        :param products:    A list of terms that are either chemicals labels (with implied stoichiometry 1),
                                or pairs (stoichiometry coefficient , chemical label).
                                If not a list, it will first get turned into one
        :param kF:          [OPTIONAL] Forward reaction rate constant
        :param kR:          [OPTIONAL] Reverse reaction rate constant
        :param delta_H:     [OPTIONAL] Change in Enthalpy (from reactants to products)
        :param delta_S:     [OPTIONAL] Change in Entropy (from reactants to products)
        :param delta_G:     [OPTIONAL] Change in Free Energy (from reactants to products), in Joules
        :param kinetic_rate_function:  [OPTIONAL] Note - the current default will be removed in later versions
                                        EXAMPLES:  ReactionKinetics.compute_rate_mass_action_kinetics  (the generalized "standard rate law")
                                                   ReactionKinetics.compute_rate_first_order (reaction is first order in all reactants and products)
        """
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        self.reactants = None       # A list of pairs (stoichiometry, chemical label)
        self.products = None        # A list of pairs (stoichiometry, chemical label)

        self.macro_enzyme = None    # The pair (macromolecule name, binding site number)
                                    #   EXAMPLE: ("M2", 5)          TODO: maybe turn into a data object

        #self.catalyst = None       # The label of a chemical that catalyzes this reaction, if applicable (at most 1)
                                    # TODO: deprecate

        self.kinetic_rate_function = kinetic_rate_function
                                    # A function used to estimate the reaction rate(aka "velocity"),
                                    # at the start of the time step.
                                    # It takes the following args:
                                    #       reactant_terms :[(int, str)] , product_terms :[(int, str)],
                                    #       kF :float, kR :float,
                                    #       conc_dict :dict
                                    #  and return a float
                                    # EXAMPLES:  ReactionKinetics.compute_rate_mass_action_kinetics  (the generalized "standard rate law")
                                    #            ReactionKinetics.compute_rate_first_order (reaction is first order in all reactants and products)

        self.reversible = reversible
        self.kF = kF                # Forward rate constant
        self.kR = kR                # Reverse rate constant
        self.delta_H = delta_H
        self.delta_S = delta_S
        self.delta_G = delta_G
        self.K = None

        if (self.kinetic_rate_function == ReactionKinetics.compute_rate_mass_action_kinetics) \
                                and (kF is not None) and (kR is not None) and (not np.allclose(self.kR, 0)):
            self.K = kF / kR    # True for reactions that follows mass-action kinetics

        if self.temp is not None:
            thermo_data = ThermoDynamics.extract_thermodynamic_data(K=self.K,
                                                                delta_H=delta_H, delta_S=delta_S, delta_G=delta_G,
                                                                temp=self.temp)
            self.K = thermo_data["K"]
            self.delta_H = thermo_data["delta_H"]
            self.delta_S = thermo_data["delta_S"]
            self.delta_G = thermo_data["delta_G"]

            if (self.kinetic_rate_function == ReactionKinetics.compute_rate_mass_action_kinetics) \
                    and (self.K is not None):
                # True for reactions that follows mass-action kinetics
                # This following is what the function set_rate_constants_from_equilibrium_constant() does for elementary reactions
                if (self.kR is None) and (kF is not None):
                    self.kR = kF / self.K

                elif (self.kF is None) and (kR is not None):
                    self.kF = self.K * kR


        assert reactants is not None, "ReactionGeneric(): the argument `reactants` is a required one; it can't be None"
        if type(reactants) != list:
            reactants = [reactants]

        assert products is not None, "ReactionGeneric(): the argument `products` is a required one; it can't be None"
        if type(products) != list:
            products = [products]

        reactant_list = [(1, r) if type(r) == str else r
                            for r in reactants]   # A list of pairs
        product_list =  [(1, p) if type(p) == str else p
                            for p in products]   # A list of pairs


        # Catch identical reaction sides, even if terms are reshuffled
        assert set(reactant_list) != set(product_list), \
            f"Reaction(): the two sides of the reaction can't be identical! " \
            f"Same reactants and products: {self._standard_form_chem_eqn(reactant_list)}"


        # Locate any catalysts that might be present - though for now a warning is issued if more than 1
        enzyme_list = []
        for reactant in reactant_list:
            if reactant in product_list:
                enzyme_list.append(self.extract_chem_label(reactant))

        number_enzymes = len(enzyme_list)

        if number_enzymes == len(reactant_list) or number_enzymes == len(product_list):
            raise Exception(f"ReactionGeneric(): all the terms in the reaction appear to be enzymes!  "
                            f"Enzymes: {enzyme_list}")


        self.reactants = reactant_list
        self.products = product_list

        #if number_enzymes >= 1:
            #self.catalyst = enzyme_list[0]  # In the irregular scenarios that there appear to be multiple catalyst, only one
                                            #   is considered, and a warning is printed out (the other catalysts
                                            #   will be treated as any other reagent/product)
        if number_enzymes > 1:
            print(f"ReactionGeneric(): WARNING - the reaction appears to have multiple enzymes:"
                  f" {enzyme_list}")




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



    def extract_intermediate(self) -> str | None:
        """
        Return the name of the reaction intermediate species,
        or None if there's no intermediate
        :return:
        """
        return None   # There are no intermediates



    def set_macro_enzyme(self, macromolecule: str, site_number: int) -> None:
        """
        Register that the given macromolecule, at the given site on it,
        catalyzes this reaction

        :param macromolecule:   Name of macromolecule acting as a catalyst
        :param site_number:     Integer to identify a binding site on the above macromolecule
        :return:                None
        """
        self.macro_enzyme = (macromolecule, site_number)


    '''
    def extract_catalyst(self) -> Union[str, None]:
        """

        :return:
        """
        #TODO: deprecate?
        return self.catalyst
    '''





    #####################################################################################################

    '''                       ~   TO READ DATA STRUCTURE of the REACTION  ~                           '''

    def ________READ_RXN_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def extract_reactants(self) -> [(int, str)]:
        """
        Return a list of pairs with details of the reactants of the given reaction,
        incl. their stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
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


    def extract_products(self) -> [(int, str)]:
        """
        Return a list of pairs with details of the products of the given reaction,
        incl. their stoichiometry and chemical label

        :return:    A list of pairs of the form (stoichiometry, chemical label)
        """
        return self.products


    def extract_products_formula(self) -> str:
        """
        Return a string with a user-friendly form of the right (products) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        products = self.extract_products()
        return self._standard_form_chem_eqn(products)



    def unpack_for_dynamics(self) -> tuple:
        """
        A convenient unpacking meant for dynamics simulations
        that need the reactants, the products, and the forward and reverse rate constants

        :return:    A 4-element tuple, containing:
                        (reactants , products , forward rate constant , reverse rate constant)
                        Note: both reactants products are lists of pairs
        """
        return (self.reactants, self.products, self.kF, self.kR)



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the labels of ALL the chemicals appearing in this reaction

        :return:    A SET of the labels of ALL the chemicals involved in this reaction
                    Note: being a set, it's NOT in any particular order
        """
        return set(self.extract_reactant_labels()) | set(self.extract_product_labels())   # Union of sets



    def extract_reactant_labels(self) -> [str]:
        """
        Return the list of the labels of ALL the reactant in this reaction,
        (including any catalysts, if applicable),
        in the order in which they appear when the reaction was first defined

        :return:    List of chemical labels
        """
        reactants = self.extract_reactants()
        reactant_names = [self.extract_chem_label(r) for r in reactants]

        return reactant_names


    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels of ALL the reaction products,
        (including any catalysts, if applicable),
        in the order in which they appear when the reaction was first defined

        :return:    List of chemical labels
        """
        products = self.extract_products()
        product_names = [self.extract_chem_label(r) for r in products]

        return product_names




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
        :return:            A string with a description of this reaction
        """
        reactants = self.extract_reactants()
        products = self.extract_products()

        left = self._standard_form_chem_eqn(reactants)       # Left side of the equation, as a user-friendly string
        right = self._standard_form_chem_eqn(products)       # Right side of the equation

        if self.reversible:
            rxn_description = f"{left} <-> {right}"
        else:
            rxn_description = f"{left} -> {right}"

        if concise:
            return rxn_description      # Minimalist description


        # If we get this far, we're looking for a more detailed description
        rxn_description += self.reaction_details(self.extract_rxn_properties())


        # If a CATALYST is involved, show it
        #if self.catalyst is not None:
            #rxn_description += f" | Enzyme: {self.catalyst}"

        if self.macro_enzyme is not None:
            rxn_description += f" | Macromolecule Enzyme: {self.macro_enzyme[0]}, at site # {self.macro_enzyme[1]}"


        # Show details about the reaction kinetics handler
        if self.kinetic_rate_function:
            rxn_description += f" | Kinetic rate function : `{self.kinetic_rate_function.__name__}()`"
        else:
            rxn_description += " | No kinetic rate function provided"

        return rxn_description




    #####################################################################################################

    '''                                     ~   ANALYSIS  ~                                           '''

    def ________ANALYSIS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def reaction_quotient(self, conc, explain=False) -> np.double | Tuple[np.double, str]:
        """
        Compute the "Reaction Quotient" (aka "Mass–action Ratio"),
        given the concentrations of chemicals involved in this reaction

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
        assert self.kinetic_rate_function == ReactionKinetics.compute_rate_mass_action_kinetics, \
            "reaction_quotient(): the reaction quotient can only be computed when the reaction follows the 'standard rate law'; " \
            "consider a call to set_rate_function(ReactionKinetics.compute_rate_mass_action_kinetics)"

        reactant_data = [(self.extract_stoichiometry(r), self.extract_chem_label(r))
                         for r in self.reactants]

        product_data = [(self.extract_stoichiometry(p), self.extract_chem_label(p) )
                               for p in self.products]

        return ThermoDynamics.compute_reaction_quotient(reactant_data=reactant_data,
                                                        product_data=product_data,
                                                        conc=conc, explain=explain)



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



    def set_thermodynamic_data(self, temp :float) -> None:
        """
        Set all the thermodynamic data derivable from the given temperature,
        and all previously passed kinetic and thermodynamic data.
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

        self.set_rate_constants_from_equilibrium_constant(K=self.K)




    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the generic reaction,
        determine its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step.

        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """
        function_to_call = self.kinetic_rate_function
        assert function_to_call is not None, \
            f"determine_reaction_rate(): no kinetic rate function was provide for the reaction `{self.describe(concise=True)}` isn't set; " \
            f"make sure to first call set_rate_function()"
        #print(f"determine_reaction_rate() - function being invoked to determine the reaction's rate: `{function_to_call.__name__}()`")
        return function_to_call(reactant_terms=self.reactants, product_terms=self.products,
                                                      kF = self.kF, kR=self.kR,
                                                      conc_dict=conc_dict)                        # Carry out the function call



    def set_rate_constants_from_equilibrium_constant(self, K) -> None:
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
        if self.K is None:
            return

        if self.kinetic_rate_function != ReactionKinetics.compute_rate_mass_action_kinetics:
            return

        if (self.kR is None) and (self.kF is not None):
            self.kR = self.kF / K
            return

        if (self.kF is None) and (self.kR is not None):
            self.kF = K * self.kR
            return

        if (self.kF is not None) and (self.kR is not None):
            assert np.allclose(K, self.kF / self.kR), \
                f"set_rate_constants_from_equilibrium_constant(): values for kR and kR already exist, " \
                f"and are inconsistent with the passed value of K ({K})"



    def step_simulation(self, delta_time, conc_dict :dict, exact=False) -> (dict, float):
        """
        Simulate the generic reaction, over the specified time interval.
        The forward Euler method is used

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentrations won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :param exact:       UNUSED - this option is unavailable
        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn     The mapping of chemical labels
                                                                    to their concentration CHANGES
                                                                    during this step
                                - rxn_rate                      The reaction rate ("velocity") for this reaction
                                EXAMPLE of increment_dict_single_rxn: {"B": -1.3, "F": 2.9, "D": -1.6}
        """

        increment_dict_single_rxn = {}      # The keys are the chemical indexes,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        delta_rxn = rxn_rate * delta_time      # forward reaction - reverse reaction


        reactants = self.extract_reactants() # A list of pairs of the form (stoichiometry, chemical label))
        products = self.extract_products()   # A list of pairs of the form (stoichiometry, chemical label))


        """
        Determine the concentration adjustments as a result of this reaction step, 
        for this individual reaction being considered
        """

        # The reactants DECREASE based on the quantity delta_rxn
        for r in reactants:
            # Unpack data from the reactant r
            species_name = self.extract_chem_label(r)
            #if species_name == self.catalyst:
                #print(f"*** SKIPPING reactant CATALYST {species_name} in reaction")
                #continue    # Skip if r is a catalyst for this reaction

            stoichiometry = self.extract_stoichiometry(r)

            delta_conc = stoichiometry * (- delta_rxn)  # Increment to this reactant from the reaction being considered

            increment_dict_single_rxn[species_name] = increment_dict_single_rxn.get(species_name,0) + delta_conc


        # The reaction products INCREASE based on the quantity delta_rxn
        for p in products:
            # Unpack data from the reactant r
            species_name = self.extract_chem_label(p)
            #if species_name == self.catalyst:
                #print(f"*** SKIPPING product CATALYST {species_name} in reaction")
                #continue    # Skip if p is a catalyst for this reaction

            stoichiometry = self.extract_stoichiometry(p)

            delta_conc = stoichiometry * delta_rxn  # Increment to this reaction product from the reaction being considered

            increment_dict_single_rxn[species_name] = increment_dict_single_rxn.get(species_name,0) + delta_conc

        '''
        # Macro-molecule related part, if applicable    TODO: implement
        if (self.macro_system_state != {}) and (rxn.macro_enzyme is not None):
            print(f"[NOT YET IMPLEMENTED] Making adjustments for macro-molecule catalysis for reaction")    #  # {rxn_index}
            print(f"    Macromolecule: {rxn.macro_enzyme[0]}, at site # {rxn.macro_enzyme[1]}")
            #print(f"    Site occupancy at the beginning of the time step:")
            #print(f"    Macromolecule count:")
        '''

        assert len(increment_dict_single_rxn) == len(self.extract_chemicals_in_reaction())  # TODO: temp test

        return (increment_dict_single_rxn, rxn_rate)



    def find_equilibrium_conc(self, conc_dict :dict) -> dict:
        """
        Determine the equilibrium concentrations that would be eventually reached
        by the chemicals of this reaction,
        given the specified concentrations,
        IN THE ABSENCE of any other reaction.

        :param conc_dict:   A dict mapping chemical labels to their initial concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"X": 4.3, "Y": 1.5, "F": 31.6, "G": 3.6}

        :return:            A dict mapping the above chemical labels to their equilibrium concentrations
        """
        assert self.kinetic_rate_function == ReactionKinetics.compute_rate_first_order, \
            "find_equilibrium_conc(): only implemented for reactions whose kinetic rate function is `ReactionKinetics.compute_rate_first_order` \n" \
            "if that's the case, make sure to first invoke:   set_rate_function(ReactionKinetics.compute_rate_first_order)"

        # To conform with ReactionKinetics.compute_equilibrium_conc,
        # we'll express the reaction in the form aA + bB <-> bC + dD

        reactants = self.extract_reactants()
        products = self.extract_products()

        assert len(reactants) == 2, \
            f"find_equilibrium_conc(): only implemented for reactions with 2 reactants (we have {len(reactants)})"

        assert len(products) == 2, \
            f"find_equilibrium_conc(): only implemented for reactions with 2 products (we have {len(products)})"

        r1, r2 = reactants
        p1, p2 = products

        A0 = conc_dict.get(self.extract_chem_label(r1))
        assert A0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the reactant `{self.extract_chem_label(r1)}` was not provided"

        B0 = conc_dict.get(self.extract_chem_label(r2))
        assert B0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the reactant `{self.extract_chem_label(r2)}` was not provided"


        C0 = conc_dict.get(self.extract_chem_label(p1))
        assert C0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the produce `{self.extract_chem_label(p1)}` was not provided"

        D0 = conc_dict.get(self.extract_chem_label(p2))
        assert D0 is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                               f"concentration of the produce `{self.extract_chem_label(p2)}` was not provided"



        eq_dict = ReactionKinetics.compute_equilibrium_conc_first_order(kF=self.kF, kR=self.kR,
                                                                        a=self.extract_stoichiometry(r1), b=self.extract_stoichiometry(r2),
                                                                        p=self.extract_stoichiometry(p1), q=self.extract_stoichiometry(p1),
                                                                        A0=A0, B0=B0, P0=C0, Q0=D0)

        # eq_dict contains the keys "A", "B", "P", "Q";
        # translate the standard names A, B, P, Q into the actual names, and also drop any missing term
        return  {self.extract_chem_label(r1): eq_dict["A"], self.extract_chem_label(r2): eq_dict["B"],
                 self.extract_chem_label(p1): eq_dict["P"], self.extract_chem_label(p2): eq_dict["Q"]}






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
            species_name = self.extract_chem_label(t)

            if stoichiometry == 1:
                term = species_name
            else:
                term = f"{stoichiometry} {species_name}"

            formula_list.append(term)

        return " + ".join(formula_list)
