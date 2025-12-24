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

    Individual reactions classes need to provide the following methods [TODO: partial list]:
        extract_reactant_labels()   Incl. enzyme, if applicable
        extract_reactants()

        extract_product_labels()    Incl. enzyme, if applicable
        extract_products()
        extract_catalyst()
        extract_chemicals_in_reaction()

        extract_rxn_properties()     A dict

        set_thermodynamic data()
        step_simulation()
        etc.
    """

    def __init__(self, active=True, temp=None):
        """
        :param active:  [NOT YET IMPLEMENTED] If False, the reaction is regarded as inactive
        :param temp:    In Kelvins
        """
        self.active = active        # TODO: not yet in use
        self.temp = temp            # In Kelvins
        self.catalyst = None        # The label of a chemical that catalyzes this reaction, if applicable (at most 1)




    def extract_stoichiometry(self, term :(int, str, int)) -> int:
        """
        Return the stoichiometry coefficient, from a reaction "TERM"

        :param term:    A triplet (int, str, int) representing a reaction term
        :return:        An integer with the stoichiometry coefficient
        """
        return term[0]


    def extract_species_name(self, term :(int, str, int)) -> str:
        """
        Return the name of the chemical species, from a reaction "TERM"

        :param term:    A triplet (int, str, int) representing a reaction term
        :return:        The name of the chemical species in the term
        """
        # TODO: rename to extract_chem_label()
        return term[1]


    def extract_rxn_order(self, term :(int, str, int)) -> int:
        """
        Return the reaction order, from a reaction "TERM"

        :param term:    A triplet (int, str, int) representing a reaction term
        :return:        An integer with the reaction order for this term
        """
        return term[2]



    def extract_catalyst(self) -> Union[str, None]:
        """

        :return:
        """
        return self.catalyst





###################################################################################################################


class ReactionOneStep(ReactionCommon):
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
                f"ReactionOneStep instantiation: irreversible reactions cannot have a value for the reverse rate constant (kR = {kR})"

        self.kF = kF                # Forward rate constant
        self.kR = kR                # Reverse rate constant
        self.delta_H = delta_H
        self.delta_S = delta_S
        self.delta_G = delta_G
        self.K = None               # Equilibrium constant

        if (kR is not None) and not np.allclose(self.kR, 0):
            self.K = kF / kR

        # Process the kinetic and thermodynamic data, and update various object attributes accordingly
        self._set_kinetic_and_thermodynamic(forward_rate=kF, reverse_rate=kR,
                                            delta_H=delta_H, delta_S=delta_S, delta_G=delta_G, temp=self.temp)



    def reaction_details(self) -> str:
        """
        Return a string with some details about the parameters of this reaction

        :return:    EXAMPLE: "  (kF = 3 / kR = 2 / Delta_G = -1,005.13 / Temp = 25 C)"
        """
        details = []
        rxn_properties = self.extract_rxn_properties()
        for k,v in rxn_properties.items():
            details.append(f"{k} = {v:,.5g}")          # EXAMPLE: "kF = 3"

        description = ""

        if self.temp:
            details.append(f"Temp = {self.temp - 273.15:,.4g} C")          # EXAMPLE: "Temp = 25 C"

        if details:
            description = "  (" + ' / '.join(details) + ")"   # EXAMPLE: "  (kF = 3 / kR = 2 / Delta_G = -1,005.13)"


        return description



    def extract_intermediates(self) -> [str]:
        """

        :return:
        """
        return []   # There are no intermediates



    def _set_kinetic_and_thermodynamic(self, forward_rate, reverse_rate,
                                       delta_H, delta_S, delta_G, temp) -> None:
        """
        Set all the kinetic and thermodynamic data derivable - directly or indirectly - from the passed arguments,
        storing it in object attributes.
        Raise an Exception if any inconsistency is detected.

        :param forward_rate:
        :param reverse_rate:
        :param delta_H:
        :param delta_S:
        :param delta_G:
        :param temp:
        :return:                None
        """
        self.kF = forward_rate
        self.kR = reverse_rate
        self.delta_H = delta_H
        self.delta_S = delta_S
        self.delta_G = delta_G


        # Process kinetic data, if available
        #       (extracting thermodynamic data when feasible)
        if (self.kF is not None) and (self.kR is not None) and not np.allclose(self.kR, 0):
            # If all the kinetic data is available...
            self.K = self.kF / self.kR    # ...compute the equilibrium constant (from kinetic data)

            if temp:
                # If the temperature is set, compute the change in Gibbs Free Energy
                delta_G_kinetic = ThermoDynamics.delta_G_from_K(K = self.K, temp = temp)
                if self.delta_G is None:
                    self.delta_G = delta_G_kinetic
                else:   # If already present (passed as argument), make sure that the two match!
                    assert np.allclose(delta_G_kinetic, self.delta_G), \
                        f"_set_kinetic_and_thermodynamic(): Kinetic data (leading to Delta_G={delta_G_kinetic}) " \
                        f"is inconsistent with the passed value of Delta_G={self.delta_G})"


        if (self.delta_H is not None) and (self.delta_S is not None) and (temp is not None):
            # If all the thermodynamic data (possibly except delta_G) is available...

            # Compute the change in Gibbs Free Energy from delta_H and delta_S, at the current temperature
            delta_G_thermo = ThermoDynamics.delta_G_from_enthalpy(delta_H = self.delta_H, delta_S = self.delta_S, temp = temp)

            if self.delta_G is None:
                self.delta_G = delta_G_thermo
            else:  # If already present (passed as argument or was set from kinetic data), make sure that the two match!
                if not np.allclose(delta_G_thermo, self.delta_G):
                    if delta_G is not None:
                        raise Exception(f"_set_kinetic_and_thermodynamic(): thermodynamic data (leading to Delta_G={delta_G_thermo}) "
                                        f"is inconsistent with the passed value of delta_G={delta_G})")
                    else:
                        raise Exception(f"_set_kinetic_and_thermodynamic(): thermodynamic data (leading to Delta_G={delta_G_thermo}) "
                                        f"is inconsistent with kinetic data (leading to Delta_G={self.delta_G})")


        if self.delta_G is not None:
            if (self.K is None) and (temp is not None):
                # If the temperature is known, compute the equilibrium constant (from the thermodynamic data)
                # Note: no need to do it if self.K is present, because we ALREADY handled that case
                self.K = ThermoDynamics.K_from_delta_G(delta_G = self.delta_G, temp = temp)

                # If only one of the Forward or Reverse rates was provided, compute the other one
                if (self.kF is None) and (self.kR is not None):
                    self.kF = self.K * self.kR
                if (self.kR is None) and (self.kF is not None):
                    self.kR = self.kF / self.K

            if temp is not None:
                # If either Enthalpy or Entropy is missing, but the other one is known, compute the missing one
                if (self.delta_H is None) and (self.delta_S is not None):
                    self.delta_H = ThermoDynamics.delta_H_from_gibbs(delta_G=self.delta_G, delta_S=self.delta_S, temp=temp)
                elif (self.delta_H is not None) and (self.delta_S is None):
                    self.delta_S = ThermoDynamics.delta_S_from_gibbs(delta_G=self.delta_G, delta_H=self.delta_H, temp=temp)



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

        if not temp:
            return      # There's nothing to do

        # Process kinetic data, if available,
        #       extracting thermodynamic data when feasible
        if self.K:
            # If the temperature is set, compute the change in Gibbs Free Energy
            delta_G_kinetic = ThermoDynamics.delta_G_from_K(K = self.K, temp = temp)
            if self.delta_G is None:
                self.delta_G = delta_G_kinetic
            else:   # If already present (passed as argument), make sure that the two match!
                assert np.allclose(delta_G_kinetic, self.delta_G), \
                    f"set_thermodynamic_data(): Kinetic data (leading to Delta_G={delta_G_kinetic}) " \
                    f"is inconsistent with the passed value of Delta_G={self.delta_G})"


        if (self.delta_H is not None) and (self.delta_S is not None):
            # If all the thermodynamic data (possibly except delta_G) is available...

            # Compute the change in Gibbs Free Energy from delta_H and delta_S, at the current temperature
            delta_G_thermo = ThermoDynamics.delta_G_from_enthalpy(delta_H = self.delta_H, delta_S = self.delta_S, temp = temp)

            if self.delta_G is None:
                self.delta_G = delta_G_thermo
            else:  # If already present (passed as argument or was set from kinetic data), make sure that the two match!
                if not np.allclose(delta_G_thermo, self.delta_G):
                    if self.delta_G is not None:
                        raise Exception(f"set_thermodynamic_data(): thermodynamic data (leading to Delta_G={delta_G_thermo}) "
                                        f"is inconsistent with the passed value of delta_G={self.delta_G})")
                    else:
                        raise Exception(f"set_thermodynamic_data(): thermodynamic data (leading to Delta_G={delta_G_thermo}) "
                                        f"is inconsistent with kinetic data (leading to Delta_G={self.delta_G})")


        if self.delta_G is not None:
            if self.K is None:
                # If the temperature is known, compute the equilibrium constant (from the thermodynamic data)
                # Note: no need to do it if self.K is present, because we ALREADY handled that case
                self.K = ThermoDynamics.K_from_delta_G(delta_G = self.delta_G, temp = temp)

                # If only one of the Forward or Reverse rates was provided, compute the other one
                if (self.kF is None) and (self.kR is not None):
                    self.kF = self.K * self.kR
                if (self.kR is None) and (self.kF is not None):
                    self.kR = self.kF / self.K


            # If either Enthalpy or Entropy is missing, but the other one is known, compute the missing one
            if (self.delta_H is None) and (self.delta_S is not None):
                self.delta_H = ThermoDynamics.delta_H_from_gibbs(delta_G=self.delta_G, delta_S=self.delta_S, temp=temp)
            elif (self.delta_H is not None) and (self.delta_S is None):
                self.delta_S = ThermoDynamics.delta_S_from_gibbs(delta_G=self.delta_G, delta_H=self.delta_H, temp=temp)



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





###################################################################################################################

class ReactionUnimolecular(ReactionOneStep):
    """
    Reactions of type A <-> P, of first order in A and P
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
        :return:            A string with a description of the specified reaction
        """
        description = f"{self.reactant} <-> {self.product}"

        if not concise:
            if self.reversible:
                description += "  (Elementary Unimolecular reaction)"
            else:
                description += "  (Elementary Unimolecular Irreversible reaction)"

            description += self.reaction_details()

        return description



    def extract_reactant_labels(self) -> [str]:
        """
        Return the list of ALL the reactant labels in this reaction

        :return:
        """
        return [self.reactant]


    def extract_reactants(self) -> [(int, str, int)]:
        """
        Return a list of triplets with details of the reactants of the given reaction,
        incl. their stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
        """
        return [(1, self.reactant, 1)]



    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels ALL the products of this reaction

        :return:
        """
        return [self.product]


    def extract_products(self) -> [(int, str, int)]:
        """
        Return a list of triplet with details of the products of the given reaction,
        incl. their stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
        """
        return [(1, self.product, 1)]



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the labels of ALL the chemicals appearing in this reaction

        :return:    A SET of the labels of the chemicals involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return {self.reactant, self.product}



    def reaction_quotient(self, conc, explain=False) -> Union[np.double, Tuple[np.double, str]]:
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
        return ReactionKinetics.compute_reaction_quotient(reactant_data=[(self.reactant, 1)], product_data=[(self.product, 1)],
                                                          conc=conc, explain=explain)



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the unimolecular reaction,
        determine its current reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate".

        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """
        return ReactionKinetics.compute_reaction_rate_first_order(reactants = [self.reactant], products=[self.product],
                                                                  kF = self.kF, kR=self.kR, reversible=self.reversible,
                                                                  conc_dict=conc_dict)



    def step_simulation(self, delta_time, conc_dict :dict) -> (dict, float):
        """
        Simulate the unimolecular reaction A <-> P, over the specified time interval,
        using the "Forward Euler" method

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6}

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical label to their concentration changes
                                                            during this step
                                                            EXAMPLE: {"B": -0.2, "F": 0.2}

                                - rxn_rate                  is the reaction rate ("velocity") for this reaction
                                                            (rate of change of the product), at the START of the simulation step
                                                            EXAMPLE: 3.5

        """
        increment_dict_single_rxn = {}      # The keys are the chemical labels,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        # In the "forward Euler" approximation, this rate is taken to remain unvaried during the entire (small) time step
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        delta_rxn = rxn_rate * delta_time   # forward reaction - reverse reaction


        # Determine the concentration adjustments as a result of this reaction step:

        # The reactant DECREASES based on the quantity delta_rxn
        r = self.reactant           # EXAMPLE: "B"
        delta_conc = - delta_rxn    # Increment to this reactant from the reaction step
        increment_dict_single_rxn[r] = delta_conc


        # The reaction product INCREASES based on the quantity delta_rxn
        p = self.product            # EXAMPLE: "F"
        delta_conc = delta_rxn      # Increment to this reaction product from the reaction step
        increment_dict_single_rxn[p] = delta_conc

        return (increment_dict_single_rxn, rxn_rate)





#######################################################################################################################

class ReactionSynthesis(ReactionOneStep):
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

        '''
        assert (r1 != r2), \
            "ReactionSynthesis instantiation: the 2 reactants cannot be the same. Use ReactionGeneric instead"
            #TODO: maybe overcome this restriction
        '''

        self.reactant_1 = r1
        self.reactant_2 = r2
        self.product = product




    def describe(self, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction

        :param concise:     If True, less detail is shown
        :return:            A string with a description of the specified reaction
        """
        description = f"{self.reactant_1} + {self.reactant_2}  <-> {self.product}"

        if not concise:
            if self.reversible:
                description += "  (Elementary Synthesis reaction)"
            else:
                description += "  (Elementary Synthesis Irreversible reaction)"

            description += self.reaction_details()

        return description



    def extract_reactant_labels(self) -> [str]:
        """
        Return the list of ALL the reactant labels in this reaction

        :return:
        """
        return [self.reactant_1, self.reactant_2]


    def extract_reactants(self) -> [(int, str, int)]:
        """
        Return a list of triplets with details of the reactants of the given reaction,
        incl. their stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
        """
        return [(1, self.reactant_1, 1) , (1, self.reactant_2, 1)]



    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels ALL the products of this reaction

        :return:
        """
        return [self.product]


    def extract_products(self) -> [(int, str, int)]:
        """
        Return a list of triplet with details of the products of the given reaction,
        incl. their stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
        """
        return [(1, self.product, 1)]



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the labels of ALL the chemicals appearing in this reaction

        :return:    A SET of the labels of the chemicals involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return {self.reactant_1, self.reactant_2, self.product}



    def reaction_quotient(self, conc, explain=False) -> Union[np.double, Tuple[np.double, str]]:
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
        return ReactionKinetics.compute_reaction_quotient(reactant_data=[(self.reactant_1, 1) , (self.reactant_2, 1)],
                                                          product_data= [(self.product, 1)],
                                                          conc=conc, explain=explain)



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the synthesis reaction,
        determine its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step.

        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """
        return ReactionKinetics.compute_reaction_rate_first_order(reactants = [self.reactant_1, self.reactant_2],
                                                                  products=[self.product],
                                                                  kF = self.kF, kR=self.kR, reversible=self.reversible,
                                                                  conc_dict=conc_dict)



    def step_simulation(self, delta_time, conc_dict :dict) -> (dict, float):
        """
        Simulate the synthesis reaction A + B <-> C, over the specified time interval,
        using the "Forward Euler" method

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"A": 1.5, "B": 31.6, "C": 19.9}

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical label to their concentration changes
                                                            during this step
                                - rxn_rate                  is the reaction rate ("velocity") for this reaction
                                                            (rate of change of the product)
                                EXAMPLE of increment_dict_single_rxn: "A": -1.3, "B": 2.9, "C": -1.6
        """

        increment_dict_single_rxn = {}      # The keys are the chemical labels,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        delta_rxn = rxn_rate * delta_time   # forward reaction - reverse reaction


        # Determine the concentration adjustments as a result of this reaction step:

        # The reactants DECREASE based on the quantity delta_rxn
        for r in [self.reactant_1, self.reactant_2]:
            # EXAMPLE of r: "A"
            # stoichiometry = 1
            delta_conc = - delta_rxn    # Increment to this reactant from the reaction step
            increment_dict_single_rxn[r] = delta_conc


        # The reaction product INCREASES based on the quantity delta_rxn
        p = self.product            # EXAMPLE: "F"
        # stoichiometry = 1
        delta_conc = delta_rxn      # Increment to this reaction product from the reaction step
        increment_dict_single_rxn[p] = delta_conc

        return (increment_dict_single_rxn, rxn_rate)





#######################################################################################################################

class ReactionDecomposition(ReactionOneStep):
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

        '''
        assert (p1 != p2), \
            "ReactionDecomposition instantiation: the 2 reaction products cannot be the same. Use ReactionGeneric instead"
            #TODO: maybe overcome this restriction
        '''

        self.reactant = reactant
        self.product_1 = p1
        self.product_2 = p2




    def describe(self, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction

        :param concise:     If True, less detail is shown
        :return:            A string with a description of the specified reaction
        """
        description = f"{self.reactant} <-> {self.product_1} + {self.product_2}"

        if not concise:
            if self.reversible:
                description += "  (Elementary Decomposition reaction)"
            else:
                description += "  (Elementary Decomposition Irreversible reaction)"

            description += self.reaction_details()

        return description



    def extract_reactant_labels(self) -> [str]:
        """
        Return the list of ALL the reactant labels in this reaction

        :return:
        """
        return [self.reactant]


    def extract_reactants(self) -> [(int, str, int)]:
        """
        Return a list of triplets with details of the reactants of the given reaction,
        incl. their stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
        """
        return [(1, self.reactant, 1)]



    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels ALL the products of this reaction

        :return:
        """
        return [self.product_1, self.product_2]


    def extract_products(self) -> [(int, str, int)]:
        """
        Return a list of triplet with details of the products of the given reaction,
        incl. their stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
        """
        return [(1, self.product_1, 1), (1, self.product_2, 1)]



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the labels of ALL the chemicals appearing in this reaction

        :return:    A SET of the labels of the chemicals involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return {self.reactant, self.product_1, self.product_2}



    def reaction_quotient(self, conc, explain=False) -> Union[np.double, Tuple[np.double, str]]:
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
        return ReactionKinetics.compute_reaction_quotient(reactant_data=[(self.reactant, 1)],
                                                          product_data=[(self.product_1, 1) , (self.product_2, 1)],
                                                          conc=conc, explain=explain)



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the decomposition reaction,
        determine its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step.

        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """
        return ReactionKinetics.compute_reaction_rate_first_order(reactants = [self.reactant],
                                                                  products=[self.product_1, self.product_2],
                                                                  kF = self.kF, kR=self.kR, reversible=self.reversible,
                                                                  conc_dict=conc_dict)



    def step_simulation(self, delta_time, conc_dict :dict) -> (dict, float):
        """
        Simulate the decomposition reaction A <-> B + C, over the specified time interval,
        using the "Forward Euler" method

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"A": 1.5, "B": 31.6, "C": 19.9}

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical label to their concentration changes
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
        for p in [self.product_1, self.product_2]:
            # EXAMPLE of p: "B"
            # stoichiometry = 1
            delta_conc = delta_rxn      # Increment to this reaction product from the reaction step
            increment_dict_single_rxn[p] = delta_conc

        return (increment_dict_single_rxn, rxn_rate)




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
        self.catalyst = enzyme      # "catalyst" is a more generic concept, shared with other types of reactions
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
        :return:            A string with a description of the specified reaction
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


    def extract_reactants(self) -> [(int, str, int)]:
        """
        Return a list of triplets with details of ALL the reactants of the given reaction
        (including the enzyme)
        Each triplet contains stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
        """
        return [(1, self.substrate, 1), (1, self.enzyme, 1)]



    def extract_intermediates(self) -> [str]:
        """

        :return:
        """
        return [self.intermediate]



    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels of ALL the products of this reaction,
        (including the enzyme)

        :return:
        """
        return [self.product, self.enzyme]


    def extract_products(self) -> [(int, str, int)]:
        """
        Return a list of triplet with details of the products of the given reaction,
        incl. their stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
        """
        return [(1, self.product, 1), (1, self.enzyme, 1)]



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
        assert self.kcat is not None, "compute_vmax(): missing value for kcat"
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



    def step_simulation(self, delta_time, conc_dict :dict) -> (dict, float):
        """
        Simulate the enzymatic reaction E + S <-> ES -> E + P, over the specified time interval,
        using the "Forward Euler" method

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"E": 1.5, "S": 31.6, "ES": 0.4, "P": 19.9}

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical label to their concentration changes
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

        rxn_rate_1 = ReactionKinetics.compute_reaction_rate_first_order(reactants = [self.enzyme, self.substrate],
                                                                  products=[self.intermediate],
                                                                  kF = self.k1_F, kR=self.k1_R, reversible=True,
                                                                  conc_dict=conc_dict)

        rxn_rate_2 = ReactionKinetics.compute_reaction_rate_first_order(reactants = [self.intermediate],
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

class ReactionGeneric(ReactionOneStep):
    """
    DEPRECATED!  May be eliminated...

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

    Internally, each Reactant and each Product is a triplet of the form: (stoichiometry, species index, reaction order).
    The "reaction order" in that triplet refers to the forward reaction for reactants, and to the reverse reaction for products.
    Note that any reactant and products might be catalysts
    """

    def __init__(self, reactants: Union[int, str, list], products: Union[int, str, list], **kwargs):
        """
        Create the structure for a new SINGLE chemical reaction,
        optionally including its kinetic and/or thermodynamic data.
        All the involved chemicals must be already registered - use add_chemical() if needed.

        NOTE: in the reactants and products, if the stoichiometry coefficients aren't specified,
              they're assumed to be 1.
              The reaction orders, if not specified, are assumed to be equal to their corresponding
              stoichiometry coefficients.

              The full structure of each term in the list of reactants and of products
              is the triplet:  (stoichiometry coefficient, name, reaction order)

              EXAMPLES of formats to use for each term in the lists of the reactants and of the products:
                "F"         is taken to mean (1, "F", 1) - default stoichiometry and reaction order
                (2, "F")    is taken to mean (2, "F", 2) - stoichiometry coefficient used as default for reaction order
                (2, "F", 1) means stoichiometry coefficient 2 and reaction order 1 - no defaults invoked
              It's equally acceptable to use LISTS in lieu of tuples for the pair or triplets

        :param reactants:   A list of triplets (stoichiometry, species name, reaction order),
                                or simplified terms in various formats; for details, see above.
                                If not a list, it will get turned into one
        :param products:    A list of triplets (stoichiometry, species name, reaction order of REVERSE reaction),
                                or simplified terms in various formats; for details, see above.
                                If not a list, it will get turned into one
        :param kF:          [OPTIONAL] Forward reaction rate constant
        :param kR:          [OPTIONAL] Reverse reaction rate constant
        :param delta_H:     [OPTIONAL] Change in Enthalpy (from reactants to products)
        :param delta_S:     [OPTIONAL] Change in Entropy (from reactants to products)
        :param delta_G:     [OPTIONAL] Change in Free Energy (from reactants to products), in Joules
        """
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        self.reactants = None       # A list of triplets (stoichiometry, species name, reaction order)
        self.products = None        # A list of triplets (stoichiometry, species name, reaction order)

        self.macro_enzyme = None    # The pair (macromolecule name, binding site number)
                                    #   EXAMPLE: ("M2", 5)          TODO: maybe turn into a data object

        assert reactants is not None, "ReactionGeneric(): the argument `reactants` is a required one; it can't be None"
        if type(reactants) != list:
            reactants = [reactants]


        assert products is not None, "ReactionGeneric(): the argument `products` is a required one; it can't be None"
        if type(products) != list:
            products = [products]


        reactant_list = [self._parse_reaction_term(r, "reactant") for r in reactants]   # A list of triplets of integers
        product_list = [self._parse_reaction_term(r, "product") for r in products]      # A list of triplets of integers

        # Catch identical reaction sides, even if terms are reshuffled
        assert set(reactant_list) != set(product_list), \
            f"Reaction(): the two sides of the reaction can't be identical! " \
            f"Same reactants and products: {self._standard_form_chem_eqn(reactant_list)}"


        # Locate any catalysts that might be present - though for now a warning is issued if more than 1
        enzyme_list = []
        for reactant in reactant_list:
            if reactant in product_list:
                enzyme_list.append(self.extract_species_name(reactant))

        number_enzymes = len(enzyme_list)

        if number_enzymes == len(reactant_list) or number_enzymes == len(product_list):
            raise Exception(f"ReactionGeneric(): all the terms in the reaction appear to be enzymes!  "
                            f"Enzymes: {enzyme_list}")


        self.reactants = reactant_list
        self.products = product_list

        if number_enzymes >= 1:
            self.catalyst = enzyme_list[0]  # In the irregular scenarios that there appear to be multiple catalyst, only one
                                            #   is considered, and a warning is printed out (the other catalysts
                                            #   will be treated as any other reagent/product)
        if number_enzymes > 1:
            print(f"ReactionGeneric(): WARNING - the reaction appears to have multiple enzymes:"
                  f" {enzyme_list}")




    def set_macro_enzyme(self, macromolecule: str, site_number: int) -> None:
        """
        Register that the given macromolecule, at the given site on it,
        catalyzes this reaction

        :param macromolecule:   Name of macromolecule acting as a catalyst
        :param site_number:     Integer to identify a binding site on the above macromolecule
        :return:                None
        """
        self.macro_enzyme = (macromolecule, site_number)




    #####################################################################################################

    '''                       ~   TO READ DATA STRUCTURE of the REACTION  ~                           '''

    def ________READ_RXN_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def extract_reactants(self) -> [(int, str, int)]:
        """
        Return a list of triplets with details of the reactants of the given reaction,
        incl. their stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
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



    def extract_products(self) -> [(int, str, int)]:
        """
        Return a list of triplet with details of the products of the given reaction,
        incl. their stoichiometry, chemical label, and reaction order

        :return:    A list of triplets of the form (stoichiometry, chemical label, reaction order)
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
                        Note: both reactants products are lists of triplets
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
        reactant_names = [self.extract_species_name(r) for r in reactants]

        return reactant_names


    def extract_product_labels(self) -> [str]:
        """
        Return the list of the labels of ALL the reaction products,
        (including any catalysts, if applicable),
        in the order in which they appear when the reaction was first defined

        :return:    List of chemical labels
        """
        products = self.extract_products()
        product_names = [self.extract_species_name(r) for r in products]

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
        :return:            A string with a description of the specified reaction
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
        rxn_description += self.reaction_details()


        # If a CATALYST is involved, show it
        if self.catalyst is not None:
            rxn_description += f" | Enzyme: {self.catalyst}"

        if self.macro_enzyme is not None:
            rxn_description += f" | Macromolecule Enzyme: {self.macro_enzyme[0]}, at site # {self.macro_enzyme[1]}"


        # Show details about the ORDER of the reaction
        high_order = False      # Is there any term whose order is greater than 1 ?
        for r in reactants:
            if r[2] > 1:
                rxn_description += f" | {self.extract_rxn_order(r)}-th order in reactant {self.extract_species_name(r)}"
                high_order = True
        for p in products:
            if p[2] > 1:
                rxn_description += f" | {self.extract_rxn_order(p)}-th order in product {self.extract_species_name(p)}"
                high_order = True

        if not high_order:
            rxn_description += " | 1st order in all reactants & products"

        return rxn_description




    #####################################################################################################

    '''                                     ~   ANALYSIS  ~                                           '''

    def ________ANALYSIS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def reaction_quotient(self, conc, explain=False) -> Union[np.double, Tuple[np.double, str]]:
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
        reactants_and_order = [(self.extract_species_name(r) , self.extract_rxn_order(r))
                                for r in self.reactants]

        products_and_order = [(self.extract_species_name(p) , self.extract_rxn_order(p))
                               for p in self.products]

        return ReactionKinetics.compute_reaction_quotient(reactant_data=reactants_and_order,
                                                          product_data=products_and_order,
                                                          conc=conc, explain=explain)



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the generic reaction,
        determine its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step.

        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """

        reactants_and_order = [(self.extract_species_name(r) , self.extract_rxn_order(r))
                                for r in self.reactants]

        products_and_order = [(self.extract_species_name(p) , self.extract_rxn_order(p))
                               for p in self.products]

        return ReactionKinetics.compute_reaction_rate(reactant_data= reactants_and_order, product_data=products_and_order,
                                                      kF = self.kF, kR=self.kR, reversible=self.reversible,
                                                      conc_dict=conc_dict)



    def step_simulation(self, delta_time, conc_dict :dict) -> (dict, float):
        """
        Simulate the generic reaction, over the specified time interval

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical label to their concentration changes
                                                            during this step
                                - rxn_rate                  is the reaction rate ("velocity") for this reaction
                                EXAMPLE of increment_dict_single_rxn: {"B": -1.3, "F": 2.9, "D": -1.6}
        """

        increment_dict_single_rxn = {}      # The keys are the chemical indexes,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        delta_rxn = rxn_rate * delta_time      # forward reaction - reverse reaction


        reactants = self.extract_reactants() # A list of triplets of the form (stoichiometry, species name, reaction order)
        products = self.extract_products()   # A list of triplets of the form (stoichiometry, species name, reaction order)


        """
        Determine the concentration adjustments as a result of this reaction step, 
        for this individual reaction being considered
        """

        # The reactants DECREASE based on the quantity delta_rxn
        for r in reactants:
            # Unpack data from the reactant r
            species_name = self.extract_species_name(r)
            if species_name == self.catalyst:
                #print(f"*** SKIPPING reactant CATALYST {species_name} in reaction")
                continue    # Skip if r is a catalyst for this reaction

            stoichiometry = self.extract_stoichiometry(r)

            delta_conc = stoichiometry * (- delta_rxn)  # Increment to this reactant from the reaction being considered

            increment_dict_single_rxn[species_name] = increment_dict_single_rxn.get(species_name,0) + delta_conc


        # The reaction products INCREASE based on the quantity delta_rxn
        for p in products:
            # Unpack data from the reactant r
            species_name = self.extract_species_name(p)
            if species_name == self.catalyst:
                #print(f"*** SKIPPING product CATALYST {species_name} in reaction")
                continue    # Skip if p is a catalyst for this reaction

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
            species_name = self.extract_species_name(t)

            if stoichiometry == 1:
                term = species_name
            else:
                term = f"{stoichiometry} {species_name}"

            formula_list.append(term)

        return " + ".join(formula_list)



    def _parse_reaction_term(self, term: Union[str, tuple, list], name="term") -> (int, str, int):
        """
        Accept various ways to specify a reaction term, and return a standardized triplet form for it.

        NOTE:   A) if the stoichiometry coefficient isn't specified, it defaults to 1
                B) if the reaction order isn't specified, it defaults to the stoichiometry coefficient

        In the passed tuples or lists:
            - required 1st entry is the stoichiometry
            - required 2nd entry is the chemical name
            - optional 3rd one is the reaction order.  If unspecified, it defaults to the stoichiometry

        If just a string is being passed, it is taken to be the chemical name,
        with stoichiometry and reaction order both 1

        EXAMPLES:
            "F"          gets turned into:  (1, "F", 1)   - defaults used for stoichiometry and reaction order
            (2, "F")                        (2, "F", 2)   - default used for reaction order
            (2, "F", 1)                     (2, "F", 1)   - no defaults invoked
            It's equally acceptable to use LISTS in lieu of tuples

        :param term:    A string (a chemical name)
                            OR  a pair (stoichiometry coeff, name)
                            OR  a triplet (stoichiometry coeff, name, reaction order)
        :param name:    An optional nickname, handy to refer to this term in error messages if needed
                            (for example, "reactant" or "product")
        :return:        A standardized triplet of the form (stoichiometry, species_name, reaction_order),
                            where stoichiometry and reaction_order are integers, while species_name is a string
        """
        if type(term) == str:
            return  (1, term, 1)    # Accept simply the chemical name as a shortcut,
                                    # for when the stoichiometry coefficient and reaction order are both 1

        if type(term) != tuple and type(term) != list:
            raise Exception(f"_parse_reaction_term(): {name} must be either a string (a chemical name), "
                            f"or a pair (stoichiometry coeff, name) or a triplet (stoichiometry coeff, name, reaction order). "
                            f"Instead, it is `{term}` (of type {type(term)})")

        # If we get thus far, term is either a tuple or a list
        assert len(term) in [2, 3],  \
            f"_parse_reaction_term(): Unexpected length for {name} tuple/list: it should be either 2 or 3. " \
            f"Instead, it is {len(term)}"

        stoichiometry = term[0]
        assert type(stoichiometry) == int, \
            f"_parse_reaction_term(): The stoichiometry coefficient must be an integer. Instead, it is {stoichiometry}"

        species_name = term[1]
        assert type(species_name) == str, \
                            f"_parse_reaction_term(): The chemical name must be a string. " \
                            f"Instead, it is `{species_name}` (of type {type(species_name)})"

        if len(term) == 2:
            return (stoichiometry, species_name, stoichiometry)
        else:   # Length is 3
            reaction_order = term[2]
            return (stoichiometry, species_name, reaction_order)
