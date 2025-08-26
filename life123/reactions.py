# 3 classes: "ReactionEnz", "ReactionGeneric" and "Reactions"

from typing import Union, Set, Tuple
import numpy as np
from life123.thermodynamics import ThermoDynamics
from life123.reaction_kinetics import ReactionKinetics
from life123.visualization.py_graph_visual import PyGraphVisual
from life123.visualization.graphic_log import GraphicLog
from life123.html_log import HtmlLog as log



###################################################################################################################

class ReactionCommon:
    """
    Base class for all individual reactions.
    Typically NOT instantiated by the user
    """
    def __init__(self, active=True, reversible=True):
        self.active = active              # TODO: not yet in use
        self.reversible = reversible



class ReactionElementary(ReactionCommon):
    """
    Base class for all elementary reactions (i.e. with no intermediaries).
    Typically NOT instantiated by the user
    """
    def __init__(self, kF=None, kR=None,
                 delta_H=None, delta_S=None, delta_G=None, **kwargs):

        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        if not self.reversible:
            assert not kR, \
                f"ReactionElementary instantiation: irreversible reactions cannot have a value for the reverse rate constant (kR = {kR})"

        self.kF = kF                # Forward rate constant
        self.kR = kR                # Reverse rate constant
        self.delta_H = delta_H
        self.delta_S = delta_S
        self.delta_G = delta_G
        self.K = None               # Equilibrium constant

        if (kR is not None) and not np.allclose(self.kR, 0):
            self.K = kF / kR



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





###################################################################################################################

class ReactionUnimolecular(ReactionElementary):
    """
    Reactions of type A <-> B

    NOT YET IN FULL USE for simulations.  Use the class "ReactionGeneric" for now
    """
    def __init__(self, reactant :str, product :str, **kwargs):
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        assert type(reactant) == str, "ReactionUnimolecular instantiation: argument `reactant` must be a string"
        assert type(product) == str, "ReactionUnimolecular instantiation: argument `product` must be a string"
        assert reactant != product, "ReactionUnimolecular instantiation: the `reactant` and the `product` cannot be identical"

        self.reactant = reactant
        self.product = product




    def describe(self, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction

        :param concise:     If True, less detail is shown
        :return:            A string with a description of the specified reaction
        """
        if self.reversible:
            return f"{self.reactant} <-> {self.product}  (Elementary Unimolecular reaction)"
        else:
            return f"{self.reactant} -> {self.product}  (Elementary Unimolecular Ir-reversible reaction)"


    def extract_reactant_names(self) -> [str]:
        # TODO: maybe rename to extract_reactant_labels()
        return [self.reactant]

    def extract_product_names(self) -> [str]:
        return [self.product]



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the chemical labels of all the chemicals appearing in this reaction.

        :return:    A SET of the labels of the chemicals involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return {self.reactant, self.product}



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the chemicals in the unimolecular reaction,
        determine its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the start of the time step.

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
        Simulate the unimolecular reaction, over the specified time interval

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in the given reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6}

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn is the mapping of chemical label to their concentration changes
                                                            during this step
                                - rxn_rate                  is the reaction rate ("velocity") for this reaction
                                EXAMPLE of increment_dict_single_rxn: {"B": -0.2, "F": 0.2}
        """

        increment_dict_single_rxn = {}      # The keys are the chemical labels,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        delta_rxn = rxn_rate * delta_time


        # Determine the concentration adjustments as a result of this reaction step:

        # The reactants DECREASE based on the quantity (forward reaction - reverse reaction)
        r = self.reactant           # EXAMPLE: "B"
        # stoichiometry = 1
        delta_conc = - delta_rxn    # Increment to this reactant from the reaction step
        increment_dict_single_rxn[r] = delta_conc


        # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
        p = self.product            # EXAMPLE: "F"
        # stoichiometry = 1
        delta_conc = delta_rxn      # Increment to this reaction product from the reaction step
        increment_dict_single_rxn[p] = delta_conc

        return (increment_dict_single_rxn, rxn_rate)




#######################################################################################################################

class ReactionSynthesis(ReactionElementary):
    """
    Reactions of type A + B <-> C

    NOT YET IN FULL USE for simulations.  Use the class "ReactionGeneric" for now
    """
    def __init__(self, reactants :(str, str), product :str, **kwargs):

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
        :return:            A string with a description of the specified reaction
        """
        if self.reversible:
            return f"{self.reactant_1} + {self.reactant_2}  <-> {self.product}  (Elementary Synthesis reaction)"
        else:
            return f"{self.reactant_1} + {self.reactant_2}  -> {self.product}  (Elementary Synthesis Ir-reversible reaction)"



    def extract_reactant_names(self) -> [str]:
        # TODO: maybe rename to extract_reactant_labels()
        return [self.reactant_1, self.reactant_2]

    def extract_product_names(self) -> [str]:
        return [self.product]



    def extract_chemicals_in_reaction(self) -> Set[str]:
        """
        Return a SET of the chemical labels of all the chemicals appearing in this reaction.

        :return:    A SET of the labels of the chemicals involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        return {self.reactant_1, self.reactant_2, self.product}



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
                                EXAMPLE of increment_dict_single_rxn: "A": -1.3, "B": 2.9, "C": -1.6
        """

        increment_dict_single_rxn = {}      # The keys are the chemical labels,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        delta_rxn = rxn_rate * delta_time


        # Determine the concentration adjustments as a result of this reaction step:

        # The reactants DECREASE based on the quantity (forward reaction - reverse reaction)
        for r in [self.reactant_1, self.reactant_2]:
            # EXAMPLE of r: "A"
            # stoichiometry = 1
            delta_conc = - delta_rxn    # Increment to this reactant from the reaction step
            increment_dict_single_rxn[r] = delta_conc


        # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
        p = self.product            # EXAMPLE: "F"
        # stoichiometry = 1
        delta_conc = delta_rxn      # Increment to this reaction product from the reaction step
        increment_dict_single_rxn[p] = delta_conc

        return (increment_dict_single_rxn, rxn_rate)





#######################################################################################################################

class ReactionDecomposition(ReactionElementary):
    """
    Reactions of type A <-> B + C

    NOT YET IN FULL USE for simulations.  Use the class "ReactionGeneric" for now
    """
    def __init__(self, reactant :str, products :(str, str), **kwargs):
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        assert type(products) == list or type(products) == tuple, \
            "ReactionSynthesis instantiation: argument `reactants` must be a list or tuple"
        assert len(products) == 2, \
            "ReactionSynthesis instantiation: argument `reactants` must be a pair"
        assert type(reactant) == str, \
            "ReactionSynthesis instantiation: argument `product` must be a string"

        (p1, p2) = products
        assert (p1 != reactant) and (p2 != reactant), \
            "ReactionSynthesis instantiation: the `reactant` cannot be identical to any of the reaction products"

        self.reactant = reactant
        self.product_1 = p1
        self.product_2 = p2




    def describe(self, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the reaction

        :param concise:     If True, less detail is shown
        :return:            A string with a description of the specified reaction
        """
        if self.reversible:
            return f"{self.reactant} <-> {self.product_1} + {self.product_2} (Elementary Decomposition reaction)"
        else:
            return f"{self.reactant} -> {self.product_1} + {self.product_2} (Elementary Decomposition Ir-reversible reaction)"



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
                                EXAMPLE of increment_dict_single_rxn: "A": -1.3, "B": 2.9, "C": -1.6
        """

        increment_dict_single_rxn = {}      # The keys are the chemical labels,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        delta_rxn = rxn_rate * delta_time


        # Determine the concentration adjustments as a result of this reaction step:

        # The reactants DECREASE based on the quantity (forward reaction - reverse reaction)
        r = self.reactant            # EXAMPLE: "F"
        # EXAMPLE of r: "A"
        # stoichiometry = 1
        delta_conc = - delta_rxn    # Increment to this reactant from the reaction step
        increment_dict_single_rxn[r] = delta_conc


        # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
        for p in [self.product_1, self.product_1]:
            # EXAMPLE of p: "B"
            # stoichiometry = 1
            delta_conc = delta_rxn      # Increment to this reaction product from the reaction step
            increment_dict_single_rxn[p] = delta_conc

        return (increment_dict_single_rxn, rxn_rate)




###################################################################################################################

class ReactionEnz(ReactionCommon):
    """
    Data about a SINGLE enzyme-catalyzed reaction that can be modeled kinetically as:

    E + S <-> ES <-> E + P

        E : Enzyme
        S : Substrate
        ES: Intermediate Enzyme-Substrate complex
        P : Product
    """

    def __init__(self, enzyme=None, substrate=None, product=None,
                 k1_F=None, k1_R=None, k2_F=None, k2_R=None,
                 kM=None, kcat=None, **kwargs):
        """
        :param enzyme:      The label for the chemical acting as enzyme
        :param substrate:
        :param product:
        :param k1_F:
        :param k1_R:
        :param k2_F:
        :param k2_R:
        :param kM:          Michaelis constant
        :param kcat:
        """
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        self.enzyme = enzyme
        self.substrate = substrate
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
                    f"Inconsistent values passed during instantiation of ReactionEnz.  " \
                    f"The passed kM value ({kM}) doesn't the value ({self.kM}) inferred from the given reaction rate constants"

        if k2_F is not None:
            self.kcat = k2_F
            if kcat is not None:
                assert np.allclose(self.kcat, kcat), \
                    f"Inconsistent values passed during instantiation of ReactionEnz.  " \
                    f"The passed kcat value ({kcat}) doesn't the value ({self.kcat}) of the given k2_F reaction rate constants"



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





###################################################################################################################

class ReactionGeneric(ReactionCommon):
    """
    Data about a generic SINGLE reaction of the most general type,
    with arbitrary number of reactants and products,
    arbitrary stoichiometry,
    and arbitrary kinetic reaction orders for each participating chemical.

    Included:
        - stoichiometry
        - kinetic data (reaction rates, reaction orders)
        - thermodynamic data (temperature, changes in enthalpy/entropy/Gibbs Free Energy)
        - list of involved enzymes


    (Note: this data will eventually be stored in a graph database)

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

    def __init__(self, reactants: Union[int, str, list], products: Union[int, str, list],
                 kF=None, kR=None,
                 delta_H=None, delta_S=None, delta_G=None, temp=None, **kwargs):
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
        :param temp:        [OPTIONAL] Temperature in Kelvins.  For now, assumed constant everywhere,
                                and unvarying (or very slowly varying)
        """
        super().__init__(**kwargs)          # Invoke the constructor of its parent class

        self.reactants = None
        self.products = None
        self.kF = kF
        self.kR = kR
        self.delta_H = delta_H
        self.delta_S = delta_S
        self.delta_G = delta_G
        self.K = None               # Equilibrium constant
        self.enzyme = None          # The INDEX of a chemical that catalyzes this reaction, if applicable
                                    #   Note: enzymes are automatically extracted from the reaction formula
        self.macro_enzyme = None    # The pair (macromolecule name, binding site number)
                                    #   EXAMPLE: ("M2", 5)          TODO: maybe turn into a data object

        assert reactants is not None, "Reaction(): the argument `reactants` is a required one; it can't be None"
        if type(reactants) != list:
            reactants = [reactants]


        assert products is not None, "Reaction(): the argument `products` is a required one; it can't be None"
        if type(products) != list:
            products = [products]


        reactant_list = [self._parse_reaction_term(r, "reactant") for r in reactants]   # A list of triples of integers
        product_list = [self._parse_reaction_term(r, "product") for r in products]      # A list of triples of integers

        # Catch identical reaction sides, even if terms are reshuffled
        assert set(reactant_list) != set(product_list), \
            f"Reaction(): the two sides of the reaction can't be identical! " \
            f"Same reactants and products: {self._standard_form_chem_eqn(reactant_list)}"


        # Locate any enzymes (catalysts) that might be present - though for now a warning is issued if more than 1
        enzyme_list = []
        for reactant in reactant_list:
            if reactant in product_list:
                enzyme_list.append(self.extract_species_name(reactant))

        number_enzymes = len(enzyme_list)

        if number_enzymes == len(reactant_list) or number_enzymes == len(product_list):
            raise Exception(f"Reaction(): all the terms in the reaction appear to be enzymes!  "
                            f"Enzymes: {enzyme_list}")


        self.reactants = reactant_list
        self.products = product_list

        if number_enzymes >= 1:
            self.enzyme = enzyme_list[0]    # In the irregular scenarios that there appear to be multiple enzymes, only one
                                            #   is considered, and a warning is printed out (the other apparent enzyme
                                            #   will be treated as any other reagent/product)
        if number_enzymes > 1:
            print(f"Reaction(): WARNING - the reaction appears to have multiple enzymes:"
                  f" {enzyme_list}")


        # Process the kinetic and thermodynamic data, and update various object attributes accordingly
        self._set_kinetic_and_thermodynamic(forward_rate=kF, reverse_rate=kR,
                                            delta_H=delta_H, delta_S=delta_S, delta_G=delta_G, temp=temp)



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
        incl. their stoichiometry, name and reaction order

        :return:    A list of triplets of the form (stoichiometry, species name, reaction order)
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
        incl. their stoichiometry, name and reaction order

        :return:    A list of triplets of the form (stoichiometry, species name, reaction order)
        """
        return self.products



    def extract_products_formula(self) -> str:
        """
        Return a string with a user-friendly form of the right (products) side of the reaction formula

        :return:    A string with one side of a chemical reaction
        """
        products = self.extract_products()
        return self._standard_form_chem_eqn(products)



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



    def unpack_for_dynamics(self) -> tuple:
        """
        A convenient unpacking meant for dynamics simulations
        that need the reactants, the products, and the forward and reverse rate constants

        :return:    A 4-element tuple, containing:
                        (reactants , products , forward rate constant , reverse rate constant)
                        Note: both reactants products are lists of triplets
        """
        return (self.reactants, self.products, self.kF, self.kR)



    def extract_stoichiometry(self, term :(int, str, int)) -> int:
        """
        Return the stoichiometry coefficient, from a reaction TERM

        :param term:    A triplet (int, str, int) representing a reaction term
        :return:        An integer with the stoichiometry coefficient
        """
        return term[0]

    def extract_species_name(self, term :(int, str, int)) -> str:
        """
        Return the name of the chemical species, from a reaction TERM

        :param term:    A triplet (int, str, int) representing a reaction term
        :return:        The name of the chemical species in the term
        """
        return term[1]

    def extract_rxn_order(self, term :(int, str, int)) -> int:
        """
        Return the reaction order, from a reaction TERM

        :param term:    A triplet (int, str, int) representing a reaction term
        :return:        An integer with the reaction order for this term
        """
        return term[2]



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



    def extract_chemicals_in_reaction(self, exclude_enzyme=False) -> Set[str]:
        """
        Return a SET of the chemical labels of all the chemicals appearing in this reaction.

        Optionally, exclude any that participate in a catalytic role
        (appearing identically on both sides of the reaction)

        :param exclude_enzyme:  If True, any enzyme, if present, won't be included
        :return:                A SET of the labels of the chemicals involved in this reaction
                                Note: being a set, it's NOT in any particular order
        """
        chem_set = set()    # Running set being built

        reactants = self.extract_reactants()
        products = self.extract_products()

        for r in reactants:
            species_name = self.extract_species_name(r)
            chem_set.add(species_name)

        for p in products:
            species_name = self.extract_species_name(p)
            chem_set.add(species_name)

        if exclude_enzyme:
            chem_set = chem_set - {self.enzyme}     # Difference between sets

        return chem_set



    def extract_reactant_names(self, exclude_enzyme=False) -> [str]:
        """
        In the order in which they appear when the reaction was first defined

        :param exclude_enzyme:  If True, any enzyme, if present, won't be included
        :return:                List of chemical names
        """
        reactants = self.extract_reactants()
        reactant_names = [self.extract_species_name(r) for r in reactants]

        if exclude_enzyme:
            reactant_names.remove(self.enzyme)

        return reactant_names


    def extract_product_names(self, exclude_enzyme=False) -> [str]:
        """
        In the order in which they appear when the reaction was first defined

        :param exclude_enzyme:  If True, any enzyme, if present, won't be included
        :return:                List of chemical names
        """
        products = self.extract_products()
        product_names = [self.extract_species_name(r) for r in products]

        if exclude_enzyme:
            product_names.remove(self.enzyme)

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
        details = []
        rxn_properties = self.extract_rxn_properties()
        for k,v in rxn_properties.items():
            details.append(f"{k} = {v:,.5g}")          # EXAMPLE: "kF = 3"

        if details:
            rxn_description += "  (" + ' / '.join(details) + ")"    # EXAMPLE: "  (kF = 3 / kR = 2 / Delta_G = -1,005.13)"


        # If an ENZYME is involved, show it
        if self.enzyme is not None:
            rxn_description += f" | Enzyme: {self.enzyme}"

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
                            The keys are the chemical names
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
        :param explain:     If True, it also returns the math formula being used for the computation
                                EXAMPLES:   "([C][D]) / ([A][B])"
                                            "[B] / [A]^2"

        :return:            If explain is False, return value for the "Reaction Quotient" (aka "Mass–action Ratio");
                                if True, return a pair with that quotient and a string with the math formula that was used.
                                Note that the reaction quotient is a Numpy scalar that might be np.inf or np.nan
        """
        numerator = np.double(1)    # The product of all the concentrations of the reaction products (adjusted for reaction order)
        denominator = np.double(1)  # The product of all the concentrations of the reactants (also adjusted for reaction order)

        numerator_text = ""      # First part of the the textual explanation
        denominator_text = ""    # Second part of the the textual explanation


        # Compute the numerator of the "Reaction Quotient"
        for p in self.products:
            # Loop over the reaction products
            species_name = self.extract_species_name(p)
            rxn_order = self.extract_rxn_order(p)

            species_conc = conc.get(species_name)
            assert species_conc is not None, f"reaction_quotient(): unable to proceed because the " \
                                             f"concentration of `{species_name}` was not provided"

            numerator *= (species_conc ** rxn_order)
            if explain:
                numerator_text += f"[{species_name}]"
                if rxn_order > 1:
                    numerator_text += f"^{rxn_order} "

        if explain and len(self.products) > 1:
            numerator_text = f"({numerator_text})"  # In case of multiple terms, enclose them in parenthesis


        # Compute the denominator of the "Reaction Quotient"
        for r in self.reactants:
            # Loop over the reactants
            species_name =  self.extract_species_name(r)
            rxn_order =  self.extract_rxn_order(r)

            species_conc = conc.get(species_name)
            assert species_conc is not None, f"reaction_quotient(): unable to proceed because the " \
                                             f"concentration of `{species_name}` was not provided"

            denominator *= (species_conc ** rxn_order)
            if explain:
                denominator_text += f"[{species_name}]"
                if rxn_order > 1:
                    denominator_text += f"^{rxn_order} "

        if explain and len(self.reactants) > 1:
            denominator_text = f"({denominator_text})"  # In case of multiple terms, enclose them in parenthesis


        with np.errstate(divide='ignore', invalid='ignore'):
            # It might be np.inf (if just the denominator is zero) or np.nan (if both are zero)
            quotient = numerator / denominator

        if explain:
            formula = f"{numerator_text} / {denominator_text}"
            return (quotient, formula)

        return quotient




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

        delta_rxn = rxn_rate * delta_time


        reactants = self.extract_reactants() # A list of triplets of the form (stoichiometry, species name, reaction order)
        products = self.extract_products()   # A list of triplets of the form (stoichiometry, species name, reaction order)


        """
        Determine the concentration adjustments as a result of this reaction step, 
        for this individual reaction being considered
        """

        # The reactants DECREASE based on the quantity (forward reaction - reverse reaction)
        for r in reactants:
            # Unpack data from the reactant r
            species_name = self.extract_species_name(r)
            if species_name == self.enzyme:
                #print(f"*** SKIPPING reactant CATALYST {species_name} in reaction")
                continue    # Skip if r is a catalyst for this reaction

            stoichiometry = self.extract_stoichiometry(r)

            delta_conc = stoichiometry * (- delta_rxn)  # Increment to this reactant from the reaction being considered

            increment_dict_single_rxn[species_name] = increment_dict_single_rxn.get(species_name,0) + delta_conc


        # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
        for p in products:
            # Unpack data from the reactant r
            species_name = self.extract_species_name(p)
            if species_name == self.enzyme:
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





###################################################################################################################

class Reactions:
    """
    Manage reaction-related data

    (this class was formerly called "AllReactions")
    """

    def __init__(self, chem_data):

        """

        :param chem_data:   Object of type "ChemData"
        """

        # TODO: consider adding arguments   =None, labels=None
        """
        assert (chem_data is not None) or (labels is not None), \
            "Reactions() instantiation: exactly one of the arguments `chem_data` or `labels` must be provided"

        assert (chem_data is None) or (labels is None), \
            "Reactions() instantiation: cannot specify both the arguments `chem_data` and `labels`"

        if labels is not None:
            chem_data = ChemData(labels=labels)
        """

        assert chem_data is not None, \
            "Reactions() instantiation: the arguments `chem_data` must be provided, and cannot be None"

        self.chem_data = chem_data

        self.reaction_list = []     # List of objects of the various individual reaction classes,
                                    # such as "ReactionGeneric" and "ReactionEnz"

        self.temp = 298.15          # Temperature in Kelvins.  (By default, the equivalent of 25 C)
                                    # For now, assumed constant everywhere, and unvarying (or very slowly varying)

        self.active_chemicals = set()   # Set of the names of the chemicals - not counting pure catalysts - involved
                                        # in any of the registered reactions
                                        # CAUTION: the concept of "active chemical" might change in future versions, where only SOME of
                                        #          the reactions are simulated

        self.active_enzymes = set()     # Set of the names of the enzymes (catalysts) involved
                                        # in any of the registered reactions
                                        # CAUTION: the concept of "active enzyme" might change in future versions, where only SOME of
                                        #          the reactions are simulated



    def number_of_reactions(self, include_inactive=False) -> int:
        """
        Return the number of registered chemical reactions

        :param include_inactive:    If True, disabled reactions are also included
        :return:                    The number of registered chemical reactions
        """
        if include_inactive:
            return len(self.reaction_list)

        count = 0
        for rxn in self.reaction_list:
            if rxn.active:
                count += 1

        return count



    def active_reaction_indices(self) -> [int]:
        """
        Return a list of the reaction index numbers of all the active reactions

        :return:    A list of integers, to identify the active reactions by their indices
        """
        l = []
        for i, rxn in enumerate(self.reaction_list):
            if rxn.active:
                l.append(i)

        return l



    def assert_valid_rxn_index(self, index :int) -> None:
        """
        Raise an Exception if the specified reaction index isn't valid

        :param index:   An integer that indexes the reaction of interest (numbering starts at 0)
        :return:        None
        """
        assert self.number_of_reactions() > 0, \
            f"ChemData.assert_valid_rxn_index(): there are no reactions defined yet.  Use add_reaction() to add them first"

        assert (type(index) == int), \
            f"ChemData.assert_valid_rxn_index(): the requested reaction index must be an integer; " \
            f"the provided value ({index}) is of type {type(index)}"

        assert 0 <= index < self.number_of_reactions(), \
            f"ChemData.assert_valid_rxn_index(): the requested reaction index is not the expected range [0 to {self.number_of_reactions() - 1}], inclusive; " \
            f"the value passed was: {index} (there is no reaction whose index is {index})"



    def get_chem_data(self):
        return self.chem_data



    def get_reaction(self, i :int):
        """
        Return the data structure of the i-th reaction,
        in the order in which reactions were added (numbering starts at 0)

        :param i:   An integer that indexes the reaction of interest (numbering starts at 0)
        :return:    An object of one of the various individual reaction classes,
                    such as "ReactionGeneric" and "ReactionEnz"
        """
        self.assert_valid_rxn_index(i)

        return self.reaction_list[i]



    def get_reactants(self, i :int) -> [(int, int, int)]:
        """
        Return a list of triplets with details of the reactants of the i-th reaction

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A list of triplets of the form (stoichiometry, species index, reaction order)
        """
        rxn = self.get_reaction(i)
        return rxn.extract_reactants()



    def get_reactants_formula(self, i :int) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A string with a user-friendly form of the left (reactants) side of the chemical reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_reactants_formula()



    def get_products(self, i :int) -> [(int, int, int)]:
        """
        Return a list of triplets with details of the products of the i-th reaction

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A list of triplets of the form (stoichiometry, species index, reaction order)
        """
        rxn = self.get_reaction(i)
        return rxn.extract_products()


    def get_products_formula(self, i :int) -> str:
        """
        Return a string with a user-friendly form of the right (products) side of the reaction formula

        :param i:   The index (0-based) to identify the reaction of interest
        :return:    A string with a user-friendly form of the right (products) side of the chemical reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_products_formula()



    def get_forward_rate(self, i :int) -> float:
        """
        Return the value of the forward rate constant of the i-th reaction

        :param i:   The integer index (0-based) to identify the reaction of interest
        :return:    The value of the forward rate constant for the above reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_forward_rate()


    def get_reverse_rate(self, i :int) -> float:
        """
        Return the value of the reverse (back) rate constant of the i-th reaction

        :param i:   The integer index (0-based) to identify the reaction of interest
        :return:    The value of the reverse (back) rate constant for the above reaction
        """
        rxn = self.get_reaction(i)
        return rxn.extract_reverse_rate()



    def get_chemicals_in_reaction(self, rxn_index :int) -> {int}:
        """
        Return a SET of indices (being a set, they're NOT in any particular order)
        of all the chemicals participating in the i-th reaction

        :param rxn_index:   An integer with the (zero-based) index to identify the reaction of interest
        :return:            A SET of indices of the chemicals involved in the above reaction
                                Note: being a set, it's NOT in any particular order
        """
        rxn = self.get_reaction(rxn_index)

        name_set = rxn.extract_chemicals_in_reaction()

        index_set = {self.chem_data.get_index(name) for name in name_set}

        return index_set


    def get_chemicals_indexes_in_reaction(self, rxn_index :int) -> [int]:
        """
        Return a sorted list of the indexes
        of all the chemicals participating in the i-th reaction

        :param rxn_index:   An integer with the (zero-based) index to identify the reaction of interest
        :return:            A sorted list of indices of the chemicals involved in the above reaction
        """
        index_set = self.get_chemicals_in_reaction(rxn_index)
        index_list = list(index_set)

        return sorted(index_list)



    def get_reactions_participating_in_TODO(self, species_index :int) -> [ReactionGeneric]:
        """
        Return a list of all the reactions that the given chemical species
        is involved in

        :param species_index:
        :return:                List of "Reaction" objects
        """
        pass        # TODO: write; also, accept a name



    def set_temp(self, temp :Union[float, int], units="K") -> None:
        """
        Specify the temperature of the environment
        (for now assumed uniform everywhere)

        :param temp:    Temperature, in Kelvins, or None
        :param units:   Not yet implemented
        :return:        None
        """
        self.temp = temp



    def add_reaction_from_object(self, rxn) -> int:
        """
        Register a new SINGLE chemical reaction,
        and set all its kinetic and/or thermodynamic data from the available information,
        including the value of the temperature (stored in object variable.)

        All the involved chemicals can be either previously registered, or not.

        :param rxn: One of the specific Reaction classes, such as
                        ReactionUnimolecular, ReactionSynthesis, ReactionDecomposition,
                        ReactionEnz, ReactionGeneric
        :return:    Integer index of the newly-added reaction
                        (in the list self.reaction_list, stored as object variable)
        """
        self.reaction_list.append(rxn)

        # Register any newly-encountered reactant not already registered
        rxn_reactants = rxn.extract_reactant_names()
        for label in rxn_reactants:
            if not self.chem_data.label_exists(label):
                self.chem_data.add_chemical(name=label)

        # Register any newly-encountered reaction product not already registered
        # (reactants are done first, because that's typically a more appealing order of appearance)
        rxn_products = rxn.extract_product_names()
        for label in rxn_products:
            if not self.chem_data.label_exists(label):
                self.chem_data.add_chemical(name=label)


        self.active_chemicals = self.active_chemicals | set(rxn_reactants) | set(rxn_products)     # Union of sets

        rxn.set_thermodynamic_data(temp=self.temp)

        return len(self.reaction_list) - 1



    def add_reaction(self, reactants :Union[int, str, list], products :Union[int, str, list],
                     forward_rate=None, reverse_rate=None,
                     delta_H=None, delta_S=None, delta_G=None) -> int:
        """
        Register a new SINGLE chemical reaction,
        optionally including its kinetic and/or thermodynamic data.
        All the involved chemicals can be either previously registered, or not.

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

        :param reactants:       A list of triplets (stoichiometry, species name, reaction order),
                                    or simplified terms in various formats; for details, see above.
                                    If not a list, it will get turned into one
        :param products:        A list of triplets (stoichiometry, species name, reaction order of REVERSE reaction),
                                    or simplified terms in various formats; for details, see above.
                                    If not a list, it will get turned into one
        :param forward_rate:    [OPTIONAL] Forward reaction rate constant
        :param reverse_rate:    [OPTIONAL] Reverse reaction rate constant
        :param delta_H:         [OPTIONAL] Change in Enthalpy (from reactants to products)
        :param delta_S:         [OPTIONAL] Change in Entropy (from reactants to products)
        :param delta_G:         [OPTIONAL] Change in Free Energy (from reactants to products)

        :return:                Integer index of the newly-added reaction
                                    (in the list self.reaction_list, stored as object variable)
        """
        rxn = ReactionGeneric(reactants, products, forward_rate, reverse_rate,
                              delta_H, delta_S, delta_G, temp=self.temp)
        self.reaction_list.append(rxn)

        # Register any newly-encountered reactant not already registered
        rxn_reactants = rxn.extract_reactant_names(exclude_enzyme=False)
        for label in rxn_reactants:
            if not self.chem_data.label_exists(label):
                self.chem_data.add_chemical(name=label)

        # Register any newly-encountered reaction product not already registered
        # Note: reactants are done first, because that's typically a more appealing order of appearance
        rxn_products = rxn.extract_product_names(exclude_enzyme=False)
        for label in rxn_products:
            if not self.chem_data.label_exists(label):
                self.chem_data.add_chemical(name=label)


        # TODO: since we already have rxn_reactants, rxn_products and rxn.enzyme, the following could be simplified!
        involved_chemicals = rxn.extract_chemicals_in_reaction(exclude_enzyme=True)

        self.active_chemicals = self.active_chemicals.union(involved_chemicals)     # Union of sets
        if rxn.enzyme is not None:
            self.active_enzymes.add(rxn.enzyme)       # Add the new entry to a set

        return len(self.reaction_list) - 1



    def clear_reactions_data(self) -> None:
        """
        Get rid of all the reactions; start again with "an empty slate" (but still with reference
        to the same data object about the chemicals and their properties)

        :return:    None
        """
        self.reaction_list = []
        self.active_chemicals = set()
        self.active_enzymes = set()



    def inactivate_reaction(self, i :int) -> None:
        """
        Mark the i-th reaction as "inactive/disabled" (essentially, "deleted", but holding its positional
        index, to avoid a change in index in other reactions)

        TODO: Not yet supported by the dynamical modules; DON'T USE YET in simulations!

        :param i:   Zero-based index of the reaction to disable
        :return:    None
        """
        rxn = self.get_reaction(i)
        rxn.active = False

        # Re-construct self.active_chemicals and self.active_enzymes
        self.active_chemicals = set()
        self.active_enzymes = set()

        for rxn in self.reaction_list:
            involved_chemicals = rxn.extract_chemicals_in_reaction(exclude_enzyme=True)
            self.active_chemicals = self.active_chemicals.union(involved_chemicals)     # Union of sets
            if rxn.enzyme is not None:
                self.active_enzymes.add(rxn.enzyme)       # Add the new entry to a set





    #####################################################################################################

    '''                             ~   TO DESCRIBE THE REACTIONS  ~                                  '''

    def ________DESCRIBE_RXNS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def describe_reactions(self, concise=False) -> None:
        """
        Print out a user-friendly plain-text form of ALL the reactions.
        If wanting to describe just 1 reaction, use single_reaction_describe()

        EXAMPLE (not concise):
            Number of reactions: 2 (at temp. 25 C)
            (0) CH4 + 2 O2 <-> CO2 + 2 H2O  (kF = 3.0 / kR = 2.0 / Delta_G = -1,005.13 / K = 1.5) | 1st order in all reactants & products
            (1) A + B <-> C  (kF = 5.0 / kR = 1.0 / Delta_G =  / K = 5.0) | 1st order in all reactants & products
            Set of chemicals involved in the above reactions: {'CH4', 'O2', 'H2O', 'A', 'B', 'C'}

        :param concise:     If True, less detail is shown
        :return:            None
        """
        print(f"Number of reactions: {self.number_of_reactions()} (at temp. {self.temp - 273.15:,.4g} C)")

        # Print a concise description of each reaction in turn
        for description in self.multiple_reactions_describe(concise=concise):
            print(description)

        chem_labels = self.labels_of_active_chemicals(sort_by_index=True)   # Set of chem labels, sorted by chemical index

        if self.chem_data.color_dict != {}:   # If plot colors were registered, show them alongside the chem labels
            chem_labels_with_colors = []
            for label in chem_labels:
                color = self.chem_data.get_plot_color(label)
                if color:
                    chem_labels_with_colors.append(f'"{label}" ({color})')
                else:
                    chem_labels_with_colors.append(f'"{label}"')

            chem_labels = "{" + ", ".join(chem_labels_with_colors) + "}"


        if self.active_enzymes == set():    # If no enzymes were involved in any reaction
            print(f"Chemicals involved in the above reactions: {chem_labels}")
        else:
            print(f"Chemicals involved in the above reactions (not counting enzymes): {chem_labels}")
            print(f"Enzymes involved in the above reactions: "
                  f"{self.names_of_enzymes()}")



    def multiple_reactions_describe(self, rxn_list=None, concise=False) -> [str]:
        """
        The counterpart of single_reaction_describe() for many reactions.
        Return a list of strings, each string being a (concise or not) user-friendly plain-text form of
        each of the reactions

        :param rxn_list:    Either a list of integers, to identify the reactions of interest,
                                or None, meaning ALL reactions
        :param concise:     If True, less detail is shown
        :return:            A list of strings; each string is the description of one of the reactions
        """
        if rxn_list is None:
            rxn_list = range(self.number_of_reactions())    # Show ALL reactions, by default

        out = []    # Output list being built

        for i in rxn_list:
            description = self.single_reaction_describe(rxn_index=i, concise=concise)
            description = f"{i}: {description}"
            out.append(description)

        return out



    def single_reaction_describe(self, rxn_index: int, concise=False) -> str:
        """
        Return as a string, a user-friendly plain-text form of the given reaction
        EXAMPLE (concise):      "CH4 + 2 O2 <-> CO2 + 2 H2O"
        EXAMPLE (not concise):  "CH4 + 2 O2 <-> CO2 + 2 H2O  (kF = 3.0 / kR = 2.0 / Delta_G = -1,005.13 / K = 1.5) | 1st order in all reactants & products"

        :param rxn_index:   Integer to identify the reaction of interest
        :param concise:     If True, less detail is shown
        :return:            A string with a description of the specified reaction
        """
        rxn = self.get_reaction(rxn_index)

        return rxn.describe(concise)



    def labels_of_active_chemicals(self, sort_by_index=False) -> [str]:
        """
        Return a list of the labels of all the chemicals
        involved in ANY of the registered reactions,
        but NOT counting chemicals that always appear
        in a catalytic role in all the reactions they participate in
        (if a chemical participates in a non-catalytic role in ANY reaction, it'll appear here)

        The list is not in any particular order, unless sort_by_index is True

        :param sort_by_index:   If True, the list is sorted by the index (order of registration)
                                    of the chemicals in it
        :return:                A set of chemical labels
        """
        if not sort_by_index:
            return list(self.active_chemicals)

        return sorted(self.active_chemicals, key=self.chem_data.get_index)



    def number_of_active_chemicals(self) -> int:
        """
        Return the number of all the chemicals
        involved in ANY of the registered reactions,
        but NOT counting chemicals that always appear
        in a catalytic role in all the reactions they participate in
        (if a chemical participates in a non-catalytic role in ANY reaction, it'll appear here)
        """
        return len(self.active_chemicals)



    def indexes_of_active_chemicals(self) -> [int]:
        """
        Return the ordered list (numerically sorted) of the INDEX numbers of all the chemicals
        involved in ANY of the registered reactions,
        but NOT counting chemicals that always appear in a catalytic role in all the reactions they
        participate in
        (if a chemical participates in a non-catalytic role in ANY reaction, it'll appear here.)

        EXAMPLE: [2, 7, 8]  if only those 3 chemicals (with indexes of, respectively, 2, 7 and 8)
                            are actively involved in ANY of the registered reactions

        CAUTION: the concept of "active chemical" might change in future versions, where only SOME of
                 the reactions are simulated
        """
        index_list = list(map(self.chem_data.get_index, self.active_chemicals))
        return sorted(index_list)



    def names_of_enzymes(self) -> Set[str]:
        """
        Return the set of the names of the enzymes (catalysts) involved
        in any of the registered reactions
        (regardless of whether they might participate in a non-enzymatic role in other reactions)
        """
        return self.active_enzymes



    #####################################################################################################

    '''                          ~   FOR CREATION OF NETWORK DIAGRAMS  ~                              '''

    def ________NETWORK_DIAGRAMS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def prepare_graph_network(self) -> dict:
        """
        Prepare and return a data structure with chemical-reaction data in a network format,
        ready to be passed to the front end, for network-diagram visualization with the Cytoscape.js library
        (in the graph module "vue_cytoscape")

        EXAMPLE of the graph structure part of the returned object for an  A <-> B reaction:
           [{'id': 'C-0', 'labels': ['Chemical'], 'name': 'A', 'diff_rate': None},
            {'id': 'C-1', 'labels': ['Chemical'], 'name': 'B', 'diff_rate': None},

            {'id': 'R-0', 'labels': ['Reaction'], 'name': 'RXN', 'kF': 3.0, 'kR': 2.0, 'K': 1.5, 'Delta_G': -1005.13},

            {'id': 'edge-1', 'name': 'produces', 'source': 'R-0', 'target': 'C-1', 'stoich': 1, 'rxn_order': 1},
            {'id': 'edge-2', 'name': 'reacts',   'source': 'C-0', 'target': 'R-0', 'stoich': 1, 'rxn_order': 1}
           ]

        :return:    A dictionary with 3 keys: 'structure', 'color_mapping', 'caption_mapping'
        """
        graph = PyGraphVisual()

        # Note: the graph nodes representing Chemicals will be given an id such as "C-123" and a label "Chemical";
        #       the graph nodes representing Reactions will be given an id such as "R-456" and a label "Reaction"

        for i, rxn in enumerate(self.reaction_list):    # Consider each REACTION in turn
            # Add a node representing the reaction
            rxn_id = f"R-{i}"               # Example: "R-456"
            node_data = {'name': 'RXN'}

            rxn_properties = rxn.extract_rxn_properties()
            for k,v in rxn_properties.items():
                node_data[k] = f"{v:,.6g}"

            graph.add_node(node_id=rxn_id, labels='Reaction', data=node_data)


            # Process all the PRODUCTS of this reaction
            products = rxn.extract_products()
            for term in products:
                species_name = rxn.extract_species_name(term)
                chemical_id = f"C-{self.chem_data.get_index(species_name)}"      # Example: "C-12"
                # Add each product to the graph as a node (if not already present)
                graph.add_node( node_id=chemical_id, labels="Chemical",
                                data={'name': species_name,
                                      'diff_rate': self.chem_data.get_diffusion_rate(name=species_name)
                                      })

                # Append edge from "reaction node" to "product node"
                graph.add_edge(from_node=rxn_id, to_node=chemical_id, name="produces",
                               data={'stoich': rxn.extract_stoichiometry(term),
                                     'rxn_order': rxn.extract_rxn_order(term)
                                     })


            # Process all the REACTANTS of this reaction
            reactants = rxn.extract_reactants()
            for term in reactants:
                species_name = rxn.extract_species_name(term)
                chemical_id = f"C-{self.chem_data.get_index(species_name)}"      # Example: "C-34"
                # Add each reactant to the graph as a node (if not already present)
                graph.add_node(node_id=chemical_id, labels="Chemical",
                               data={'name': species_name,
                                     'diff_rate': self.chem_data.get_diffusion_rate(name=species_name)
                                     })

                # Append edge from "reactant node" to "reaction node"
                graph.add_edge(from_node=chemical_id, to_node=rxn_id, name="reacts",
                               data={'stoich': rxn.extract_stoichiometry(term),
                                     'rxn_order': rxn.extract_rxn_order(term)
                                     })


        graph.assign_color_mapping(label='Chemical', color='graph_green')
        graph.assign_color_mapping(label='Reaction', color='graph_lightbrown')

        graph.assign_caption(label='Chemical', caption='name')
        graph.assign_caption(label='Reaction', caption='name')

        #print(graph)

        return graph.serialize()



    def plot_reaction_network(self, graphic_component :str, unpack=False) -> None:
        """
        Send a plot of the network of reactions to the HTML log file,
        also including a brief summary of all the reactions

        EXAMPLE of usage:  plot_reaction_network("vue_cytoscape_2")

        :param graphic_component:   The name of a Vue component that accepts a "graph_data" argument,
                                        an object with the following keys
                                        'structure', 'color_mapping' and 'caption_mapping'
                                        For more details, see ChemData.prepare_graph_network()
        :param unpack:              Use True for Vue components that require their data unpacked into individual arguments;
                                        False for that accept a single data argument, named "graph_data"
        :return:                    None
        """
        # Send a brief summary of all the reactions to the HTML log file
        log.write("List of reactions:", style=log.h3, newline=False, also_print=False)

        rxn_descriptions = self.multiple_reactions_describe(concise=True)
        for desc in rxn_descriptions:
            log.write(desc, indent=4, also_print=False)

        #log.blank_line()

        graph_data = self.prepare_graph_network()
        # A dictionary with 3 keys: 'structure', 'color_mapping' and 'caption_mapping'

        # Send a plot of the network of reactions to the HTML log file
        GraphicLog.export_plot(graph_data, graphic_component, unpack=unpack)
