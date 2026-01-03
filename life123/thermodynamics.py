import math
import numpy as np
from typing import Tuple


class ThermoDynamics:
    """
    Manage the Thermodynamics aspects of reactions:
    changes in Gibbs Free Energy, Enthalpy, Entropy - and how
    they relate to equilibrium constant, at a given temperature.

    This class does NOT get instantiated.

            "K"       (equilibrium constant - from either kinetic or thermodynamic data;
                       if both present, they must match up!)
            "delta_H" (change in Enthalpy: Enthalpy of Products - Enthalpy of Reactants)
            "delta_S" (change in Entropy)
            "delta_G" (change in Gibbs Free Energy)

            Note - at constant temperature T :
                Delta_G = Delta_H - T * Delta_S
                Equilibrium constant = exp(-Delta_G / RT)
    """

    # Class Attribute
    R = 8.31446261815324    # Ideal gas constant, in units of Joules / (Kelvin * Mole)



    @classmethod
    def K_from_delta_G(cls, delta_G, temp) -> float:
        """
        Compute a reaction's equilibrium constant from the thermodynamic data

        :param delta_G: Change in Gibbs Free Energy (from reactants to products), in Joules
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's equilibrium constant
        """
        return math.exp(- delta_G / (cls.R * temp))



    @classmethod
    def delta_G_from_K(cls, K, temp) -> float:
        """
        Compute a reaction's change in its Gibbs Free Energy from its equilibrium constant,
        at the specified temperature

        :param K:       The reaction's equilibrium constant
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Gibbs Free Energy (from reactants to products),
                            in Joules
        """
        return -cls.R * temp * math.log(K)      # Natural log



    @classmethod
    def delta_G_from_enthalpy(cls, delta_H, delta_S, temp) -> float:
        """
        Compute the change in Gibbs Free Energy, from Enthalpy and Entropy changes

        :param delta_H: The reaction's change in Enthalpy (from reactants to products)
        :param delta_S: The reaction's change in Entropy (from reactants to products)
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Free Energy (from reactants to products)
        """
        return delta_H - temp * delta_S


    @classmethod
    def delta_H_from_gibbs(cls, delta_G, delta_S, temp) -> float:
        """
        Compute the change in Enthalpy, from changes in Gibbs Free Energy and in Entropy

        :param delta_G: The reaction's change in Gibbs Free Energy (from reactants to products)
        :param delta_S: The reaction's change in Entropy (from reactants to products)
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Enthalpy (from reactants to products)
        """
        return delta_G + temp * delta_S


    @classmethod
    def delta_S_from_gibbs(cls, delta_G, delta_H, temp) -> float:
        """
        Compute the change in Entropy, from  changes in the Gibbs Free Energy and in Enthalpy

        :param delta_G: The reaction's change in Gibbs Free Energy (from reactants to products)
        :param delta_H: The reaction's change in Enthalpy (from reactants to products)
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Entropy (from reactants to products)
        """
        return (delta_H - delta_G) / temp



    @classmethod
    def compute_reaction_quotient(cls, reactant_data :[(int, str)], product_data :[(int, str)],
                                  conc :dict, explain=False) -> np.double | Tuple[np.double, str]:
        """
        Compute the "Reaction Quotient" Q (aka "Mass–action Ratio"),
        for the reaction with the specified parameters,
        given the concentrations of chemicals involved in the reaction.

        EXAMPLE: use reactant_data=[(2, "A")] and product_data=[(1, "B")] ,
                 for a reaction of the form 2 A <-> B , 
                 alongside a dictionary with the concentrations (activities) of A and B

        Note: in a heterogeneous mixture, solids, pure liquids and solvents have an activity that has a fixed value of 1,
              and should be omitted from the parameters passed to this function.
              We're using the term "concentrations" instead of "chemical activities";
              concentrations approximate the activities of ideal dilute solutions

        :param reactant_data:   List of PAIRS of the form (label of reactant , order of reaction with respect to it);
                                    in elementary reactions, the orders will be equal to their respective stoichiometry coefficients
        :param product_data:    List of PAIRS of the form (label of reaction product , order of reaction with respect to it);
                                    in elementary reactions, the orders will be equal to their respective stoichiometry coefficients
        :param conc:            Dictionary with the concentrations (activities) of the species involved in the reaction.
                                The keys are the chemical labels
                                    EXAMPLE: {'A': 23.9, 'B': 36.1}
        :param explain:         If True, it also returns the math formula being used for the computation
                                    EXAMPLES:   "([C][D]) / ([A][B])"
                                                "[B] /  [A]^2 "

        :return:                If explain is False, return a value for the "Reaction Quotient" (aka "Mass–action Ratio");
                                    if True, return a pair with that quotient and a string with the math formula that was used.
                                    Note that the reaction quotient is a Numpy scalar that might be np.inf or np.nan
        """
        # TODO: maybe also accept strings, in lieu of pairs (for cases where the stoichiometry coefficient is 1)
        # TODO: could be tidier in avoiding unnecessary blanks in the explanations
        numerator = np.double(1)    # The product of all the concentrations of the reaction products (adjusted for reaction order)
        denominator = np.double(1)  # The product of all the concentrations of the reactants (also adjusted for reaction order)

        numerator_text = ""      # First part of the textual explanation
        denominator_text = ""    # Second part of the textual explanation


        # Compute the numerator of the "Reaction Quotient"
        for (stoich_coeff, p) in product_data:
            # Loop over the reaction products
            assert type(stoich_coeff) == int, f"reaction_quotient(): the argument `product_data` " \
                                              f"must be a list of pairs (integer and string).  `{stoich_coeff}` is not an integer"
            assert type(p) == str, f"reaction_quotient(): the argument `product_data` " \
                                   f"must be a list of pairs (integer and string).  {p} is not a string"

            species_name = p
            species_conc = conc.get(p)
            assert species_conc is not None, f"reaction_quotient(): unable to proceed because the " \
                                             f"concentration of product `{species_name}` was not provided"

            numerator *= (species_conc ** stoich_coeff)
            if explain:
                if stoich_coeff > 1:
                    numerator_text += f" [{species_name}]^{stoich_coeff} "
                else:
                    numerator_text += f"[{species_name}]"

        if explain and len(product_data) > 1:
            numerator_text = f"({numerator_text})"  # In case of multiple terms, enclose them in parenthesis


        # Compute the denominator of the "Reaction Quotient"
        for (stoich_coeff, r) in reactant_data:
            # Loop over the reactants
            assert type(stoich_coeff) == int, f"reaction_quotient(): the argument `reactant_data` " \
                                              f"must be a list of pairs (integer and string).  `{stoich_coeff}` is not an integer"
            assert type(r) == str, f"reaction_quotient(): the argument `reactant_data` " \
                                   f"must be a list of pairs (integer and string).  {r} is not a string"
            species_name =  r
            species_conc = conc.get(species_name)
            assert species_conc is not None, f"reaction_quotient(): unable to proceed because the " \
                                             f"concentration of reactant `{species_name}` was not provided"

            denominator *= (species_conc ** stoich_coeff)
            if explain:
                if stoich_coeff > 1:
                    denominator_text += f" [{species_name}]^{stoich_coeff} "
                else:
                    denominator_text += f"[{species_name}]"

        if explain and len(reactant_data) > 1:
            denominator_text = f"({denominator_text})"  # In case of multiple terms, enclose them in parenthesis


        with np.errstate(divide='ignore', invalid='ignore'):
            # It might be np.inf (if just the denominator is zero) or np.nan (if both are zero)
            quotient = numerator / denominator

        if explain:
            formula = f"{numerator_text} / {denominator_text}"
            return (quotient, formula)

        return quotient



    @classmethod
    def extract_thermodynamic_data(cls, temp, K=None, delta_H=None, delta_S=None, delta_G=None) -> dict:
        """

        :param temp:                System's temperature, in degree Kelvins
        :param K:       [OPTIONAL]
        :param delta_H: [OPTIONAL]
        :param delta_S: [OPTIONAL]
        :param delta_G: [OPTIONAL]

        :return:        A dict with the 4 keys  "K", "delta_H", "delta_S", "delta_G"
        """
        #print(f"In extract_thermodynamic_data(): temp={temp}, K={K}, delta_H={delta_H}, delta_S={delta_S}, delta_G={delta_G}")

        assert temp is not None, \
            "extract_thermodynamic_data(): a temperature value (in K) must be passed to argument `temp`"


        if (K is not None):
            # If the temperature is set, compute the change in Gibbs Free Energy
            delta_G_derived = cls.delta_G_from_K(K = K, temp = temp)
            if delta_G is None:
                delta_G = delta_G_derived
            else:   # If already present (passed as argument), make sure that the two match!
                assert np.allclose(delta_G_derived, delta_G), \
                    f"extract_thermodynamic_data(): inconsistency between the derived ({delta_G_derived}) " \
                    f"and the passed ({delta_G}) values of Delta_G"


        if (delta_H is not None) and (delta_S is not None):
            # If all the thermodynamic data (possibly except delta_G) is available...

            # Compute the change in Gibbs Free Energy from delta_H and delta_S, at the current temperature
            delta_G_derived = cls.delta_G_from_enthalpy(delta_H = delta_H, delta_S = delta_S, temp = temp)

            if delta_G is None:
                delta_G = delta_G_derived
            else:  # If already present (passed as argument), make sure that the two match!
                assert np.allclose(delta_G_derived, delta_G), \
                    f"extract_thermodynamic_data(): inconsistency between the value of Delta_G ({delta_G_derived}) " \
                    f"derived from enthalpy/entropy, and the its value ({delta_G}) passed as argument or derived from K"


        if delta_G is not None:
            if K is None:
                # Compute the equilibrium constant (from the thermodynamic data)
                # Note: no need to do it if K is present, because we ALREADY checked for consistency
                K = ThermoDynamics.K_from_delta_G(delta_G = delta_G, temp = temp)

            # If either Enthalpy or Entropy is missing, but the other one is known, compute the missing one
            if (delta_H is None) and (delta_S is not None):
                delta_H = ThermoDynamics.delta_H_from_gibbs(delta_G=delta_G, delta_S=delta_S, temp=temp)
            elif (delta_H is not None) and (delta_S is None):
                delta_S = ThermoDynamics.delta_S_from_gibbs(delta_G=delta_G, delta_H=delta_H, temp=temp)


        return {"K": K, "delta_H": delta_H, "delta_S": delta_S, "delta_G": delta_G}
