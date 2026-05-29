import math
import numpy as np


class ThermoDynamics:
    """
    Manage the Thermodynamics aspects of reactions:
    changes in Gibbs Free Energy, Enthalpy, Entropy - and how
    they relate to equilibrium constant, at a given temperature.

    This class does NOT get instantiated.

            "K"       (equilibrium constant - from thermodynamic data)
            "delta_H" (change in Enthalpy: Enthalpy of Products - Enthalpy of Reactants)
            "delta_S" (change in Entropy)
            "delta_G" (change in Gibbs Free Energy)

            No correction is made for the temperature dependency of delta_H

            Note - at constant temperature T :
                Delta_G = Delta_H - T * Delta_S
                Equilibrium constant = exp(-Delta_G / RT)
    """

    # Class Attribute
    R = 8.31446261815324    # Ideal gas constant, in units of Joules / (Kelvin * Mole)



    @classmethod
    def K_from_gibbs(cls, delta_G, temp) -> float:
        """
        Compute a reaction's equilibrium constant from its change in Gibbs Free Energy,
        at the specified temperature

        :param delta_G: Change in Gibbs Free Energy (from reactants to products), in Joules
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's equilibrium constant
        """
        return math.exp(- delta_G / (cls.R * temp))



    @classmethod
    def gibbs_from_K(cls, K, temp) -> float:
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
    def gibbs_from_enthalpy(cls, delta_H, delta_S, temp) -> float:
        """
        Compute the change in Gibbs Free Energy, from Enthalpy and Entropy changes,
        at the specified temperature

        :param delta_H: The reaction's change in Enthalpy (from reactants to products)
        :param delta_S: The reaction's change in Entropy (from reactants to products)
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Free Energy (from reactants to products)
        """
        return delta_H - temp * delta_S


    @classmethod
    def enthalpy_from_gibbs(cls, delta_G, delta_S, temp) -> float:
        """
        Compute the change in Enthalpy, from changes in Gibbs Free Energy and in Entropy,
        at the specified temperature

        :param delta_G: The reaction's change in Gibbs Free Energy (from reactants to products)
        :param delta_S: The reaction's change in Entropy (from reactants to products)
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Enthalpy (from reactants to products)
        """
        return delta_G + temp * delta_S


    @classmethod
    def entropy_from_gibbs(cls, delta_G, delta_H, temp) -> float:
        """
        Compute the change in Entropy, from  changes in the Gibbs Free Energy and in Enthalpy,
        at the specified temperature

        :param delta_G: The reaction's change in Gibbs Free Energy (from reactants to products)
        :param delta_H: The reaction's change in Enthalpy (from reactants to products)
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Entropy (from reactants to products)
        """
        return (delta_H - delta_G) / temp



    @classmethod
    def extract_thermodynamic_data(cls, temp, K=None, delta_H=None, delta_S=None, delta_G=None) -> dict:
        """
        Given the system temperature, and any combination of thermodynamic data,
        verify the consistency of the given data,
        derive whatever is missing and can be derived

        :param temp:                System's temperature, in degree Kelvins
        :param K:       [OPTIONAL]  The reaction's equilibrium constant
        :param delta_H: [OPTIONAL]  The reaction's change in Enthalpy (from reactants to products)
        :param delta_S: [OPTIONAL]  The reaction's change in Entropy (from reactants to products)
        :param delta_G: [OPTIONAL]  The reaction's change in Gibbs Free Energy (from reactants to products)

        :return:        A dict with the 4 keys  "K", "delta_H", "delta_S", "delta_G";
                            any value that wasn't passed, or derivable from the others, will be None
        """
        #print(f"In extract_thermodynamic_data(): temp={temp}, K={K}, delta_H={delta_H}, delta_S={delta_S}, delta_G={delta_G}")

        assert temp is not None, \
            "extract_thermodynamic_data(): a temperature value (in K) must be passed to argument `temp`"


        if (K is not None):
            # If the temperature is set, compute the change in Gibbs Free Energy
            delta_G_derived = cls.gibbs_from_K(K = K, temp = temp)
            if delta_G is None:
                delta_G = delta_G_derived
            else:   # If already present (passed as argument), make sure that the two match!
                assert np.allclose(delta_G_derived, delta_G), \
                    f"extract_thermodynamic_data(): inconsistency between the derived ({delta_G_derived}) " \
                    f"and the passed ({delta_G}) values of Delta_G"


        if (delta_H is not None) and (delta_S is not None):
            # If all the thermodynamic data (possibly except delta_G) is available...

            # Compute the change in Gibbs Free Energy from delta_H and delta_S, at the current temperature
            delta_G_derived = cls.gibbs_from_enthalpy(delta_H = delta_H, delta_S = delta_S, temp = temp)

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
                K = ThermoDynamics.K_from_gibbs(delta_G = delta_G, temp = temp)

            # If either Enthalpy or Entropy is missing, but the other one is known, compute the missing one
            if (delta_H is None) and (delta_S is not None):
                delta_H = ThermoDynamics.enthalpy_from_gibbs(delta_G=delta_G, delta_S=delta_S, temp=temp)
            elif (delta_H is not None) and (delta_S is None):
                delta_S = ThermoDynamics.entropy_from_gibbs(delta_G=delta_G, delta_H=delta_H, temp=temp)


        return {"K": K, "delta_H": delta_H, "delta_S": delta_S, "delta_G": delta_G}
