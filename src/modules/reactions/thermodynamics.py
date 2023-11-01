import math


class ThermoDynamics:
    """
    Manage the Thermodynamics aspects of reactions:
    changes in Gibbs Free Energy, Enthalpy, Entropy - and how
    they relate to equilibrium constant, at a given temperature.

    This class does NOT get instantiated
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
