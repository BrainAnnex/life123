import math
import numpy as np
from life123.units import standardize_units


class ThermoDynamics:
    """
    Manage the Thermodynamics aspects of reactions:
    changes in Gibbs Free Energy, Enthalpy, Entropy - and how
    they relate to equilibrium constant, at a given temperature.

    This class does NOT get instantiated.

    UNITS internally used (same as the DEFAULT units for the function calls):
        Temperature:    degree K
        Energy:         kJ              <- NOTE the divergence from SI units, for familiarity of convention
        Molar energy:   kJ/mol
        Molar entropy:  Joules/(mol·K)  <- NOTE : Joules, NOT kJ, for the entropy!

    USING ALTERNATE UNITS:
        We're *beginning* to phase in the units.py module.
        Some - but not yet all - functions now accept arguments such as temp=(25, C)

    ENTITIES:
            "K"       (equilibrium constant - from thermodynamic data)
            "delta_H" (change in Enthalpy: Enthalpy of Products - Enthalpy of Reactants)
            "delta_S" (change in Entropy)
            "delta_G" (change in Gibbs Free Energy)

            No correction is made for the temperature dependency of delta_H

            Note - at constant temperature T :
                Delta_G = Delta_H - T * Delta_S
                Equilibrium constant = exp(-Delta_G / RT)
    """

    # Class Attributes
    R_SI_units = 8.31446261815324   # Ideal gas constant, in units of Joules / (Kelvin * mol)
    R = 0.00831446261815324         #  Ideal gas constant, in units of kJ mol−1 K−1



    @classmethod
    def equilibrium_constant_from_gibbs_energy(cls, delta_G, temp :float) -> float:
        """
        Compute a reaction's equilibrium constant from its change in Gibbs Free Energy,
        at the specified temperature

        :param delta_G: Change in Gibbs Free Energy (from reactants to products), in kJ/mol
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's equilibrium constant
        """
        return math.exp(- delta_G / (cls.R * temp))



    @classmethod
    def gibbs_energy_from_equilibrium_constant(cls, K, temp) -> float:
        """
        Compute a reaction's change in its Gibbs Free Energy from its equilibrium constant,
        at the specified temperature

        :param K:       The reaction's equilibrium constant
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Gibbs Free Energy (from reactants to products),
                            in kJ/mol
        """
        return -cls.R * temp * math.log(K)      # Natural log



    @classmethod
    def gibbs_energy_from_enthalpy_entropy(cls, delta_H, delta_S, temp) -> float:
        """
        Compute the change in Gibbs Free Energy, from Enthalpy and Entropy changes,
        at the specified temperature

        :param delta_H: The reaction's change in Enthalpy (from reactants to products),
                            in kJ/mol           NOTE it's KILO-Joules
        :param delta_S: The reaction's change in Entropy (from reactants to products),
                            in Joules/(mol·K)   NOTE it's Joules, NOT kJ
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Free Energy (from reactants to products),
                            in kJ/mol
        """
        delta_S_internal = delta_S/1000     # Convert to kJ/(mol·K)
        return delta_H - temp * delta_S_internal


    @classmethod
    def enthalpy_from_gibbs_energy(cls, delta_G, delta_S, temp) -> float:
        """
        Compute the change in Enthalpy, from changes in Gibbs Free Energy and in Entropy,
        at the specified temperature

        :param delta_G: The reaction's change in Gibbs Free Energy (from reactants to products),
                            in kJ/mol           NOTE it's KILO-Joules
        :param delta_S: The reaction's change in Entropy (from reactants to products),
                            in Joules/(mol·K)   NOTE it's Joules, NOT kJ
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Enthalpy (from reactants to products),
                            in kJ/mol
        """
        delta_S_internal = delta_S/1000     # Convert to kJ/(mol·K)
        return delta_G + temp * delta_S_internal


    @classmethod
    def entropy_from_gibbs_energy(cls, delta_G, delta_H, temp) -> float:
        """
        Compute the change in Entropy, from  changes in the Gibbs Free Energy and in Enthalpy,
        at the specified temperature

        :param delta_G: The reaction's change in Gibbs Free Energy (from reactants to products),
                            in kJ/mol           NOTE it's KILO-Joules
        :param delta_H: The reaction's change in Enthalpy (from reactants to products),
                            in kJ/mol           NOTE it's KILO-Joules
        :param temp:    System's temperature, in degree Kelvins
        :return:        The reaction's change in Entropy (from reactants to products),
                            in J mol-1 K-1      NOTE it's Joules, NOT kJ
        """
        delta_S_internal = (delta_H - delta_G) / temp
        return delta_S_internal * 1000



    @classmethod
    def extract_thermodynamic_data(cls, temp, K=None, delta_H=None, delta_S=None, delta_G=None) -> dict:
        """
        Given the system temperature, and any combination of thermodynamic data,
        verify the consistency of the given data,
        derive whatever is missing and can be derived

        :param temp:    System's temperature, in degree Kelvins
        :param K:       [OPTIONAL] The reaction's equilibrium constant
        :param delta_H: [OPTIONAL] The reaction's change in Enthalpy (from reactants to products),
                            in kJ/mol           NOTE it's KILO-Joules
        :param delta_S: [OPTIONAL] The reaction's change in Entropy (from reactants to products),
                            in J mol-1 K-1      NOTE it's Joules, NOT kJ
        :param delta_G: [OPTIONAL] The reaction's change in Gibbs Free Energy (from reactants to products),
                            in kJ/mol           NOTE it's KILO-Joules

        :return:        A dict with the 4 keys  "K", "delta_H", "delta_S", "delta_G";
                            any value that wasn't passed, or derivable from the others, will be None
        """
        #print(f"In extract_thermodynamic_data(): temp={temp}, K={K}, delta_H={delta_H}, delta_S={delta_S}, delta_G={delta_G}")

        assert temp is not None, \
            "extract_thermodynamic_data(): a temperature value (in K) must be passed to argument `temp`"


        if (K is not None):
            # If the temperature is set, compute the change in Gibbs Free Energy
            delta_G_derived = cls.gibbs_energy_from_equilibrium_constant(K = K, temp = temp)
            if delta_G is None:
                delta_G = delta_G_derived
            else:   # If already present (passed as argument), make sure that the two match!
                assert np.allclose(delta_G_derived, delta_G), \
                    f"extract_thermodynamic_data(): inconsistency between the derived ({delta_G_derived}) " \
                    f"and the passed ({delta_G}) values of Delta_G"


        if (delta_H is not None) and (delta_S is not None):
            # If all the thermodynamic data (possibly except delta_G) is available...

            # Compute the change in Gibbs Free Energy from delta_H and delta_S, at the current temperature
            delta_G_derived = cls.gibbs_energy_from_enthalpy_entropy(delta_H = delta_H, delta_S = delta_S, temp = temp)

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
                K = ThermoDynamics.equilibrium_constant_from_gibbs_energy(delta_G = delta_G, temp = temp)

            # If either Enthalpy or Entropy is missing, but the other one is known, compute the missing one
            if (delta_H is None) and (delta_S is not None):
                delta_H = ThermoDynamics.enthalpy_from_gibbs_energy(delta_G=delta_G, delta_S=delta_S, temp=temp)
            elif (delta_H is not None) and (delta_S is None):
                delta_S = ThermoDynamics.entropy_from_gibbs_energy(delta_G=delta_G, delta_H=delta_H, temp=temp)


        return {"K": K, "delta_H": delta_H, "delta_S": delta_S, "delta_G": delta_G}



    @classmethod
    def relative_population_states(cls, delta_molar_energy :float|tuple, temp :float|tuple):
        """
        In the presence of two states with respective energies εi and εj (given our knowledge of εj-εi),
        use the Boltzmann distribution to determine the ratio of the populations of those state, Ni/Nj

        EXAMPLE: the second state is 6 kJ/mol lower than the first state;
                 at temp=25 C, the ratio of the populations of the 1st state and the 2nd one is
                    relative_population_states(delta_molar_energy=-6, temp=25)

        :param delta_molar_energy:  Molar energy change from state 1 to state 2 (e.g. E2 - E1)
                                        In kJ/mol, or in pair format such as (5000, J_PER_MOL)
        :param temp:                In degrees Kelvin, or in pair format such as (25, C)
        :return:                    Ratio of the populations of the first state over that of the second one
        """
        delta_molar_energy = standardize_units(delta_molar_energy, dimension="molar energy")
        temp = standardize_units(temp, dimension="temperature")
        return math.exp( delta_molar_energy / (cls.R * temp) )
