from __future__ import annotations      # To facilitate type annotations
import numpy as np
import math
from typing import Set, Mapping
from dataclasses import dataclass, field, asdict
from life123.thermodynamics import ThermoDynamics
from life123.reaction_kinetics import ReactionKinetics
from life123.species_registry import SpeciesRegistry
from life123.units import show_standard_units, convert, K, C




class MichaelisMentenLaw:
    name = "Michaelis-Menten"

    def set_parameters(self, parameters :dict) -> None:
        """
        Validate and set the passed kinetic parameters,
        as well as any others derivable from them

        :param parameters:
        :return:            None
        """
       # Validate that at most only the allowed key were passed
        ALLOWED_KEYS = {"k1_F", "k1_R", "k2_F", "kM", "kcat"}
        """
        :param k1_F:    [OPTIONAL] The forward reaction rate of the 1st part of the reaction
        :param k1_R:    [OPTIONAL] The reverse reaction rate of the 1st part of the reaction
        :param k2_F:    [OPTIONAL] The forward reaction rate of the 2nd part of the reaction
        :param kM:      [OPTIONAL] "Michaelis constant"
        :param kcat:    [OPTIONAL] "Catalytic rate constant" aka "Turnover number" aka "Collective rate constant"
                            (equal to k2_F)
        """
        unexpected_keys = set(parameters.keys()) - ALLOWED_KEYS
        if unexpected_keys:
            raise TypeError(f"set_parameters(): Unexpected parameter keys:  {sorted(unexpected_keys)} ")

        k1_F = parameters.get("k1_F")
        k1_R = parameters.get("k1_R")
        k2_F = parameters.get("k2_F")

        kM = parameters.get("kM")
        kcat = parameters.get("kcat")

        if all(v is not None for v in [k1_F, k1_R, k2_F]):
            kM_derived = (k2_F + k1_R) / k1_F
            if kM is not None:
                assert np.allclose(kM, kM_derived), \
                    f"set_parameters(): inconsistent arguments.  " \
                    f"The passed `kM` value ({kM}) doesn't match the value ({kM_derived}) inferred from the given reaction rate constants"
            else:
                kM = kM_derived

        if k2_F is not None:
            kcat_derived = k2_F
            if kcat is not None:
                assert np.allclose(kcat, kcat_derived), \
                    f"set_parameters(): inconsistent arguments.  " \
                    f"The passed `kcat` value ({kcat}) doesn't match the value ({kcat_derived}) of the given `k2_F` reaction rate constant"
            else:
               kcat = kcat_derived

        self.parameters = {"k1_F": k1_F, "k1_R": k1_R, "k2_F": k2_F, "kM": kM, "kcat": kcat}

    def rate(self, concentrations):
        pass



#########################

class MassActionLaw:
    name = "mass action"

    def __init__(self):
        self.kF: float | None = None
        self.kR: float | None = None
        self.K: float | None = None
        self.reversible: bool = False



    def get_parameters(self):
        return {"kF": self.kF, "kR": self.kR, "K": self.K, "reversible": self.reversible}



    def set_parameters(self, parameters :dict) -> None:
        """
        Validate and set the passed kinetic parameters,
        as well as any others derivable from them

        :param parameters:
        :return:
        """
        # Validate that at most only the allowed key were passed
        ALLOWED_KEYS = {"kR", "kF", "K"}
        unexpected_keys = set(parameters.keys()) - ALLOWED_KEYS

        if unexpected_keys:
            raise TypeError(f"set_parameters(): Unexpected parameter keys:  {sorted(unexpected_keys)} ")



        # --------------------------------------------------------------
        # 1. Merge current state with supplied values
        # --------------------------------------------------------------
        # (TODO: make sure to coordinate with thermodynamic


        values = {
            "kF": self.kF,
            "kR": self.kR,
            "K": self.K,
        }

        values.update(parameters)


        # --------------------------------------------------------------
        # 2. Validate values
        # --------------------------------------------------------------

        for name in ("kF", "kR"):
            value = values[name]

            if value is None:
                continue

            if not isinstance(value, (int, float)):
                raise TypeError(f"set_parameters(): `{name}` must be a number or None; value passed was of type {type(value)}")

            if value < 0:
                raise ValueError(f"set_parameters(): `{name}` must be non-negative; value passed was {value}")


        K = values["K"]

        if K is not None:
            if not isinstance(K, (int, float)) or math.isnan(K):
                raise TypeError(f"set_parameters(): `K` must be a number or None or math.inf; value passed was of type {type(K)}")

            if K < 0:
                raise ValueError(f"set_parameters(): `K` must be non-negative; value passed was {K}")      # # +inf is explicitly allowed


        # --------------------------------------------------------------
        # 3. Resolve what can be resolved
        # --------------------------------------------------------------

        kF = values["kF"]
        kR = values["kR"]


        # ---- kF + kR known: derive K (and check consistency with the given K, if applicable) -------------------

        if (kF is not None) and (kR is not None):

            derived_K = None

            if kF > 0 and kR > 0:
                derived_K = kF / kR

            elif kF > 0 and kR == 0:
                derived_K = math.inf

            if K is None:
                K = derived_K

            elif not (math.isclose(K, derived_K)) or (math.isinf(K) and math.isinf(derived_K)):
                raise ValueError(
                    f"set_parameters(): Inconsistent kinetic parameters: "
                    f"kF={kF}, kR={kR}, K={K}"
                )


        # ---- kF + K known: derive kR (and check consistency with the given kR, if applicable) -------------------

        elif (kF is not None) and (K is not None):

            if K > 0 and math.isfinite(K):
                kR = kF / K

            elif K == 0:
                if kF > 0:
                    raise ValueError(
                        f"set_parameters(): Inconsistent kinetic parameters: "
                        f"kF={kF} cannot coexist with K=0"
                    )

                # kF == 0, K == 0:
                # kR is not uniquely determined.
                kR = None

            else:
                # K == + infinity
                if kF > 0:
                    kR = 0.0
                else:
                    # kF == 0, K == infinity is impossible:
                    # 0 / 0 is undefined.
                    raise ValueError("set_parameters(): K=inf is inconsistent with kF=0")



        # ---- kR + K known: derive kF -------------------------------

        elif (kR is not None) and (K is not None):

            if K > 0 and math.isfinite(K):
                kF = kR * K

            elif K == 0:
                if kR > 0:
                    kF = 0.0
                else:
                    # kR == 0, K == 0:
                    # 0 / 0 is undefined.
                    raise ValueError("set_parameters(): K=0 is inconsistent with kR=0")

            else:
                # K == + infinity
                if kR > 0:
                    raise ValueError(f"set_parameters(): K=inf is inconsistent with kR={kR}")

                # kR == 0, K == +inf:
                # kF can be any positive value, so it cannot
                # be uniquely derived.
                kF = None


        # --------------------------------------------------------------
        # Commit only after everything succeeds
        # --------------------------------------------------------------

        self.kF = kF
        self.kR = kR
        self.K = K

        self.reversible = True if (self.kR is not None and self.kR > 0) else False



    def rate(self, concentrations):
        pass


############################################################################

class Kinetics:
    """

    """

    def __init__(self, law :str, parameters=None):
        """
        :param law:
            EXAMPLES of `law`:
                "mass action"
                "MM"                [Michaelis-Menten]
                "custom"            [user-supplied Python function]
                "Hill"              [not yet supported]
                "enzyme inhibition" [hypothetical future extension; not available]

        :param parameters:  A dict of data that is specific to the given rate law
        """
        # TODO: maybe have a separate class for each rate law

        AVAILABLE_RATE_LAWS = ["mass action", "MM", "custom"]

        self.KINETIC_LAWS = {
            "mass action": MassActionLaw(),
            "MM": MichaelisMentenLaw()
        }

        if law is not None:
            # Validate that it's a known rate law
            assert law in AVAILABLE_RATE_LAWS, \
                f"Kinetics instantiation: the value passed to the `law` argument (\"{law}\") " \
                f"is not one the allowed values: {AVAILABLE_RATE_LAWS}"


        self.law = law
        self.kinetic_rate_function = None
        #self.parameters = {}

        self.set_parameters(parameters)


    def get_parameters(self) -> dict:
        """

        :return:
        """
        law_module = self.KINETIC_LAWS[self.law]
        return law_module.get_parameters()



    def rate(self, concentrations):
        law_module = self.KINETIC_LAWS[self.law]
        return law_module.rate(self.parameters, concentrations)



    def set_parameters(self, parameters :dict) -> None:
        """
        Validate and set the passed kinetic parameters,
        as well as any others derivable from them

        :param parameters:
        :return:            None
        """
        if parameters is None:
            parameters = {}

        law_module = self.KINETIC_LAWS[self.law]
        law_module.set_parameters(parameters)

        return

        if self.law == "MM":
            # Validate that at most only the allowed key were passed
            ALLOWED_KEYS = {"k1_F", "k1_R", "k2_F", "kM", "kcat"}
            """
            :param k1_F:    [OPTIONAL] The forward reaction rate of the 1st part of the reaction
            :param k1_R:    [OPTIONAL] The reverse reaction rate of the 1st part of the reaction
            :param k2_F:    [OPTIONAL] The forward reaction rate of the 2nd part of the reaction
            :param kM:      [OPTIONAL] "Michaelis constant"
            :param kcat:    [OPTIONAL] "Catalytic rate constant" aka "Turnover number" aka "Collective rate constant"
                                (equal to k2_F)
            """
            unexpected_keys = set(parameters.keys()) - ALLOWED_KEYS
            if unexpected_keys:
                raise ValueError(f"set_parameters(): Unexpected parameter keys:  {unexpected_keys} ")

            k1_F = parameters.get("k1_F")
            k1_R = parameters.get("k1_R")
            k2_F = parameters.get("k2_F")

            kM = parameters.get("kM")
            kcat = parameters.get("kcat")

            if all(v is not None for v in [k1_F, k1_R, k2_F]):
                kM_derived = (k2_F + k1_R) / k1_F
                if kM is not None:
                    assert np.allclose(kM, kM_derived), \
                        f"set_parameters(): inconsistent arguments.  " \
                        f"The passed `kM` value ({kM}) doesn't match the value ({kM_derived}) inferred from the given reaction rate constants"
                else:
                    kM = kM_derived

            if k2_F is not None:
                kcat_derived = k2_F
                if kcat is not None:
                    assert np.allclose(kcat, kcat_derived), \
                        f"set_parameters(): inconsistent arguments.  " \
                        f"The passed `kcat` value ({kcat}) doesn't match the value ({kcat_derived}) of the given `k2_F` reaction rate constant"
                else:
                   kcat = kcat_derived

            self.parameters = {"k1_F": k1_F, "k1_R": k1_R, "k2_F": k2_F, "kM": kM, "kcat": kcat}



    def to_dict(self) -> dict:
        """
        Return a dictionary form of the dataclass.
        Unset fields are omitted

        :return:    A dictionary populated with the public fields of this data class
        """
        properties = {"kinetics_type": self.law}

        parameters = self.get_parameters()

        # Only include the fields that were set
        for k,v in parameters.items():
            if v is not None:
                properties[k] = v

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




