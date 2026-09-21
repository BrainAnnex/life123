from __future__ import annotations      # To facilitate type annotations
import math
from life123.reaction_kinetics import ReactionKinetics



class OrderedBiBi_Model:
    pass

class PingPongBiBi_Model:
    pass


################################################################################

class MassAction_Model:
    name = "mass action"
    supports_equilibrium_constant = True


    def __init__(self, stoichiometry):
        self.kF: float | None = None
        self.kR: float | None = None
        self.K: float | None = None
        self.reversible: bool = False   # Model metadata/state

        self.derived_pars : set[str] = set()    # Set of names of parameters that were DERIVED - i.e.
                                                # not directly supplied by the user
        self.stoichiometry = stoichiometry



    def __str__(self):
        return f"MassAction_Model object.  Parameters: {self.get_parameters()}"



    def get_parameters(self) -> dict:
        """
        Return the model's parameters
        :return:
        """
        # TODO: maybe take out "reversible", and implement is_reversible() instead??
        return {"kF": self.kF, "kR": self.kR, "K": self.K, "reversible": self.reversible}



    def set_parameters(self, parameters :dict, derived_pars=None) -> None:
        """
        Validate and set the passed kinetic parameters,
        as well as any others derivable from them.

        Any existing values will be updated as applicable

        :param parameters:
        :param derived_pars:
        :return:                None
        """
        if derived_pars:
            self.derived_pars |= derived_pars   # Set union



        # --------------------------------------------------------------
        # 1. Merge current state with supplied values
        # --------------------------------------------------------------

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
                self.derived_pars.add("K")  # add an element to a set

            elif not (math.isclose(K, derived_K)) or (math.isinf(K) and math.isinf(derived_K)):
                raise ValueError(
                    f"set_parameters(): Inconsistent kinetic parameters: "
                    f"kF={kF}, kR={kR}, K={K}"
                )


        # ---- kF + K known: derive kR (and check consistency with the given kR, if applicable) -------------------

        elif (kF is not None) and (K is not None):

            if K > 0 and math.isfinite(K):
                kR = kF / K
                self.derived_pars.add("kR")  # add an element to a set

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
                self.derived_pars.add("kF")  # add an element to a set

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
        self.derived_pars.add("reversible")  # add an element to a set



    def rate(self, conc_dict :dict) -> float:
        """
        For the specified reaction and its species concentrations,
        determine the instantaneous reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the current system state.

        :param conc_dict:
        :return:
        """
        reactants = self.stoichiometry.get_reactant_ids()
        products  = self.stoichiometry.get_product_ids()
        return ReactionKinetics.compute_rate_elementary(reactants=reactants, products=products,
                                                        kF=self.kF, kR=self.kR, reversible=self.reversible,
                                                        conc_dict=conc_dict)




################################################################################

class MichaelisMenten_Model:
    name = "Michaelis-Menten"
    supports_equilibrium_constant = False


    def __init__(self, stoichiometry):
        self.kM: float|None = None      # "Michaelis constant"
        self.kcat: float|None = None    # "Catalytic rate constant" aka "Turnover number" aka "Collective rate constant"

        self.k1_F: float|None = None
        self.k1_R: float|None = None
        self.k2_F: float|None = None

        self.derived_pars : set[str] = set()    # Set of names of parameters that were DERIVED - i.e.
                                                # not directly supplied by the user
        self.stoichiometry = stoichiometry

        if len(stoichiometry.catalysts) != 1:
            if stoichiometry.catalysts == []:
                raise Exception(f"MichaelisMenten_Model instantiation: "
                                f"Missing enzyme in reaction {stoichiometry.standard_chemical_formula()}")
            else:
                raise Exception(f"MichaelisMenten_Model instantiation: Too many enzymes ({len(stoichiometry.catalysts)}) "
                                f"in reaction {stoichiometry.standard_chemical_formula()}")

        self.E = stoichiometry.catalysts[0]

        reactants = stoichiometry.get_reactant_ids(exclude_catalysts=True)
        assert len(reactants) == 1, \
            f"MichaelisMenten_Model instantiation: Incorrect number of reactants " \
            f"in reaction {stoichiometry.standard_chemical_formula()}"

        products = stoichiometry.get_product_ids(exclude_catalysts=True)
        assert len(products) == 1, \
            f"MichaelisMenten_Model instantiation: Incorrect number of products " \
            f"in reaction {stoichiometry.standard_chemical_formula()}"

        (self.S, ) = reactants  # Unpack
        (self.P, ) = products  # Unpack



    def __str__(self):
        return f"MichaelisMenten_Model object.  Parameters: {self.get_parameters()}"



    def get_parameters(self) -> dict:
        """
        Return the model's parameters
        :return:
        """
        return {"kM": self.kM, "kcat": self.kcat, "Substrate": self.S, "Enzyme": self.E, "Product": self.P}



    def set_parameters(self, parameters :dict, derived_pars=None) -> None:
        """
        Validate and set the passed kinetic parameters,
        as well as any others derivable from them

        TODO: Any existing values will be updated as applicable.  Do it as for "mass action"

        :param parameters:
        :param derived_pars:    [OPTIONAL] Set of names of parameters that were derived
                                    (i.e. not passed by the user)
        :return:                None
        """
        if derived_pars:
            self.derived_pars |= derived_pars   # Set union

        # Validate that at most only the allowed key were passed
        #ALLOWED_KEYS = {"kM", "kcat"}
        """

        """
        #unexpected_keys = set(parameters.keys()) - ALLOWED_KEYS
        #if unexpected_keys:
        #    raise TypeError(f"set_parameters(): Unexpected parameter keys:  {sorted(unexpected_keys)} ")


        #self.parameters = {"k1_F": k1_F, "k1_R": k1_R, "k2_F": k2_F, "kM": kM, "kcat": kcat}
        #self.kM = parameters.get("kM")
        #self.kcat = parameters.get("kcat")

        if parameters is None:
            parameters = {}

        # Validate values
        for name, value in parameters.items():
            if value is None:
                continue

            if not isinstance(value, (int, float)):
                raise TypeError(f"set_parameters(): `{name}` must be a number or None; value passed was of type {type(value)}")

            if value < 0:
                raise ValueError(f"set_parameters(): `{name}` must be non-negative; value passed was {value}")


        # Resolve what can be resolved

        kM = parameters.get("kM")
        kcat = parameters.get("kcat")
        k1_F = parameters.get("k1_F")
        k1_R = parameters.get("k1_R")
        k2_F = parameters.get("k2_F")

        if k2_F is not None:
            derived_kcat = k2_F

            if kcat is None:
                kcat = derived_kcat
                self.derived_pars.add("kcat")
            else:
                if not math.isclose(kcat, derived_kcat):
                    raise ValueError(
                        f"MichaelisMenten_Model.set_parameters(): Inconsistent kinetic parameters: "
                        f"passed kcat={kcat}, derived kcat={derived_kcat}"
                    )


            if (k1_F is not None) and (k1_R is not None):   # still inside the earlier clause (k2_F is not None)
                # Issue advisory
                print("INFO: values for k1_F, k1_R and k2_F were all provided.  "
                      "Consider using the more accurate reaction model 'single substrate mechanism'")

                assert not math.isclose(k1_F, 0), \
                        f"set_parameters(): Cannot use the reaction model 'michaelis menten' when k1_F is zero"

                derived_kM = (k2_F + k1_R) / k1_F

                if kM is None:
                    kM = derived_kM
                    self.derived_pars.add("kM")
                else:
                    if not math.isclose(kM, derived_kM):
                        raise ValueError(
                            f"MichaelisMenten_Model.set_parameters(): Inconsistent kinetic parameters: "
                            f"passed kM={kM}, derived kM={derived_kM}"
                        )


        # Commit only after everything succeeds
        self.kM = kM
        self.kcat = kcat
        self.k1_F = k1_F
        self.k1_R = k1_R
        self.k2_F = k2_F



    def rate(self, conc_dict :dict) -> float:
        """
        For the specified reaction and its species concentrations,
        determine the instantaneous reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the current system state.

        :param conc_dict:
        :return:
        """
        V_max = self.kcat * conc_dict[self.E]

        S_conc = conc_dict[self.S]

        return (V_max * S_conc) / (self.kM + S_conc)





################################################################################


class Custom_Model:
    name = "custom"
    supports_equilibrium_constant = True


    def __init__(self, stoichiometry):
        self.kF: float | None = None
        self.kR: float | None = None
        self.K: float | None = None
        self.reversible: bool = False   # Model metadata/state
        self.rate_function = None

        self.derived_pars : set[str] = set()    # Set of names of parameters that were DERIVED - i.e.
                                                # not directly supplied by the user
        self.stoichiometry = stoichiometry



    def __str__(self):
        return f"Custom_Model object.  Parameters: {self.get_parameters()}"



    def get_parameters(self) -> dict:
        """
        Return the model's parameters
        :return:
        """
        # TODO: maybe take out "reversible", and implement is_reversible() instead??
        return {"kF": self.kF, "kR": self.kR, "K": self.K, "rate_function": self.rate_function,
                "reversible": self.reversible}



    def set_parameters(self, parameters :dict, derived_pars=None) -> None:
        """
        Validate and set the passed kinetic parameters,
        as well as any others derivable from them.

        Any existing values will be updated as applicable

        :param parameters:
        :param derived_pars:
        :return:                None
        """
        if derived_pars:
            self.derived_pars |= derived_pars   # Set union



        # --------------------------------------------------------------
        # 1. Merge current state with supplied values
        # --------------------------------------------------------------

        values = {
            "kF": self.kF,
            "kR": self.kR,
            "K": self.K,
            "rate_function": self.rate_function
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
                self.derived_pars.add("K")  # add an element to a set

            elif not (math.isclose(K, derived_K)) or (math.isinf(K) and math.isinf(derived_K)):
                raise ValueError(
                    f"set_parameters(): Inconsistent kinetic parameters: "
                    f"kF={kF}, kR={kR}, K={K}"
                )


        # ---- kF + K known: derive kR (and check consistency with the given kR, if applicable) -------------------

        elif (kF is not None) and (K is not None):

            if K > 0 and math.isfinite(K):
                kR = kF / K
                self.derived_pars.add("kR")  # add an element to a set

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
                self.derived_pars.add("kF")  # add an element to a set

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
        self.rate_function = values["rate_function"]

        self.reversible = True if (self.kR is not None and self.kR > 0) else False
        self.derived_pars.add("reversible")  # add an element to a set



    def rate(self, conc_dict :dict) -> float:
        """
        For the specified reaction and its species concentrations,
        determine the instantaneous reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the current system state.

        :param conc_dict:
        :return:
        """
        function_to_call = self.rate_function
        assert function_to_call is not None, \
            "Custom_Model.rate(): no kinetic rate function was provide for this reaction.  " \
            "Use set_parameters({'rate_function': YOUR_FUNCTION_NAME}) to specify one"

        #print(f"Custom_Model.rate() - function being invoked to determine the reaction's rate: `{function_to_call.__name__}()`")

        return function_to_call(stoichiometry = self.stoichiometry,
                                kinetic_parameters = {"kF": self.kF, "kR": self.kR},
                                conc_dict = conc_dict)             # Carry out the invocation of the custom function call
