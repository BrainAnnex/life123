from __future__ import annotations      # To facilitate type annotations
import numpy as np
import math
from typing import Set, Tuple,  Mapping
from dataclasses import dataclass, field, asdict
from life123.thermodynamics import ThermoDynamics
from life123.reaction_kinetics import ReactionKinetics
from life123.species_registry import SpeciesRegistry
from life123.kinetics import MassAction_Model, MichaelisMenten_Model, Custom_Model
from life123.units import show_standard_units, convert, K, C



@dataclass(frozen=True)
class Stoichiometry:
    """
    Managing the SIGNED stoichiometric coefficients,
    for all species in the reaction.

    Reaction reagents have negative signs, and products have positive signs.

    Catalysts (which would have a net coefficient of zero) are stored separately
    """

    # EXAMPLES below are for the reaction A + E -> B + E

    vector: Mapping[str, int]
    # Signed stoichiometric reaction vector.
    # Maps species id -> signed stoichiometric coefficient.
    # Reactants are negative; products are positive.
    # Catalysts are not represented here, since their net coefficient is zero.
    # EXAMPLE: {"A": -1, "B": 1}

    catalysts : list[str] = field(default_factory=list)
    # EXAMPLE: ["E"]    (list of species that have a net coefficient of zero in the reaction vector)


    def __post_init__(self) -> None:
        """
        Automatically invoked by the constructor, just before it terminates.

        Enforce non zero values in argument `coefficients`

        :return: None
        """
        for sp, coeff in self.vector.items():
            assert coeff != 0, f"Stoichiometry instantiation: coefficients cannot be zero (species \"{sp}\") - " \
            f"if you want to specify an enzyme/catalyst, pass it as a list to the argument `catalysts`"



    def to_dict(self) -> dict:
        """
        Return a dictionary form of the stoichiometry,
        mapping all species id's to their SIGNED stoichiometric coefficients.
        IMPORTANT: catalysts, if any, carry a net signed stoichiometric coefficient of zero

        :return:    A dictionary with the mapping of the species id's
                    to their signed stoichiometric coefficients
        """
        d = dict(self.vector).copy()  # Clone the dictionary
        for cat in self.catalysts:
            d[cat] = 0

        return d



    def get_reaction_vector(self) -> dict:
        """
        Following Martin Feinberg's "Foundations of Chemical Reaction Network Theory",
        we define the "reaction vector" of a reaction y -> y' (where y any y' are vectors)
        as:  y' - y

        IMPORTANT: Note that catalysts, if any, are NOT included.

        The component of (y′ − y) corresponding to species s is  y′_s − y_s,
        i.e the difference between the stoichiometric coefficient of s in the product complex y′ (the right-hand side of the reaction)
        and its stoichiometric coefficient in the reactant complex y (the left-hand side of the equation).
        This difference is the net number of molecules of s
        produced with each occurrence of the reaction y → y′.

        Other terms for reaction vectors: "reaction difference vector", or "reaction increment"
        ~~~
        EXAMPLE: for reaction  A + E -> 2P + Q + E
                 the reactant complex y is:     A + E
                 while product complex y′ is:   2P + Q + E
                 and the corresponding reaction vector y′ - y is:  2P + Q - A
                 The non-zero components of the reaction vector,
                 written as a mapping, are: {"A": -1, "P": 2, "Q": 1}
        ~~~

        :return:    The non-zero components of the reaction vector,
                        written as a dict mapping of species id to its component value
        """
        # Note that catalysts, if any, are NOT included
        return dict(self.vector)



    def get_reactant_list(self) -> list:
        """
        Return all reactants as a list of pairs
        of the form (stoichiometry coefficient, species id).
        Catalysts, if any, are also included.
        ~~~
        EXAMPLE: {"A": -2, "B":- 2, "P": 3, "E": 0}
                 would return to [(2, "A"), (2, "B"), (1, "E")]
        :return:
        """
        # TODO: consider returning a set instead of a list
        return [  (-v, k) for k,v in self.vector.items() if v < 0]  \
                + [(1, c) for c in self.catalysts]



    def get_reactant_ids(self, exclude_catalysts=False) -> Set[str]:
        """
        Return the set of all reactant species id's

        :param exclude_catalysts:
        :return:
        """
        s = { k for k,v in self.vector.items()  if v < 0 }    #  Set construction

        if not exclude_catalysts:
            s |= { c for c in self.catalysts}   # Set union

        return s


    def get_product_list(self) -> list:
        """
        Return all products as a list of pairs
        of the form (stoichiometry coefficient, species id).
        Catalysts, if any, are also included.
        ~~~
        EXAMPLE: {"A": -2, "B":- 2, "P": 3, "E": 0}
                 would return to [(3, "P"), (1, "E")]
        :return:
        """
        # TODO: consider returning a set instead of a list
        return [   (v, k) for k,v in self.vector.items() if v > 0] \
                + [(1, c) for c in self.catalysts]



    def get_product_ids(self, exclude_catalysts=False) -> Set[str]:
        """
        Return the set of all product species id's

        :param exclude_catalysts:
        :return:
        """
        s = { k for k,v in self.vector.items()  if v > 0 }    #  Set construction

        if not exclude_catalysts:
            s |= { c for c in self.catalysts}   # Set union

        return s



    def get_all_species_ids(self) -> Set[str]:
        """
        Return the SET of the id's of ALL the species appearing in this reaction

        :return:    A SET of the id's of the species involved in this reaction
                        Note: being a set, it's NOT in any particular order
        """
        # Use set construction and set union
        return    { k for k,_ in self.vector.items() }  \
                | { c for c in self.catalysts}



    def get_reaction_complexes(self) -> tuple:
        """
        It separates the reactants and products in a pair;
        basically (the left-hand side, the right-hand side).
        Following Martin Feinberg's "Foundations of Chemical Reaction Network Theory",
        a reaction is:  reactant complex -> product complex
        ~~~
        EXAMPLE: for reaction  A + E -> 2P + Q + E
                 the reactant complex is:     A + E
                 while product complex is:   2P + Q + E
                 So, this function will return  ( {"A": 1, "E": 1} ,  {"P": 2, "Q": 1, "E": 1} )
        ~~
        :return:    The pair (reactants complex , products complex)
        """
        d = self.vector
        reactants = {k:-v  for k,v in d.items() if v < 0}
        for cat in self.catalysts:
            reactants[cat] = 1

        products = {k:v    for k,v in d.items() if v > 0}
        for cat in self.catalysts:
            products[cat] = 1

        return (reactants, products)



    def _standard_form_complex(self, complex :dict[str]) -> str:
        """
        Return a user-friendly form of a "complex" (a side of a chemical equation)

        EXAMPLE:  turn {"Fe": 1,  "Cl": 2]  into  "Fe + 2 Cl"

        :param complex:     A dictionary encoding either side of a chemical equation
        :return:            A string with a user-friendly form of a side of a chemical equation
        """
        formula_list = []
        for species_name, stoichiometry in complex.items():

            if stoichiometry == 1:
                term = species_name
            else:
                term = f"{stoichiometry} {species_name}"

            formula_list.append(term)

        return " + ".join(formula_list)



    def standard_chemical_formula(self, reversible=False) -> str:
        """
        Return,  as a string, a user-friendly plain-text form of the reaction

        :return:
        """
        reactants, products = self.get_reaction_complexes()

        arrow = " <-> " if reversible else " -> "

        return self._standard_form_complex(reactants) + arrow + self._standard_form_complex(products)



    def reaction_pattern(self) -> tuple[int, int, int]:
        """
        Return the triplet (n_reactants, n_products, n_catalysts)
        Catalysts, if present, are counted separated, and NOT included under either "reactants" or "products"
        ~~~
        EXAMPLES:
            The reaction A + B -> C      gives (2, 1, 0)
            The reaction A + E -> B + E  gives (1, 1, 1)
        ~~~

        :return:    A triplet of integers (n_reactants, n_products, n_catalysts)
        """
        values = self.vector.values()     # EXAMPLE: dict_values([-1, 1])

        negative_count = sum(v < 0 for v in values)
        positive_count = sum(v > 0 for v in values)
        zero_count     = len(self.catalysts)

        return (negative_count, positive_count, zero_count)



    def consistency_checker(self, conc_before :dict, conc_after :dict) -> None:
        """
        Investigate the change in the concentration of the species involved in the reaction,
        to ascertain whether the change is consistent with the reaction's stoichiometry.

        In other words, a "stoichiometric-ratio consistency check"

        More formally stated, for a single reaction,
        with a reaction vector ν = (ν_A, ν_B, ...)  , where the coefficients are signed numbers,
        and the observed concentration change Delta c = conc_after - conc_before ,
        the question is whether there exists a single scalar "reaction extent" ξ such that
        Δc = ξ ν

        :param conc_before: A dict that maps species `id` to its initial concentration
        :param conc_after:  A dict that maps species `id` to its final concentration
        :return:            None.  Raise an Exception if the change in reactant/product concentrations
                                is not consistent with the reaction's stoichiometry
        """
        assert len(conc_before) == len(self.vector), \
            f"consistency_checker(): the argument `conc_before` must contain exactly the same keys " \
            f"as the species in this reaction: {list(self.vector.keys())}"

        assert len(conc_after) == len(self.vector), \
            f"consistency_checker(): the argument `conc_after` must contain exactly the same keys " \
            f"as the species in this reaction: {list(self.vector.keys())}"

        delta_conc = {}
        for sp,conc in conc_before.items():
            assert sp in self.vector, \
                f"consistency_checker(): the species \"{sp}\" appearing in the argument `conc_before` " \
                f"is not part of the reaction: {self.vector}"

            assert sp in conc_after, \
                f"consistency_checker(): the species \"{sp}\" appearing in the argument `conc_before` " \
                f"is not found among the species provided to the argument `conc_after`: {conc_after}"

            delta_conc[sp] = conc_after[sp] - conc_before[sp]

        #print("delta_conc: ", delta_conc)

        # Example: if the reaction's signed stoichiometric coefficient are
        #          {"A": -1, "B": 1, "E": 0}
        #          then a delta_conc of {"A": -10, "B": 10, "E": 0}
        #          is consistent because of multiple of the reaction vector
        ratio = None
        # Loop over all the values of delta_conc,
        # and make sure they are all the same multiple of the corresponding reaction coefficients
        for sp,conc in delta_conc.items():
            if self.vector[sp] == 0:
                continue
            if ratio is None:
                ratio = delta_conc[sp] / self.vector[sp]
                #print("ratio: ", ratio)
            else:
                assert np.allclose(delta_conc[sp], ratio * self.vector[sp]), \
                    f"consistency_checker(): the delta concentration {delta_conc} " \
                    f"is incompatible with the reaction's stoichiometry of {self.vector}"






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
    K_eq: float | None = None
    temp: float | None = None

    derived_pars : set = field(default_factory=set) # Set of names of parameters that were DERIVED - i.e.
                                                    # not directly supplied by the user

    def __post_init__(self) -> None:
        """
        Automatically invoked by the constructor, just before it terminates

        :return: None
        """
        self.set_temperature(self.temp)



    def to_dict(self) -> dict:
        """
        Return a dictionary form of the dataclass.
        Unset (missing) fields are omitted

        :return:    A dictionary populated with the public fields of this data class
        """
        d = asdict(self)
        result = {}
        for k, v in d.items():
            if k == "derived_pars":
                continue
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
        #TODO: maybe don't pass `temp` as arg, since available as self.temp
        #TODO: maybe rename to "apply_temperature()" or "derive_parameters()"
        # Process the thermodynamic data, and update various object attributes accordingly

        if temp is None:
            return      # Can't do anything!

        thermo_data = ThermoDynamics.extract_thermodynamic_data(K=self.K_eq,
                                                                delta_H=self.delta_H, delta_S=self.delta_S, delta_G=self.delta_G,
                                                                temp=temp)

        #print(f"thermo_data : {thermo_data}")
        if self.K_eq is None and thermo_data["K"] is not None:
            self.derived_pars.add("K_eq")

        if thermo_data["K"] is not None:
            if self.K_eq is None:
                self.derived_pars.add("K_eq")

            self.K_eq = thermo_data["K"]

        if thermo_data["delta_H"] is not None:
            if self.delta_H is None:
                self.derived_pars.add("delta_H")

            self.delta_H = thermo_data["delta_H"]

        if thermo_data["delta_S"] is not None:
            if self.delta_S is None:
                self.derived_pars.add("delta_S")

        self.delta_S = thermo_data["delta_S"]

        if thermo_data["delta_G"] is not None:
            if self.delta_G is None:
                self.derived_pars.add("delta_G")

        self.delta_G = thermo_data["delta_G"]

        self.temp = temp




#############################################################################################

class SimulationReaction:
    """
    The reaction object actually handed to the simulation engine.
    Simulation level.
    The actual chemical reaction, as seen by the simulator.

    It can handle a variety of kinetic laws.

    Note: at a future date, the simulation engine might have things that aren't strictly "kinetic reactions";
    for example, transport events, membrane events, diffusion operators, binding events, etc.
    """
    def __init__(self, model, stoichiometry, source_object :ReactionDefinition, derivation=None,
                 analytic_solution_family=None):
        """

        :param model:           Object of type such as "MassAction_Model" or "MichaelisMenten_Model"
        :param stoichiometry:   Object of type "Stoichiometry"
        :param source_object:   Provenance info: the "ReactionDefinition" source object
        :param derivation:      [OPTIONAL] To explain how this object came about:
                                    either "direct" (essentially just an expansion of the original reaction definition)
                                    or "generated" (one of multiple sub-reactions used to model the original reaction)
        :param analytic_solution_family:   [OPTIONAL]
        """
        self.model = model
        self.stoichiometry: Stoichiometry|None = stoichiometry

        self.source_object = source_object          # Provenance object
        self.derivation : str | None = derivation

        self.analytic_solution_family : str|None = analytic_solution_family

        # TODO: maybe add another variable "role", such as "binding", "catalysis", "ES formation"
        #       "ES breakdown" ("what role this particular generated reaction plays")



    def get_parameters(self) -> dict:
        """
        Get the kinetic parameters

        :return:
        """
        return self.model.get_parameters()



    def get_thermodynamics(self) -> ReactionThermodynamics:
        """
        Get the "ReactionThermodynamics" object
        associated to the source reaction

        :return:
        """
        return self.source_object.thermodynamics



    def set_parameters(self, parameters :dict, derived_pars=None) -> None:
        """
        Set the kinetic parameters (won't affect other existing values not passed here)

        :param parameters:
        :param derived_pars:
        :return:
        """
        self.model.set_parameters(parameters=parameters, derived_pars=derived_pars)



    def describe(self, concise=False) -> str:
        """
        This is the "SimulationReaction" version of ReactionDefinition.describe()

        :param concise:
        :return:
        """
        reversible = getattr(self.model, 'reversible', False)    # Note: the "reversible" attribute may or may not be present

        rxn_description = self.stoichiometry.standard_chemical_formula(reversible=reversible)

        if not concise:
            rxn_description += f'Type: "{self.model.name}" {self.source_object.format_reaction_details(self.model.get_parameters())}'

        return rxn_description



    def standard_chemical_formula(self):
        """

        :return:
        """
        reversible = getattr(self.model, 'reversible', False)    # Note: the "reversible" attribute may or may not be present
        return self.stoichiometry.standard_chemical_formula(reversible=reversible)



    def extract_forward_rate_constant(self) -> float | None:
        """

        :return:    The value of the forward rate constant for this reaction,
                        IF it exists for this reaction type, and is set
        """
        return getattr(self.model, "kF", None)  # Note: the "kF" attribute may or may not be present,
                                                #       depending on the reaction type


    def extract_reverse_rate_constant(self) -> float:
        """

        :return:    The value of the reverse (back) rate constant for this reaction,
                        IF it exists for this reaction type, and is set
        """
        return getattr(self.model, "kR", None)  # Note: the "kR" attribute may or may not be present,
                                                #       depending on the reaction type



    def reaction_quotient(self, conc, explain=False) -> np.double | tuple[np.double, str]:
        """
        Compute the "Reaction Quotient" (aka "Mass–action Ratio"),
        given the concentrations of chemicals involved in this reaction.

        Note: this implementation only covers reactions that have "mass action" kinetics

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
        assert self.model.name == "mass action", \
            "reaction_quotient(): only 'mass action' reaction models are currently supported"

        return ReactionKinetics.compute_reaction_quotient(reactant_data=self.stoichiometry.get_reactant_list(),
                                                          product_data=self.stoichiometry.get_product_list(),
                                                          conc=conc, explain=explain)



    def determine_reaction_rate(self, conc_dict :dict) -> float:
        """
        For the specified concentrations of the species in the reaction,
        determine its initial reaction's "rate" (aka "velocity"),
        i.e. its "forward rate" minus its "reverse rate",
        at the current system state.

        :param conc_dict:   A dict mapping species id's to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :return:            The differences between the reaction's forward and reverse rates
        """
        return self.model.rate(conc_dict=conc_dict)



    def step_simulation(self, delta_time, conc_dict :dict, exact=False) -> Tuple[dict, float]:
        """
        Simulate the generic reaction, over the specified time interval.
        If exact=False, the forward Euler method is used.

        :param delta_time:  The time duration of this individual reaction step - assumed to be small enough that the
                                concentrations won't vary significantly during this span
        :param conc_dict:   A dict mapping chemical labels to their concentrations,
                                for all the chemicals involved in this reaction
                                EXAMPLE:  {"B": 1.5, "F": 31.6, "D": 19.9}
        :param exact:       Only available if this reaction type has a known analytical solution

        :return:            The pair (increment_dict_single_rxn, rxn_rate)
                                - increment_dict_single_rxn     The mapping of chemical labels
                                                                    to their concentration CHANGES
                                                                    during this step
                                - rxn_rate                      The reaction rate ("velocity") for this reaction
                                EXAMPLE of increment_dict_single_rxn: {"B": -1.3, "F": 2.9, "D": -1.6}
        """
        increment_dict_single_rxn = {}      # The keys are the species id's,
                                            # and the values are their respective concentration changes as a result of this reaction

        # Compute the reaction rate ("velocity"), at the current system chemical concentrations, for this reaction
        rxn_rate = self.determine_reaction_rate(conc_dict=conc_dict)

        reactants = self.stoichiometry.get_reactant_list()     # A list of pairs of the form (stoichiometry coefficient, species id))
        products = self.stoichiometry.get_product_list()       # A list of pairs of the form (stoichiometry coefficient, species id))


        if exact:
            if self.analytic_solution_family == "ONE_TO_ONE":
                r = reactants[0][1]           # EXAMPLE: "R"
                p = products[0][1]            # EXAMPLE: "P"

                R0 = conc_dict[r]
                P0 = conc_dict[p]
                # Compute the respective increments of R0 and P0
                if self.model.reversible:
                    delta_p = ReactionKinetics.exact_advance_unimolecular_reversible(kF=self.model.kF, kR=self.model.kR,
                                                                                     A0=R0, P0=P0, t=delta_time, incremental=True)
                else:
                    delta_p = ReactionKinetics.exact_advance_unimolecular_irreversible(kF=self.model.kF,
                                                                                       A0=R0, P0=P0, t=delta_time, incremental=True)

                # Work out the stoichiometry for all the species
                increment_dict_single_rxn = {r: -delta_p, p: delta_p}
                return (increment_dict_single_rxn, rxn_rate)
            else:
                raise Exception("step_simulation(): no exact analytical solution is available for this reaction type")



        # If we get thus far, exact=False

        # In the "forward Euler" approximation, the following rate is taken to remain unvaried during the entire (small) time step
        delta_rxn = rxn_rate * delta_time      # forward reaction - reverse reaction


        """
        Determine the concentration adjustments as a result of this reaction step, 
        for this individual reaction being considered
        """

        # The reactants DECREASE based on the quantity delta_rxn
        for stoichiometry, species_id in reactants:         # Unpack data from each reactant
            delta_conc = stoichiometry * (- delta_rxn)      # Increment to this reactant from the reaction being considered

            increment_dict_single_rxn[species_id] = increment_dict_single_rxn.get(species_id,0) + delta_conc


        # The reaction products INCREASE based on the quantity delta_rxn
        for stoichiometry, species_id in products:      # Unpack data from each product
            delta_conc = stoichiometry * delta_rxn      # Increment to this reaction product from the reaction being considered

            increment_dict_single_rxn[species_id] = increment_dict_single_rxn.get(species_id,0) + delta_conc


        assert len(increment_dict_single_rxn) == len(self.stoichiometry.get_all_species_ids())  # TODO: temporary check to eventually drop

        return (increment_dict_single_rxn, rxn_rate)






#############################################################################################

class Reconciler:
    # Experimental

    @staticmethod
    def reconcile(thermodynamics_data, model):
        """
        If available and feasible, propagate the thermodynamic value K_eq to the kinetic model

        :param thermodynamics_data:
        :param model:
        :return:
        """
        #TODO: maybe just pass `K_eq` as arg, in lieu of `thermodynamics_data`

        if thermodynamics_data.K_eq is None:
            return      # Nothing to do

        if model.supports_equilibrium_constant:
            if model.K is not None:
                # Check for inconsistency
                assert math.isclose(model.K, thermodynamics_data.K_eq), \
                       "conflict between thermodynamic value K_eq and kinetic value K"
            else:
                # Enrich the kinetic data, from thermodynamic values
                model.set_parameters(parameters={"K":thermodynamics_data.K_eq}, derived_pars={"K"})




#############################################################################################


class ReactionCompiler_MassAction:
    @staticmethod
    def compile(stoichiometry, kinetic_parameters, thermodynamics_data, source_object,
               species_registry=None, analytic_solution_family=None) -> tuple:
        """

        :param stoichiometry:
        :param kinetic_parameters:
        :param thermodynamics_data:
        :param source_object:
        :param species_registry:
        :return:
        """
        #print("In compile() method of class 'ReactionCompiler_MassAction'")
        #print("analytic_solution_family: ", analytic_solution_family)

        ALLOWED_KEYS = {"kR", "kF", "K"}
        unexpected_keys = set(kinetic_parameters.keys()) - ALLOWED_KEYS

        if unexpected_keys:
            raise TypeError(f"ReactionCompiler_MassAction.compile(): Unexpected parameter keys:  {sorted(unexpected_keys)} ")

        m = MassAction_Model(stoichiometry)

        pars = kinetic_parameters.copy()    # Clone the dictionary
        derived_pars = set()

        # If available and feasible, propagate the thermodynamic value K_eq to the kinetic model
        if (thermodynamics_data.K_eq is not None) and (m.supports_equilibrium_constant):
            if pars.get("K") is not None:
                # Check for inconsistency
                assert np.allclose(pars["K"], thermodynamics_data.K_eq), \
                       "ReactionCompiler_MassAction(): conflict between thermodynamic value K_eq and kinetic value K"
            else:
                # Enrich the kinetic data, from thermodynamic values
                pars["K"] =  thermodynamics_data.K_eq
                derived_pars = {"K"}


        m.set_parameters(parameters=pars, derived_pars=derived_pars)     # Pass thru the parameters (possibly enhanced by additional parameters derived from the thermodynamics)
        #print(f"    name of model being used: {m.name!r}")

        sim_rxn = SimulationReaction(model=m, stoichiometry=stoichiometry,
                                     source_object=source_object,
                                     derivation="direct",
                                     analytic_solution_family=analytic_solution_family)
                                     # "Directly modeled as specified by the user"
                                     # "this simulation reaction corresponds directly to the reaction defined by the user"

        # TODO: try out:
        #Reconciler.reconcile(thermodynamics_data=thermodynamics_data, model=m)

        # If available and feasible, propagate the kinetic value K to the thermodynamic K_eq
        if m.supports_equilibrium_constant:
            if (m.K is not None) and (m.K != math.inf):     # If the kinetic value K is available
                #print(m.K)
                if thermodynamics_data.K_eq:
                    assert np.allclose(m.K, thermodynamics_data.K_eq), \
                       "ReactionCompiler_MassAction(): conflict between the kinetic value K and the thermodynamic value K_eq"
                else:
                    thermodynamics_data.K_eq = m.K      # Thermodynamic K_eq can be set from kinetic K
                    thermodynamics_data.derived_pars.add("K_eq")
                    thermodynamics_data.set_temperature(thermodynamics_data.temp)
        # TODO: END of part to replace

        return (sim_rxn,)



class ReactionCompiler_MichaelisMenten:
    """
    Coarse-grained model based on the Michaelis-Menten mechanism:
        E + S <-> ES* -> E + P

    -- reduced model
    -- kM + kcat
    -- ES* intermediate eliminated (i.e. not represented)

    More details under "class MichaelisMenten_Model"

    SEE ALSO:  class ReactionCompiler_SingleSubstrateMechanism
    """
    @staticmethod
    def compile(stoichiometry, kinetic_parameters, thermodynamics_data, source_object,
               species_registry=None, analytic_solution_family=None):
        """

        :param stoichiometry:
        :param kinetic_parameters:
        :param thermodynamics_data:
        :param source_object:
        :param species_registry:
        :param analytic_solution_family:
        :return:
        """
        #print("In compile() method of class 'ReactionCompiler_MichaelisMenten'")

        ALLOWED_KEYS = {"kM", "kcat", "k1_F", "k1_R", "k2_F"}
        unexpected_keys = set(kinetic_parameters.keys()) - ALLOWED_KEYS
        if unexpected_keys:
            raise TypeError(f"set_parameters(): Unexpected parameter keys:  {sorted(unexpected_keys)} ")

        if (kinetic_parameters.get("k1_F") is not None) \
            and (kinetic_parameters.get("k1_R") is not None) \
            and (kinetic_parameters.get("k2_F") is not None):
            # Issue advisory
            advisory = """
            The supplied parameters (k1_F, k1_R, k2_F) define an explicit single-substrate mechanism. 
            The requested "michaelis menten" model uses only the derived kM and kcat. 
            Consider using the model "single substrate mechanism" if transient ES dynamics or enzyme sequestration are important.
            """
            print("INFO: ", advisory)

        m = MichaelisMenten_Model(stoichiometry)

        m.set_parameters(parameters=kinetic_parameters)     # Pass thru the parameters
        #print(f"    name of model being used: {m.name!r}")

        r1 = SimulationReaction(model=m, stoichiometry=stoichiometry, source_object=source_object)
        return (r1,)



class ReactionCompiler_SingleSubstrateMechanism:
    """
    Finer-grained model based on the Michaelis-Menten mechanism:
        E + S <-> ES* -> E + P

    It uses an explicit ES* intermediate.

        -- explicit mechanism
        -- k1_F + k1_R + k2_F
        -- ES* represented
        -- compiles to 2 mass-action reactions

    SEE ALSO:  class ReactionCompiler_MichaelisMenten
    """
    @staticmethod
    def compile(stoichiometry, kinetic_parameters, thermodynamics_data,  source_object,
               species_registry, analytic_solution_family=None):
        """

        :param stoichiometry:
        :param kinetic_parameters:
        :param thermodynamics_data:
        :param source_object:
        :param species_registry:
        :param analytic_solution_family:
        :return:
        """
        #print("In compile() method of class 'ReactionCompiler_SingleSubstrateMechanism'")

        ALLOWED_KEYS = {"k1_F", "k1_R", "k2_F", "kM", "kcat"}
        unexpected_keys = set(kinetic_parameters.keys()) - ALLOWED_KEYS

        if unexpected_keys:
            raise TypeError(f"ReactionCompiler_SingleSubstrateMechanism.compile(): Unexpected parameter keys:  {sorted(unexpected_keys)} ")

        assert stoichiometry.reaction_pattern() == (1, 1, 1), \
            "ReactionCompiler_SingleSubstrateMechanism.compile(): reaction stoichiometry " \
            "isn't compatible with the requested model"

        reactants, products = stoichiometry.get_reaction_complexes()

        S, coeff = next(iter(reactants.items()))    # Unpack the single-element dictionary
        assert coeff == 1

        P, coeff = next(iter(products.items()))     # Unpack the single-element dictionary
        assert coeff == 1

        E = stoichiometry.catalysts[0]

        ES = E + S + "*"
        species_registry.add_species(id= ES,
                                     annotation="reaction intermediary from single substrate mechanism")

        # Reaction 1: S + E <-> ES
        st = Stoichiometry(vector={S: -1, E: -1, ES: 1})
        m1 = MassAction_Model(stoichiometry=st)
        m1.set_parameters(parameters={"kF": kinetic_parameters.get("k1_F"),
                                      "kR": kinetic_parameters.get("k1_R")})
        r1 = SimulationReaction(model=m1, stoichiometry=st,
                                source_object=source_object)

        # Reaction 2: ES -> P + E
        st = Stoichiometry(vector={ES: -1, P: 1, E: 1})
        m2 = MassAction_Model(stoichiometry=st)
        m2.set_parameters(parameters={"kF": kinetic_parameters.get("k2_F")})
        r2 = SimulationReaction(model=m2, stoichiometry=st,
                                source_object=source_object)

        return (r1, r2)




class ReactionCompiler_Custom:
    @staticmethod
    def compile(stoichiometry, kinetic_parameters, thermodynamics_data, source_object,
               species_registry=None, analytic_solution_family=None) -> tuple:
        """

        :param stoichiometry:
        :param kinetic_parameters:
        :param thermodynamics_data:
        :param source_object:
        :param species_registry:
        :return:
        """
        #print("In compile() method of class 'ReactionCompiler_Custom'")
        #print("analytic_solution_family: ", analytic_solution_family)

        ALLOWED_KEYS = {"kR", "kF", "K", "rate_function"}
        unexpected_keys = set(kinetic_parameters.keys()) - ALLOWED_KEYS

        if unexpected_keys:
            raise TypeError(f"ReactionCompiler_Custom.compile(): Unexpected parameter keys:  {sorted(unexpected_keys)} ")

        m = Custom_Model(stoichiometry)

        pars = kinetic_parameters.copy()    # Clone the dictionary
        derived_pars = set()

        # If available and feasible, propagate the thermodynamic value K_eq to the kinetic model
        if (thermodynamics_data.K_eq is not None) and (m.supports_equilibrium_constant):
            if pars.get("K") is not None:
                # Check for inconsistency
                assert np.allclose(pars["K"], thermodynamics_data.K_eq), \
                       "ReactionCompiler_MassAction(): conflict between thermodynamic value K_eq and kinetic value K"
            else:
                # Enrich the kinetic data, from thermodynamic values
                pars["K"] =  thermodynamics_data.K_eq
                derived_pars = {"K"}


        m.set_parameters(parameters=pars, derived_pars=derived_pars)     # Pass thru the parameters (possibly enhanced by additional parameters derived from the thermodynamics)
        #print(f"    name of model being used: {m.name!r}")

        sim_rxn = SimulationReaction(model=m, stoichiometry=stoichiometry,
                                     source_object=source_object,
                                     derivation="direct",
                                     analytic_solution_family=analytic_solution_family)
                                     # "Directly modeled as specified by the user"
                                     # "this simulation reaction corresponds directly to the reaction defined by the user"

        # TODO: try out:
        #Reconciler.reconcile(thermodynamics_data=thermodynamics_data, model=m)

        # If available and feasible, propagate the kinetic value K to the thermodynamic K_eq
        if m.supports_equilibrium_constant:
            if (m.K is not None) and (m.K != math.inf):     # If the kinetic value K is available
                #print(m.K)
                if thermodynamics_data.K_eq:
                    assert np.allclose(m.K, thermodynamics_data.K_eq), \
                       "ReactionCompiler_MassAction(): conflict between the kinetic value K and the thermodynamic value K_eq"
                else:
                    thermodynamics_data.K_eq = m.K      # Thermodynamic K_eq can be set from kinetic K
                    thermodynamics_data.derived_pars.add("K_eq")
                    thermodynamics_data.set_temperature(thermodynamics_data.temp)
        # TODO: END of part to replace

        return (sim_rxn,)





###########################################################################

class ReactionModelRegistry:
    """
    Registry for all the Reaction Models.

    Invoked by ReactionDefinition,
    to pick and utilize the appropriate reaction "compiler" (model)
    """
    REGISTERED_REACTION_MODELS = \
        {
            "mass action": ReactionCompiler_MassAction,
            "michaelis menten": ReactionCompiler_MichaelisMenten,
            "single substrate mechanism": ReactionCompiler_SingleSubstrateMechanism,
            "custom": ReactionCompiler_Custom
        }

        # "single substrate mechanism" means  S + E <-> SE -> P + E  , with 3 parameters
        # Note: "single substrate mechanism" isn't merely a kinetic law;
        #        it is a model specification that happens to compile into kinetic laws.
        #        That's why we're using the term "reaction_model" rather than "kinetic_law"


    @classmethod
    def get_compiler_class(cls, model_name :str):
        """
        EXAMPLE:  get_compiler_class("mass action")

        :param model_name:
        :return:            A python class
        """
        cl = cls.REGISTERED_REACTION_MODELS.get(model_name)
        assert cl is not None, \
            f"get_compiler_class(): no python handler class register for the model name \"{model_name}\""

        return cl





###################################################################################################################


class ReactionDefinition:
    """
    The user-facing object; what the user specifies as an overall reaction (simple or complex).
    This is the authoritative biological object at the user/model level.

    The simulation engine never sees this.

    A reaction definition is what the modeler specifies;
    a simulation reaction is what the numerical engine executes.

    A ReactionDefinition object may expand into one or more SimulationReaction objects.

    In other words, some reaction definitions are compound models and expand into multiple simulation reactions.

            user-level reaction model --->  simulation-level reaction model

    EXAMPLES - specific:
        1) An ordinary mass-action reaction can simply compile as:

            ReactionDefinition
                    ↓
                   [itself / equivalent]
                    ↓
            SimulationReaction

        2) A mechanistic enzyme reaction (such as S + E <-> SE -> P + E , when we are given kF_1, kR_1 and kF_2:

            ReactionDefinition
                    ↓
              expansion
                    ↓
            SimulationReaction A
            SimulationReaction B


    EXAMPLES - tabulation:

        | User enters                                | ReactionDefinition | SimulationRepresentation |
        | -------------------------------------------| ------------------ | ------------------------ |
        | A + B -> C,   mass action                  | one                | one                      |
        | S -> P, MM mechanism, with kcat, kM        | one                | one                      |
        | E+S ⇌ ES -> E+P , with kF_1, kR_1 and kF_2 | one                | two                      |
        | ordered Bi-Bi                              | one                | several                  |
        | ping-pong Bi-Bi (2 substrates, 2 products) | one                | several                  |

    """

    def __init__(self, reactants :str|list, products :str|list,
                 species_registry :SpeciesRegistry,
                 reaction_model=None,
                 autoregister_species=False,
                 name=None, id=0,
                 thermodynamic_parameters=None,
                 kinetic_parameters=None):
        """

        :param reactants:   A list/tuple of terms that are either species id's (with implied stoichiometry 1),
                                or pairs (stoichiometry coefficient , species id).
                                If not a list, it will first get turned into one
        :param products:    A list/tuple of terms that are either chemicals labels (with implied stoichiometry 1),
                                or pairs (stoichiometry coefficient , chemical label).
                                If not a list, it will first get turned into one
        :param species_registry:
        :param auto_register_species:
        :param reaction_model:[OPTIONAL] Primarily meant for the kinetics.
                                Allowed values are "mass action", "michaelis menten",
                                "single substrate mechanism", "custom" - as detailed
                                in class ReactionModelRegistry

        :param name:        [OPTIONAL]
        :param id:

        :param delta_H:     [OPTIONAL] Change in Enthalpy (from reactants to products), in kJ/mol
        :param delta_S:     [OPTIONAL] Change in Entropy (from reactants to products), in Joules/(mol·K)
        :param temp:        [OPTIONAL]


        """
        if thermodynamic_parameters is None:
            thermodynamic_parameters = {"delta_H": None, "delta_S": None, "delta_G": None, "K_eq": None, "temp": None}

        delta_H=thermodynamic_parameters.get("delta_H")
        delta_S=thermodynamic_parameters.get("delta_S")
        delta_G=thermodynamic_parameters.get("delta_G")
        K_eq=thermodynamic_parameters.get("K_eq")
        temp=thermodynamic_parameters.get("temp")


        self.name = name
        self.id = id

        self.thermodynamics: ReactionThermodynamics | None = None

        self.source_kinetic_parameters : dict|None = kinetic_parameters if kinetic_parameters is not None else {}
        self.source_thermodynamic_parameters = thermodynamic_parameters
        #self.kinetics: Kinetics | None = None

        self.stoichiometry = None   # A "Stoichiometry" object
                                    #   managing all the stoichiometric coefficients
                                    #   (incl. for catalysts, if applicable)
                                    #   for all species in the reaction

        self.analytic_solution_family = None    # Available values: "ONE_TO_ONE", "ONE_TO_TWO", "TWO_TO_ONE"
        self.reaction_category = None

        self.species_registry = species_registry

        self.reaction_model = reaction_model
        self.sim_reactions :tuple|None = ()     # Tuple of "ReactionSimulation" objects

        self.annotations :str|None = None       # Not in current use


        self._parse(reactants=reactants, products=products, autoregister_species=autoregister_species)

        # if self._detect_elementary_reaction(reaction_model):
        #    reaction_model = "mass action"


        # Process the given thermodynamic data
        self.thermodynamics = ReactionThermodynamics(delta_H=delta_H, delta_S=delta_S, delta_G=delta_G,
                                                     K_eq=K_eq, temp=temp)


        #self.kinetics = Kinetics(law=reaction_model, parameters=kinetic_parameters)

        #self.reaction_type = self._determine_reaction_type()
        #print(f"detected reaction type `{self.reaction_type}`")

        self.reaction_category = self._determine_reaction_category()
        #print(f"detected reaction category `{self.reaction_category}`")

        self.analytic_solution_family = self._determine_analytic_solution_family()

        if reaction_model is not None:
            self._build_model()



    def _parse(self, reactants :list, products :list, autoregister_species :bool):
        """
        Parse the reactants and products,
        and set the object variable self.stoichiometry accordingly.
        Possibly modify self.species_registry as needed

        :param reactants:           A list of pairs (stoichiometry, species id)
        :param products:            A list of pairs (stoichiometry, species id)
        :param autoregister_species:
        :return:                    None
        """
        #TODO: unit test
        assert reactants is not None, \
            "ReactionDefinition() instantiation: the argument `reactants` is a required one"
        if type(reactants) == str:
            reactants = [reactants]
        else:
            assert type(reactants) is list, \
                "ReactionDefinition() instantiation: the argument `reactants` must be a list or a string"

        assert products is not None, \
            "ReactionDefinition() instantiation: the argument `products` is a required one"
        if type(products) == str:
            products = [products]
        else:
            assert type(products) is list, \
                "ReactionDefinition() instantiation: the argument `products` must be a list or a string"


        # Normalize the elements of each list to be (int, str) pairs; i.e. turn any single string "X" into the pair (1, "X")
        reactant_list = [(1, r) if type(r) == str else r
                            for r in reactants]   # A list of pairs
        product_list =  [(1, p) if type(p) == str else p
                            for p in products]   # A list of pairs

        # Catch identical reaction sides, even if terms are reshuffled
        assert set(reactant_list) != set(product_list), \
            f"ReactionDefinition(): the two sides of the reaction can't be identical! " \
            f"Same reactant and product complexes: \"{self._standard_form_chem_eqn(reactant_list)}\""


        # Check whether all the species in the reaction are registered ones
        for _, s_id in reactant_list:
            if not self.species_registry.species_exists(s_id):
                if autoregister_species:
                    self.species_registry.add_species(id=s_id)
                else:
                    raise Exception(f'No species with id "{s_id}" exists in the species registry')

        for _, s_id in product_list:
            if not self.species_registry.species_exists(s_id):
                if autoregister_species:
                    self.species_registry.add_species(id=s_id)
                else:
                    raise Exception(f'No species with id "{s_id}" exists in the species registry')


        c = self.get_signed_stoichiometric_coefficients(reactants=reactant_list, products=product_list)
        self.stoichiometry = Stoichiometry(vector= {k: v for k,v in c.items() if v != 0},
                                           catalysts =    [k  for k,v in c.items() if v == 0])



    def _build_model(self) -> None:
        """
        Process the kinetic data

        :return:                    None
        """
        # Look up the appropriate "reaction compiler"
        #print(f"reaction_model: {self.reaction_model!r}")    # EXAMPLE: 'mass action'
        reaction_compiler = ReactionModelRegistry.get_compiler_class(model_name=self.reaction_model)
        # EXAMPLE: the `ReactionCompiler_MassAction` class
        #print("reaction_compiler: ", reaction_compiler.__name__)    # EXAMPLE: ReactionCompiler_MassAction

        # Invoke the appropriate member of the "reaction compiler" family of classes;
        # a tuple of "SimulationReaction" objects is returned
        #print("self.analytic_solution_family: ", self.analytic_solution_family)
        sr_tuple = reaction_compiler.compile(stoichiometry=self.stoichiometry,
                                             kinetic_parameters=self.source_kinetic_parameters,
                                             thermodynamics_data=self.thermodynamics,
                                             source_object=self,
                                             species_registry=self.species_registry,
                                             analytic_solution_family=self.analytic_solution_family)

        self.sim_reactions = sr_tuple
        #print("self.sim_reactions: ", self.sim_reactions)



    def _detect_elementary_reaction_NOT_IN_USE(self, kinetics_type) -> bool:
        """
        TODO: maybe turn into a generator for default model type
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

        # TODO: switch to using signed terms
        reactants = self.stoichiometry.get_reactant_list()
        products = self.stoichiometry.get_product_list()

        if (r == 1 and reactants[0][0] == 1) and (p == 1 and products[0][0] == 1):
            # Reaction is of the type A <-> B               {"A": -1, "B": 1}
            return "Unimolecular rearrangement/isomerization"

        if (r == 1 and reactants[0][0] == 1) \
                and (p == 2 and products[0][0] == 1  and products[1][0] == 1):
            # Reaction is of the type A <-> B + C           {"A": -1, "B": 1, "C": 1}
            return "Unimolecular decomposition"

        if (r == 1 and reactants[0][0] == 1) and (p == 1 and products[0][0] == 2):
            # Reaction is of the type A <-> 2 B             {"A": -1, "B": 2}
            return "Unimolecular decomposition"

        if (r == 2 and reactants[0][0] == 1 and reactants[1][0] == 1) \
            and (p == 1 and products[0][0] == 1):
            # Reaction is of the type A + B <-> C           {"A": -1, "B": -1, "C": 1}
            return "Bimolecular synthesis"

        if (r == 1 and reactants[0][0] == 2) and (p == 1 and products[0][0] == 1):
            # Reaction is of the type 2 A <-> C             {"A": -2, "C": 1}
            return "Bimolecular synthesis"

        return "General one-step"



    def _determine_reaction_type(self):
        """

        :return:
        """
        reactant_list = self.stoichiometry.get_reactant_list()
        product_list = self.stoichiometry.get_product_list()

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



    def get_signed_stoichiometric_coefficients(self, reactants :list[tuple], products :list[tuple]) -> dict:
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
        # TODO: maybe move to class Stoichiometry (and turn it from dataclass to regular class, to allow multiple ways to initialize)
        coeffs = {}

        for c, species in reactants:        # Example: (2, "A")
            coeffs[species] = coeffs.get(species, 0) - c    # Accumulate the sum of the stoichiometric coefficients for this species

        for c, species in products:         # Example: (1, "P")
            coeffs[species] = coeffs.get(species, 0) + c    # Accumulate the sum of the stoichiometric coefficients for this species

        return coeffs



    def extract_rxn_properties(self) -> dict:
        """
        Create a dictionary with the numerical properties of the given reaction
        (skipping any None values)
        Possible values include:
            - forward and reverse reaction rates (kR and kR, respectively)
            - ΔH, ΔS, ΔG,
            - K (equilibrium constant)

        :return:    EXAMPLE: {'reaction_model': 'mass action', 'delta_H': -30, 'K_eq': 5.0, 'K': 5, 'kF': 10, 'kR': 2, 'reversible': True}
        """
        thermo_properties = self.thermodynamics.to_dict()
        sim_rxn_tuple = self.sim_reactions
        if len(sim_rxn_tuple) == 1:
            sim_rxn = sim_rxn_tuple[0]
            kinetic_properties = sim_rxn.model.get_parameters()
            kinetic_properties["reaction_model"] = sim_rxn.model.name
        else:
            kinetic_properties = {}
            for i, sim_rxn in enumerate(sim_rxn_tuple):
                kinetic_properties[f"Rxn{i}: reaction_model"] = sim_rxn.model.name
                for k, v in sim_rxn.model.get_parameters().items():
                    kinetic_properties[f"Rxn{i}: {k}"] = v


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

        if self.thermodynamics.K_eq is not None:
            for sr in self.sim_reactions:
                Reconciler.reconcile(thermodynamics_data=self.thermodynamics, model=sr.model)
            #self.kinetics.set_rate_constants_from_equilibrium_constant(K=self.thermodynamics.K_eq)



    def extract_intermediate(self) -> str|None:
        """
        Return the name of the reaction intermediate species,
        or None if there's no intermediate.

        If more than 1 intermediate is present, raise an Exception

        :return:    The species ID of the reaction intermediate, if present;
                        or None if not present
        """
        sim_rxn_tuple = self.sim_reactions
        #print(len(sim_rxn_tuple))
        if len(sim_rxn_tuple) < 2:
            # If at most one SimulationReaction
            return None

        assert len(sim_rxn_tuple) < 3, \
            "extract_intermediate(): currently not implemented for cases when the reaction definition " \
            "compiles into more than 2 simulation reactions"

        # If we get thus far, we have exactly 2 elements in the tuple
        st_0 = sim_rxn_tuple[0].stoichiometry.to_dict()
        st_1 = sim_rxn_tuple[1].stoichiometry.to_dict()
        overlap = set(st_0) & set(st_1)     # Set intersection
        #print("overlap: ", overlap)        # EXAMPLE: {'ES*', 'E'}

        overlap -= sim_rxn_tuple[0].stoichiometry.get_reactant_ids()     # Set difference
        overlap -= sim_rxn_tuple[1].stoichiometry.get_product_ids()      # Set difference

        assert len(overlap) < 2, \
            "extract_intermediate(): currently not implemented for cases when " \
            "there is more than 1 intermediary"

        if overlap == set():
            return None     # No overlap

        return overlap.pop()    # Extract one element from the set

        #TODO: generalize



    def describe(self, concise=False) -> str:
        """
        Return, as a string, a user-friendly plain-text form of the reaction,
        plus a fair deal of optional information

        :param concise:     If True, less detail is shown
        :return:            A string with a description of this reaction
        """
        # TODO: put to good use the new describe() method of the "SimulationReaction" class
        sim_rxn_tuple = self.sim_reactions

        reversible = False
        if len(sim_rxn_tuple) == 1:
            # If there's exactly 1 derived reaction, extract the "reversible" attribute from it, if possible
            sim_rxn = sim_rxn_tuple[0]
            reversible = getattr(sim_rxn.model, 'reversible', False)    # Note: the "reversible" attribute may or may not be present

        rxn_description = self.stoichiometry.standard_chemical_formula(reversible=reversible)

        if concise:
            return rxn_description      # Minimalist description


        # If we get this far, we're looking for a more detailed description

        INDENT = "        "

        rxn_description += "\n" + INDENT + self.reaction_category + " reaction"
        if self.reaction_model:
            rxn_description += f', with Reaction Model: "{self.reaction_model}"'
        else:
            rxn_description += ', with no Reaction Model specified'

        if self.id:
             rxn_description += f"   (Reaction ID {self.id})"


        rxn_description += f"\n{INDENT}Thermodynamics - passed:  "
        s = self.format_reaction_details(self.source_thermodynamic_parameters)
        rxn_description += s if s else "None"

        rxn_description += f"\n{INDENT}Thermodynamics - derived: "
        s = self.format_reaction_details(self.thermodynamics.to_dict())
        rxn_description += s if s else "None"

        rxn_description += f"\n{INDENT}Kinetics - passed: "
        s = self.format_reaction_details(self.source_kinetic_parameters)
        rxn_description += s if s else "None"

        rxn_description += f"\n{INDENT}Kinetics - derived:  "
        s = f"{len(sim_rxn_tuple)} derived reaction"
        rxn_description += s if len(sim_rxn_tuple) > 0 else "None"

        if len(sim_rxn_tuple) > 1:
            rxn_description += "s"  # The plural form

        for i, sim_rxn in enumerate(sim_rxn_tuple):
            rxn_description += f"\n{INDENT}    "
            #if len(sim_rxn_tuple) > 1:
            rxn_description += f"({i+1}) "     # Show the numbering, if more than one
            rxn_description += f'Type: "{sim_rxn.model.name}" {self.format_reaction_details(sim_rxn.model.get_parameters())}'

        return rxn_description



    def extract_reactant_ids(self) -> set[str]:
        """
        Return the set of ALL the reactant species id's in this reaction
        (including any catalysts, if applicable)

        :return:    A set of species id's
        """
        return self.stoichiometry.get_reactant_ids()



    def extract_reactants_formula(self) -> str:
        """
        Return a string with a user-friendly form of the left (reactants) side of the reaction formula
        (aka the reactant "complex")

        :return:    A string with the left (reactant) side of the reaction formula
        """
        return self._standard_form_chem_eqn(self.stoichiometry.get_reactant_list())



    def extract_product_ids(self) -> set[str]:
        """
        Return the set of ALL the product id's in this reaction
        (including any catalysts, if applicable)

        :return:    A set of species id's
        """
        return self.stoichiometry.get_product_ids()



    def extract_products_formula(self) -> str:
        """
        Return a string with a user-friendly form of the right (products) side of the reaction formula
        (aka the product "complex")

        :return:    A string with the right (product) side of the reaction formula
        """
        return self._standard_form_chem_eqn(self.stoichiometry.get_product_list())



    def extract_species_in_reaction(self, include_intermediaries=True) -> Set[str]:
        """
        Return a SET of the id's of ALL the species appearing in this reaction

        :return:    A SET of the id's of the species involved in this reaction.
                        Note: being a set, it's NOT in any particular order
        """
        species = self.stoichiometry.get_all_species_ids()

        if not include_intermediaries:
            return species

        sim_rxn_list = self.sim_reactions
        for sim_rxn in sim_rxn_list:
            species |= sim_rxn.stoichiometry.get_all_species_ids()         # set union

        return species



    def reaction_quotient(self, conc, explain=False) -> np.double | tuple[np.double, str]:
        """
        Compute the "Reaction Quotient" (aka "Mass–action Ratio"),
        given the concentrations of chemicals involved in this reaction.

        Note: this implementation only covers reactions that have "mass action" kinetics

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
        assert self.reaction_model == "mass action", \
            "reaction_quotient(): only 'mass action' reaction models are currently supported"

        sim_rxn = self.sim_reactions[0]
        return sim_rxn.reaction_quotient(conc=conc, explain=explain)



    def find_equilibrium_conc(self, conc_dict :dict) -> dict:
        """
        Determine the equilibrium concentrations that would be eventually reached
        by the species in this reaction,
        given their initial concentrations,
        IN THE ABSENCE of any other reaction.

        :param conc_dict:   A dict mapping species id's to their initial concentrations,
                                for all the species involved in this reaction
                                EXAMPLE:  {"X": 4.3, "Y": 1.5, "F": 31.6, "G": 3.6}

        :return:            A dict mapping the above chemical id's to their equilibrium concentrations
        """
        # TODO: move to "SimulationReaction" object
        reactants = self.stoichiometry.get_reactant_list()
        products = self.stoichiometry.get_product_list()

        if self.reaction_model != "mass action":
            raise Exception("find_equilibrium_conc(): only 'mass action' reaction models are currently supported")


        """
        # If we get thus far, we have MASS-ACTION kinetics
        """
        assert len(reactants) <= 2, \
                f"find_equilibrium_conc(): for reactions that exhibit mass-action kinetics, " \
                f"it's only implemented when there are no more than 2 reactants (number of reactants is {len(reactants)})"
        assert len(products) <= 2, \
                f"find_equilibrium_conc(): for reactions that exhibit mass-action kinetics, " \
                f"it's only implemented when there are no more than 2 products (number of products is {len(products)})"


        # To conform with functions available in ReactionKinetics,
        # we'll express the reaction in the form   A + B <-> P + Q   or  2 A <-> P  or  A <-> 2 P   (no other solver is available)

        standard_names = ["A", "B", "P", "Q"]
        coeffs = [0, 0, 0, 0]       # For general reaction A + B <-> P + Q
        concs  = [0., 0., 0., 0.]   # For A0, B0, P0, Q0

        reaction_vector = self.stoichiometry.get_reaction_vector()  # Note: catalysts will NOT appear
        #print("reaction vectors (reactants) (products): ", vec_r, vec_p)

        name_map = {}   # To map standard names into actual species names


        # Process reactants (which will set the first element or first two, in coeffs and concs)
        index = 0
        for sp, c in reaction_vector.items():
            if c > 0:
                continue    # Skip the products

            coeffs[index] = -c      # Negative of stoichiometric coefficient because it's a reactant
            conc = conc_dict.get(sp)
            assert conc is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                       f"concentration of the reactant `{sp}` was not provided"
            concs[index] = conc
            std_name = standard_names[index]
            name_map[std_name] = sp
            index += 1

        # Process products (which will set the next element or two, in coeffs and concs)
        index = 2
        for sp, c in reaction_vector.items():
            if c < 0:
                continue    # Skip the reactants

            coeffs[index] = c
            conc = conc_dict.get(sp)
            assert conc is not None, f"find_equilibrium_conc(): unable to proceed because the " \
                                       f"concentration of the product `{sp}` was not provided"
            concs[index] = conc
            std_name = standard_names[index]
            name_map[std_name] = sp
            index += 1


        # Unpack
        a, b, p, q = coeffs
        A0, B0, P0, Q0 = concs

        sim_rxm = self.sim_reactions[0]

        """
        print(f"coeffs: {coeffs} | concs: {concs} | name_map: {name_map}")
        print(f"a: {a} | b: {b} | p: {p} | q: {q}")
        print(f"A0: {A0} | B0: {B0} | P0: {P0} | Q0: {Q0}")
        print(f"kF: {sim_rxm.model.kF} | kR: {sim_rxm.model.kR}")
        print(self.analytic_solution_family)
        """

        if (self.analytic_solution_family == "ONE_TO_ONE"):
            # Reaction is of the form A <-> P
            eq_dict = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=sim_rxm.model.kF,
                                                                            kR=sim_rxm.model.kR,
                                                                            A0=A0, P0=P0)

        elif (self.analytic_solution_family == "TWO_TO_ONE") and (a == 1):
            # Reaction is of the form A + B <-> P
            eq_dict = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=sim_rxm.model.kF,
                                                                            kR=sim_rxm.model.kR,
                                                                            A0=A0, B0=B0, P0=P0)

        elif (self.analytic_solution_family == "TWO_TO_ONE") and (a == 2):
            # Reaction is of the form 2 A <-> P
            eq_dict = ReactionKinetics.compute_equilibrium_conc_elementary_synthesis(kF=sim_rxm.model.kF,
                                                                                     kR=sim_rxm.model.kR,
                                                                                     A0=A0, P0=P0)

        elif (self.analytic_solution_family == "ONE_TO_TWO") and (p == 1):
            # Reaction is of the form A <-> P + Q
            eq_dict = ReactionKinetics.compute_equilibrium_conc_mass_action(kF=sim_rxm.model.kF,
                                                                            kR=sim_rxm.model.kR,
                                                                            A0=A0, P0=P0, Q0=Q0)

        elif (self.analytic_solution_family == "ONE_TO_TWO") and (p == 2):
            # Reaction is of the form A <-> 2 P
            eq_dict = ReactionKinetics.compute_equilibrium_conc_elementary_decomposition(kF=sim_rxm.model.kF,
                                                                                         kR=sim_rxm.model.kR,
                                                                                         A0=A0, P0=P0)

        else:
            raise Exception(f"find_equilibrium_conc(): Not implemented for this reaction type ({self.analytic_solution_family})")

        """                                                                      
        eq_dict = ReactionKinetics._compute_equilibrium_conc_first_order(kF=sim_rxm.model.kF, kR=sim_rxm.model.kR,
                                                                         a=a, b=b, p=c, q=d,
                                                                         A0=A0, B0=B0, P0=C0, Q0=D0)
        """
        #print("eq_dict", eq_dict)

        # eq_dict contains the keys "A", "B", "P", "Q";
        # Translate the standard names A, B, P, Q into the actual species names, and also drop any missing term
        converted_dict = {}
        for k, v in eq_dict.items():
            converted_name = name_map[k]
            converted_dict[converted_name] = v

        return converted_dict





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
        # TODO: probably switch to using the new Stoichiometry dataclass
        #       -> See the new Stoichiometry._standard_form_complex()

        assert type(eqn_side) == list, \
            f"Reaction._standard_form_chem_eqn(): the argument must be a list (it was of type {type(eqn_side)})"

        formula_list = []
        for t in eqn_side:
            stoichiometry, species_name = t

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
        #if self.kinetics.law != "mass action":
        if (self.reaction_model is not None) and (self.reaction_model != "mass action"):
            return None

        r, p, c = self.stoichiometry.reaction_pattern()     # number of reactants, products, catalysts

        if c > 0:      # If enzymes were involved
            return None

        # TODO: switch to using signed terms
        reactants = self.stoichiometry.get_reactant_list()
        products = self.stoichiometry.get_product_list()
        if (r == 1 and reactants[0][0] == 1) and (p == 1 and products[0][0] == 1):
            # Reaction is of the type A <-> B               {"A": -1, "B": 1}
            return "ONE_TO_ONE"

        if (r == 1 and reactants[0][0] == 1) \
                and (p == 2 and products[0][0] == 1  and products[1][0] == 1):
            # Reaction is of the type A <-> B + C           {"A": -1, "B": 1, "C": 1}
            return "ONE_TO_TWO"

        if (r == 1 and reactants[0][0] == 1) and (p == 1 and products[0][0] == 2):
            # Reaction is of the type A <-> 2 B             {"A": -1, "B": 2}
            return "ONE_TO_TWO"

        if (r == 2 and reactants[0][0] == 1 and reactants[1][0] == 1) \
            and (p == 1 and products[0][0] == 1):
            # Reaction is of the type A + B <-> C           {"A": -1, "B": -1, "C": 1}
            return "TWO_TO_ONE"

        if (r == 1 and reactants[0][0] == 2) and (p == 1 and products[0][0] == 1):
            # Reaction is of the type 2 A <-> C             {"A": -2, "C": 1}
            return "TWO_TO_ONE"

        return None



    def format_reaction_details(self, rxn_properties :dict) -> str:
        """
        Format and return a string with some details about the parameters of this reaction,
        contained in the passed dictionary.
        Also, add units of measurements; any property named "temp" gets converted from degree K to C.

        :param rxn_properties:  A dictionary with numerical properties of interest for the reaction.
                                    Any None value will be dropped
                                    EXAMPLE: {'kF': 3.0, 'kR': None, 'delta_G': 1.2345, 'K': 1.5}

        :return:                A string with some details about the parameters of this reaction
                                    EXAMPLE: "  (kF = 3 | delta_G = 1.2345 kJ/mol | Temp = 25 C)"
        """
        # TODO: probably move to SimulationReaction
        #print("rxn_properties: ", rxn_properties)
        details = []    # Running list of strings with each of the individual details

        for k,v in rxn_properties.items():
            if v is None:
                continue

            if k == "temp":
                single_detail = f"Temp = {convert(v, from_unit=K, to_unit=C):,.4g} C"
                # EXAMPLE: "Temp = 25 C"
                details.append(single_detail)
                continue

            if type(v) is str:
                single_detail = f"{k} = '{v}'"
            elif type(v) is bool:
                single_detail = f"{k} = {v}"
            elif callable(v):
                #single_detail = f'{k} = function "{v.__name__}"'
                single_detail = f'{k} = {v.__name__}()'
            else:   # Numeric
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