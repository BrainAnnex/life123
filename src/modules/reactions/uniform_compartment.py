import math
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from typing import Union
from src.modules.chemicals.chem_data import ChemData
from src.modules.movies.movies import MovieTabular
from src.modules.numerical.numerical import Numerical as num
from src.modules.visualization.plotly_helper import PlotlyHelper



class RxnDynamics:
    # TODO: rename "ReactionDynamics", and move to separate file

    @staticmethod
    def norm_A(delta_conc :np.array) -> float:
        """
        Return a measure of system change, based on the average concentration changes
        of ALL the specified chemicals across a time step, adjusted for the number of chemicals.
        A square-of-sums computation (the square of an L2 norm) is used.

        :param delta_conc:  A Numpy array with the concentration changes
                                of the chemicals of interest across a time step
        :return:            A measure of change in the concentrations across the simulation step
        """
        n_active_chemicals = len(delta_conc)

        # The following are normalized by the number of chemicals
        #L2_rate = np.linalg.norm(delta_concentrations) / n_chems
        #L2_rate = np.sqrt(np.sum(delta_concentrations * delta_concentrations)) / n_chems
        #print("    L_inf norm:   ", np.linalg.norm(delta_concentrations, ord=np.inf) / delta_time)
        #print("    Adjusted L1 norm:   ", np.linalg.norm(delta_concentrations, ord=1) / n_chems)

        adjusted_L2_rate = np.sum(delta_conc * delta_conc) / (n_active_chemicals * n_active_chemicals)
        return adjusted_L2_rate


    @staticmethod
    def norm_B(baseline_conc: np.array, delta_conc: np.array) -> float:
        """
        Return a measure of system change, based on the max absolute relative concentration
        change of all the chemicals across a time step (based on an L infinity norm - but disregarding
        any baseline concentration that is very close to zero)

        :param baseline_conc:   A Numpy array with the concentration of the chemicals of interest
                                    at the start of a simulation time step
        :param delta_conc:      A Numpy array with the concentration changes
                                    of the chemicals of interest across a time step
        :return:                A measure of change in the concentrations across the simulation step
        """
        arr_size = len(baseline_conc)
        assert len(delta_conc) == arr_size, "norm_B(): mismatch in the sizes of the 2 passed array arguments"

        to_keep = ~ np.isclose(baseline_conc, 0)    # Element-wise negation; this will be an array of Booleans
        # with True for all the elements of baseline_conc that aren't too close to 0

        ratios = delta_conc[to_keep] / baseline_conc[to_keep]   # Using boolean indexing to only select some of the elements :
        # the non-zero denominators, and their corresponding numerators
        return max(abs(ratios))


    @staticmethod
    def norm_C(prev_conc: np.array, baseline_conc: np.array, delta_conc: np.array) -> float:
        """
        Return a measure of system short-period oscillation; larger values might be heralding
        onset of simulation instability

        :param prev_conc:       A numpy array with the concentration of the chemicals of interest,
                                    in the step prior to the current one (i.e. an "archive" value)
        :param baseline_conc:   A Numpy array with the concentration of the chemicals of interest
                                    at the start of a simulation time step
        :param delta_conc:      A Numpy array with the concentration changes
                                    of the chemicals of interest across a time step
        :return:
        """
        if prev_conc is None:
            return 0            # Unable to compute a norm; 0 represents "perfect"

        D1 = baseline_conc - prev_conc
        D2 = delta_conc

        #print("****** D1: ", D1)
        #print("****** D2: ", D2)

        sign_flip = ((D1 >= 0) & (D2 < 0)) | ((D1 < 0) & (D2 >= 0))

        criterion_met = sign_flip & (abs(D2) > abs(D1)) & ~np.isclose(D1, 0) & (abs(D2) < 50 * abs(D1))
        #criterion_met = sign_flip & ~np.isclose(D1, 0) & (abs(D2) < 50 * abs(D1))
        # Note: values with D1 very close to zero are ignored
        #       likewise, values where |D2| dwarfs |D1| are ignored
        #print("criterion_met: ", criterion_met)

        # Use boolean indexing to select elements from delta1 where the criterion is True
        D1_selected = D1[criterion_met]
        D2_selected = D2[criterion_met]

        #print("**** delta1_selected: ", D1_selected)
        #print("*** delta2_selected: ", D2_selected)

        ratios = abs(D2_selected / D1_selected)
        #print(ratios)

        if len(ratios) == 0:
            return 0
        else:
            return np.max(ratios)    # An argument might be made for taking the avg SUM instead


    @staticmethod
    def norm_D(prev_conc: np.array, baseline_conc: np.array, delta_conc: np.array) -> float:
        """
        Return a measure of curvature in the concentration vs. time curves; larger values might be heralding
        onset of simulation instability

        :param prev_conc:       A numpy array with the concentration of the chemicals of interest,
                                    in the step prior to the current one (i.e. an "archive" value)
        :param baseline_conc:   A Numpy array with the concentration of the chemicals of interest
                                    at the start of a simulation time step
        :param delta_conc:      A Numpy array with the concentration changes
                                    of the chemicals of interest across a time step
        :return:
        """
        if prev_conc is None:
            return 0            # Unable to compute a norm; 0 represents "perfect"

        D1 = baseline_conc - prev_conc
        D2 = delta_conc

        #print("\nnorm_D ****** D1: ", D1)
        #print("norm_D ****** D2: ", D2)

        criterion_met = ~np.isclose(D1, 0) & (abs(D2) < 100 * abs(D1))
        # Note: values with D1 very close to zero are ignored
        #       likewise, values where |D2| dwarfs |D1| are ignored
        #print("criterion_met: ", criterion_met)

        # Use boolean indexing to select elements from delta1 where the criterion is True
        D1_selected = D1[criterion_met]
        D2_selected = D2[criterion_met]

        #print("norm_D **** delta1_selected: ", D1_selected)
        #print("norm_D *** delta2_selected: ", D2_selected)

        ratios = abs(D2_selected / D1_selected)
        #print(ratios)

        res = np.sum(ratios)    # An argument might be made for taking the MAX instead

        arr_size = len(baseline_conc)

        normalized_res = res / arr_size
        #print("norm_D *** : ", normalized_res)
        return normalized_res





#############################################################################################

class ExcessiveTimeStepHard(Exception):
    """
    Used to raise Exceptions arising from excessively large time steps
    (that lead to negative concentration values, i.e. "HARD" errors)
    """
    pass

class ExcessiveTimeStepSoft(Exception):
    """
    Used to raise Exceptions arising from excessively large time steps
    (that lead to norms regarded as excessive because of user-specified values, i.e. "SOFT" errors)
    """
    pass


#############################################################################################


class UniformCompartment:
    """
    Used to simulate the dynamics of reactions (in a single compartment.)
    In the context of Life123, this may be thought of as a "zero-dimensional system"

    (This class was formerly named "ReactionDynamics")
    """

    #TODO: maybe split off part of this class into a separate "ReactionDynamics" class


    def __init__(self, chem_data=None, names=None, preset="mid"):
        """
        Note: AT MOST 1 of the following 2 arguments can be passed
        :param chem_data:   [OPTIONAL 1] Object of type "ChemData" (with data about the chemicals and their reactions)
                                       NOTE: This is the recommended argument; the other 2 options are, or might be, deprecated
        :param names:       [OPTIONAL 2 - MIGHT GET DEPRECATED] A single name, or list or tuple of names, of the chemicals
                                       the reactions can be added later, with calls to add_reaction()

        :param preset:  String with code that can be adjusted make the time resolution finer or coarser;
                        it will stay in effect from now on, unless explicitly changed later
        """
        number_args = 0
        if chem_data:
            number_args += 1
        if names:
            number_args += 1


        assert number_args <= 1, \
            f"UniformCompartment instantiation: Can only pass at most one of the arguments " \
            f"`chem_data`, and `names` ({number_args} were passed)"

        self.chem_data = chem_data  # Object of type "ChemData" (with data about the chemicals and their reactions,
                                    #                            incl. macromolecules)

        self.system_time = 0.       # Global time of the system, from initialization on

        # TODO: maybe rename "system" to "system_state", and use "system" to store a list or dict of the chemicals
        #       actually involved in this dynamic simulation
        self.system = None  # Concentration data in the single compartment we're simulating, for all the (bulk) chemicals
                            # A Numpy array of the concentrations of all the chemical species, in their index order
                            # 1-dimensional NumPy array of floats, whose size is the number of chemical species.
                            # Each entry is the concentration of the species with that index (in the "ChemData" object)
                            # Note that this is the counterpart - with 1 less dimension - of the array by the same name
                            #       in the class BioSim1D

        self.previous_system = None # Concentration data of all the chemicals at the previous simulation step


        self.macro_system = {}      # A dict mapping macromolecule names to their counts
                                    # EXAMPLE:  {"M1": 1, "M2": 3, "M3": 1}

        self.macro_system_state = {}  # The counterpart of the system data for macro-molecules, if present
                                    # Binding fractions of the applicable transcription factors, for all the macro-molecules,
                                    # over the previous time step, indexed by macromolecule and by binding site number
                                    # EXAMPLE:   {"M1": {1: ("A", 0.2), 2: ("F", 0.93), ("P", 0.)},
                                    #             "M2": {1: ("C", 0.5), 2: ("F", 0.1)},
                                    #             "M3": {} }

                                    # For background, see: https://www.annualreviews.org/doi/10.1146/annurev-cellbio-100617-062719

        self.history = MovieTabular()   # To store user-selected snapshots of (some of) the chemical concentrations,
                                        #   whenever requested by the user.


        # ***  FOR AUTOMATED ADAPTIVE TIME STEP SIZES  ***
        # Note: The "aborts" below are "elective" aborts - i.e. not aborts from hard errors (further below)
        #       The default values get packed into a "preset", specified by a name, and optionally passed when
        #       instantiating this object.
        #       Users in need to change these values post-instantiation will generally use use_adaptive_preset(), below;
        #       or, for more control, set_thresholds() and set_step_factors()
        self.thresholds = []    # A list of "rules"
        # EXAMPLE:                  [
        #                             {"norm": "norm_A", "low": 0.5, "high": 0.8, "abort": 1.44},
        #                             {"norm": "norm_B", "low": 0.08, "high": 0.5, "abort": 1.5}
        #                           ]
        #       low:    The value below which the norm is considered to be good (and the variable steps are being made larger);
        #                if None, it doesn't get used
        #       high:   The value above which the norm is considered to be excessive (and the variable steps get reduced);
        #               if None, it doesn't get used
        #       abort:  The value above which the norm is considered to be dangerously large (and the last variable step gets
        #               discarded and back-tracked with a smaller size)
        #               if None, it doesn't get used

        self.step_factors = {}
                            # EXAMPLE: {"upshift": 1.2, "downshift": 0.5, "abort": 0.4, "error": 0.2}
                            # "upshift" must be > 1 ; all the other values must be < 1
                            # Generally, error <= abort <= downshift
                            # "Error" value: Factor by which to multiply the time step
                            #   in case of negative-concentration error from excessive step size
                            #   NOTE: this is from ERROR aborts,
                            #   not to be confused with general aborts based on reaching high threshold



        self.use_adaptive_preset(preset)

        self.number_neg_concs = 0
        self.number_soft_aborts = 0

        self.norm_usage = {"norm_A": 0, "norm_B": 0, "norm_C": 0, "norm_D": 0}  # Number of times that each norm
                                                                                #   got involved in step-size decision



        # ***  FOR DIAGNOSTICS  ***     TODO: maybe turn all diagnostic data/methods into a separate object

        self.verbose_list = []          # A list of integers or strings with the codes of the desired verbose checkpoints
                                        #   EXAMPLE: [1, "my_ad_hoc_tag"] to invoke sections of code marked as 1 or 3
                                        #   Those sections will have entry points such as:  if "my_ad_hoc_tag" in self.verbose_list


        self.diagnostics = False        # Overall flag about whether using diagnostics

        self.diagnostic_rxn_data = {}   # "Diagnostic reaction data", PER REACTION: a dict with as many entries as reactions.
                                        #   The keys are the reaction indices; the values are objects of type "MovieTabular",
                                        #   which contain Pandas dataframes with the following columns
                                        #   (referring to 1 reaction):
                                        #           'START TIME' , 'Delta A' , 'Delta B' ...'time_step' , 'caption'
                                        #
                                        #   Notes:  - entries are always added, even if an interval run is aborted
                                        #           - the various 'Delta concentrations' are for ALL the chemicals (in the reaction or not),
                                        #                   over the time interval that *STARTS* at the value in the "TIME" column

        self.diagnostic_conc_data = MovieTabular(parameter_name="TIME")
                                        # An expanded version of the normal System History.
                                        #   Columns of the dataframes:
                                        #       'TIME' 	'A' 'B' ...  'caption'
                                        #
                                        #   Note: if an interval run is aborted, NO entry is created here
                                        #         (this approach DIFFERS from that of other diagnostic data)

        self.diagnostic_decisions_data = MovieTabular(parameter_name="START_TIME")
                                        #   Columns of the dataframes:
                                        #       'START_TIME' 	'Delta A' 'Delta B' ...
                                        #               [plus, if applicable, other fields such as 'action', 'norm_A', 'norm_B', 'step_factors']
                                        #
                                        #   Note: entries are always added, even if an interval run is aborted

        if names:
            self.chem_data = ChemData(names=names)






    #####################################################################################################

    '''                             ~   TO SET/READ SYSTEM DATA   ~                                   '''

    def ________TO_SET_AND_READ_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def set_conc(self, conc: Union[list, tuple, dict], snapshot=True) -> None:
        """
        Set the concentrations of ALL the chemicals at once

        :param conc:    EITHER
                            (1) a list or tuple of concentration values for ALL the registered chemicals,
                                in their index order
                            OR
                            (2) a dict indexed by the chemical names, for some or all of the chemicals
                                EXAMPLE:  {"A": 12.4, "B": 0.23, "E": 2.6}  Anything not specified will be set to zero

                        Note: any previous values will get over-written

        :param snapshot:[OPTIONAL] If True (default), add to the history
                            a snapshot of this state being set
        :return:        None
        """
        # TODO: rename set_conc_all()
        # TODO: more validations, incl. of individual values being wrong data type
        # TODO: also allow a Numpy array; make sure to do a copy() to it!

        if type(conc) == list or type(conc) == tuple:
            assert len(conc) == self.chem_data.number_of_chemicals(), \
                f"UniformCompartment.set_conc(): when a list or tuple is passed as the argument 'conc', " \
                f"its size must match the number of declared chemicals ({self.chem_data.number_of_chemicals()})"

            assert min(conc) >= 0, \
                f"UniformCompartment.set_conc(): values meant to be chemical concentrations cannot be negative " \
                f"(such as the passed value {min(conc)})"

            self.system = np.array(conc, dtype='d')      # float64      TODO: allow users to specify the type

        elif type(conc) == dict:
            for name, conc_value in conc.items():
                self.set_single_conc(conc=conc_value, species_name=name, snapshot=False)

        if snapshot:
            self.add_snapshot(caption="Initialized state")



    def set_single_conc(self, conc, species_index=None, species_name=None, snapshot=True) -> None:
        """
        Set the concentrations of 1 chemical
        Note: if both species_index and species_name are provided, species_name is used     TODO: generate an error instead; and pytest this part

        :param conc:            A non-negative number with the desired concentration value
                                    of the chemical specified below.  (Any previous value will get over-written)
        :param species_index:   (OPTIONAL) An integer that indexes the chemical of interest (numbering starts at 0)
        :param species_name:    (OPTIONAL) A name for the chemical of interest.
                                    If both species_index and species_name are provided, species_name is used
                                    At least one of "species_index" and "species_name" must be specified
        :param snapshot:        (OPTIONAL) boolean: if True, add to the history
                                    a snapshot of this state being set.  Default: True
        :return:                None
        """
        # Validate the arguments
        assert conc >= 0, \
            f"UniformCompartment.set_single_conc(): chemical concentrations cannot be negative (value passed: {conc})"

        if species_name is not None:
            species_index = self.chem_data.get_index(species_name)
        elif species_index is not None:
            self.chem_data.assert_valid_species_index(species_index)
        else:
            raise Exception("UniformCompartment.set_single_conc(): at least one "
                            "of the arguments `species_index` or `species_name` must be provided")


        if self.system is None:
            self.system = np.zeros(self.chem_data.number_of_chemicals(), dtype='d')      # float64      TODO: allow users to specify the type

        # TODO: if setting concentrations of a newly-added chemical, needs to first expand self.system
        self.system[species_index] = conc

        if snapshot:
            self.add_snapshot(caption=f"Set concentration of `{self.chem_data.get_name(species_index)}`")



    def get_system_conc(self) -> np.array:
        """
        Retrieve the concentrations of ALL the chemicals as a Numpy array

        :return:        A Numpy array with the concentrations of ALL the chemicals,
                        in their index order
                            EXAMPLE:  array([12.3, 4.56, 0.12])
                                      The 0-th chemical has concentration 12.3, and so on...
        """
        return self.system



    def get_chem_conc(self, name: str) -> float:
        """
        Return the current system concentration of the given chemical, specified by its name.
        If no chemical by that name exists, an Exception is raised

        :param name:    The name of a chemical species
        :return:        The current system concentration of the above chemical
        """
        species_index = self.chem_data.get_index(name)
        return self.system[species_index]



    def get_conc_dict(self, species=None, system_data=None) -> dict:
        """
        Retrieve the concentrations of the requested chemicals (by default all),
        as a dictionary indexed by the chemical's name

        :param species:     (OPTIONAL) list or tuple of names of the chemical species; by default, return all
        :param system_data: (OPTIONAL) a Numpy array of concentration values, in the same order as the
                                index of the chemical species; by default, use the SYSTEM DATA
                                (which is set and managed by various functions)
        :return:            A dictionary, indexed by chemical name, of the concentration values;
                                EXAMPLE: {"A": 1.2, "D": 4.67}
        """
        if system_data is None:
            system_data = self.system
        else:
            assert system_data.size == self.chem_data.number_of_chemicals(), \
                f"UniformCompartment.get_conc_dict(): the argument `system_data` must be a 1-D Numpy array with as many entries " \
                f"as the declared number of chemicals ({self.chem_data.number_of_chemicals()})"


        if species is None:
            if system_data is None:
                return {}
            else:
                return {self.chem_data.get_name(index): system_data[index]
                        for index, conc in enumerate(system_data)}
        else:
            assert type(species) == list or  type(species) == tuple, \
                f"UniformCompartment.get_conc_dict(): the argument `species` must be a list or tuple" \
                f" (it was of type {type(species)})"

            conc_dict = {}
            for name in species:
                species_index = self.chem_data.get_index(name)
                conc_dict[name] = system_data[species_index]

            return conc_dict




    '''
    Management of reactions
    '''


    def clear_reactions(self) -> None:
        """
        Get rid of all reactions; start again with "an empty slate" (but still with reference
        to the same data object about the chemicals)

        # TODO: maybe offer an option to clear just one reaction, or a list of them
        # TODO: provide support for "inactivating" reactions

        :return:    None
        """
        self.chem_data.clear_reactions_data()





    #####################################################################################################

    '''                                 ~  TO VISUALIZE SYSTEM  ~                                     '''

    def ________TO_VISUALIZE_SYSTEM________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def describe_state(self) -> None:
        """
        Print out various data on the current state of the system
        :return:        None
        """
        print(f"SYSTEM STATE at Time t = {self.system_time:,.8g}:")

        n_species = self.chem_data.number_of_chemicals()
        print(f"{n_species} species:")

        # Show a line of line of data for each chemical species in turn
        for species_index, name in enumerate(self.chem_data.get_all_names()):
            if name:    # If a name was provided, show it
                name = f" ({name})"
            else:
                name = ""

            if self.system is None:
                print(f"  Species {species_index}{name}. No concentrations set yet")
            else:
                print(f"  Species {species_index}{name}. Conc: {self.system[species_index]}")

        if self.macro_system != {}:
            print("Macro-molecules, with their counts: ", self.macro_system)

        if self.macro_system_state != {}:
            print("Fractional Occupancy at the various binding sites for each macro-molecule:")
            for mm, state_dict in self.macro_system_state.items():
                state_list = [f"{a}: {b[1]} ({b[0]})" for a, b in state_dict.items()]   # EXAMPLE: ["3: 0.1 (A)", "8: 0.6 (B)"]
                state_str = " | ".join(state_list)                                      # EXAMPLE: "3: 0.1 (A) | "8: 0.6 (B)"
                print(f"     {mm} || {state_str}")

        if self.chem_data.active_enzymes == set():    # If no enzymes were involved in any reaction
            print(f"Set of chemicals involved in reactions: {self.chem_data.names_of_active_chemicals()}")
        else:
            print(f"Set of chemicals involved in reactions (not counting enzymes): {self.chem_data.names_of_active_chemicals()}")
            print(f"Set of enzymes involved in reactions: {self.chem_data.names_of_enzymes()}")




    #####################################################################################################

    '''                            ~  PARAMETERS FOR ADAPTIVE STEPS  ~                                '''

    def ________PARAMETERS_FOR_ADAPTIVE_STEPS________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def set_thresholds(self, norm :str, low=None, high=None, abort=None) -> None:
        """
        Over-ride default values for the rules (simulation parameters) that affect the adaptive variable step sizes.
        The default None values can be used to eliminate some of the threshold rules, for the specified norm.
        If the specified norm was not previously set, it will be added.

        :param norm:    The name of the "norm" (criterion) being used - either a previously-set one or not
        :param low:     The value below which the norm is considered to be good (and the variable steps are being made larger);
                            if None, any existing value gets taken out
        :param high:    The value above which the norm is considered to be excessive (and the variable steps get reduced);
                            if None, any existing value gets taken out
        :param abort:   The value above which the norm is considered to be dangerously large (and the last variable step gets
                            discarded and back-tracked with a smaller size) ;
                            if None, any existing value gets taken out
        :return:        None
        """
        assert type(norm) == str and norm != "", \
            "set_thresholds(): the `norm` argument must be a non-empty string"

        if (abort is not None) and (high is not None):
            assert abort >= high, \
                f"set_thresholds(): `abort` value ({abort}) must be >= `high' value ({high})"

        if (high is not None) and (low is not None):
            assert high >= low, \
                f"set_thresholds(): `high` value ({high}) must be >= `low' value ({low})"

        for i, t in enumerate(self.thresholds):
            if t.get("norm") == norm:
                if low is None:
                    del t["low"]
                else:
                    t["low"] = low

                if high is None:
                    del t["high"]
                else:
                    t["high"] = high

                if abort is None:
                    del t["abort"]
                else:
                    t["abort"] = abort

                if len(t) == 1:
                    del self.thresholds[i]      # Completely eliminate this un-used norm
                return

        # If we get here, it means that we're handling a norm
        # not currently present in the list self.thresholds
        new_t = {"norm": norm}
        if low is not None:
            new_t["low"] = low
        if high is not None:
            new_t["high"] = high
        if abort is not None:
            new_t["abort"] = abort

        if self.thresholds is None:
            self.thresholds = [new_t]
        else:
            self.thresholds.append(new_t)



    def update_thresholds(self, norm :str, low=None, high=None, abort=None) -> None:
        """
        Change current values for the rules (simulation parameters) that affect the adaptive variable step sizes.
        The default None values can be used to eliminate some of the threshold rules, for the specified norm.
        If the specified norm was not previously defined, an Exception will be raised

        :param norm:    The name of the "norm" (criterion) whose threshold(s) to update.  It must have previously defined,
                            or an Exception will get raised
        :param low:     The value below which the norm is considered to be good (and the variable steps are being made larger);
                            if None, it doesn't get changed
        :param high:    The value above which the norm is considered to be excessive (and the variable steps get reduced);
                            if None, it doesn't get changed
        :param abort:   The value above which the norm is considered to be dangerously large (and the last variable step gets
                            discarded and back-tracked with a smaller size) ;
                            if None, it doesn't get changed
        :return:        None
        """
        #TODO: unit test
        assert type(norm) == str and norm != "", \
            "update_thresholds(): the `norm` argument must be a non-empty string"

        if (abort is not None) and (high is not None):
            assert abort >= high, \
                f"update_thresholds(): `abort` value ({abort}) must be >= `high' value ({high})"

        if (high is not None) and (low is not None):
            assert high >= low, \
                f"update_thresholds(): `high` value ({high}) must be >= `low' value ({low})"

        for rule in self.thresholds:
            if rule.get("norm") == norm:
                if low is not None:
                    rule["low"] = low

                if high is not None:
                    rule["high"] = high

                if abort is not None:
                    rule["abort"] = abort

                return      # A rule with the given norm was found

        # If we get here, it means that the specified norm
        # is not currently present in the list of rules in self.thresholds
        raise Exception(f"update_thresholds(): no norm named '{norm}' was found.  "
                        f"To create a new norm, use set_thresholds()")



    def set_step_factors(self, upshift=None, downshift=None, abort=None, error=None) -> None:
        """
        Over-ride current values for simulation parameters that affect the adaptive variable step sizes.
        Values not explicitly passed will remain the same.

        :param upshift:     Fraction by which to increase the variable step size
        :param downshift:   Fraction by which to decrease the variable step size
        :param abort:       Fraction by which to decrease the variable step size in cases of step re-do
        :param error:       Fraction by which to decrease the variable step size in cases of error
        :return:            None
        """
        if upshift is not None:
            assert upshift > 1, "set_step_factors(): `upshift` value must be > 1"
            self.step_factors["upshift"] = upshift

        if downshift is not None:
            assert 0 < downshift < 1, "set_step_factors(): `downshift` value must be a positive number < 1"
            self.step_factors["downshift"] = downshift

        if abort is not None:
            assert 0 < abort < 1, "set_step_factors(): `abort` value must be a positive number < 1"
            self.step_factors["abort"] = abort

        if error is not None:
            assert 0 < error < 1, "set_step_factors(): `error` value must be a positive number < 1"
            self.step_factors["error"] = error



    def set_error_step_factor(self, value) -> None:
        """
        Over-ride the default value for the simulation parameter error_abort_step_factor

        :param value:   Fraction by which to decrease the variable step size in cases of
                            value errors (such as negative concentrations) during a variable step
        :return:        None
        """
        # TODO: OBSOLETE

        print("\n\n**************  set_error_step_factor() is DEPRECATED : use set_step_factors() instead  **************\n\n")

        self.set_step_factors(error=value)



    def show_adaptive_parameters(self) -> None:
        """
        Print out the current values for the adaptive time-step parameters

        :return:    None
        """
        print("Parameters used for the automated adaptive time step sizes -")
        print("    THRESHOLDS: ", self.thresholds)
        print("    STEP FACTORS: ", self.step_factors)



    def use_adaptive_preset(self, preset :str) -> None:
        """
        Lets the user choose a preset to use from now on, unless later explicitly changed,
        for use in ALL reaction simulations involving adaptive time steps.
        The preset will affect the degree to which the simulation will be "risk-taker" vs. "risk-averse" about
        taking larger steps.

        Note:   for more control, use set_thresholds() and set_step_factors()

                For example, using the "mid" preset is the same as issuing:
                    dynamics.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.44)
                    dynamics.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=1.5)
                    dynamics.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4, error=0.25)

        :param preset:  String with one of the available preset names;
                            allowed values are (in generally-increasing speed):
                            'heavy_brakes', 'slower', 'slow', 'mid', 'fast'
        :return:        None
        """
        if preset == "heavy_brakes":   # It slams on the "brakes" hard in case of abort or errors
            self.thresholds = [{"norm": "norm_A", "low": 0.02, "high": 0.025, "abort": 0.03},
                               {"norm": "norm_B", "low": 0.05, "high": 1.0, "abort": 2.0}]
            self.step_factors = {"upshift": 1.6, "downshift": 0.15, "abort": 0.08, "error": 0.05}

        elif preset == "slower":   # Very conservative about taking larger steps
            self.thresholds = [{"norm": "norm_A", "low": 0.2, "high": 0.5, "abort": 0.8},
                               {"norm": "norm_B", "low": 0.03, "high": 0.05, "abort": 0.5}]
            self.step_factors = {"upshift": 1.01, "downshift": 0.5, "abort": 0.1, "error": 0.1}

        elif preset == "slow":   # Conservative about taking larger steps
            self.thresholds = [{"norm": "norm_A", "low": 0.2, "high": 0.5, "abort": 0.8},
                               {"norm": "norm_B", "low": 0.05, "high": 0.4, "abort": 1.3}]
            self.step_factors = {"upshift": 1.1, "downshift": 0.3, "abort": 0.2, "error": 0.1}

        elif preset == "mid":     # A "middle-of-the road" heuristic: somewhat "conservative" but not overly so
            self.thresholds = [{"norm": "norm_A", "low": 0.5, "high": 0.8, "abort": 1.44},
                               {"norm": "norm_B", "low": 0.08, "high": 0.5, "abort": 1.5}]
            self.step_factors = {"upshift": 1.2, "downshift": 0.5, "abort": 0.4, "error": 0.25}

        elif preset == "fast":   # Less conservative (more "risk-taker") about taking larger steps
            self.thresholds = [{"norm": "norm_A", "low": 0.8, "high": 1.2, "abort": 1.7},
                               {"norm": "norm_B", "low": 0.15, "high": 0.8, "abort": 1.8}]
            self.step_factors = {"upshift": 1.5, "downshift": 0.8, "abort": 0.6, "error": 0.5}

        elif preset == "mid_inclusive": # A "middle-of-the road" heuristic that makes use of more norms
            self.thresholds = [{'norm': 'norm_A', 'low': 0.2, 'high': 0.8, 'abort': 1.44},
                               {'norm': 'norm_B', 'low': 0.08, 'high': 0.5, 'abort': 1.5},
                               {'norm': 'norm_C', 'low': 0.5, 'high': 1.2, 'abort': 1.6},
                               {'norm': 'norm_D', 'low': 1.3, 'high': 1.7, 'abort': 1.8}]
            self.step_factors = {'upshift': 1.1, 'downshift': 0.5, 'abort': 0.4, 'error': 0.25}

        else:
            raise Exception(f"set_adaptive_parameters(): unknown value for the `preset` argument ({preset}); "
                            f"allowed values are 'heavy_brakes', 'slower', 'slow', 'mid', 'fast', 'mid_inclusive'")




    #####################################################################################################

    '''                          ~  TO SET/DESCRIBE THE REACTIONS  ~                                  '''

    def ________TO_SET_AND_DESCRIBE_REACTIONS________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def add_reaction(self, **kwargs) -> int:
        """
        Register a new SINGLE chemical reaction,
        optionally including its kinetic and/or thermodynamic data.
        All the involved chemicals must be already registered - use add_chemical() if needed.

        For details, see ChemData.add_reaction()

        :param kwargs:  Any arbitrary named arguments
        :return:        Integer index of the newly-added reaction
        """
        assert self.chem_data, \
            "add_reaction(): must first register the names of the chemicals.  You may use add_chemical()"

        return self.chem_data.add_reaction(**kwargs)



    def describe_reactions(self, **kwargs) -> None:
        """
        Print out a user-friendly plain-text form of ALL the reactions.
        For details, see ChemData.describe_reactions()

        :param kwargs:  Any arbitrary named arguments
        :return:        None
        """
        self.chem_data.describe_reactions(**kwargs)



    def number_of_reactions(self) -> int:
        """
        Return the number of registered chemical reactions

        :return:    The number of registered chemical reactions
        """
        return self.chem_data.number_of_reactions()



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
        self.chem_data.plot_reaction_network(graphic_component=graphic_component, unpack=False)





    #####################################################################################################

    '''                            ~  TO PERFORM THE REACTIONS  ~                                     '''

    def ________TO_PERFORM_THE_REACTIONS________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def specify_steps(self, total_duration=None, time_step=None, n_steps=None) -> (float, int):
        """
        If either the time_step or n_steps is not provided (but at least 1 of them must be present),
        determine the other one from total_duration

        Their desired relationship is: total_duration = time_step * n_steps

        :param total_duration:  Float with the overall time advance (i.e. time_step * n_steps)
        :param time_step:       Float with the size of each time step
        :param n_steps:         Integer with the desired number of steps
        :return:                The pair (time_step, n_steps)
        """
        assert (not total_duration or not time_step or not n_steps), \
            "UniformCompartment.specify_steps(): cannot specify all 3 arguments: `total_duration`, `time_step`, `n_steps` (specify only any 2 of them)"

        assert (total_duration and time_step) or (total_duration and n_steps) or (time_step and n_steps), \
            "UniformCompartment.specify_steps(): must provide exactly 2 arguments from:  `total_duration`, `time_step`, `n_steps`"

        if not time_step:
            time_step = total_duration / n_steps

        if not n_steps:
            n_steps = math.ceil(total_duration / time_step)

        return (time_step, n_steps)     # Note: could opt to also return total_duration if there's a need for it



    def single_compartment_react(self, duration=None, target_end_time=None, stop=None,
                                 initial_step=None, n_steps=None, max_steps=None,
                                 snapshots=None, silent=False,
                                 variable_steps=True, explain_variable_steps=None,
                                 reaction_duration=None) -> None:
        """
        Perform ALL the reactions in the single compartment -
        based on the INITIAL concentrations,
        which are used as the basis for all the reactions.

        Update the system state and the system time accordingly
        (object attributes self.system and self.system_time)

        :param duration:        The overall time advance for the reactions (it may be exceeded in case of variable steps)
        :param reaction_duration: [OBSOLETE OLD NAME for "duration"; being phased out]
        :param target_end_time: The final time at which to stop the reaction; it may be exceeded in case of variable steps
                                    If both `target_end_time` and `duration` are specified, an error will result

        :param initial_step:    The suggested size of the first step (it might be reduced automatically,
                                    in case of "hard" errors resulting from overly-large steps)

        :param stop:            Pair of the form (termination_keyword, termination_parameter), to indicate
                                    the criterion to use to stop the reaction
                                    EXAMPLES:
                                        ("conc_below", (chem_name, conc))  Stop when conc first dips below
                                        ("conc_above", (chem_name, conc))  Stop when conc first rises above

                                        TODO: add more options, such as
                                        ("before_time", t)                  Stop just before the given target time
                                        ("after_time", t)                   Stop just after the given target time
                                        ("equilibrium", tolerance)          Stop when equilibrium reached

        :param n_steps:         The desired number of steps

        :param max_steps:       (OPTIONAL) Max numbers of steps; if reached, it'll terminate regardless of any other criteria

        :param snapshots:       (OPTIONAL) Dict that may contain any the following keys:
                                        -"frequency" (default 1)
                                        -"show_intermediates" (default True)
                                        -"species" (default None, meaning all species)
                                        -"initial_caption" (default blank)
                                        -"final_caption" (default blank)
                                    If provided, take a system snapshot after running a multiple
                                    of "frequency" reaction steps (default 1, i.e. at every step.)
                                    EXAMPLE: snapshots={"frequency": 2, "species": ["A", "H"]}

        :param silent:              If True, less output is generated

        :param variable_steps:      If True, the steps sizes will get automatically adjusted, based on thresholds
        :param explain_variable_steps:  If not None, a brief explanation is printed about how the variable step sizes were chosen,
                                            when the System time inside that range;
                                            only applicable if variable_steps is True

        :return:                None.   The object attributes self.system and self.system_time get updated
        """

        if reaction_duration and not duration:
            print("single_compartment_react(): the argument `reaction_duration` is deprecated; use `duration` instead")
            duration = reaction_duration    # For backward compatibility.   TODO: phase out

        # Validation
        assert self.system is not None, "UniformCompartment.single_compartment_react(): " \
                                        "the concentration values of the various chemicals must be set first"

        if variable_steps and (n_steps is not None):
            raise Exception("UniformCompartment.single_compartment_react(): if `variable_steps` is True, cannot specify `n_steps` "
                            "(because the number of steps will vary); specify `duration` or `target_end_time` instead")

        if stop is not None:
            assert type(stop) == tuple and len(stop) == 2, \
                f"UniformCompartment.single_compartment_react(): the argument `stop`, if passed, must be a pair of values"


        """
        Determine all the various time parameters that were not explicitly provided
        """

        if stop is not None:
            assert initial_step > 0, \
                "single_compartment_react(): when using the `stop` argument, an `initial_step` argument must be provided"
            assert max_steps is not None, \
                "single_compartment_react(): when using the `stop` argument, a `max_steps` argument must be provided"
            time_step = initial_step


        else:
            if target_end_time is not None:
                if duration is not None:
                    raise Exception("single_compartment_react(): cannot provide values for BOTH `target_end_time` and `duration`")
                else:
                    assert target_end_time > self.system_time, \
                        f"single_compartment_react(): `target_end_time` must be larger than the current System Time ({self.system_time})"
                    duration = target_end_time - self.system_time

            # Determine the time step,
            # as well as the required number of such steps
            time_step, n_steps = self.specify_steps(total_duration=duration,
                                                    time_step=initial_step,
                                                    n_steps=n_steps)
            # Note: if variable steps are requested then n_steps stops being particularly meaningful; it becomes a
            #       hypothetical value, in the (unlikely) event that the step size were never changed

            if target_end_time is None:
                if variable_steps:
                    target_end_time = self.system_time + duration
                else:
                    target_end_time = self.system_time + time_step * n_steps


        if snapshots:
            frequency = snapshots.get("frequency", 1)   # If not present, it will be 1
            species = snapshots.get("species")          # If not present, it will be None (meaning show all)
            first_snapshot = True
            if not snapshots.get("show_intermediates"):
                snapshots["show_intermediates"] = True
        else:
            snapshots = {"frequency": 1, "show_intermediates": True, "species": None, "first_snapshot": True}
            frequency = 1
            species = None
            first_snapshot = True


        if self.diagnostics:
            # Save up the current System State, with some extra info, as "diagnostic 'concentration' data"
            system_data = self.get_conc_dict(system_data=self.system)   # The current System State, as a dict
            self.save_diagnostic_conc_data(system_data)


        step_count = 0

        while True:
            # Check various criteria for termination
            if (max_steps is not None) and (step_count >= max_steps):
                print(f"single_compartment_react(): computation stopped because max # of steps ({max_steps}) reached")
                break       # We have reached the max allowable number of steps

            if (target_end_time is not None) and (self.system_time >= target_end_time):
                break       # The system time has reached the target endtime

            if (stop is not None):
                (termination_keyword, termination_parameter) = stop
                if termination_keyword == "conc_below":
                    chem_name, conc_threshold = termination_parameter
                    if self.get_chem_conc(chem_name) < conc_threshold:
                        break   # The concentration of the specified chemical has dropped the requested threshold
                elif termination_keyword == "conc_above":
                    chem_name, conc_threshold = termination_parameter
                    if self.get_chem_conc(chem_name) > conc_threshold:
                        break   # The concentration of the specified chemical has risen above the requested threshold

            if (not variable_steps) and (step_count == n_steps)\
                    and (target_end_time is not None) and np.allclose(self.system_time, target_end_time):
                break       # When dealing with fixed steps, catch scenarios where after performing n_steps,
                            #   the System Time is below the target_end_time because of roundoff error


            # ----------  CORE OPERATION OF MAIN LOOP  ----------
            delta_concentrations, step_actually_taken, recommended_next_step = \
                    self.reaction_step_common(delta_time=time_step,
                                              variable_steps=variable_steps, explain_variable_steps=explain_variable_steps,
                                              step_counter=step_count)

            # Update the System State
            self.previous_system = self.system.copy()
            self.system += delta_concentrations
            if min(self.system) < 0:    # Check for negative concentrations. TODO: redundant, since reaction_step_common() now does that
                print(f"+++++++++++ SYSTEM STATE ERROR: FAILED TO CATCH negative concentration upon advancing reactions from system time t={self.system_time:,.5g}")

            self.system_time += step_actually_taken

            # Preserve some of the data, as requested
            if snapshots and ((step_count+1)%frequency == 0):
                if first_snapshot and "initial_caption" in snapshots:
                    self.add_snapshot(species=species, caption=snapshots["initial_caption"])
                    first_snapshot = False
                else:
                    self.add_snapshot(species=species)

            step_count += 1

            if (n_steps is not None) and (step_count > 1000 * n_steps):  # Another approach to catch infinite loops
                raise Exception("single_compartment_react(): "
                                "the computation is taking a very large number of steps, probably from automatically trying to correct instability;"
                                " trying reducing the time_step")   # TODO: is the explanation correctly phrased?

            if self.diagnostics:
                # Save up the current System State, with some extra info, as "diagnostic 'concentration' data"
                system_data = self.get_conc_dict(system_data=self.system)   # The current System State, as a dict
                self.save_diagnostic_conc_data(system_data)

            if variable_steps:
                time_step = recommended_next_step   # Follow the recommendation of the ODE solver for the next time step to take
        # END while


        # Report whether extra steps were automatically added, as well as the total # taken
        n_steps_taken = step_count

        if n_steps is not None:
            extra_steps = n_steps_taken - n_steps
            if extra_steps > 0:
                if variable_steps:
                    print(f"Some steps were backtracked and re-done, "
                          f"to prevent negative concentrations or excessively large concentration changes")
                else:
                    print(f"The computation took {extra_steps} extra step(s) - "
                          f"automatically added to prevent negative concentrations")


        if not silent:
            print(f"{n_steps_taken} total step(s) taken")
            if variable_steps:
                print(f"Number of step re-do's because of negative concentrations: {self.number_neg_concs}")
                print(f"Number of step re-do's because of elective soft aborts: {self.number_soft_aborts}")
                print(f"Norm usage:", self.norm_usage)


        if snapshots and "final_caption" in snapshots:
            self.history.set_caption_last_snapshot(snapshots["final_caption"])



    def reaction_step_common(self, delta_time: float, conc_array=None,
                             variable_steps=False, explain_variable_steps=None, step_counter=1) -> (np.array, float, float):
        """
        This is the common entry point for both single-compartment reactions,
        and the reaction part of reaction-diffusions in 1D, 2D and 3D.

        "Compartments" may or may not correspond to the "bins" of the higher layers;
        the calling code might have opted to merge some bins into a single "compartment".

        Using the given concentration data for all the applicable species in a single compartment,
        do a single reaction time step for ALL the reactions -
        based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        Return the increment vector for all the chemical species concentrations in the compartment

        TODO: no longer pass conc_array .  Use the object variable self.system instead

        NOTES:  * the actual system concentrations are NOT changed
                * this method doesn't decide on step sizes - except in case of ("hard" or "soft") aborts, which are
                    followed by repeats with a smaller step.  Also, it makes suggestions
                    to the calling module about the next step to best take (whether as a result of an abort,
                    or for other considerations)

        :param delta_time:      The requested time duration of the reaction step
        :param conc_array:      [OPTIONAL]All initial concentrations at the start of the reaction step,
                                    as a Numpy array for ALL the chemical species, in their index order.
                                    If not provided, self.system is used instead
        :param variable_steps:  If True, the step sizes will get automatically adjusted with an adaptive algorithm
        :param explain_variable_steps:  If not None, a brief explanation is printed about how the variable step sizes were chosen,
                                            when the System time inside that range;
                                            only applicable if variable_steps is True
        :param step_counter:

        :return:                The triplet:
                                    1) increment vector for the concentrations of ALL the chemical species,
                                        in their index order, as a Numpy array
                                        EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):
                                            array([7. , -21.])  TODO: is this really necessary?  Maybe update self.system here?
                                    2) time step size actually taken - which might be smaller than the requested one
                                        because of reducing the step to avoid negative-concentration errors
                                    3) recommended_next_step : a suggestions to the calling module
                                       about the next step to best take
        """

        if conc_array is not None:
            self.system = conc_array    # For historical reasons, as a convenience to Bio1D, etc.

        # Validate arguments
        assert self.system is not None, "UniformCompartment.reaction_step_common(): " \
                                        "the concentration values of the various chemicals must be set first"


        if explain_variable_steps is not None:
            assert (type(explain_variable_steps) == list) and (len(explain_variable_steps) == 2), \
                "reaction_step_common(): the argument `explain_variable_steps`, if provided, must be a pair of numbers [t_start, t_end]"


        #print(f"************ At SYSTEM TIME: {self.system_time:,.4g}, calling reaction_step_common() with:")
        #print(f"             delta_time={delta_time}, system={self.system}, ")


        recommended_next_step = delta_time     # Baseline value; no reason yet to suggest a change in step size


        delta_concentrations = None

        SMALLEST_VALUE_TO_TRY = delta_time / 2000.       # Used to prevent infinite loops

        normal_exit = False

        while delta_time > SMALLEST_VALUE_TO_TRY:       # TODO: consider moving the inside of the WHILE loop into a separate function
            try:    # We want to catch Exceptions that can arise from excessively large time steps
                    #   that lead to negative concentrations or violation of user-set thresholds ("HARD" or "SOFT" aborts)
                (delta_concentrations, recommended_next_step) = \
                        self._attempt_reaction_step(delta_time, variable_steps, explain_variable_steps, step_counter)
                normal_exit = True

                break       # IMPORTANT: this is needed because, in the absence of errors, we need to go thru the WHILE loop only once!


            # CATCH any 'ExcessiveTimeStepHard' exception raised in the loop  (i.e. a HARD ABORT)
            except ExcessiveTimeStepHard as ex:
                # Single reactions steps can fail with this error condition if the attempted time step was too large,
                # under the following scenarios:
                #       1. negative concentrations from any one reaction - caught by  validate_increment()
                #       2. negative concentration from the combined effect of multiple reactions - caught in this function
                #print("*** CAUGHT a HARD ABORT")
                self.number_neg_concs += 1
                if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                    print(ex)

                delta_time *= self.step_factors["error"]       # Reduce the excessive time step by a pre-set factor
                recommended_next_step = delta_time
                # At this point, the loop will generally try the simulation again, with a smaller step (a revised delta_time)


            # CATCH any 'ExcessiveTimeStepSoft' exception raised in the loop  (i.e. a SOFT ABORT)
            except ExcessiveTimeStepSoft as ex:
                # Single reactions steps can fail with this error condition if the attempted time step was too large,
                # under the following scenario:
                #       3. excessive norm(s) measures in the overall step - caught in this function (currently only checked in case of the variable-steps option)
                #print("*** CAUGHT a soft ABORT")
                self.number_soft_aborts += 1
                if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                    print(ex)
                delta_time *= self.step_factors["abort"]       # Reduce the excessive time step by a pre-set factor
                recommended_next_step = delta_time
                # At this point, the loop will generally try the simulation again, with a smaller step (a revised delta_time)

        # END while


        if not normal_exit:         # i.e., if no reaction simulation took place in the WHILE loop, above
            raise Exception(f"reaction_step_common(): unable to complete the reaction step.  "
                            f"In spite of numerous automated reductions of the time step, it continues to lead to concentration changes that are considered excessive; "
                            f"consider reducing the original time step, and/or increasing the 'abort' thresholds with set_thresholds(). Current values: {self.thresholds}")

        # If we get thus far, it's the normal exit of the reaction step


        return  (delta_concentrations, delta_time, recommended_next_step)     # TODO: consider returning tentative_updated_system , since we already computed it




    def _attempt_reaction_step(self, delta_time, variable_steps, explain_variable_steps, step_counter) -> (np.array, float):
        """
        Attempt to perform the core reaction step, and then raise an Exception if it needs to be aborted,
        based on various criteria.
        If variable_steps is True, determine a new value for the "recommended next step"

        :param delta_time:              The requested time duration of the reaction step
        :param variable_steps:          If True, the step sizes will get automatically adjusted with an adaptive algorithm
        :param explain_variable_steps:  If not None, a brief explanation is printed about how the variable step sizes were chosen,
                                            when the System time inside that range;
                                            only applicable if variable_steps is True
        :param step_counter:            A pair with a time range inside which to show in the explanations about the variable step sizes;
                                            only applicable if explain_variable_steps is True

        :return:                The pair (delta_concentrations, recommended_next_step)
        """
        # *****  CORE OPERATION  *****
        delta_concentrations = self._reaction_elemental_step(delta_time=delta_time, rxn_list=None)


        if self.diagnostics:
            diagnostic_data_snapshot = self._delta_conc_dict(delta_concentrations)      # A dict


        recommended_next_step = delta_time       # Baseline value; no reason yet to suggest a change in step size

        if variable_steps:
            decision_data = self.adjust_timestep(delta_conc=delta_concentrations, baseline_conc=self.system, prev_conc=self.previous_system)
            step_factor = decision_data['step_factor']
            action = decision_data['action']
            all_norms = decision_data['norms']
            applicable_norms = decision_data['applicable_norms']

            if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                if action == "abort":
                    step_status = "aborted"
                else:
                    step_status = "completed"

                print(f"\n(STEP {step_counter} {step_status}) SYSTEM TIME {self.system_time:.5g} : Examining Conc. changes "
                      f"due to tentative Δt={delta_time:.5g} ...")
                print("    Previous: ", self.previous_system)
                print("    Baseline: ", self.system)
                print("    Deltas:   ", delta_concentrations)

                if len(self.chem_data.active_chemicals) < self.chem_data.number_of_chemicals():
                    print(f"    Restricting adaptive time step analysis to {len(self.chem_data.active_chemicals)} "
                    f"chemicals only: {self.chem_data.names_of_active_chemicals()} , with indexes: {self.chem_data.active_chemicals}")

                print("    Norms:    ", all_norms)
                print("    Thresholds:    ")
                for rule in self.thresholds:
                    print(self.display_thresholds(rule, all_norms.get(rule['norm'])))

                if action != "stay":    # The step is trivially 1 when the action is "stay"
                    print("    Step Factors:    ", self.step_factors)

                print(f"    => Action: '{action.upper()}'  (with step size factor of {step_factor})")


            # Abort the current step if some rate of change is deemed excessive.
            # TODO: maybe ALWAYS check this, regardless of variable-steps option
            if action == "abort":       # NOTE: this is a "strategic" abort, not a hard one from error
                msg =   f"* INFO: the tentative time step ({delta_time:.5g}) " \
                        f"leads to a value of {applicable_norms} > its ABORT threshold:\n" \
                        f"      -> will backtrack, and re-do step with a SMALLER Δt, x{step_factor:.5g} (now set to {delta_time * step_factor:.5g}) " \
                        f"[Step started at t={self.system_time:.5g}, and will rewind there]"
                #print("WARNING: ", msg)
                if self.diagnostics:
                    # Expand the dict diagnostic_data_snapshot
                    diagnostic_data_snapshot['norm_A'] = all_norms.get('norm_A')
                    diagnostic_data_snapshot['norm_B'] = all_norms.get('norm_B')        # TODO: add 'norm_C'
                    diagnostic_data_snapshot['action'] = "ABORT"
                    diagnostic_data_snapshot['step_factor'] = step_factor
                    diagnostic_data_snapshot['time_step'] = delta_time
                    self.save_diagnostic_decisions_data(data=diagnostic_data_snapshot, caption="excessive norm value(s)")

                    # Make a note of the abort action in all the reaction-specific diagnostics
                    self.comment_diagnostic_rxn_data("aborted: excessive norm value(s)")

                raise ExcessiveTimeStepSoft(msg)    # ABORT THE CURRENT STEP


            # Put together a recommendation to the higher-level functions, about the next best step size
            recommended_next_step = delta_time * step_factor


            if self.diagnostics:
                # Expand the dict diagnostic_data_snapshot
                diagnostic_data_snapshot['norm_A'] = all_norms.get('norm_A')
                diagnostic_data_snapshot['norm_B'] = all_norms.get('norm_B')
                diagnostic_data_snapshot['action'] = f"OK ({action})"
                diagnostic_data_snapshot['step_factor'] = step_factor
                diagnostic_data_snapshot['time_step'] = delta_time


            if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                #msg =   f"the tentative time step ({delta_time:.5g}) results in norm values that lead to the following:\n"
                msg = "    "

                if step_factor > 1:         # "INCREASE
                    msg +=  f"ACTION: COMPLETE STEP NORMALLY and MAKE THE INTERVAL LARGER, " \
                            f"multiplied by {step_factor} (set to {recommended_next_step:.5g}) at the next round, because all norms are low"
                elif step_factor < 1:     # "DECREASE"
                    msg +=  f"ACTION: COMPLETE STEP NORMALLY and MAKE THE INTERVAL SMALLER, " \
                            f"multiplied by {step_factor} (set to {recommended_next_step:.5g}) at the next round, because {applicable_norms} is high"
                else:   # "STAY THE COURSE"
                    msg +=  f"ACTION: COMPLETE NORMALLY - we're inside the target range of all norms.  No change to step size."

                msg += f"\n    [The current step started at System Time: {self.system_time:.5g}, and will continue to {self.system_time + delta_time:.5g}]"
                print(msg)
        # END if variable_steps


        if self.diagnostics:
            self.save_diagnostic_decisions_data(data=diagnostic_data_snapshot)


        # Check whether the COMBINED delta_concentrations will make any conc negative;
        # if so, raised an "ExcessiveTimeStepHard" exception (a custom exception)
        tentative_updated_system = self.system + delta_concentrations
        if min(tentative_updated_system) < 0:
            if explain_variable_steps and (explain_variable_steps[0] <= self.system_time <= explain_variable_steps[1]):
                print(f"*** CAUTION: negative concentration resulting from the combined effect of all reactions, "
                      f"upon advancing reactions from system time t={self.system_time:,.5g}\n"
                      f"         It'll be AUTOMATICALLY CORRECTED with a reduction in time step size")

            # A type of HARD ABORT is detected (a negative concentration resulting from the combined effect of all reactions)
            if self.diagnostics:
                self.save_diagnostic_decisions_data(data={"action": "ABORT",
                                                          "step_factor": self.step_factors["error"],
                                                          "caption": "neg. conc. from combined effect of all rxns",
                                                          "time_step": delta_time})
                for rxn_index in range(self.chem_data.number_of_reactions()):
                    self.save_diagnostic_rxn_data(rxn_index=rxn_index, time_step=delta_time,
                                                  increment_vector_single_rxn=None,
                                                  caption=f"aborted: neg. conc. from combined multiple rxns")

            raise ExcessiveTimeStepHard(f"INFO: the tentative time step ({delta_time:.5g}) "
                                        f"leads to a NEGATIVE concentration of one of the chemicals: "
                                        f"\n      -> will backtrack, and re-do step with a SMALLER delta time, "
                                        f"multiplied by {self.step_factors['error']} (set to {delta_time * self.step_factors['error']:.5g}) "
                                        f"[Step started at t={self.system_time:.5g}, and will rewind there.  Baseline values: {self.system} ; delta conc's: {delta_concentrations}]")


        return  (delta_concentrations, recommended_next_step)       # Maybe also return tentative_updated_system



    def adjust_timestep(self, delta_conc: np.array, baseline_conc=None, prev_conc=None) -> dict:
        """
        Computes some measures of the change of concentrations, from the last step, in the context of the
        baseline initial concentrations of that same step, and the concentrations in the step before that.
        Based on the magnitude of the measures, propose a course of action about what to do for the next step.

        :param delta_conc:      A numpy array of changes in concentrations for the chemicals of interest,
                                    across a simulation time step (typically, the current step a run in progress)
        :param baseline_conc:   A numpy array of baseline concentration values for those same chemicals,
                                    prior to the above change, at the start of a simulation time step
        :param prev_conc:       A numpy array of concentration values for those same chemicals,
                                    in the step prior to the current one (i.e. an "archive" value)
        :return:                A dict:
                                    "action"           - String with the name of the computed recommended action:
                                                            either "low", "stay", "high" or "abort"
                                    "step_factor"      - A factor by which to multiply the time step at the next iteration round;
                                                            if no change is deemed necessary, 1
                                    "norms"            - A dict of all the computed norm name/values (any of the norms, except norm_A,
                                                            may be missing)
                                    "applicable_norms" - The name of the norm that triggered the decision; if all norms were involved,
                                                            it will be "ALL"
        """
        # TODO: move to separate static class
        # TODO: pass the value of self.chem_data.indexes_of_active_chemicals() as argument
        #       to this method

        n_chems = self.chem_data.number_of_chemicals()

        if baseline_conc is not None:
            assert len(baseline_conc) == len(delta_conc), \
                f"adjust_speed(): the number of entries in the passed array `delta_conc` ({len(delta_conc)}) " \
                f"does not match the number of entries in the passed array `baseline_conc` ({len(baseline_conc)})"

        assert n_chems == len(delta_conc), \
            f"adjust_speed(): the number of entries in the passed array `delta_conc` ({len(delta_conc)}) " \
            f"does not match the number of registered chemicals ({n_chems})"


        # If some chemicals are not dynamically involved in the reactions
        # (i.e. if they don't occur in any reaction, or occur as enzyme),
        # restrict our consideration to only the dynamically involved ones
        # CAUTION: the concept of "active chemical" might change in future versions, where only SOME of
        #          the reactions are simulated
        if len(self.chem_data.names_of_enzymes()) > 0:
            delta_conc = delta_conc[self.chem_data.indexes_of_active_chemicals()]
            #print(f"\nadjust_speed(): restricting adaptive time step analysis to {n_chems} chemicals only; their delta_conc is {delta_conc}")
            if baseline_conc is not None:
                baseline_conc = baseline_conc[self.chem_data.indexes_of_active_chemicals()]
                prev_conc = prev_conc[self.chem_data.indexes_of_active_chemicals()]
        # Note: setting delta_conc, etc, only affects local variables, and won't mess up the arrays passed as arguments

        all_norms = {}

        all_small = True            # Provisional answer to the question: "do ALL the rule yield a 'low'?"
                                    # Any rule failing to yield a 'low' will flip this status

        high_seen_at = []           # List of rule names at which a "high" is encountered, if applicable

        for rule in self.thresholds:
            norm_name = rule["norm"]

            if norm_name == "norm_A":
                result = RxnDynamics.norm_A(delta_conc)
            elif norm_name == "norm_B":
                result = RxnDynamics.norm_B(baseline_conc, delta_conc)
            elif norm_name == "norm_C":
                result = RxnDynamics.norm_C(prev_conc, baseline_conc, delta_conc)
            else:
                result = RxnDynamics.norm_D(prev_conc, baseline_conc, delta_conc)

            all_norms[norm_name] = result

            if ("abort" in rule) and (result > rule["abort"]):
                # If any rules declares an abort, no need to proceed further: it's an abort
                self.norm_usage[norm_name] += 1
                return {"action": "abort", "step_factor": self.step_factors["abort"], "norms": all_norms, "applicable_norms": [norm_name]}

            if ("high" in rule) and (result > rule["high"]):
                # If any rules declares a "high", still need to consider the other rules - in case any of them over-rides
                # the "high" with an "abort"
                high_seen_at.append(norm_name)
                all_small = False           # No longer the case of all 'low` yields

            if all_small and ("low" in rule) and (result > rule["low"]):
                all_small = False           # No longer the case of all 'low` yields
        # END for


        if high_seen_at:
            for n in high_seen_at:
                self.norm_usage[n] += 1
            return {"action": "high", "step_factor": self.step_factors["downshift"], "norms": all_norms, "applicable_norms": high_seen_at}


        if all_small:
            for i in self.norm_usage:
                self.norm_usage[i] += 1     # All the norms were used
            return {"action": "low", "step_factor": self.step_factors["upshift"], "norms": all_norms, "applicable_norms": "ALL"}


        # If we get thus far, none of the rules were found applicable
        return {"action": "stay", "step_factor": 1, "norms": all_norms, "applicable_norms": "ALL"}



    def display_thresholds(self, rule :dict, value) -> str:
        """
        Examine how the specified value fits
        relatively to the 'low', 'high', and 'abort' stored in the given rule

        :param rule:    A dict that must contain the key 'norm',
                            and may contain the keys: 'low', 'high', and 'abort'
                            (referring to 3 increasingly-high thresholds)
        :param value:   Either None, or a number to compare to the thresholds
        :return:        A string that visually highlights the relative position of the value
                            relatively to the given thresholds
        """
        s = f"                   {rule['norm']} : "     # Name of the norm being used

        if value is None:
            return s + " (skipped; not needed)"

        # Extract the 3 thresholds (some might be missing)
        low = rule.get('low')
        high = rule.get('high')
        abort = rule.get('abort')

        if low is not None and value <= low:
            # The value is below the `low` threshold
            return f"{s}(VALUE {value:.5g}) | low {low} | high {high} | abort {abort}"

        if high is not None and value < high:
            # The value is between the `low` and `high` thresholds
            return f"{s}low {low} | (VALUE {value:.5g}) | high {high} | abort {abort}"

        if abort is not None and value < abort:
            # The value is above the `high` threshold
            return f"{s}low {low} | high {high} | (VALUE {value:.5g}) | abort {abort}"

        # If we get thus far, the value is above the `abort` threshold
        return f"{s}low {low} | high {high} | abort {abort} | (VALUE {value:.5g})"



    def _reaction_elemental_step(self, delta_time: float, rxn_list=None) -> np.array:
        """
        Using the system concentration data of ALL the chemical species,
        do the specified SINGLE TIME STEP for ONLY the requested reactions (by default all).

        All computations are based on the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions (in "forward Euler" approach.)

        Return the Numpy increment vector for ALL the chemical species concentrations, in their index order
        (whether involved in these reactions or not)

        NOTES:  - the actual System Concentrations and the System Time (stored in object variables) are NOT changed
                - if any of the concentrations go negative, an Exception is raised

        :param delta_time:      The time duration of this individual reaction step - assumed to be small enough that the
                                    concentration won't vary significantly during this span.
        :param rxn_list:        OPTIONAL list of reactions (specified by their indices) to include in this simulation step ;
                                    EXAMPLE: [1, 3, 7]
                                    If None, do all the reactions

        :return:            The increment vector for the concentrations of ALL the chemical species
                                (whether involved in the reactions or not),
                                as a Numpy array for all the chemical species, in their index order
                            EXAMPLE (for a single-reaction reactant and product with a 3:1 stoichiometry):   array([7. , -21.])
        """

        # The increment vector is cumulative for ALL the requested reactions
        increment_vector = np.zeros(self.chem_data.number_of_chemicals(), dtype=float)       # One element per chemical species

        # Compute the forward and reverse "conversions" of all the applicable reactions
        delta_dict = self.compute_all_reaction_deltas(delta_time=delta_time, rxn_list=rxn_list)
        if 3 in self.verbose_list:
            print(f"      delta_list: {delta_dict}")


        if rxn_list is None:    # Meaning ALL (active) reactions
            # A list of the reaction indices of all the active reactions
            rxn_list = self.chem_data.active_reaction_indices()


        number_chemicals = self.chem_data.number_of_chemicals()

        # For each applicable reaction, find the needed adjustments ("deltas")
        #   to the concentrations of the reactants and products,
        #   based on the forward and reverse rates of the reaction
        for rxn_index in rxn_list:      # Consider each reaction in turn
            # TODO: maybe switch to a call to the experimental _reaction_elemental_step_SINGLE_REACTION()

            # One element per chemical species; notice that this array is being RESET for EACH REACTION
            # TODO: instead of using a potentially very large array (mostly of zeros) for each rxn, consider a dict instead
            #       then combine then at end
            increment_vector_single_rxn = np.zeros(number_chemicals, dtype=float)

            rxn = self.chem_data.get_reaction(rxn_index)
            reactants = rxn.extract_reactants()
            products = rxn.extract_products()

            """
            Determine the concentration adjustments as a result of this reaction step, 
            for the reaction being considered
            """

            # The reactants DECREASE based on the quantity (forward reaction - reverse reaction)
            for r in reactants:
                # Unpack data from the reactant r
                species_name = rxn.extract_species_name(r)
                species_index = self.chem_data.get_index(species_name)
                if species_name == rxn.enzyme:
                    #print(f"*** SKIPPING reactant ENZYME {species_index} in reaction {rxn_index}")
                    continue    # Skip if r is an enzyme for this reaction

                stoichiometry = rxn.extract_stoichiometry(r)

                delta_conc = stoichiometry * (- delta_dict[rxn_index])  # Increment to this reactant from the reaction being considered
                # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
                # Note: not enough to detect conc going negative from combined changes from multiple reactions!
                #       Further testing done upstream
                self.validate_increment(delta_conc=delta_conc, baseline_conc=self.system[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                increment_vector_single_rxn[species_index] += delta_conc


            # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
            for p in products:
                # Unpack data from the reactant r
                species_name = rxn.extract_species_name(p)
                species_index = self.chem_data.get_index(species_name)
                if species_name == rxn.enzyme:
                    #print(f"*** SKIPPING product ENZYME {species_index} in reaction {rxn_index}")
                    continue    # Skip if p is an enzyme for this reaction

                stoichiometry = rxn.extract_stoichiometry(p)

                delta_conc = stoichiometry * delta_dict[rxn_index]  # Increment to this reaction product from the reaction being considered
                # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
                # Note: not enough to detect conc going negative from combined changes from multiple reactions!
                #       Further testing done upstream
                self.validate_increment(delta_conc=delta_conc, baseline_conc=self.system[species_index],
                                        rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

                increment_vector_single_rxn[species_index] += delta_conc

            # Macro-molecule related part, if applicable    TODO: implement
            if (self.macro_system_state != {}) and (rxn.macro_enzyme is not None):
                print(f"Making adjustments for macro-molecule catalysis for reaction # {rxn_index}")
                print(f"    Macromolecule: {rxn.macro_enzyme[0]}, at site # {rxn.macro_enzyme[1]}")
                print(f"    Site occupancy at the beginning of the time step:")
                print(f"    Macromolecule count:")


            increment_vector += increment_vector_single_rxn     # Accumulate the increment vector across ALL reactions
                                                                # TODO: consider using a small dict in lieu of increment_vector_single_rxn

            if self.diagnostics:
                self.save_diagnostic_rxn_data(rxn_index=rxn_index, time_step=delta_time,
                                              increment_vector_single_rxn=increment_vector_single_rxn)
        # END for (over rxn_list)

        return increment_vector



    def _reaction_elemental_step_SINGLE_REACTION(self, delta_time: float, increment_vector,
                                                 rxn_index :int,
                                                 delta_dict):     # TODO: EXPERIMENTAL; not yet in use
        """
        :param delta_time:
        :param increment_vector:
        :param rxn_index:       The integer index (0-based) to identify the reaction of interest
        :param delta_dict:
        :return:
        """
        if 2 in self.verbose_list:
            print(f"      Determining the conc.'s changes as a result of rxn # {rxn_index}")

        # NOTE: instead of using a potentially very large array (mostly of zeros) for each rxn,
        #       we use a dict instead; then combine all at end
        #increment_vector_single_rxn = np.zeros(self.reaction_data.number_of_chemicals(), dtype=float)
        increment_dict_single_rxn = {}      # The key will be the elements in the reaction

        # TODO: turn into a more efficient single step, as as:
        #(reactants, products) = cls.all_reactions.unpack_terms(rxn_index)
        reactants = self.chem_data.get_reactants(rxn_index)
        products = self.chem_data.get_products(rxn_index)


        """
        Determine the concentration adjustments as a result of this reaction step, 
        for the reaction being considered
        """

        # The reactants DECREASE based on the quantity (forward reaction - reverse reaction)
        for r in reactants:
            stoichiometry, species_index, order = r                 # Unpack
            delta_conc = stoichiometry * (- delta_dict[rxn_index])  # Increment to this reactant from the reaction being considered
            # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
            # Note: not enough to detect conc going negative from combined changes from multiple reactions!  Further testing done upstream
            # TODO: it might be more efficient to check this in bulk than with many function calls
            self.validate_increment(delta_conc=delta_conc, baseline_conc=self.system[species_index],
                                    rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

            increment_dict_single_rxn[species_index] = increment_dict_single_rxn.get(species_index, 0) + delta_conc
            #increment_vector_single_rxn[species_index] += delta_conc


        # The reaction products INCREASE based on the quantity (forward reaction - reverse reaction)
        for p in products:
            stoichiometry, species_index, order = p             # Unpack
            delta_conc = stoichiometry * delta_dict[rxn_index]  # Increment to this reaction product from the reaction being considered
            # Do a validation check to avoid negative concentrations; an Exception will get raised if that's the case
            # Note: not enough to detect conc going negative from combined changes from multiple reactions!  Further testing done upstream
            # TODO: it might be more efficient to check this in bulk than with many function calls
            self.validate_increment(delta_conc=delta_conc, baseline_conc=self.system[species_index],
                                    rxn_index=rxn_index, species_index=species_index, delta_time=delta_time)

            increment_dict_single_rxn[species_index] = increment_dict_single_rxn.get(species_index, 0) + delta_conc
            #increment_vector_single_rxn[species_index] += delta_conc


        for k, v in increment_dict_single_rxn.items():
            increment_vector[k] += v                # Accumulate the increment vector across all the increments from this reaction
        #increment_vector += increment_vector_single_rxn     # Accumulate the increment vector across all reactions


        # If requested to evaluate the reaction "speeds"
        # TODO: need a version of examine_increment_array() and save_diagnostic_data() that accept increment_dict_single_rxn
        #       instead of increment_vector_single_rxn [also, maybe move one or both of them to calling function]
        '''  
        if self.diagnostics:
            self.save_diagnostic_data(delta_time, increment_vector_single_rxn, rxn_index)
        '''



    def validate_increment(self,  delta_conc, baseline_conc: float,
                           rxn_index :int, species_index: int, delta_time) -> None:
        """
        Examine the requested concentration change given by delta_conc
        (typically, as computed by an ODE solver),
        relative to the baseline (pre-reaction) value baseline_conc,
        for the given SINGLE chemical species and SINGLE reaction.

        If the concentration change would render the concentration negative,
        raise an Exception (of custom type "ExcessiveTimeStepHard")

        :param delta_conc:              The change in concentration computed by the ode solver
                                            (for the specified chemical, in the given reaction)
        :param baseline_conc:           The initial concentration

        [The remaining arguments are ONLY USED for error printing]
        :param rxn_index:               The index (0-based) to identify the reaction of interest (ONLY USED for error printing)
        :param species_index:           The index (0-based) to identify the chemical species of interest (ONLY USED for error printing)
        :param delta_time:              The time duration of the reaction step (ONLY USED for error printing)

        :return:                        None (an Exception is raised if a negative concentration is detected)
        """
        if (baseline_conc + delta_conc) < 0:
            #print(f"\n*** CAUTION: negative concentration in chemical `{self.chem_data.get_name(species_index)}` in step starting at t={self.system_time:.4g}.  "
            #      f"It will be AUTOMATICALLY CORRECTED with a reduction in time step size, as follows:")

            # A type of HARD ABORT is detected (a single reaction that, by itself, would lead to a negative concentration;
            #   while it's possible that other coupled reactions might counterbalance this - nonetheless, it's taken as a sign of excessive step size)
            if self.diagnostics:
                self.save_diagnostic_decisions_data(data={"action": "ABORT",
                                                          "step_factor": self.step_factors['error'],
                                                          "caption": f"neg. conc. in {self.chem_data.get_name(species_index)} from rxn # {rxn_index}",
                                                          "time_step": delta_time})
                self.save_diagnostic_rxn_data(rxn_index=rxn_index, time_step=delta_time,
                                              increment_vector_single_rxn=None,
                                              caption=f"aborted: neg. conc. in `{self.chem_data.get_name(species_index)}`")

            raise ExcessiveTimeStepHard(f"    INFO: the tentative time step ({delta_time:.5g}) "
                                    f"leads to a NEGATIVE concentration of `{self.chem_data.get_name(species_index)}` "
                                    f"from reaction {self.chem_data.single_reaction_describe(rxn_index=rxn_index, concise=True)} (rxn # {rxn_index}): "
                                    f"\n      Baseline value: {baseline_conc:.5g} ; delta conc: {delta_conc:.5g}"
                                    f"\n      -> will backtrack, and re-do step with a SMALLER delta time, "
                                    f"multiplied by {self.step_factors['error']} (set to {delta_time * self.step_factors['error']:.5g}) "
                                    f"[Step started at t={self.system_time:.5g}, and will rewind there]")



    def compute_all_reaction_deltas(self, delta_time: float, rxn_list=None) -> dict:
        """
        For an explanation of the "reaction delta", see compute_reaction_delta().
        Compute the "reaction delta" for all the specified reaction (by default, all.)
        Return a list with an entry for each reaction, in their index order.

        For background info: https://life123.science/reactions

        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                                concentration won't vary significantly during this span
        :param rxn_list:    OPTIONAL list of reactions (specified by their integer index);
                                if None, do all the reactions.  EXAMPLE: [1, 3, 7]

        :return:            A dict of the differences between forward and reverse "conversions" -
                                for explanation, see compute_reaction_delta().
                                The dict is indexed by the reaction number, and contains as many entries as the
                                number of reactions being investigated
        """
        delta_dict = {}

        if rxn_list is None:    # Meaning ALL reactions
            rxn_list = range(self.chem_data.number_of_reactions())

        # Process the requested reactions
        for i in rxn_list:      # Consider each reaction in turn
            rxn = self.chem_data.get_reaction(i)
            delta = delta_time * self.compute_reaction_delta_rate(rxn=rxn)
            delta_dict[i] = delta

        return delta_dict



    def compute_reaction_delta_rate(self, rxn) -> float:
        """
        For the SINGLE given reaction, and the current concentrations of chemicals in the system,
        compute the reaction's "forward rate" minus its "reverse rate",
        as defined in https://life123.science/reactions

        :param rxn:         An object of type "Reaction"
        :return:            The differences between the reaction's forward and reverse rates
        """
        reactants, products, fwd_rate_constant, rev_rate_constant = rxn.unpack_for_dynamics()

        forward_rate = fwd_rate_constant
        for r in reactants:
            # Unpack data from the reactant r
            species_name = rxn.extract_species_name(r)
            order = rxn.extract_rxn_order(r)
            species_index = self.chem_data.get_index(species_name)
            conc = self.system[species_index]
            #assert conc is not None, \
            #   f"UniformCompartment.compute_reaction_delta_rate(): lacking the value for the concentration of the chemical species `{self.reaction_data.get_name(species_index)}`"
            forward_rate *= conc ** order      # Raise to power

        reverse_rate = rev_rate_constant
        for p in products:
            # Unpack data from the reaction product p
            species_name = rxn.extract_species_name(p)
            order = rxn.extract_rxn_order(p)
            species_index = self.chem_data.get_index(species_name)
            conc = self.system[species_index]
            #assert conc is not None, \
            #   f"UniformCompartment.compute_reaction_delta_rate(): lacking the concentration value for the species `{self.reaction_data.get_name(species_index)}`"
            reverse_rate *= conc ** order     # Raise to power

        return forward_rate - reverse_rate



    def solve_exactly(self, rxn_index :int, A0 :float, B0 :float, t_arr) -> (np.array, np.array):
        """
        Return the exact solution of the reaction with the requested index,
        PROVIDED that it is a 1st Order Reaction of the type A <=> B.

        Use the given initial conditions,
        and return the solutions sampled at the specified times.

        For details, see https://life123.science/reactions

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest
        :param A0:
        :param B0:
        :param t_arr:       A Numpy array with the desired times at which the solutions are desired
        :return:            A pair of Numpy arrays
        """
        rxn = self.chem_data.get_reaction(rxn_index)
        reactants, products, kF, kR = rxn.unpack_for_dynamics()

        assert len(reactants) == 1, "Currently only works for `A <-> B` reactions"
        assert len(products) == 1, "Currently only works for `A <-> B` reactions"
        assert rxn.extract_stoichiometry(reactants[0]) == 1, \
            "Currently only works for `A <-> B` reactions"
        assert rxn.extract_stoichiometry(products[0]) == 1, \
            "Currently only works for `A <-> B` reactions"
        # TODO: should also verify the reaction orders to be 1

        #A = rxn.extract_species_index(reactants[0])     # The index of the "A" chemical in the expected reaction
        #B = rxn.extract_species_index(products[0])      # The index of the "B" chemical in the expected reaction

        TOT = A0 + B0
        #print(kF, kR, A0, TOT)

        return self._exact_solution(kF, kR, A0, TOT, t_arr)


        
    def _exact_solution(self, kF, kR, A0, TOT, t_arr) -> (np.array, np.array):
        """
        Return the exact solution of the 1st Order Reaction A <=> B,
        with the specified parameters, 
        sampled at the given times.
        
        For details, see https://life123.science/reactions

        :param kF:
        :param kR:
        :param A0:
        :param TOT:
        :param t_arr:   A Numpy array with the desired times at which the solutions are desired
        :return:        A pair of Numpy arrays
        """
        # (A0 - (kR TOT) / (kF + kR)) Exp[-(kF + kR) t] + kR TOT / (kF + kR)
        A_arr = (A0 - (kR * TOT) / (kF + kR)) * np.exp(-(kF + kR) * t_arr) + kR * TOT / (kF + kR)
        B_arr = TOT - A_arr
        return (A_arr, B_arr)




    #####################################################################################################

    '''                                 ~  MACROMOLECULE DYNAMICS ~                                   '''

    def ________MACROMOLECULE_DYNAMICS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def set_macromolecules(self, data=None) -> None:
        """
        Specify the macromolecules, and their counts, to be included in the system.
        The fractional occupancy is set to 0 at all binding sites of all the specified macromolecules.
        Any previous data gets over-written.

        Note: to set a single fractional occupancy value, use set_occupancy()

        :param data:    A dict mapping macromolecule names to their counts
                            EXAMPLE:  {"M1": 1, "M2": 3, "M3": 1}
                        If any of the requested macromolecules isn't registered, an Exception will be raised
                        If data=None, then the set of registered macromolecules is used,
                            and all their counts are set to 1
        :return:        None.
                        The object variables self.macro_system and self.macro_system_state get set
        """
        if data is None:
            # Use the registered macromolecules, and set all counts to 1
            data = {}
            for mm in self.chem_data.get_macromolecules():
                data[mm] = 1        # EXAMPLE, after this operation: data = {"M1": 1}

        self.macro_system = data


        self.macro_system_state = {}    # Reset
        for mm in data.keys():          # For each macromolecule in our system
            binding_sites_and_ligands = self.chem_data.get_binding_sites_and_ligands(mm)     # EXAMPLE: {1: "A", 2: "C"}
            d = {}
            for (site_number, ligand) in binding_sites_and_ligands.items():
                d[site_number] = (ligand, 0.)           # All "binding occupancy fractions" are set to 0.

            self.macro_system_state[mm] = d



    def set_occupancy(self, macromolecule, site_number: int, fractional_occupancy: float) -> None:
        """
        Set the fractional occupancy at the given binding site of the specified macromolecule,
        using the requested value.
        If the specified macromolecule hasn't yet been added to the dynamical system state,
        automatically add it with count 1

        :param macromolecule:           Name of a previously-registered macromolecule
        :param site_number:             Integer to identify a binding site on the macromolecule
        :param fractional_occupancy:    A number between 0. and 1., inclusive
        :return:                        None
        """
        assert 0. <= fractional_occupancy <= 1., \
            f"set_occupancy(): the value for the fractional occupancy must be a number between 0. and 1. " \
            f"Value given: {fractional_occupancy}"

        ligand = self.chem_data.get_ligand_name(macromolecule=macromolecule, site_number=site_number)

        # If the specified macromolecule hasn't yet been added to the dynamic system, automatically add it
        # with count 1
        if self.macro_system_state == {}:
            self.set_macromolecules({macromolecule: 1})

        self.macro_system_state[macromolecule][site_number] = (ligand, fractional_occupancy)



    def get_occupancy(self, macromolecule, site_number) -> float:
        """
        Get the fractional occupancy at the given binding site of the specified macromolecule.

        :param macromolecule:           Name of a previously-registered macromolecule
        :param site_number:             Integer to identify a binding site on the macromolecule
        :return:                        A number between 0. and 1., representing the fractional occupancy
        """
        assert self.macro_system_state != {}, \
            "get_occupancy(): The system state for macromolecules has not been set yet;  " \
            "use set_macromolecules() or set_occupancy()"

        assert macromolecule in self.macro_system_state, \
            f"get_occupancy(): No occupancy data yet set for macromolecule `{macromolecule}`"

        assert site_number in self.macro_system_state[macromolecule], \
            f"get_occupancy(): No occupancy data yet set for site number {site_number} " \
            f"of macromolecule `{macromolecule}`"

        (ligand, fractional_occupancy) = self.macro_system_state[macromolecule][site_number]
        return fractional_occupancy



    def update_occupancy(self) -> None:
        """
        Update the fractional occupancy at all binding sites,
        based on the current system concentrations of the relevant ligands

        :return:    None
        """
        for mm in self.chem_data.get_macromolecules():
            # For each macromolecule
            d = self.chem_data.get_binding_sites_and_ligands(mm)    # EXAMPLE: {1: "A", 2: "C"}
            for (site_number, ligand) in d.items():
                aff_data = self.chem_data.get_binding_site_affinity(mm, site_number)
                conc = self.get_chem_conc(ligand)
                fractional_occupancy = self.sigmoid(conc=conc, Kd=aff_data.Kd)

                self.set_occupancy(macromolecule=mm, site_number=site_number, fractional_occupancy=fractional_occupancy)



    def sigmoid(self, conc: float, Kd: float) -> float:
        """
        Return an estimate of fractional occupancy (between 0 and 1)
        on a particular binding site on a particular macromolecule,
        from the concentration of the ligand (such as a Transcription Factor)
        and its affinity to that binding site.

        A sigmoid curve is expected.

        Based on fig. 3A of the 2019 paper "Low-Affinity Binding Sites and the
        Transcription Factor Specificity Paradox in Eukaryotes"
        (https://doi.org/10.1146/annurev-cellbio-100617-062719):

            - at extremely low concentration, the occupancy is 0
            - when the concentration is 10% of Kd, the occupancy is about 0.1
            - when the concentration matches Kd, the occupancy is 1/2 by definition
            - when the concentration is 10 times Kd, the occupancy is about 0.9
            - at concentrations beyond that, the occupancy saturates to 1.0

        :param conc:    Concentration of the ligand (such as a Transcription Factor), in microMolars
        :param Kd:      Binding-side Affinity, in microMolars
        :return:        Estimated binding-site fractional occupancy : a value between
                            0. (no occupancy at all during the previous time step) and 1. (continuous occupancy)
        """
        if conc == 0:
            conc = 1e-15     # To avoid taking log of 0

        return self.logistic(x = math.log10(conc), x0 = math.log10(Kd), k = 2.1972245)



    def logistic(self, x: float, x0 = 0., k = 1.) -> float:
        """
        Compute the value of the Logistic function, in the range (0, 1), at the given point
        See: https://en.wikipedia.org/wiki/Logistic_function

        :param x:
        :param x0:
        :param k:
        :return:    The value of the Logistic function at the given point x
        """
        return 1. / ( 1 + math.exp( -k * (x-x0) ) )




    #####################################################################################################

    '''                                  ~   FOR DIAGNOSTICS   ~                                      '''

    def ________FOR_DIAGNOSTICS________(DIVIDER):
        pass         # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def stoichiometry_checker(self, rxn_index :int, conc_arr_before: np.array, conc_arr_after: np.array) -> bool:
        """
        For the indicated reaction, investigate the change in the concentration of the involved chemicals,
        to ascertain whether the change is consistent with the reaction's stoichiometry.
        See https://life123.science/reactions

        IMPORTANT: this function is currently meant for simulations involving only 1 reaction (TODO: generalize)

        NOTE: the concentration changes in chemicals not involved in the specified reaction are ignored

        :param rxn_index:       The integer index (0-based) to identify the reaction of interest
        :param conc_arr_before: Numpy array with the concentrations of ALL the chemicals (whether involved
                                    in the reaction or not), in their index order, BEFORE the reaction step
        :param conc_arr_after:  Same as above, but after the reaction
                                TODO: maybe also accept a Panda's dataframe row

        :return:                True if the change in reactant/product concentrations is consistent with the
                                    reaction's stoichiometry, or False otherwise
        """
        assert self.chem_data.number_of_reactions() == 1, \
            f"stoichiometry_checker() currently only works for 1-reaction simulations, but " \
            f"{self.chem_data.number_of_reactions()} reactions are present."

        assert len(conc_arr_before) == self.chem_data.number_of_chemicals(), \
            f"stoichiometry_checker() : the argument `conc_arr_before` must be a Numpy array " \
            f"with as many elements as the number of registered chemicals ({self.chem_data.number_of_chemicals()}); " \
            f"instead, it has {len(conc_arr_before)}"

        assert len(conc_arr_after) == self.chem_data.number_of_chemicals(), \
            f"stoichiometry_checker() : the argument `conc_arr_after` must be a Numpy array " \
            f"with as many elements as the number of registered chemicals ({self.chem_data.number_of_chemicals()}); " \
            f"instead, it has {len(conc_arr_after)}"

        self.chem_data.assert_valid_rxn_index(rxn_index)

        return self._stoichiometry_checker_from_deltas(rxn_index = rxn_index,
                                                       delta_arr = conc_arr_after-conc_arr_before)



    def _stoichiometry_checker_from_deltas(self, rxn_index :int, delta_arr: np.array) -> bool:
        """
        Helper function.

        For the indicated reaction, investigate the change in the concentration of the involved chemicals,
        to ascertain whether the change is consistent with the reaction's stoichiometry.
        See https://life123.science/reactions

        IMPORTANT: this function is currently meant for simulations involving only 1 reaction (TODO: generalize)

        NOTE: the concentration changes in chemicals not involved in the specified reaction are ignored

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest
        :param delta_arr:   Numpy array of numbers, with the concentrations changes of ALL the chemicals (whether involved
                                in the reaction or not), in their index order,
                                as a result of JUST the reaction of interest
        :return:            True if the change in reactant/product concentrations is consistent with the
                                reaction's stoichiometry, or False otherwise
                                Note: if any of the elements of the passed Numpy array is NaN, then True is returned
                                      (because NaN values are indicative of aborted steps; can't invalidate the stoichiometry
                                      check because of that)
        """
        #print("delta_arr: ", delta_arr)

        if np.isnan(delta_arr).any():
            return True         # The presence of a NaN, anywhere in delta_arr, is indicative of an aborted step

        rxn = self.chem_data.get_reaction(rxn_index)
        reactants = self.chem_data.get_reactants(rxn_index)
        products = self.chem_data.get_products(rxn_index)

        # Pick (arbitrarily) the first reactant,
        # to establish a baseline change in concentration relative to its stoichiometric coefficient
        baseline_term = reactants[0]
        baseline_species_name = rxn.extract_species_name(baseline_term)
        baseline_species_index = self.chem_data.get_index(baseline_species_name)
        baseline_stoichiometry = rxn.extract_stoichiometry(baseline_term)
        baseline_ratio =  (delta_arr[baseline_species_index]) / baseline_stoichiometry
        #print("\nbaseline_ratio: ", baseline_ratio)

        for i, term in enumerate(reactants):
            if i != 0:
                species_name = rxn.extract_species_name(term)
                species_index = self.chem_data.get_index(species_name)
                stoichiometry = rxn.extract_stoichiometry(term)
                ratio =  (delta_arr[species_index]) / stoichiometry
                #print(f"ratio for `{self.chem_data.get_name(species)}`: {ratio}")
                if not np.allclose(ratio, baseline_ratio):
                    return False

        for term in products:
            species_name = rxn.extract_species_name(term)
            species_index = self.chem_data.get_index(species_name)
            stoichiometry = rxn.extract_stoichiometry(term)
            ratio =  - (delta_arr[species_index]) / stoichiometry     # The minus in front is b/c we're on the other side of the eqn
            #print(f"ratio for `{self.chem_data.get_name(species)}`: {ratio}")
            if not np.allclose(ratio, baseline_ratio):
                return False

        return True



    def stoichiometry_checker_entire_run(self) -> bool:
        """
        Verify that the stoichiometry is satisfied in all the reaction steps,
        using the diagnostic data from an earlier run

        IMPORTANT: this function is currently meant for simulations involving only 1 reaction (TODO: generalize)

        :return:    True if everything checks out, or False otherwise
        """
        if self.diagnostic_rxn_data == {}:
            print("WARNING *** In order to run stoichiometry_checker_entire_run(), "
                  "the diagnostics must be turned on, with set_diagnostics(), prior to the simulation run!")
            return False

        number_rxns = self.chem_data.number_of_reactions()
        assert number_rxns == 1, \
                f"stoichiometry_checker_entire_run(): this function is currently designed for just 1 reaction " \
                f"(whereas {number_rxns} are present)"

        for rxn_index in range(number_rxns):
            diagnostic_data = self.get_diagnostic_rxn_data(rxn_index=rxn_index, print_reaction=False)
            for row_index in range(len(diagnostic_data)):
                df_row = diagnostic_data.loc[row_index]     # Row in the Panda's data frame of diagnostic data
                chemical_delta_list = self._delta_names()           # EXAMPLE: ["Delta A", "Delta B", "Delta C"]
                delta = df_row[chemical_delta_list].to_numpy(dtype='float32')      # Extract select columns from the data frame row, and turn into Numpy array
                status = self._stoichiometry_checker_from_deltas(rxn_index=rxn_index, delta_arr=delta)
                if not status:
                    print(f"Stoichiometry NOT satisfied by reaction # {rxn_index}: "
                          f"see row # {row_index} in the diagnostic data for that reaction, from get_diagnostic_rxn_data()")
                    print(df_row )
                    return False

        return True



    def _delta_names(self) -> [str]:
        """
        Return a list of strings, with the names of ALL the registered chemicals,
        in their index order, each prefixed by the string "Delta "
        EXAMPLE: ["Delta A", "Delta B", "Delta X"]

        :return:    A list of strings
        """
        chemical_list = self.chem_data.get_all_names()      # EXAMPLE: ["A", "B", "X"]
        chemical_delta_list = ["Delta " + name
                               for name in chemical_list]
        return chemical_delta_list



    def _delta_conc_dict(self, delta_conc_arr: np.ndarray) -> dict:
        """
        Convert a Numpy array into a dict, based on all the registered chemicals.
        The keys are the chemical names, prefixed by "Delta "

        :param delta_conc_arr:  A Numpy array of "delta concentrations".  EXAMPLE: array[1.23, 52.2]
        :return:                A dictionary such as {"Delta A": 1.23, "Delta X": 52.2}
        """
        chemical_delta_list = self._delta_names()    # EXAMPLE: ["Delta A", "Delta X"]
        assert len(chemical_delta_list) == len(delta_conc_arr), \
            f"_delta_conc_dict(): mismatch in number of chemicals ({len(chemical_delta_list)} vs. {len(delta_conc_arr)})"

        delta_conc_dict = {}
        for i, delta_name in enumerate(chemical_delta_list):
            delta_conc_dict[delta_name] = delta_conc_arr[i]

        return delta_conc_dict




    #####  DIAGNOSTICS  #####

    def set_diagnostics(self):
        # Turn on the overall diagnostics mode
        self.diagnostics = True

    def unset_diagnostics(self):
        # Turn off the overall diagnostics mode; existing diagnostics data, if any, is left untouched
        self.diagnostics = False



    #####  1. diagnostic_rxn_data  #####

    def save_diagnostic_rxn_data(self, time_step, increment_vector_single_rxn: Union[np.array, None],
                                 rxn_index :int, caption="") -> None:
        """
        Save up diagnostic data for 1 reaction, for a simulation step
        (by convention, regardless of whether the step is completed or aborted)

        :param time_step:                   The duration of the current simulation step
        :param increment_vector_single_rxn: A Numpy array of size equal to the total number of chemical species,
                                                containing the "delta concentrations" for
                                                ALL the chemicals (whether involved in the reaction or not)
        :param rxn_index:                   The integer index (0-based) to identify the reaction of interest
        :param caption:                     OPTIONAL string to describe the snapshot
        :return:                            None
        """
        # Validate the reaction index
        self.chem_data.assert_valid_rxn_index(rxn_index)

        # Validate increment_vector_single_rxn
        if increment_vector_single_rxn is None:     # If the values aren't available (as is the case in aborts)
            increment_vector_single_rxn = np.full(self.chem_data.number_of_chemicals(), np.nan)
        else:
            assert len(increment_vector_single_rxn) == self.chem_data.number_of_chemicals(), \
                f"save_diagnostic_rxn_data(): the length of the Numpy array increment_vector_single_rxn " \
                f"({len(increment_vector_single_rxn)}) doesn't match the total numbers of chemicals ({self.chem_data.number_of_chemicals()})"


        if self.diagnostic_rxn_data == {}:    # INITIALIZE the dictionary self.diagnostic_rxn_data, if needed
            for i in range(self.chem_data.number_of_reactions()):
                self.diagnostic_rxn_data[i] = MovieTabular(parameter_name="START_TIME")       # One per reaction

        data_snapshot = {}
        # Add entries to the above dictionary, starting with the Delta conc. for all the chemicals
        for index, conc in enumerate(increment_vector_single_rxn):
            data_snapshot["Delta " + self.chem_data.get_name(index)] = conc

        #data_snapshot["reaction"] = rxn_index           # TODO: now redundant because factored out into separate data frames
        data_snapshot["time_step"] = time_step
        self.diagnostic_rxn_data[rxn_index].store(par=self.system_time,
                                                  data_snapshot=data_snapshot, caption=caption)



    def comment_diagnostic_rxn_data(self, msg: str) -> None:
        """
        Set the comment field of the last record for EACH of the reaction-specific dataframe

        :param msg: Value to set the comment field to
        :return:    None
        """
        for rxn_index in range(self.chem_data.number_of_reactions()):
            self.diagnostic_rxn_data[rxn_index].set_caption_last_snapshot(msg)



    def get_diagnostic_rxn_data(self, rxn_index :int, head=None, tail=None,
                                t=None, print_reaction=True) -> pd.DataFrame:
        """
        Return a Pandas dataframe with the diagnostic run data of the requested SINGLE reaction,
        from the time that the diagnostics were activated by a call to set_diagnostics().

        In particular, the dataframe contains the "Delta" values for each of the chemicals
        involved in the reaction - i.e. the change in their concentrations
        over the time interval that *STARTS* at the value in the "TIME" column.
        (So, there'll be no row with the final current System Time)

        Note: entries are always added, even if an interval run is aborted, and automatically re-done.

        Optionally, print out a brief description of the reaction.

        Optionally, limit the dataframe to a specified numbers of rows at the end,
        or just return one entry corresponding to a specific time
        (the row with the CLOSEST time to the requested one, which will appear in an extra column
        called "search_value")

        :param rxn_index:       The integer index (0-based) to identify the reaction of interest
                                TODO: if not specified, show all reactions in turn

        :param head:            (OPTIONAL) Number of records to return,
                                    from the start of the diagnostic dataframe.
        :param tail:            (OPTIONAL) Number of records to return,
                                    from the end of the diagnostic dataframe.
                                    If either the "head" arguments is passed, this argument will get ignored

        :param t:               (OPTIONAL) Individual time to pluck out from the dataframe;
                                    the row with closest time will be returned.
                                    If this parameter is specified, an extra column - called "search_value" -
                                    is inserted at the beginning of the dataframe
                                    If either the "head" or the "tail" arguments are passed, this argument will get ignored

        :param print_reaction:  (OPTIONAL) If True (default), concisely print out the requested reaction

        :return:                A Pandas data frame with (all or some of)
                                    the diagnostic data of the specified reaction.
                                    Columns of the dataframes:
                                    'START_TIME' 'Delta A' 'Delta B'... 'time_step' 'caption'
        """
        # Validate the reaction index
        self.chem_data.assert_valid_rxn_index(rxn_index)

        if print_reaction:
            print("Reaction: ", self.chem_data.single_reaction_describe(rxn_index=rxn_index, concise=True))

        movie_obj = self.diagnostic_rxn_data.get(rxn_index)    # Object of type "MovieTabular"

        if movie_obj is None:
            raise Exception(f"get_diagnostic_rxn_data(): no diagnostics data exists for reaction index {rxn_index} ;"
                            f" did you call set_diagnostics() prior to the simulation run?")

        return movie_obj.get_dataframe(head=head, tail=tail, search_col="START_TIME", search_val=t)



    #####  2. diagnostic_conc_data  #####

    def save_diagnostic_conc_data(self, system_data) -> None:
        """
        To save the diagnostic concentration data during the run, indexed by the current System Time.
        Note: if an interval run is aborted, by convention NO entry is created here

        :return: None
        """
        self.diagnostic_conc_data.store(par=self.system_time,
                                        data_snapshot=system_data)



    def get_diagnostic_conc_data(self) -> pd.DataFrame:
        """
        Return the diagnostic concentration data saved during the run.
        This will be a complete set of simulation steps,
        even if we only saved part of the history during the run

        Note: if an interval run is aborted, by convention NO entry is created here

        :return: A Pandas dataframe, with the columns:
                    'TIME' 	'A' 'B' ... 'caption'
                 where 'A', 'B', ... are all the chemicals
        """
        return self.diagnostic_conc_data.get_dataframe()



    #####  3. save_diagnostic_decisions_data  #####

    def save_diagnostic_decisions_data(self, data, caption="") -> None:
        """
        Used to save the diagnostic concentration data during the run, indexed by the current System Time.
        Note: if an interval run is aborted, by convention an entry is STILL created here

        :return: None
        """
        self.diagnostic_decisions_data.store(par=self.system_time,
                                             data_snapshot=data, caption=caption)


    def get_diagnostic_decisions_data(self) -> pd.DataFrame:
        """
        Determine and return the diagnostic data about concentration changes at every step - EVEN aborted ones

        :return:    A Pandas dataframe with a "TIME" column, and columns for all the "Delta concentration" values
        """
        return self.diagnostic_decisions_data.get_dataframe()



    def get_diagnostic_decisions_data_ALT(self) -> pd.DataFrame:
        """
        Determine and return the diagnostic data about concentration changes at every step - EVEN aborted ones
        TODO: OBSOLETE - BEING PHASED OUT

        :return:    A Pandas dataframe with a "TIME" column, and columns for all the "Delta concentration" values
        """
        fields = self._delta_names()    # EXAMPLE: ["Delta A", "Delta B", "Delta C"]

        rxn_0 = self.get_diagnostic_rxn_data(rxn_index=0, print_reaction=False)

        df_sum = rxn_0[fields]

        time_col = rxn_0["START_TIME"]

        n_rows = len(df_sum)

        for rxn_index in range(1, self.chem_data.number_of_reactions()):
            rxn = self.get_diagnostic_rxn_data(rxn_index=rxn_index, print_reaction=False)
            assert len(rxn)  == n_rows, "get_diagnostic_delta_data(): this function cannot be used because different dataframes" \
                                        "of diagnostic reaction data have different numbers of rows"
            df_sum += rxn[fields]

        df_sum.insert(0, "START_TIME", time_col)

        return df_sum



    def explain_reactions(self) -> bool:
        """
        Provide a detailed explanation of all the steps of the reactions,
        from the saved diagnostic data

        WARNING: Currently designed only for exactly 2 reactions!  TODO: generalize to any number of reactions

        TODO: test and validate usefulness, now that substeps got eliminated
        TODO: allow arguments to specify the min and max reaction times during which to display the explanatory data

        :return:    True if the diagnostic data is consistent for all the steps of all the reactions,
                    or False otherwise
        """
        if not self.diagnostics:
            print("No diagnostic information is available:\nIn order to run explain_reactions(), "
                  "call set_diagnostics() prior to running single_compartment_react()")
            return False

        number_reactions = self.chem_data.number_of_reactions()

        assert number_reactions == 2, \
            "explain_reactions() currently ONLY works when exactly 2 reactions are present. " \
            "Future versions will lift this restriction"

        row_baseline = 0
        row_list = [0, 0]       # TODO: generalize
        active_list = [0, 1]    # ALL the reactions.  TODO: generalize

        self._explain_reactions_helper(active_list=active_list,
                                       row_baseline=row_baseline, row_list=row_list)


        while row_baseline < len(self.get_diagnostic_conc_data()) - 2:
            row_baseline += 1

            print("Advance a single-step across all tables")

            for i in active_list:
                row_list[i] += 1

            status = self._explain_reactions_helper(active_list=active_list, row_baseline=row_baseline, row_list=row_list) # TODO: generalize


            if not status:
                return False    # Error termination

        # END while

        return True             # Successful termination



    def _explain_reactions_helper(self, active_list, row_baseline, row_list) -> bool:
        """
        Helper function for explain_reactions()

        :param active_list:
        :param row_baseline:
        :param row_list:
        :return:            True is the diagnostic data is consistent for this step, or False otherwise
        """
        print("-----------")
        print("ROW of baseline data: ", row_baseline)

        current_time = self.get_diagnostic_conc_data().loc[row_baseline]['TIME']
        print(f"TIME = {current_time:.5g}")

        print("row_list: ", row_list)
        print("active_list: ", active_list)

        chemical_list = self.chem_data.get_all_names()
        chemical_delta_list = self._delta_names()

        conc_arr_before = self.get_diagnostic_conc_data().loc[row_baseline][chemical_list].to_numpy(dtype='float16')
        print("baseline concentrations: ", conc_arr_before)

        delta_cumulative = np.zeros(self.chem_data.number_of_chemicals(),
                                    dtype=float)  # One element per chemical species

        # For each reaction
        for rxn_index in range(self.chem_data.number_of_reactions()):
            if (rxn_index in active_list):  # If reaction is tagged as "fast"
                row = row_list[rxn_index]   # Row in the data frame for the diagnostic data on this reaction
                delta_rxn = self.get_diagnostic_rxn_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy(dtype='float16')
                print(f"From fast rxn {rxn_index}: delta_rxn = {delta_rxn}")
            else:                           # If reaction is tagged as "slow"
                row = row_list[rxn_index]   # Row in the data frame for the diagnostic data on this reaction
                delta_rxn = self.get_diagnostic_rxn_data(rxn_index=rxn_index).loc[row][chemical_delta_list].to_numpy(dtype='float16')
                #delta_rxn = np.zeros(self.reaction_data.number_of_chemicals(),
                #                     dtype=float)   # This is the way it was in release beta 19
                print(f"From slow rxn {rxn_index}: delta_rxn = {delta_rxn}")


            delta_cumulative += delta_rxn
        # END for

        print("delta_cumulative: ", delta_cumulative)

        conc_after = conc_arr_before + delta_cumulative
        print("updated concentrations: ", conc_after)

        next_system_state = self.get_diagnostic_conc_data().loc[row_baseline + 1][chemical_list].to_numpy(dtype='float16')
        print(f"concentrations from the next row ({row_baseline + 1}) of the system state: ", next_system_state)

        status = np.allclose(conc_after.astype(float), next_system_state.astype(float))
        if status:
            print("Match OK")
        else:
            print("****************************************   MISMATCH!!!  ****************************************")

        print("-----------")

        return status



    def explain_time_advance(self, return_times=False, silent=False, use_history=False) -> Union[None, tuple]:
        """
        Use the saved-up diagnostic data, to print out details of the timescales of the reaction run

        If diagnostics weren't enabled ahead of calling this function, an Exception is raised

        EXAMPLE of output:
            From time 0 to 0.0304, in 17 FULL steps of 0.0008
            (for a grand total of 38 FULL steps)

        :param return_times:    If True, all the critical times (times where the interval steps change)
                                    are saved and returned as a list
        :param silent:          If True, nothing gets printed out
        :param use_history:     If True, use the system history in lieu of the diagnostic data;
                                    to keep in mind is the fact that the user might only have asked
                                    for PART of the history to be saved
        :return:                Depending on the argument return_times, either None, or a pair with 2 lists:
                                        1 - list of time values
                                        2 - list of step sizes  (will have one less element than the first list)
        """
        if use_history:
            df = self.get_history()
            assert not df.empty , \
                "UniformCompartment.explain_time_advance(): no history data found.  " \
                "Did you run the reaction simulation prior to calling this function?"
        else:
            df = self.get_diagnostic_conc_data()
            assert not df.empty , \
                "UniformCompartment.explain_time_advance(): no diagnostic data found.  " \
                "Did you call set_diagnostics() prior to the reaction run?"

        n_entries = len(df)
        if use_history:
            t = list(df["SYSTEM TIME"])    # List of times (simulation step points)
        else:
            t = list(df["TIME"])    # List of times (simulation step points)

        if n_entries == 1:
            print(f"UniformCompartment.explain_time_advance(): no time advance found. "
                  f"Diagnostics only contain an initial System Time of {t[0]:.3g} ;  "
                  f"did you run the reaction simulation?")
            return


        # If we get thus far, we have at least 2 entries in the list "t" of times
        grad = np.diff(t)
        grad_shifted = np.insert(grad, 0, 0)  # Insert a zero at the top (shifting down the array)
        #print(grad)
        #print(grad_shifted)
        start_interval = t[0]

        critical_times = [start_interval]  # List of times to report on (start time, plus times where the steps change)
        step_sizes = []

        total_number_full_steps = 0
        for i in range(1, n_entries-1):
            #print(i)
            if not np.allclose(grad[i], grad_shifted[i]):
                #print (f"   Detection at element {i} : time={t[i]}")
                #print (f"   From time {start_interval} to {t[i]}, in steps of {grad_shifted[i]}")
                n_full_steps_taken = self._explain_time_advance_helper(t_start=start_interval, t_end=t[i],
                                                                       delta_baseline=grad_shifted[i], silent=silent)

                start_interval = t[i]
                if n_full_steps_taken > 0:
                    total_number_full_steps += n_full_steps_taken
                    critical_times.append(t[i])
                    step_sizes.append(grad_shifted[i])


        # Final wrap-up of the interval's endpoint
        n_full_steps_taken = self._explain_time_advance_helper(t_start=start_interval, t_end=t[-1],
                                                               delta_baseline=grad_shifted[-1], silent=silent)

        if n_full_steps_taken > 0:
            total_number_full_steps += n_full_steps_taken
            critical_times.append(t[-1])
            step_sizes.append(grad_shifted[-1])

        if not silent:
            print(f"({total_number_full_steps} steps total)")

        if return_times:
            assert len(critical_times) == len(step_sizes) + 1, \
                "explain_time_advance(): validation error in the values to return"
            return (critical_times, step_sizes)



    def _explain_time_advance_helper(self, t_start, t_end, delta_baseline, silent: bool) -> Union[int, float]:
        """
        Using the provided data, about a group of same-size steps, create and print a description of it for the user

        :param t_start:
        :param t_end:
        :param delta_baseline:
        :param silent:          If True, nothing gets printed; otherwise, a line is printed out
        :return:                The corresponding number of FULL steps taken
        """
        if np.allclose(t_start, t_end):
            #print(f"   [Ignoring interval starting and ending at same time {t_start:.3g}]")
            return 0

        n_steps = round((t_end - t_start) / delta_baseline)
        step_s = "step" if n_steps == 1 else "steps"  # singular vs. plural
        if not silent:
            print(f"From time {t_start:.4g} to {t_end:.4g}, in {n_steps} {step_s} of {delta_baseline:.3g}")
        return n_steps




    #####################################################################################################

    '''                                    ~   GRAPHICS   ~                                           '''

    def ________GRAPHICS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def plot_history(self, chemicals=None, colors=None, title=None, title_prefix=None, xrange=None,
                     ylabel=None,
                     vertical_lines=None, show_intervals=False, show=False) -> go.Figure:
        """
        Using plotly, draw the plots of concentration values over time, based on history data that gets
        automatically saved when running reactions.

        Note: if this plot is to be later combined with others, use PlotlyHelper.combine_plots()

        :param chemicals:   (OPTIONAL) List of the names of the chemicals whose concentration changes are to be plotted,
                                or a string with just one name;
                                if None, then display all
        :param colors:      (OPTIONAL) Either a single color (string with standard plotly name, such as "red"),
                                or list of names to use, in order; if None, then use the hardwired defaults
        :param title:       (OPTIONAL) Title for the plot;
                                if None, use default titles that will vary based on the # of reactions; EXAMPLES:
                                    "Changes in concentrations for 5 reactions"
                                    "Reaction `A <-> 2 B` .  Changes in concentrations with time"
                                    "Changes in concentration for `2 S <-> U` and `S <-> X`"
        :param title_prefix: (OPTIONAL) If present, it gets prefixed (followed by ".  ") to the title,
                                    whether the title is specified by the user or automatically generated
        :param xrange:          (OPTIONAL) list of the form [t_start, t_end], to initially show only a part of the timeline.
                                    Note: it's still possible to zoom out, and see the excluded portion
        :param ylabel:          (OPTIONAL) Caption to use for the y-axis.
                                    By default, the name in `the chemicals` argument, in square brackets, if only 1 chemical,
                                    or "Concentration" if more than 1 (a legend also shown)
        :param vertical_lines:  (OPTIONAL) Ignored if the argument `show_intervals` is specified.
                                    TODO: rename add_vertical_lines
                                    List or tuple or Numpy array or Pandas series
                                    of x-coordinates at which to draw thin vertical dotted gray lines.
                                    If the number of vertical line is so large as to overwhelm the plot,
                                    only a sample of them is shown.
                                    Note that vertical lines, if requested, go into the plot's "layout";
                                    as a result they might not appear if this plot is later combined with another one.
        :param show_intervals:  (OPTIONAL) If True, it over-rides any value passed to the `vertical_lines` argument,
                                    and draws thin vertical dotted gray lines at all the x-coords
                                    of the data points in the saved history data;
                                    also, it adds a comment to the title.
        :param show:        If True, the plot will be shown
                                Note: on JupyterLab, simply returning a plot object (without assigning it to a variable)
                                      leads to it being automatically shown

        :return:            A plotly "Figure" object
        """
        # TODO: allow alternate label for x-axis
        # TODO: allow specifying a yrange

        MAX_NUMBER_VERTICAL_LINES = 150     # Used to avoid extreme clutter in the plot, in case
                                            # a very large number of vertical lines is requested;
                                            # if this value is exceeded, then the vertical lines are sampled
                                            # infrequently enough to bring the total number below this value

        df = self.get_history()     # The expected columns are "SYSTEM TIME",
                                    # followed by concentrations for the various chemicals,
                                    # and a final "caption" column

        if chemicals is None:
            chemicals = self.chem_data.get_all_names()      # List of the chemical names.  EXAMPLE: ["A", "B", "H"]

        number_of_curves = len(chemicals)

        if colors is None:
            colors = PlotlyHelper.get_default_colors(number_of_curves)
        elif type(colors) == str:
            colors = [colors]


        if title is None:   # If no title was specified, create one based on how many reactions are present
            number_of_rxns = self.chem_data.number_of_reactions()
            if number_of_rxns > 2:
                title = f"Changes in concentrations for {number_of_rxns} reactions"
            elif number_of_rxns == 1:
                rxn_text = self.chem_data.single_reaction_describe(rxn_index=0, concise=True)   # The only reaction
                title = f"Reaction `{rxn_text}` .  Changes in concentrations with time"
            else:   # Exactly 2 reactions
                rxn_text_0 = self.chem_data.single_reaction_describe(rxn_index=0, concise=True)
                rxn_text_1 = self.chem_data.single_reaction_describe(rxn_index=1, concise=True)
                title = f"Changes in concentration for `{rxn_text_0}` and `{rxn_text_1}`"

        if title_prefix is not None:
            title = f"{title_prefix}  <br>{title}"

        if show_intervals:
            vertical_lines = df["SYSTEM TIME"]  # Make use of the simulation times
            title += " (time steps shown in dashed lines)"


        # Create the main plot
        fig = px.line(data_frame=df, x="SYSTEM TIME", y=chemicals,
                      title=title, range_x=xrange,
                      color_discrete_sequence = colors,
                      labels={"value": "Concentration", "variable": "Chemical"})

        if type(chemicals) == str:  # Somehow, the `labels` argument in px.line is ignored when chemicals is just a string
            if ylabel is None:
                fig.update_layout(yaxis_title=f"[{chemicals}]")     # EXAMPLE:  "[A]"
            else:
                fig.update_layout(yaxis_title=ylabel)


        if vertical_lines is not None:
            assert (type(vertical_lines) == list) or (type(vertical_lines) == tuple) \
                or (type(vertical_lines) == np.ndarray) or (type(vertical_lines) == pd.core.series.Series), \
                    "plot_curves(): the argument `vertical_lines`, " \
                    "if not None, must be a list or tuple or Numpy array or Pandas series of numbers (x-axis coords)"

            vline_list = []
            if xrange:
                step = 1    # Always show all vertical lines if a range on the x-axis was specified
            else:
                # Possibly limit the number of vertical lines shown
                step = 1 + len(vertical_lines) // MAX_NUMBER_VERTICAL_LINES
                if step > 1:
                    print(f"plot_curves() WARNING: Excessive number of vertical lines ({len(vertical_lines)}) - only showing 1 every {step} lines")

            for xi in vertical_lines[::step] :  # Notice that if step > 1 then we're sampling a subset of the elements
                # The following is the internal data structure used by Plotly Express,
                # for each of the vertical lines
                vline = {  'line': {'color': 'gray', 'dash': 'dot', 'width': 1},
                           'type': 'line',
                           'x0': xi,
                           'x1': xi,
                           'xref': 'x',
                           'y0': 0,
                           'y1': 1,
                           'yref': 'y domain'
                        }
                vline_list.append(vline)
                # Strangely, a direct call to fig.add_vline(), as done below, dramatically slows things down in case
                # of a large number of vertical lines; so, we'll be directly modifying the data structure of the "fig" dictionary
                #fig.add_vline(x=xi, line_width=1, line_dash="dot", line_color="gray")
            # END for
            fig['layout']['shapes'] = vline_list    # The vertical lines are regarded by Plotly Express as "shapes"
                                                    # that are stored in the figure's "layout"
        if show:
            fig.show()  # Actually display the plot

        return fig



    def plot_step_sizes(self, show_intervals=False) -> None:
        """
        Using plotly, draw the plot of the step sizes vs. time
        (only meaningful when the variable-step option was used).
        The same scale as plot_curves() will be used.
        This function requires the diagnostics option to be turned on, prior to running the simulation

        :param show_intervals:  If True, will add to the plot thin vertical dotted gray lines
                                    at the time steps
        :return:                None
        """
        (transition_times, step_sizes) = self.explain_time_advance(return_times=True, silent=True)

        x=transition_times
        y=step_sizes

        # Create a step plot (TODO: there might be a way to directly do this in plotly)
        new_x = [x[0]]
        new_y = [y[0]]
        for i, xi in enumerate(x[1:-1]) :   # Drop the first and last elements
            new_x.append(xi)
            new_y.append(y[i])

            new_x.append(xi)
            new_y.append(y[i+1])

        new_x.append(x[-1])
        new_y.append(y[-1])


        df = self.get_diagnostic_conc_data()    # Pandas dataframe with a column called "TIME"

        # Note: the step size at the final end time isn't a defined quantity - so, we'll just repeat
        #       the last value, to maintain the full x-axis size
        fig = px.line(x=new_x, y=new_y)

        if show_intervals:
            for xi in df["TIME"]:
                fig.add_vline(x=xi, line_width=1, line_dash="dot", line_color="gray")


        fig.update_layout(title='Simulation step sizes',
                          xaxis_title='SYSTEM TIME',
                          yaxis_title='Step size')

        fig.show()




    #####################################################################################################

    '''                                      ~   HISTORY   ~                                          '''

    def ________HISTORY________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def add_snapshot(self, species=None, caption="", time=None, system_data=None) -> None:
        """
        Preserve some or all the chemical concentrations into the history,
        linked to the passed time (by default the current System Time),
        with an optional caption.

        EXAMPLES:  add_snapshot()
                    add_snapshot(species=['A', 'B']), caption="Just prior to infusion")

        :param species:     (OPTIONAL) list of name of the chemical species whose concentrations we want to preserve for later use.
                                If not specified, save all
        :param caption:     (OPTIONAL) caption to attach to this preserved data
        :param time:        (OPTIONAL) time value to attach to the snapshot (default: current System Time)
        :param system_data: (OPTIONAL) a Numpy array of concentration values, in the same order as the
                                index of the chemical species; by default, use the SYSTEM DATA
                                (which is set and managed by various functions)
        :return:            None
        """

        data_snapshot = self.get_conc_dict(species=species, system_data=system_data)    # This will be a dict

        if time is None:
            time = self.system_time     # By default, use the system time

        self.history.store(par=time,
                           data_snapshot=data_snapshot, caption=caption)



    def get_history(self, t_start=None, t_end=None, head=None, tail=None, t=None, columns=None) -> pd.DataFrame:
        """
        Retrieve and return a Pandas dataframe with the system history that had been saved
        using add_snapshot().
        Optionally, restrict the result with a start and/or end times,
        or by limiting to a specified numbers of rows at the end

        :param t_start: (OPTIONAL) Start time in the "SYSTEM TIME" column.  Watch out for roundoff errors!
        :param t_end:   (OPTIONAL) End time.  Watch out for roundoff errors!
        :param head:    (OPTIONAL) Number of records to return,
                                   from the start of the history dataframe.
        :param tail:    (OPTIONAL) Number of records to consider, from the end of the history dataframe
        :param t:       (OPTIONAL) Individual time to pluck out from the dataframe;
                                   the row with closest time will be returned.
                                   If this parameter is specified, an extra column - called "search_value" -
                                   is inserted at the beginning of the dataframe.
                                   If either the "head" or the "tail" arguments are passed, this argument will get ignored
        :param columns: (OPTIONAL) List of columns to return; if not specified, all are returned.
                                   Make sure to include "SYSTEM TIME" in the list, if the time variable needs to be included

        :return:        A Pandas dataframe
        """
        #TODO: allow searches also for columns other than "SYSTEM TIME"

        # Note self.history is an object of class MovieTabular
        df = self.history.get_dataframe(head=head, tail=tail, search_val=t,
                                        search_col="SYSTEM TIME", val_start=t_start, val_end=t_end)

        if columns:
            assert type(columns) == list, \
                "get_history(): the argument `columns`, if specified, must be a list"
            return df[columns]

        return df



    def get_historical_concentrations(self, row=None, t=None, df=None) -> np.array:
        """
        Typically used to retrieve a system snapshot from its history.

        Return a Numpy array with ALL the chemical concentrations (in their index order)
        from one row in the given Pandas data frame (by default, the system history);
        the row can be identified either by it row number, or by the system time.

        :param row: (OPTIONAL) Integer with the zero-based row number of the system history
                        (which is a Pandas data frame)
        :param t:   (OPTIONAL) Individual time to pluck out from the dataframe;
                        the row with closest time will be returned.
                        Exactly one of "t" and "row" must be specified
        :param df:  (OPTIONAL) A Pandas data frame with concentration information in columns that have
                        the names of the chemicals (if None, the system history is used)
        :return:    A Numpy array of floats.  EXAMPLE: array([200., 40.5])
        """
        assert row is None or t is None, "Cannot specify both arguments `row` and `t`"
        assert row is not None or t is not None, "Must specify either argument `row` and `t`"

        if df is None:
            if t:
                df = self.get_history(t=t)
            else:
                df = self.get_history()

        chem_list = self.chem_data.get_all_names()  # List of all the chemicals' names

        if row:
            return df.loc[row][chem_list].to_numpy(dtype='float32')
        else:
            return df.iloc[0][chem_list].to_numpy(dtype='float32')




    #####################################################################################################

    '''                                ~   RESULT ANALYSIS   ~                                        '''

    def ________RESULT_ANALYSIS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def is_in_equilibrium(self, rxn_index=None, conc=None, tolerance=1, explain=True, verbose=None) -> Union[bool, dict]:
        """
        Ascertain whether the given concentrations (by default the current System concentrations)
        are in equilibrium for the specified reactions (by default, for all reactions)

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest;
                                if None, then check all the reactions
        :param conc:        Dict with the concentrations of the species involved in the reaction(s).
                            The keys are the chemical names
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
                            If None, then use the current System concentrations instead
        :param tolerance:   Allowable relative tolerance, as a PERCENTAGE,
                                to establish satisfactory match with expected values
        :param explain:     If True, print out details about the analysis,
                                incl. the formula(s) being used to check the equilibrium
                                EXAMPLES:   "([C][D]) / ([A][B])"
                                            "[B] / [A]^2"
        :param verbose:     Alternate name for argument `explain`

        :return:            Return True if ALL the reactions are close enough to an equilibrium,
                                as allowed by the requested tolerance;
                                otherwise, return a dict of the form {False: [list of reaction indexes]}
                                for all the reactions that failed the criterion
                                EXAMPLE:  {False: [3, 6]}
        """
        # TODO: optionally display last lines in diagnostic data, if available
        # TODO: phase out "explain" in favour of "verbose"

        if verbose is not None:
            explain = verbose

        if conc is None:
            conc=self.get_conc_dict()   # Use the current System concentrations, as a dict.
                                        # EXAMPLE: {'A': 23.9, 'B': 36.1}

        failures_dict = {False: []}     # 1-element dict whose value is
                                        # a list of reactions that fail to meet the criterion for equilibrium

        if rxn_index is not None:
            # Check the 1 reaction that was requested
            if explain:
                description = self.chem_data.single_reaction_describe(rxn_index=rxn_index, concise=True)
                print(description)

            status = self.reaction_in_equilibrium(rxn_index=rxn_index, conc=conc, tolerance=tolerance, explain=explain)
            if not status:
                failures_dict = {False: [rxn_index]}    # Make a note that this reaction failed the equilibrium test

        else:
            # Check ALL the reactions
            status = True       # Overall status
            description_list = self.chem_data.multiple_reactions_describe(concise=True)
            for rxn_index in range(self.chem_data.number_of_reactions()):
                # For each reaction
                if explain:
                    print(description_list[rxn_index])

                single_status = self.reaction_in_equilibrium(rxn_index=rxn_index, conc=conc, tolerance=tolerance, explain=explain)

                if not single_status:
                    status = False
                    failures_dict[False].append(rxn_index)      # Make a note that this reaction failed the equilibrium test

        if status:
            return True
        else:
            return failures_dict



    def reaction_in_equilibrium(self, rxn_index :int, conc, tolerance, explain :bool) -> bool:
        """
        Ascertain whether the given concentrations are in equilibrium for the specified SINGLE reaction;
        return True or False, accordingly.

        Pathological case: if at least one of the reactants AND at least one of the products have zero
        concentration, then the reaction is "stuck" - and thus regarded in "equilibrium"

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest
        :param conc:        Dictionary with the concentrations of the species involved in the reaction.
                            The keys are the chemical names
                                EXAMPLE: {'A': 23.9, 'B': 36.1}
        :param tolerance:   Allowable relative tolerance, as a PERCENTAGE,
                                to establish satisfactory match with expected values
        :param explain:     If True, print out the formula being used, as well as the concentrations (reactants, then products)
                                and some other info.
                                EXAMPLES of formulas:   "([C][D]) / ([A][B])"
                                                        "[B] / [A]^2"
        :return:            True if the given reaction is close enough to an equilibrium,
                                as allowed by the requested tolerance
        """
        rxn = self.chem_data.get_reaction(rxn_index)    # Look up the requested reaction

        rate_ratio = rxn.extract_forward_rate() / rxn.extract_reverse_rate()  # Ratio of forward/reverse reaction rates

        result = rxn.reaction_quotient(conc=conc, explain=explain)

        # Unpack the result
        if explain:
            rxn_quotient, formula = result
        else:
            rxn_quotient = result

        if explain:
            # Prepare a concise listing from the given concentrations,
            # only including the concentrations that are applicable to this reaction
            all_applicable_concs = []
            '''
            for species_name in rxn.extract_chemicals_in_reaction(exclude_enzyme=False):
                s = f"[{species_name}] = {conc[species_name]:,.4g}"         # EXAMPLE: "[A] = 20.3"
                all_applicable_concs.append(s)
            '''
            reactants = rxn.extract_reactant_names(exclude_enzyme=False)
            for species_name in reactants:
                s = f"[{species_name}] = {conc[species_name]:,.4g}"         # EXAMPLE: "[A] = 20.3"
                all_applicable_concs.append(s)

            products = rxn.extract_product_names(exclude_enzyme=False)
            for species_name in products:
                if species_name not in reactants:           # Don't report the same concentration twice!
                    s = f"[{species_name}] = {conc[species_name]:,.4g}"         # EXAMPLE: "[B] = 0.3"
                    all_applicable_concs.append(s)

            all_applicable_concs_str = " ; ".join(all_applicable_concs)     # EXAMPLE: "[A] = 20.3 ; [B] = 0.3"


        # Handle the special case of zero concentration in some of the reactants AND in some of the products
        if np.isnan(rxn_quotient):
            if explain:
                print(f"Final concentrations: {all_applicable_concs_str}")
                print("Reaction IS in equilibrium because it can't proceed in either direction "
                      "due to zero concentrations in both some reactants and products!\n")
            return True

        if np.isinf(rxn_quotient):      # If only the denominator is 0
            if explain:
                print(f"Final concentrations: {all_applicable_concs_str}")
                print("Reaction is NOT in equilibrium because either the products contain a zero concentration\n")
            return False


        # Normal case
        status = np.allclose(rxn_quotient, rate_ratio, rtol=tolerance/100., atol=0)

        if explain:
            print(f"Final concentrations: {all_applicable_concs_str}")
            print(f"1. Ratio of reactant/product concentrations, adjusted for reaction orders: {rxn_quotient:,.6g}")
            print(f"    Formula used:  {formula}")
            print(f"2. Ratio of forward/reverse reaction rates: {rate_ratio:,.6g}")
            print(f"Discrepancy between the two values: {100 * abs(rxn_quotient - rate_ratio)/rate_ratio :,.4g} %")
            if status:
                print(f"Reaction IS in equilibrium (within {tolerance:.2g}% tolerance)\n")
            else:
                print(f"Reaction is NOT in equilibrium (not within {tolerance:.2g}% tolerance)\n")

        return status



    def find_equilibrium_conc(self, rxn_index :int) -> dict:
        """
        Determine the equilibrium concentrations that would be reached by the chemicals
        participating in the specified reaction, given their current concentrations,
        IN THE ABSENCE of any other reaction.

        IMPORTANT: currently limited to just aA + bB <-> cC + dD reactions, first-order in all chemicals,
                   (some of terms can be missing)
                   An Exception will be raised in all other cases

        :param rxn_index:   The integer index (0-based) to identify the reaction of interest
        :return:            A dictionary of the equilibrium concentrations of the
                                chemicals involved in the specified reaction
                            EXAMPLE:  {'A': 24.0, 'B': 36.0, 'C': 1.8}
        """
        #TODO: generalize to reactions with more terms and higher order

        rxn = self.chem_data.get_reaction(rxn_index)    # Look up the requested reaction

        reactants = rxn.extract_reactants()
        products = rxn.extract_products()
        K = rxn.extract_equilibrium_constant()


        assert len(reactants) <= 2, \
                "find_equilibrium_conc(): Currently only implemented " \
                "for reactions with at most 2 reactants and 2 products"
        assert len(products) <= 2, \
                "find_equilibrium_conc(): Currently only implemented " \
                "for reactions with at most 2 reactants and 2 products"

        '''
        For reactions of the form aA + bB <-> cC + dD  that are first-order in all chemicals,
        the equilibrium equation is:  
                [(C0 + c*m) (D0 + d*m)] / [(A0 - a*m) (B0 - b*m)]  =  K
            
        where K is the Equilibrium constant,
        and the unknown m (to be solved for) is the number of "moles of forward reaction", 
        from the starting point to the equilibrium point
        
        For reaction terms that aren't present (for example the "D" part in the reaction A + B <-> C),
        we'll use 0 for the "stoichiometry coefficient" and 1 for the "initial concentration";
        that's a hack to make the multiplicative terms of the form X0 + x*m  (where the X's are the A,B,C,D and the x's are the a,b,c,d)
        to be identical equal to 1, and thus have no effect on the solution of the above equation
        '''
        a = 0
        b = 0
        c = 0
        d = 0
        A0 = 1
        B0 = 1
        C0 = 1
        D0 = 1


        name_map = {}   # To transform the A, B, C, D into the actual names

        # Look at the denominator of the "Reaction Quotient"
        for i, r in enumerate(reactants):
            # Loop over the reactants
            species_name =  rxn.extract_species_name(r)
            rxn_order =  rxn.extract_rxn_order(r)
            coefficient = rxn.extract_stoichiometry(r)

            assert rxn_order == 1, \
                "find_equilibrium_conc(): Currently only implemented for 1st order reactions"

            if i == 0:
                name_map["A"] = species_name
                a = coefficient
                A0 = self.get_chem_conc(species_name)
                assert A0 is not None, f"equilibrium_concentration(): unable to proceed because the " \
                                       f"concentration of `{species_name}` was not provided"
            else:
                name_map["B"] = species_name
                b = coefficient
                B0 = self.get_chem_conc(species_name)
                assert B0 is not None, f"equilibrium_concentration(): unable to proceed because the " \
                                       f"concentration of `{species_name}` was not provided"


        # Look at the numerator of the "Reaction Quotient"
        for i, p in enumerate(products):
            # Loop over the reaction products
            species_name =  rxn.extract_species_name(p)
            rxn_order =  rxn.extract_rxn_order(p)
            coefficient = rxn.extract_stoichiometry(p)

            assert rxn_order == 1, "find_equilibrium_conc(): Currently only implemented for 1st order reactions"

            if i == 0:
                name_map["C"] = species_name
                c = coefficient
                C0 = self.get_chem_conc(species_name)
                assert C0 is not None, f"equilibrium_concentration(): unable to proceed because the " \
                                       f"concentration of `{species_name}` was not provided"
            else:
                name_map["D"] = species_name
                d = coefficient
                D0 = self.get_chem_conc(species_name)
                assert D0 is not None, f"equilibrium_concentration(): unable to proceed because the " \
                                       f"concentration of `{species_name}` was not provided"


        #print("Initial values: ", A0, B0, C0, D0)
        #print("Stoichiometry coefficients: ", a, b, c, d)

        '''
        The equation we saw earlier,
                [(C0 + c*m) (D0 + d*m)] / [(A0 - a*m) (B0 - b*m)]  =  K
                
        can be expanded into a standard quadratic form for the unknown m :
        
                alpha * m**2 + beta * m + gamma = 0
                
        and then solved for m    
        '''
        alpha = (c * d - K * a * b)
        beta = d * C0 + c * D0 + K * (A0 * b + B0 * a)
        gamma = C0 * D0 - K * A0 * B0

        #print(name_map)
        #print("alpha, beta, gamma : ", alpha, beta, gamma)

        if alpha == 0:
            # The quadratic reduces to the linear equation:  beta * m + gamma = 0
            m = -gamma / beta
        else:
            sqrt_discriminant = math.sqrt(beta**2 - 4 * alpha * gamma)
            m1 = (-beta + sqrt_discriminant) / (2 * alpha)
            m2 = (-beta - sqrt_discriminant) / (2 * alpha)
            #print("m1, m2 : ", m1, m2)
            m = m1  # Let's start with one of the 2 possible solutions of the quadratic equation


        # After m "moles of forward reaction", the concentration of the reactant "A"
        # in aA + bB <-> cC + dD gets reduced by a*m . Likewise for the other terms.
        # Reaction products get increased.  Values for missing terms will be meaningless
        std_result = {"A" : A0 - a*m, "B" : B0 - b*m, "C" : C0 + c*m, "D" : D0 + d*m}

        if min(std_result.values()) < 0:    # If there's any negative value in the concentrations...
            # ...then repeat the computation using the other solution to the quadratic
            m = m2
            std_result = {"A" : A0 - a*m, "B" : B0 - b*m, "C" : C0 + c*m, "D" : D0 + d*m}


        # Let's translate our standard names A, B, C, D into the actual names,
        # and also drop any missing term
        result = {}
        for k, v in std_result.items():
            actual_name = name_map.get(k)
            if actual_name:     # Missing terms will get dropped out
                result[actual_name] = v

        return result




    def estimate_rate_constants(self, t :np.array, reactant_conc :np.array, product_conc :np.array,
                                reactant_name="Reactant", product_name="Product"):
        """
        Estimate the rate constants for a 1-st order reaction of the type A <-> B,
        given time evolution of [A] and [B] on a grid of time points (don't need to be equally spaced)

        IMPORTANT : Currently restricted to reactions with a 1:1 stoichiometry between the given reactant and product

        :param t:               A numpy array of time grid points where the other functions are specified
        :param reactant_conc:
        :param product_conc:
        :param reactant_name:
        :param product_name:
        :return:                A plotly "Figure" object.  The estimated rate constants are printed out
        """
        total_conc_arr = reactant_conc + product_conc
        total_conc = np.median(total_conc_arr)    # TODO: give warning or abort if there's too much variance
        sd = np.std(total_conc_arr)

        print(f"Total REACTANT + PRODUCT has a median of {total_conc:,.4g}, "
              f"\n    with standard deviation {sd:,.4g} (ideally should be zero)")


        # The rate of change of [reactant] with time
        deriv_reactant_conc = np.gradient(reactant_conc, t, edge_order=2)
        # The rate of change of [product] with time
        deriv_product_conc = np.gradient(product_conc, t, edge_order=2)


        median_sum_derivs = np.median(deriv_reactant_conc + deriv_product_conc)
        print(f"The sum of the time derivatives of reactant and product "
              f"\n    has a median of {median_sum_derivs:,.4g} (ideally should be zero)")


        # Do a least-square fit of the time-gradient of the PRODUCT concentration,
        # as a function of the REACTANT's concentration
        Y = deriv_product_conc  # The dependent variable:   B'(t)
        X = reactant_conc       # The independent variable: A(t)

        M = np.vstack([np.ones(len(Y)), X]).T
        # M is an nx2 matrix , where n is the number of data points.
        # The 1st column is all 1's, and the 2nd column contains the values of X

        a, b = np.linalg.lstsq(M, Y, rcond=None)[0]  # Carry out the least-square fit as: Y = a + b X
        print(f"Least square fit: Y = {a:,.4g} + {b:,.4g} X"
              f"\n    where X is the array [{reactant_name}] and Y is the time gradient of {product_name}")

        # Plot both Y and its least-square fit, as functions of X
        fig = PlotlyHelper.plot_curves(x=X, y=[Y , a + b*X],
                                       title=f"d/dt {product_name}(t) as a function of {reactant_name}(t), alongside its least-square fit",
                                       xlabel=f"{reactant_name}(t)", ylabel=f"{product_name}'(t)",
                                       curve_labels=[f"{product_name}'(t)", "Linear Fit"], legend_title="Curve vs Fit:",
                                       colors=['green', 'red'])

        kR = -a / total_conc

        kF = b - kR

        print(f"\n-> ESTIMATED RATE CONSTANTS: kF = {kF:,.4g} , kR = {kR:,.4g}")

        return fig




    def estimate_rate_constants_OLD(self, t :np.array, reactant_conc :np.array, product_conc :np.array,
                                product_name="Product"):
        """
        IMPORTANT : Currently restricted to reactions with a 1:1 stoichiometry between the given reactant and product

        :param t:               A numpy array of time grid points where the other functions are specified
        :param reactant_conc:
        :param product_conc:
        :param product_name:
        :return:                A plotly "Figure" object
        """
        total_conc_arr = reactant_conc + product_conc
        total_conc = np.median(total_conc_arr)    # TODO: give warning or abort if there's too much variance
        sd = np.std(total_conc_arr)

        print(f"Total REACTANT + PRODUCT has a median of {total_conc:,.4g}, "
              f"\n    with standard deviation {sd:,.4g} (expected: zero)")


        # The rate of change of [reactant] with time
        deriv_reactant_conc = np.gradient(reactant_conc, t, edge_order=2)
        # The rate of change of [product] with time
        deriv_product_conc = np.gradient(product_conc, t, edge_order=2)


        median_sum_derivs = np.median(deriv_reactant_conc + deriv_product_conc)
        print(f"The sum of the derivatives has a median of {median_sum_derivs:,.4g} (expected: zero)")


        # Do a least-square fit of the time-gradient of the product, as a function of the product's concentration
        Y = deriv_product_conc  # The dependent variable
        X = product_conc        # The independent variable

        M = np.vstack([np.ones(len(Y)), X]).T
        # M is an nx2 matrix , where n is the number of data points.  The 2nd column contains the values of X

        a, b = np.linalg.lstsq(M, Y, rcond=None)[0]  # Carry out the least-square fit as: Y = a + b X
        print(f"Least square fit: Y = {a:,.4g} + {b:,.4g} X"
              f"\n    where X is the array [{product_name}] and Y is its gradient")

        fig = PlotlyHelper.plot_curves(x=X, y=[Y , a + b*X],
                                       title=f"d/dt {product_name}(t) as a function of {product_name}(t), alongside its least-square fit",
                                       xlabel=f"{product_name}(t)", ylabel=f"{product_name}'(t)",
                                       curve_labels=[f"{product_name}'(t)", "Linear Fit"], legend_title="Curve vs Fit:",
                                       colors=['green', 'red'])

        kF = a / total_conc

        kR = -(a / total_conc) - b

        print(f"\n-> ESTIMATED RATE CONSTANTS: kF = {kF:,.4g}, kR = {kR:,.4g}")

        return fig



    def reach_threshold(self, chem :str, threshold) -> Union[float, None]:
        """

        :param chem:        The name of the chemical of interest
        :param threshold:   A number with the concentration value being sought
        :return:            The time at which the linearly-interpolated concentration of the specified chemical
                                first reaches the given threshold;
                                or None if it never does (in which case a message is printed)
        """
        # Prepare a Pandas dataframe with 2 columns
        df = self.get_history(columns=["SYSTEM TIME", chem])

        x_intersection = num.reach_threshold(df, x="SYSTEM TIME", y=chem, y_threshold=threshold)
        if x_intersection is None:
            print(f"reach_threshold(): the concentrations of `{chem}` never reaches the specified threshold of {threshold}")

        return x_intersection



    def curve_intersect(self, chem1 :str, chem2 :str, t_start, t_end, explain=False) -> (float, float):
        """
        Find and return the intersection of the concentrations of the chemicals chem1 and chem2,
        in the time interval [t_start, t_end]
        If more than one is present, only the first (smallest value of x-coord) one will be returned;
        so, the specified time interval should be narrow enough to bracket the intersection of interest

        :param chem1:   The name of the 1st chemical of interest
        :param chem2:   The name of the 2nd chemical of interest
        :param t_start: The start of the time interval being considered
        :param t_end:   The end of the time interval being considered
        :param explain: [OPTIONAL] If True, print out some details of the computation
        :return:        The pair (time of intersection, common value) ;
                            if not found, a message is printed, and None is returned
        """
        # Prepare a Pandas dataframe with 3 columns
        df = self.get_history(t_start=t_start, t_end=t_end, columns=["SYSTEM TIME", chem1, chem2])

        intersection = num.curve_intersect(df, x="SYSTEM TIME", var1=chem1, var2=chem2, explain=explain)
        if intersection is None:
            print(f"curve_intersect(): No intersection detected between the concentrations of `{chem1}` and `{chem2}` "
                  f"in the time interval [{t_start} - {t_end}]")

        return intersection



    def curve_intersection(self, chem1, chem2, t_start, t_end) -> (float, float):
        """
        TODO: OBSOLETE
        Find and return the intersection of the 2 curves in the columns var1 and var2,
        in the time interval [t_start, t_end]
        If there's more than one intersection, only one - in an unpredictable choice - is returned;
        so, the specified time interval should be narrow enough to avoid that

        :param chem1:   The name of the 1st chemical of interest
        :param chem2:   The name of the 2nd chemical of interest
        :param t_start: The start of the time interval being considered
        :param t_end:   The end of the time interval being considered
        :return:        The pair (time of intersection, common value)
        """
        print("\n*** The function curve_intersection() is now OBSOLETE : use curve_intersect() instead\n")

        return self.curve_intersect(chem1, chem2, t_start, t_end, explain=False)



    def extract_delta_concentrations(self, df, row_from :int, row_to :int, chem_list: [str]) -> np.array:
        """
        Extract the concentration changes of the specified chemical species from a Pandas dataframe
        of concentration values

        EXAMPLE:  extract_delta_concentrations(my_dataframe, 7, 8, ['A', 'B'])

        :param df:          A Pandas dataframe of concentration values (it MUST contain columns
                                with the names given in chem_list
        :param row_from:    Row number of the first row of data we're interested in
        :param row_to:      Row number of the last row of data we're interested in
        :param chem_list:   A list of names of chemicals

        :return:            A Numpy array of floats
        """
        # TODO: add validations
        from_values = df.loc[row_from][chem_list]
        to_values = df.loc[row_to][chem_list]
        return (to_values - from_values).astype(float).to_numpy(dtype='float32')
