import numpy as np
import pandas as pd
import math
from scipy.fft import rfft, rfftfreq    # Fast Fourier Transforms to extract frequency components
from scipy.stats import norm
from typing import Union, List, Tuple
from life123.collections import CollectionTabular
from life123.uniform_compartment import UniformCompartment
from life123.history import HistoryBinConcentration
import plotly.express as px
import plotly.graph_objects as pgo
from life123.html_log import HtmlLog as log
from life123.visualization.graphic_log import GraphicLog
from life123.visualization.plotly_helper import PlotlyHelper
from life123.visualization.colors import Colors



class BioSim1D:
    """
    1D simulations of diffusion and reactions,
    with an early partial implementation of membranes
    """

    def __init__(self, n_bins :int, chem_data=None, reaction_handler=None):
        """
        :param n_bins:          The number of compartments (bins) to use in the simulation
        :param chem_data:       [OPTIONAL] Object of class "ChemData";
                                    if not specified, it will get extracted
                                    from the "UniformCompartment" class (if passed to the next argument)
        :param reaction_handler:[OPTIONAL] Object of class "UniformCompartment";
                                    if not specified, it'll get instantiated here
        """
        self.debug = False

        self.n_bins = 0         # Number of spatial compartments (bins) used in the simulation

        self.n_species = 1      # The number of (non-water) chemical species   TODO: phase out?

        self.chem_data = None   # Object of type "ChemData", with info on the individual chemicals

        self.reactions = None   # Object of type "Reactions", with info on all the reactions

        self.reaction_dynamics = None   # Object of class "UniformCompartment"
                                        # TODO: for now just 1 object is instantiated;
                                        #       in the future, it might be 1 per bin (or bin cluster)

        self.system_length = None       # The linear extension of the system,
                                        # from the middle of the leftmost bin to the middle of the rightmost one.
                                        #  The de-facto default value, though not used, is (n_bins-1)

        self.global_Dx = 1              # Used in cases when not using ad-hoc local changes in x-scale resolution

        self.system = None      # Concentration data in the System we're simulating, for all the chemicals:
                                #   NumPy array of floats, of dimension: (n_species) x (n_bins)
                                #   Each row represents a species;
                                #   e.g, self.system[0] is the linear sequence for the 0-th chemical

        self.system_earlier = None  # NOT IN CURRENT USE.  Envisioned for simulations where the past 2 time states are used
                                    # to compute the state at the next time step

        self.membranes = None   # NumPy boolean array of dimension n_bins;
                                #   each boolean value marks the presence of a membrane in that bin.
                                #   Any given bin is expected to either contain, or not contain, a single membrane passing thru it
                                #   TODO - explore possible alternative: list of bin addresses the contains a membrane,
                                #          or other sparse-matrix representation

        self.A_fraction = None  # TODO - DEPRECATED
                                # NumPy array of floats, of dimension n_bins;
                                #   each value records the fraction (between 0. and 1.) of the LEFT one
                                #   of the 2 parts into which the membrane splits the bin.
                                #   0. if N/A

        self.system_B = None    # TODO - DEPRECATED
                                # Just like "system", but only for the "A" fractions of the bin, where applicable;
                                #   0. if N/A

        self.delta_diffusion = None  # Buffer for the concentration changes from diffusion step (n_species x n_bins)
        self.delta_reactions = None  # Buffer for the concentration changes from reactions step (n_species x n_bins)

        self.delta_reactions_B = None   # Same as above, but for the "other sides" of bins with membranes
                                        # TODO: explore sparse-matrix representations

        self.sealed = True           # If True, no exchange with the outside; if False, immersed in a "bath"

        # Only applicable if "sealed" is False:
        self.bath_concentrations = None      # A NumPy array for each species
        self.container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container

        self.time_step_threshold = 0.33333  # This is used to set an Upper Bound on the single time steps
                                            #   in the diffusion process.
                                            #   See explanation in file overly_large_single_timesteps.py

        self.history = CollectionTabular()  # TODO: phase out in favor of the new self.conc_history
                                            # To store user-selected snapshots of (parts of) the system,
                                            #   whenever requested by the user.
                                            #   Note that we're using the "tabular" format - friendly to Pandas

        self.conc_history = HistoryBinConcentration(active=False)

        self.system_time = None             # Global time of the system, from initialization on


        self._initialize_system(n_bins=n_bins, chem_data=chem_data, reaction_handler=reaction_handler)




    #########################################################################
    #                                                                       #
    #                           SYSTEM-WIDE                                 #
    #                                                                       #
    #########################################################################

    def _initialize_system(self, n_bins: int, chem_data=None, reaction_handler=None) -> None:
        """
        Initialize all concentrations to zero.
        Membranes, if present, need to be set later.

        :param n_bins:          The number of compartments (bins) to use in the simulation
        :param chem_data:       (OPTIONAL) Object of class "ReactionData";
                                    if not specified, it will get extracted from the "UniformCompartment" class
        :param reaction_handler:(OPTIONAL) Object of class "UniformCompartment";
                                    if not specified, it'll get instantiated here
        :return:                None
        """
        #TODO?: maybe allow optionally passing n_species in lieu of chem_data,
        #       and let it create and return the "Chemicals" object in that case

        assert type(n_bins) == int, "BioSim1D() instantiation: the argument `n_bins` must be an integer"
        assert n_bins >= 1, "BioSim1D() instantiation: the number of bins must be at least 1"

        assert chem_data is not None or reaction_handler is not None, \
            "BioSim1D() instantiation: at least one of the arguments `chem_data` or `reaction_handler` must be set"
            # TODO: maybe drop this requirement?  And then set the system matrix later on?

        if chem_data:
            self.chem_data = chem_data
        else:
            self.chem_data = reaction_handler.chem_data

        if reaction_handler:
            self.reaction_dynamics = reaction_handler
        else:
            self.reaction_dynamics = UniformCompartment(chem_data=self.chem_data)

        self.reactions = self.reaction_dynamics.get_reactions()

        self.n_bins = n_bins

        self.n_species = self.chem_data.number_of_chemicals()

        assert self.n_species >= 1, \
            "BioSim1D() instantiation: At least 1 chemical species must be declared prior to instantiating class"

        # Initialize all bin concentrations to zero
        self.system = np.zeros((self.n_species, n_bins), dtype=float)

        self.system_time = 0             # "Start the clock"





    #####################################################################################################

    '''                                    ~   VIEW/UPDATE SYSTEM   ~                                           '''

    def ________VIEW_UPDATE_SYSTEM________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def system_size(self) -> int:
        """
        Return the number of bins in the system
        Note: the bin numbers will range between 0 and (system_size - 1)

        :return:    The number of bins in the system
        """
        return self.n_bins



    def get_chem_data(self):
        """
        Return all the associated chemical data,
        incl. diffusion rate constants (but EXCLUSIVE of reactions)

        :return:    An Object of type "ChemData"
        """
        return self.chem_data



    def get_reactions(self):
        """
        Return all the associated reactions

        :return:    Object ot type "Reactions" (with data about all the reactions)
        """
        return self.reactions



    def get_reaction_handler(self, bin_address=0):
        """
        Return the object that manages the reactions.
        Note that for now just 1 object is present;
        in the future, it might be 1 per bin (or per bin cluster)

        :param bin_address: CURRENTLY NOT USED
        :return:            Object ot type "UniformCompartment"
        """
        return self.reaction_dynamics



    def reset_system(self) -> None:
        """
        WARNING - THIS IS VERY PARTIAL.

        :return:    None
        """
        # TODO: expand, or drop (not sure if really needed anymore)
        self.system = None
        self.system_earlier = None
        self.system_B = None
        self.membranes = None
        self.A_fraction = None
        self.chem_data = None
        self.reaction_dynamics = None
        self.history = None

        self.n_bins = 0
        self.system_length = None
        self.global_Dx = 1
        self.system_time = 0

        self.delta_reactions = None
        self.delta_reactions_B = None



    def save_system(self) -> dict:
        """
        For now, just return a copy of self.system, with a "frozen" snapshot of the current system state

        :return:    A dict of (for now a part of) the elements needed to later restore the complete system state
        """
        # TODO: deal with membranes, and anything else needed to later restore the complete system state

        return {"system": self.system.copy(), "system_time": self.system_time}



    def restore_system(self, new_state: dict) -> None:
        """
        Replace (some, for now, of) the various parts the System's internal state.
        For details of the data structure, see the class variable "system"

        :param new_state:   Numpy array containing the desired new System's internal state
        :return:
        """
        #TODO: membranes aren't yet managed. System length and global_Dx are currently not modified

        self.system = new_state["system"]
        self.n_species, self.n_bins = self.system.shape     # Extract from the new state
        assert self.n_species == self.chem_data.number_of_chemicals(), \
            f"restore_system(): inconsistency in the number of chemical species in the specified new state ({self.n_species}) " \
            f"vs. what's stored in the `Chemicals` object ({self.chem_data.number_of_chemicals()})"

        self.system_time = new_state["system_time"]




    def replace_system(self, new_state: np.array) -> None:
        """
        Replace the System's internal state.
        For details of the data structure, see the class variable "system"
        IMPORTANT: membranes aren't handled. System length and global_Dx are currently not modified

        :param new_state:   Numpy array containing the desired new System's internal state
        :return:
        """
        self.system = new_state
        self.n_species, self.n_bins = new_state.shape     # Extract from the new state
        assert self.n_species == self.chem_data.number_of_chemicals(), \
            "replace_system(): inconsistency in the number of chemical species vs. what's stored in the `Chemicals` object"



    def system_snapshot_arr(self, chem_label=None, chem_index=None) -> np.ndarray:
        """
        Return a snapshot of all the concentrations of the given chemical species,
        across ALL BINS, as a 1D Numpy array.
        If a Pandas dataframe is desired, use system_snapshot()

        :param chem_label:  String with the label to identify the chemical of interest
        :param chem_index:  Integer to identify the chemical of interest.  Cannot specify both chem_label and chem_index
        :return:            A 1-D Numpy array of concentration values along bin coordinates
        """
        assert (chem_label is not None) or (chem_index is not None), \
            "system_snapshot_arr(): at least one of the args `chem_label` or `chem_index` must be provided"

        if chem_label is not None:
            assert chem_index is None, "system_snapshot_arr(): cannot pass both arguments `chem_label` and `chem_index`"
            chem_index = self.chem_data.get_index(chem_label)
        else:
            assert chem_index is not None, "system_snapshot_arr(): must pass one of the arguments `chem_label` or `chem_index`"
            self.chem_data.assert_valid_species_index(chem_index)

        arr = self.system[chem_index]    # A 1-D Numpy array with the chemical data along bin coordinates

        return arr



    def system_snapshot(self) -> pd.DataFrame:
        """
        Return a snapshot of all the concentrations of all the species, across all bins
        as a Pandas dataframe

        :return:    A Pandas dataframe: each row is a bin,
                        and each column a chemical species
        """
        #TODO: make allowance for membranes
        all_chem_names = self.chem_data.get_all_labels()
        if self.system is None:
            return pd.DataFrame(columns = all_chem_names)   # Empty dataframe

        matrix = self.system.T

        df = pd.DataFrame(matrix, columns = all_chem_names)

        return df




    #########################################################################
    #                                                                       #
    #                SET/MODIFY CONCENTRATIONS (or membranes)               #
    #                                                                       #
    #########################################################################

    def set_uniform_concentration(self, conc: float, species_index=None, species_name=None) -> None:
        """
        Assign the given concentration to all the bins of the specified species (identified by its index or name.)
        Any previous values get over-written

        :param conc:            The desired value of chemical concentration for the above species
        :param species_index:   Zero-based index to identify a specific chemical species
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :return:                None
        """
        if species_name is not None:
            species_index = self.chem_data.get_index(species_name)
        else:
            self.chem_data.assert_valid_species_index(species_index)

        assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"

        self.system[species_index] = np.full(self.n_bins, conc, dtype=float)

        if self.uses_membranes():
            # If membranes are present, also set "side B" of the bins with membranes, to the same concentration
            # TODO: see if it can be done as vector operation
            for bin_address in self.bins_with_membranes():
                if self.membranes[bin_address]:
                    self.system_B[species_index, bin_address] = conc




    def set_all_uniform_concentrations(self, conc_list: Union[list, tuple], snapshot=True) -> None:
        """
        Set the concentrations of all chemical species at once, uniformly across all bins

        :param conc_list:   List or tuple of concentration values for each of the chemical species
        :param snapshot:    [OPTIONAL] If True (default), add to the history
                                a snapshot of this state being set
        :return:            None
        """
        assert len(conc_list) == self.chem_data.number_of_chemicals(), \
            f"set_all_uniform_concentrations(): the argument must be a list or tuple of size {self.chem_data.number_of_chemicals()}"

        for i, conc in enumerate(conc_list):
            self.set_uniform_concentration(species_index=i, conc=conc)

        if snapshot:
            self.capture_snapshot(caption="Set concentration")  # Save this operation in the history (if enabled)




    def set_bin_conc(self, bin_address: int, conc: float, species_index=None, species_name=None,
                     across_membrane=False, both_sides=False) -> None:
        """
        Assign the requested concentration value to the given bin, for the specified chemical species.
        Optionally, set the value for the "alternate bin" ("other side" of the membrane)

        :param bin_address:     The zero-based bin number of the desired compartment
        :param conc:            The desired concentration value to assign to the specified location
        :param species_index:   Zero-based index to identify a specific chemical species
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :param across_membrane: It True, consider the "other side" of the bin, i.e. the portion across the membrane
        :param both_sides:      If True, set the "regular" bin and the "other side" as well
        :return:                None
        """
        self.assert_valid_bin(bin_address)

        assert conc >= 0., \
            f"set_bin_conc(): the concentration must be a positive number or zero (the requested value was {conc})"


        if species_name is not None:
            species_index = self.chem_data.get_index(species_name)
        else:
            self.chem_data.assert_valid_species_index(species_index)


        if across_membrane or both_sides:
            assert self.system_B is not None, \
                "set_bin_conc(): the `other_side` option cannot be used unless membranes are first set"
            self.system_B[species_index, bin_address] = conc

        if both_sides or (not across_membrane):
            self.system[species_index, bin_address] = conc



    def set_species_conc(self, conc_list: Union[list, tuple, np.ndarray], species_index=None, species_name=None) -> None:
        """
        Assign the requested list of concentration values to all the bins, in bin order, for the single specified species.

        :param conc_list:       A list, tuple or Numpy array with the desired concentration values
                                    to assign to all the bins.
                                    The dimensions must match the system's dimensions.
        :param species_index:   Zero-based index to identify a specific chemical species
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :return:                None
        """
        if species_name is not None:
            # If the chemical is being identified by name, look up its index
            species_index = self.chem_data.get_index(species_name)
        elif species_index is None:
            raise Exception("BioSim1D.set_species_conc(): must provide a `species_name` or `species_index`")
        else:
            self.chem_data.assert_valid_species_index(species_index)

        assert (type(conc_list) == list) or (type(conc_list) == tuple) or (type(conc_list) == np.ndarray), \
            f"BioSim1D.set_species_conc(): the argument `conc_list` must be a list, tuple or Numpy array; " \
            f"the passed value was of type {type(conc_list)})"

        assert len(conc_list) == self.n_bins, \
            f"BioSim1D.set_species_conc(): the argument `conc_list` must be a list of concentration values for ALL the bins " \
            f"(the length should be {self.n_bins}, rather than {len(conc_list)})"

        # Verify that none of the concentrations are negative
        assert min(conc_list) >= 0, \
            f"BioSim1D.set_species_conc(): concentrations cannot be negative (values like {min(conc_list)} aren't permissible)"

        # Update the system state
        self.system[species_index] = conc_list



    def inject_conc_to_bin(self, bin_address: int, delta_conc: float, chem_label=None, chem_index=None, zero_clip = False) -> None:
        """
        Add the requested concentration to the cell with the given address, for the specified chem species

        :param bin_address:     The zero-based bin number of the desired cell
        :param chem_label:
        :param chem_index:      Zero-based index to identify a specific chemical species
        :param delta_conc:      The concentration to add to the specified location
        :param zero_clip:       If True, any requested increment causing a concentration dip below zero, will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        self.assert_valid_bin(bin_address)

        # TODO: turn into a utility private method
        assert (chem_label is not None) or (chem_index is not None), \
            "inject_conc_to_bin(): at least one of the args `chem_label` or `chem_index` must be provided"
        if chem_label is not None:
            assert chem_index is None, "inject_conc_to_bin(): cannot pass both arguments `chem_label` and `chem_index`"
            chem_index = self.chem_data.get_index(chem_label)
        else:
            assert chem_index is not None, "inject_conc_to_bin(): must pass one of the arguments `chem_label` or `chem_index`"
            self.chem_data.assert_valid_species_index(chem_index)


        if (self.system[chem_index, bin_address] + delta_conc) < 0. :
            # Take special action if the requested change would make the bin concentration negative
            if zero_clip:
                self.system[chem_index, bin_address] = 0
                return
            else:
                raise Exception("inject_conc_to_bin(): The requested concentration change would result in a negative final value")

        # Normal scenario, not leading to negative values for the final concentration
        self.system[chem_index, bin_address] += delta_conc



    def inject_gradient(self, species_name, conc_left = 0., conc_right = 0.) -> None:
        """
        Add to the concentrations of the specified chemical species a linear gradient spanning across all bins,
        with the indicated values at the endpoints of the system.

        :param species_name:    The name of the chemical species whose concentration we're modifying
        :param conc_left:       The desired amount of concentration to add to the leftmost bin (the start of the gradient)
        :param conc_right:      The desired amount of concentration to add to the rightmost bin (the end of the gradient)
        :return:                None
        """
        assert conc_left >= 0. and conc_right >= 0., \
                    f"BioSim1D.inject_gradient(): the concentration values cannot be negative"

        assert self.n_bins > 1, \
                    f"BioSim1D.inject_gradient(): minimum system size must be 2 bins"

        species_index = self.chem_data.get_index(species_name)

        # Create an array of equally-spaced values from conc_left to conc_right
        # Size of array is same as the number of bins in the system
        increments_arr = np.linspace(conc_left, conc_right, self.n_bins)

        # Update the concentrations across all bins, for the requested chemical species
        # Note: all the values we're adding are non-negative; so, no danger of creating negative concentrations
        self.system[species_index] += increments_arr



    def inject_sine_conc(self, species_name, frequency, amplitude, bias=0, phase=0, zero_clip = False) -> None:
        """
        Add to the concentrations of the specified chemical species a sinusoidal signal across all bins.

        Note:   A sine wave of the form  f(x) = A sin(B x - C)
                has an amplitude of A, a period of 2Pi/B and a right phase shift of C (in radians)

        In Mathematica:  Plot[Sin[B x - C] /. {B -> 2 Pi, C -> 0} , {x, 0, 1}, GridLines -> Automatic]

        :param species_name:    The name of the chemical species whose concentration we're modifying
        :param frequency:       Number of waves along the length of the system
        :param amplitude:       Amplitude of the Sine wave.  Note that peak-to-peak values are double the amplitude
        :param bias:            Amount to be added to all values (akin to "DC bias" in electrical circuits)
        :param phase:           In degrees: phase shift to the RIGHT.  EXAMPLE: 180 to flip the Sine curve
        :param zero_clip:       If True, any requested change causing a concentration dip below zero,
                                will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        species_index = self.chem_data.get_index(species_name)

        period = self.n_bins / frequency
        #print("period: ", period)

        phase_radians = phase * math.pi / 180

        B = 2*math.pi / period
        #print("B: ", B)

        for x in range(self.n_bins):
            # Update each bin concentration in turn
            conc = amplitude * math.sin(B*x - phase_radians) + bias
            self.inject_conc_to_bin(bin_address = x, chem_index = species_index,
                                   delta_conc = conc, zero_clip = zero_clip)



    def inject_bell_curve(self, species_name, mean=0.5, sd=0.15, amplitude=1., bias=0) -> None:
        """
        Add to the concentrations of the specified chemical species a signal across all bins in the shape of a Bell curve.
        The default values provide bell shape centered in the middle of the system, and fairly spread out
        (but pretty close to zero at the endpoints)

        :param species_name:    The name of the chemical species whose concentration we're modifying
        :param mean:            A value, generally between 0 and 1, indication the position of the mean relative to the system;
                                    if less than 0 or greater than 1, only one tail of the curve will be seen
        :param sd:              Standard deviation, in units of the system length
        :param amplitude:       Amount by which to multiply the signal
        :param bias:            Positive amount to be added to all values (akin to "DC bias" in electrical circuits)
        :return:                None
        """
        assert bias >= 0, \
            f"BioSim1D.inject_bell_curve(): the value for the `bias` ({bias}) cannot be negative"

        assert amplitude >= 0, \
            f"BioSim1D.inject_bell_curve(): the value for the `amplitude` ({amplitude}) cannot be negative"

        species_index = self.chem_data.get_index(species_name)

        # Create an array of equally-spaced values from 0. to 1.
        # Size of array is same as the number of bins in the system
        x = np.linspace(0, 1, self.n_bins)

        # Create a range of y-values that correspond to normal pdf with the given mean and SD
        increments_arr = norm.pdf(x, mean, sd)

        # Stretch and shift the normal pdf, as requested;
        # use those values to update the concentrations across all bins,
        # for the requested chemical species
        # Note: all the values we're adding are non-negative; so, no danger of creating negative concentrations
        self.system[species_index] += amplitude * increments_arr + bias




    ########  MEMBRANE-RELATED  ################


    def uses_membranes(self) -> bool:
        """
        Return True if membranes are part of the system

        :return:
        """
        return self.membranes is not None



    def bins_with_membranes(self) -> [int]:
        """

        :return:
        """
        # Create and return a list of bin numbers whose entry in self.membranes is True
        return [bin_number for bin_number in range(self.n_bins)
                if self.membranes[bin_number]]



    def set_membranes(self, membrane_pos: Union[List, Tuple]) -> None:
        """
        DEPRECATED!

        Set the presence of all membranes in the system,
        and optionally specify the fraction of the "A" part of the bin (by default 0.5)

        Initialize the class variables "membranes" and "A_fraction"

        EXAMPLES:   set_membranes([4, 20, 23])
                    set_membranes([4, (20, 0.7), 23])

        IMPORTANT: any previously-set membrane information is lost.

        :param membrane_pos:    A list or tuple of:
                                    1) EITHER indexes of bins that contain membranes
                                    2) OR pairs of bins numbers and fractional values
        :return:                None
        """
        print("**** DEPRECATED, and expected to be removed soon")

        self.membranes = np.zeros(self.n_bins, dtype=bool)
        self.A_fraction = np.zeros(self.n_bins, dtype=float)
        if self.system_B is None:
            self.system_B = np.zeros((self.n_species, self.n_bins), dtype=float)


        # Loop over all the values passed as a list or tuple
        for bin_info in membrane_pos:
            if type(bin_info) == int:
                bin_number = bin_info
                fraction = 0.5
            elif (type(bin_info) == tuple) or (type(bin_info) == list):
                assert len(bin_info) == 2, \
                    "set_membranes(): tuples or lists must contain exactly 2 elements"
                bin_number, fraction = bin_info

                assert (type(fraction) == float) or (type(fraction) == int), \
                    "set_membranes(): requested fractions must be floats or integers"
                assert 0. <= fraction <= 1., \
                    "set_membranes(): requested fractions must be in the [0.-1.0] range, inclusive"
            else:
                raise Exception("set_membranes(): the elements of the argument `membrane_pos` must be tuples or lists")

            self.assert_valid_bin(bin_number)
            self.membranes[bin_number] = True
            self.A_fraction[bin_number] = fraction

            # Equalize all the concentrations across the compartments on both sides of the membrane
            if self.chem_data:
                for chem_index in range(self.chem_data.number_of_chemicals()):
                    # TODO: do as vector operation
                    conc = self.bin_concentration(bin_address=bin_number, species_index=chem_index)
                    self.set_bin_conc(bin_address=bin_number, species_index=chem_index, conc=conc, across_membrane=True)



    def assert_valid_bin(self, bin_address: int) -> None:
        """
        Raise an Exception if the given bin number isn't valid

        :param bin_address: An integer that ought to be between 0 and (self.n_bins-1), inclusive
        :return:            None
        """
        assert type(bin_address) == int, \
            f"BioSim1D: the requested bin address ({bin_address}) is not an integer; its type is {type(bin_address)}"

        if bin_address < 0 or bin_address >= self.n_bins:
            raise Exception(f"BioSim1D: the requested bin address ({bin_address}) is out of bounds for the system size; "
                            f"allowed range is [0-{self.n_bins-1}], inclusive")



    def lookup_species(self, species_index=None, species_name=None, trans_membrane=False, copy=False) -> np.array:
        """
        Return the NumPy array of concentration values across the all bins (from left to right)
        for the single specified chemical species.
        NOTE: what is being returned NOT a copy, unless specifically requested

        :param species_index:   The index order of the chemical species of interest
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :param trans_membrane:  If True, consider only the "other side" of the bins, i.e. the portion across the membrane
                                    (it will be zero for bins without membrane)
        :param copy:            If True, an independent numpy array will be returned: a *copy* rather than a view
        :return:                A NumPy 1-D array of concentration values across the bins (from left to right);
                                    the size of the array is the number of bins
        """
        if species_name is not None:
            species_index = self.chem_data.get_index(species_name)
        else:
            self.chem_data.assert_valid_species_index(species_index)

        if trans_membrane:
            species_conc =  self.system_B[species_index]
        else:
            species_conc = self.system[species_index]

        if copy:
            return species_conc.copy()
        else:
            return species_conc



    def bin_concentration(self, bin_address: int, species_index=None, species_label=None, trans_membrane=False) -> float:
        """
        Return the concentration at the requested bin of the specified chemical species

        :param bin_address:     The bin number
        :param species_index:   The index order of the chemical species of interest
        :param species_label:   [OPTIONAL] If provided, it over-rides the value for species_index
        :param trans_membrane:  If True, consider the "other side" of the bin, i.e. the portion across the membrane
        :return:                A concentration value at the indicated bin, for the requested species
        """
        if species_label is not None:
            species_index = self.chem_data.get_index(species_label)

        self.chem_data.assert_valid_species_index(species_index)

        if trans_membrane:
            return self.system_B[species_index, bin_address]
        else:
            return self.system[species_index, bin_address]



    def bin_snapshot(self, bin_address: int) -> dict:
        """
        Extract the concentrations of all the chemical species at the specified bin,
        as a dict whose keys are the names of the species
        EXAMPLE:  {'A': 10.0, 'B': 50.0}

        :param bin_address: An integer with the bin number
        :return:            A dict of concentration values; the keys are the names of the species
        """
        self.assert_valid_bin(bin_address)

        d = {}
        for species_index in range(self.n_species):
            name = self.chem_data.get_label(species_index)
            conc = self.bin_concentration(bin_address, species_index)
            d[name] = conc

        return d



    def bin_snapshot_array(self, bin_address: int, trans_membrane=False) -> np.array:
        """
        Extract the concentrations of all the chemical species at the specified bin,
        as a Numpy array in the index order of the species
        EXAMPLE: np.array([10., 50.)]

        :param bin_address:     An integer with the bin number
        :param trans_membrane:  If True, consider the "other side" of the bin, i.e. the portion across the membrane
        :return:                A Numpy array  of concentration values, in the index order of the species
        """
        self.assert_valid_bin(bin_address)

        if trans_membrane:
            return self.system_B[:, bin_address]
        else:
            return self.system[: , bin_address]



    def show_system_snapshot(self) -> pd.DataFrame:
        """
        Print a header, and return a dataframe

        :return:    A Pandas dataframe
        """
        print(f"SYSTEM SNAPSHOT at time {self.system_time:,.8g}:")
        return self.system_snapshot()



    def selected_concentrations(self, bins=None, chem_labels=None) -> dict:
        """
        Extract and return the concentration values of one or more (use None for all) chemicals,
        in one or more bins.
        The value is returned as a dictionary where the keys are bin addresses, and the values are dicts of
        concentration values for the various chemicals (identified by their labels)
            EXAMPLE:
                    {   5: {"A": 1.3, "B": 3.9},
                        8: {"A": 4.6, "B": 2.7}
                    }

            TODO: alternate returned structure - a Pandas dataframe such as
                    BIN ADDRESS      A       B
                        5           1.3     3.9
                        8           4.6     2.7

        :param bins:        Bin address (integer), or list of bin addresses. Use None to indicate all
        :param chem_labels: Chemical label, or list of labels. Use None to indicate all
        :return:            A dict indexed by bin address
        """
        if bins is None:
            bins = list(range(self.n_bins))  # All bins
        elif type(bins) != list:
            bins = [bins]

        if chem_labels is None:
            chem_labels = self.chem_data.get_all_labels()
        elif type(chem_labels) != list:
            chem_labels = [chem_labels]

        result = {}
        for bin_address in bins:
            self.assert_valid_bin(bin_address)
            bin_values = {}
            for chem_label in chem_labels:
                conc = self.bin_concentration(bin_address=bin_address, species_label=chem_label)
                bin_values[chem_label] = conc

            result[bin_address] = bin_values

        return result



    def describe_state(self, concise=False) -> Union[pd.DataFrame, None]:
        """
        A simple printout of the state of the system, for now useful only for small systems.

        EXAMPLE (concise):
            SYSTEM STATE at Time t = 0:
            [[0. 0. 0. 0.]
             [0. 0. 0. 0.]]

        EXAMPLE (not concise):
            SYSTEM STATE at Time t = 0:
            4 bins and 2 chemical species:
               <PANDAS data frame returned>

        TODO: The goal for ASCII-based printouts involving membranes is something like.
           ____________________
        A: |20|18|12|8()  2| 6|    Diff rate: 0.2 (M: 0.01)
        B: |30| 2|56|4()3.5|12|    Diff rate: 0.8 (M: 0.3)
           --------------------
             0  1  2      3  4

        :param concise: If True, only produce a minimalist printout with just the concentration values
        :return:        None, if membranes are present, or a Pandas dataframe otherwise
        """
        print(f"SYSTEM STATE at Time t = {self.system_time:,.8g}:")

        if concise:             # A minimalist printout...
            print(self.system)   # ...only showing the concentration data (a Numpy array)
            if self.uses_membranes():
                print("Membranes:")
                print(self.system_B)
            return

        # If we get thus far, it's a FULL printout

        print(f"{self.n_bins} bins and {self.n_species} chemical species:")

        if self.uses_membranes():
            # Show a line of line of data for each chemical species in turn
            for species_index in range(self.n_species):
                name = self.chem_data.get_label(species_index)
                if name:    # If a name was provided, show it
                    name = f" ({name})"
                else:
                    name = ""

                if self.membranes is None:
                    all_conc = self.system[species_index]
                else:
                    all_conc = "|"
                    for bin_no in range(self.n_bins):
                        all_conc += str(self.bin_concentration(bin_no, species_index=species_index))
                        if self.membranes[bin_no]:
                            # Add a symbol for the membrane, and the additional membrane concentration data (on the "other side")
                            all_conc += "()"    # To indicate a membrane
                            all_conc += str(self.bin_concentration(bin_no, species_index=species_index, trans_membrane=True))
                        all_conc += "|"

                if not self.chem_data.get_all_diffusion_rates():
                    print(f"  Species {species_index}{name}. Diff rate: NOT SET. Conc: {all_conc}")
                else:
                    print(f"  Species {species_index}{name}. Diff rate: {self.chem_data.get_diffusion_rate(species_index=species_index)}. Conc: {all_conc}")

        else:   # This alternate view, using Pandas, currently only usable when there are no membranes
            df = pd.DataFrame(self.system, columns=[f"Bin {i}" for i in range(self.n_bins)])

            df.insert(0, "Species", self.chem_data.get_all_labels())
            df.insert(1, "Diff rate", self.chem_data.get_all_diffusion_rates()) # Unset values will show up as None

            return df





    def show_membranes(self, n_decimals=1) -> str:
        """
        DEPRECATED!

        A simple-minded early method to visualize where the membranes are.
        Print, and return, a string with a diagram to visualize membranes and the fractions
        of their "left" sides

        EXAMPLE (with 2 membranes on the right part of a 5-bin system):
                _____________________
                |   |   |0.8|   |0.3|
                ---------------------

        :param n_decimals:  Number of decimal places to show in the fractions
        :return:            A string with the character-based diagram;
                            if no membranes were defined, return an empty string
        """
        print("**** DEPRECATED; a completely-new system is expected")

        if self.membranes is None:
            print("No membranes present.  Call set_membranes() to set them")
            return ""

        # Prepare the middle line
        box_contents = "|"
        for bin_no, val in enumerate(self.membranes):
            if val:
                fraction = round(self.A_fraction[bin_no], n_decimals)
                box_contents += str(fraction) + "|"
            else:
                box_contents += "   |"

        box_width = len(box_contents)

        box = "\n"
        box += "_" * box_width + "\n"   # The top of the box
        box += box_contents + "\n"
        box += "-" * box_width          # The bottom of the box

        print(box)
        return box




    #####################################################################################################

    '''                                    ~   SPATIAL ELEMENTS   ~                                           '''

    def ________SPATIAL_ELEMENTS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def set_dimensions(self, length) -> None:
        """
        Set the overall length of the system.
        Doing so, will permit to convert bin numbers to positional values

        :param length:
        :return:
        """
        assert (type(length) == float) or (type(length) == int), "set_dimensions(): length must be a number"
        assert length > 0, "set_dimensions(): length must be positive"

        self.system_length = length
        self.global_Dx = length / (self.n_bins-1)



    def x_coord(self, bin_address):
        """
        Return the x coordinate of the middle of the specified bin.
        By convention, for the leftmost bin, it's zero,
        and for the rightmost, it's the overall length of the system
        """
        assert self.system_length, "x_coord(): must first call set_dimensions()"

        return bin_address * self.global_Dx



    def increase_spatial_resolution(self, factor:int) -> None:
        """
        Increase the spatial resolution of the system by cloning and repeating
        each bin, by the specified number of times.
        Replace the System's internal state
        (note that the number of bins will increase by the requested factor)

        EXAMPLE: if the (2-chemicals) system is
                        [[11. 12. 13.]
                         [ 5. 15. 25.]]
                and factor=2, then the result will be
                        [[11. 11. 12. 12. 13. 13.]
                         [ 5.  5. 15. 15. 25. 25.]]

        :param factor:  Number of bins into which to split each bin (replicating their concentration values)
        :return:        None
        """
        assert type(factor) == int, "The argument `factor` must be an integer"
        new_state = np.repeat(self.system, factor, axis=1)
        self.replace_system(new_state)



    def double_spatial_resolution_linear(self) -> None:
        """
        Increase the spatial resolution of the system by inserting between bins
        their average value (of concentrations, for every chemical species.)
        Replace the System's internal state
        (note that if the number of bins is initially N, it'll become 2N-1).
        If the system has fewer than 2 bins, an Exception will be raised.

        EXAMPLE: if the (2-chemicals) system is
                        [[11. 12. 13.]
                         [ 5. 15. 25.]]
                then the result will be
                        [[11.  11.5 12.  12.5 13. ]
                         [ 5.  10.  15.  20.  25. ]]

        :return:    None
        """
        assert self.n_bins >= 2, \
            "double_spatial_resolution_linear(): function can only be used if the system contains at least 2 bins"

        new_state = np.zeros((self.n_species, self.n_bins * 2 - 1), dtype=float)
        #print(new_state)

        for i in range(self.n_bins-1):
            two_col = self.system[ : , i:i+2 ]
            avg_col = two_col.mean(axis=1)

            #print("two_col: ", two_col)
            #print("avg_col: ", avg_col)

            new_state[ : , 2*i] = two_col[ : , 0]
            new_state[ : , 2*i+1] = avg_col

        new_state[ : , -1] = self.system[ : , -1]

        self.replace_system(new_state)



    def decrease_spatial_resolution(self, factor:int) -> None:
        """

        TODO: eliminate the restriction that the number of bins must be a multiple of factor

        EXAMPLE: if the system is
                        [[10., 20., 30., 40., 50., 60.]
                         [ 2., 8.,   5., 15., 4.,   2.]]
                and factor=2, then the result will be
                        [[15., 35., 55.]
                         [ 5., 10.,  3.]]

        :param factor:
        :return:
        """
        assert type(factor) == int, "The argument `factor` must be an integer"
        assert self.n_bins % factor == 0, f"The number of bins (currently {self.n_bins}) must be a multiple of the requested scaling factor"

        reduced_n_bins = int(self.n_bins / factor)
        # The result matrix will have the same number of chemical species, but fewer bins
        new_state = np.zeros((self.n_species, reduced_n_bins), dtype=float)

        for i in range(reduced_n_bins):
            start_col = factor * i      # The start column will initially be 0, and will get incremented by factor
            col_group = self.system[ : , start_col:start_col+factor] # Extract a submatrix containing the number of columns specified by "factor",
                                                                    # starting with the column specified by start_col
            compressed_col_group = np.sum(col_group, axis=1, keepdims=True) / factor    # Create a single column that is the average of the columns in the group
            new_state[:, [i]] = compressed_col_group                                       # Store the newly-computed column of averages in the appropriate place
                                                                                        # in the result matrix
        self.replace_system(new_state)



    def smooth_spatial_resolution(self) -> None:
        """
        EXAMPLE: if the system is
                        [[10., 20., 30.]
                         [ 2., 8.,   4.]]
                then the result will be
                        [[10., 15., 20., 25., 30.]
                         [ 2.,  5.,  8.,  6.,  4.]]

        :return:
        """
        n_bins = self.n_bins
        new_n_bins = n_bins * 2 - 1     # The final number of bins
        new_state = np.zeros((self.n_species, new_n_bins), dtype=float)

        for start_col in range(n_bins-1):                           # The start column will be between 0 and (self.n_bins-2), inclusive
            col_group = self.system[ : , start_col:start_col+2]      # Extract a submatrix containing 2 columns,
                                                                    # starting with the one in position start_col
            avg_col = np.sum(col_group, axis=1, keepdims=True) / 2. # Create a single column that is the average of the columns in the group

            new_state[:, [2*start_col]] = self.system[ : , start_col:start_col+1]   # Set one column, from the first column in the group
            new_state[:, [2*start_col+1]] = avg_col                                # Set the next column with the newly-computed column of averages


        new_state[ : , -1:] = self.system[ : , -1:]                 # Set the last column of result to the last column of self.system

        self.replace_system(new_state)





    #####################################################################################################

    '''                                    ~   UTILITIES   ~                                          '''

    def ________UTILITIES________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def check_mass_conservation(self, expected :float, chem_label=None, chem_index=None) -> bool:
        """
        Check whether the sum of all the concentrations of the specified chemical,
        across all bins, adds up to the passed value

        :param expected:
        :param chem_label:  String with the label to identify the chemical of interest
        :param chem_index:  Integer to identify the chemical of interest.  Cannot specify both chem_label and chem_index

        :return:
        """
        arr = self.system_snapshot_arr(chem_label=chem_label, chem_index=chem_index)
        total = np.sum(arr)
        return np.allclose(expected, total)



    def reaction_in_equilibrium(self, bin_address: int, rxn_index :int, tolerance=1, explain=True) -> bool:
        """
        Ascertain whether the given system concentrations are in equilibrium at the given bin,
        for the specified reactions (by default, check all reactions)

        :param bin_address: The zero-based bin number of the desired compartment
        :param rxn_index:   The integer index (0-based) to identify the reaction of interest;
                                if None, then check all the reactions
        :param tolerance:   Allowable relative tolerance, as a PERCENTAGE,
                                to establish satisfactory match with expected values
        :param explain:     If True, print out details about the analysis,
                                incl. the formula(s) being used to check the equilibrium
                                EXAMPLES:   "([C][D]) / ([A][B])"
                                            "[B] / [A]^2"

        :return:            Return True if ALL the reactions are close enough to an equilibrium,
                                as allowed by the requested tolerance;
                                otherwise, return a dict of the form {False: [list of reaction indexes]}
                                for all the reactions that failed the criterion
                                EXAMPLE:  {False: [3, 6]}
        """
        #TODO: put together a version for 2D
        return self.get_reaction_handler().is_in_equilibrium(conc=self.bin_snapshot(bin_address),
                                                             rxn_index=rxn_index, tolerance=tolerance, explain=explain)






    #####################################################################################################

    '''                                      ~   HISTORY   ~                                          '''

    def ________HISTORY________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def enable_history(self, bins=None, frequency=1, chem_labels=None, take_snapshot=False) -> None:
        """
        Request history capture, with the specified parameters.
        If history was already enabled, this function can be used to alter its capture parameters.

        :param bins:            Bin address (integer), or list of bin addresses. Use None to indicate all
        :param frequency:
        :param chem_labels:     [OPTIONAL] List of chemicals to include in the history;
                                    if None (default), include them all.
        :param take_snapshot:   If True, a snapshot of the system's current configuration is added to the history

        :return:                None
        """
        # Make sure that all bin addresses, if specified, are valid
        if bins:
            for b in bins:
                self.assert_valid_bin(b)

        self.conc_history.enable_history(frequency=frequency, chem_labels=chem_labels, bins=bins)
        if take_snapshot:
            self.capture_snapshot()

        print(f"History enabled for bins {bins} and chemicals {chem_labels} (None means 'all')")



    def capture_snapshot(self, step_count=None, caption="") -> None:
        """

        :param step_count:
        :param caption:
        :return:            None
        """
        if not self.conc_history.to_capture(step_count):
            return

        data_snapshot = self.selected_concentrations(bins=self.conc_history.restrict_bins,
                                                     chem_labels=self.conc_history.restrict_chemicals)
        '''
           EXAMPLE of data_snapshot:
                { 5: {"A": 1.3, "B": 3.9},
                  8: {"A": 4.6, "B": 2.7}
                }        
        '''
        self.conc_history.save_snapshot(step_count=step_count, system_time=self.system_time,
                                        data_snapshot=data_snapshot,
                                        caption=caption)



    def get_bin_history(self, bin_address :int) -> pd.DataFrame:
        """
        Get the concentration history at the given bin(s) of all the chemicals
        whose history was requested by a call of enable_history()

        :param bin_address: A single bin address (an integer)
        :return:            A Pandas data frame
        """
        return self.conc_history.bin_history(bin_address = bin_address)



    def add_snapshot(self, data_snapshot: dict, caption ="") -> None:
        """
        TODO: being phased out
        Preserve some data value (passed as dictionary) in the history, linked to the
        current System Time.

        EXAMPLE:  add_snapshot(data_snapshot = {"concentration_A": 12.5, "concentration_B": 3.7},
                                caption="Just prior to infusion")

        IMPORTANT: if the data is not immutable, then it ought to be cloned first,
                   to preserve it from possible later modifications

        :param data_snapshot:   A dictionary of data to preserve for later use
        :param caption:         Optional caption to attach to this preserved data
        :return:                None
        """
        self.history.store(par=self.system_time,
                           data_snapshot = data_snapshot, caption=caption)



    def get_history(self) -> pd.DataFrame:
        """
        TODO: being phased out
        Retrieve and return a Pandas dataframe with the system history that had been saved
        using add_snapshot()

        :return:        a Pandas dataframe
        """
        return self.history.get_dataframe()





    #####################################################################################################

    '''                                    ~   SIMULATIONS   ~                                           '''

    def ________SIMULATIONS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def diffuse(self, total_duration=None, time_step=None, n_steps=None, delta_x=1, algorithm=None) -> dict:
        """
        Uniform-step diffusion, with 2 out of 3 criteria specified:
            1) until reaching, or just exceeding, the desired time duration
            2) using the given time step
            3) carrying out the specified number of steps

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param delta_x:         Distance between consecutive bins
        :param algorithm:          (Optional) code specifying the method to use to solve the diffusion equation.
                                    Currently available options: "5_1_explicit"
        :return:                A dictionary with data about the status of the operation
                                    (for now, just the number of steps run; key: "steps")
        """
        time_step, n_steps = self.reaction_dynamics.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)
        for i in range(n_steps):
            if self.debug:
                if (i < 2) or (i >= n_steps-2):
                    print(f"    Performing diffusion step {i}...")
                elif i == 2:
                    print("    ...")

            self.diffuse_step(time_step, delta_x=delta_x, algorithm=algorithm)
            self.system += self.delta_diffusion     # Matrix operation to update all the concentrations
            self.system_time += time_step
            self.capture_snapshot(step_count=i)     # Save historical values (if enabled)

        if self.debug:
            print(f"\nSystem after Delta time {total_duration}, at end of {n_steps} steps of size {time_step}:")
            self.describe_state(concise=True)
            print()

        status = {"steps": n_steps}
        return status



    def diffuse_step(self, time_step, delta_x=1, algorithm=None) -> None:
        """
        Diffuse all the species for the given time step, across all bins;
        clear the delta_diffusion array, and then re-compute it from all the species.

        IMPORTANT: the actual system concentrations are NOT changed.

        :param time_step:   Time step over which to carry out the diffusion
                            If too large - as determined by the method is_excessive() - an Exception will be raised
        :param delta_x:     Distance between consecutive bins
        :param algorithm:      (Optional) code specifying the method to use to solve the diffusion equation.
                                Currently available options: "5_1_explicit"
        :return:            None (the array in the class variable "delta_diffusion" gets set)
        """
        # TODO: parallelize the independent computations

        # 2-D array of incremental changes at every bin, for each chemical species
        self.delta_diffusion = np.zeros((self.n_species, self.n_bins), dtype=float)

        for species_index in range(self.n_species):

            if algorithm is None:
                increment_vector = self.diffuse_step_single_species(time_step, species_index=species_index, delta_x=delta_x)
            elif algorithm == "5_1_explicit":
                increment_vector = self.diffuse_step_single_species_5_1_stencils(time_step, species_index=species_index, delta_x=delta_x)
            else:
                raise Exception(f"diffuse_step(): unknown method: `{algorithm}`")

            #print("Increment vector is: ", increment_vector)

            # For each bin, update the concentrations from the buffered increments
            self.delta_diffusion[species_index] = increment_vector      # Vector operation to a row of the matrix delta_diffusion



    def diffuse_step_single_species(self, time_step: float, species_index=0, delta_x=1) -> np.array:
        """
        Diffuse the specified single chemical species, for the given time step, across all bins,
        and return a 1-D array of the changes in concentration ("Delta concentration")
        for the given species across all bins.

        IMPORTANT: the actual system concentrations are NOT changed.

        We're assuming an isolated environment, with nothing diffusing thru the "walls"

        This approach is based on a "3+1 stencil", aka "Explicit Forward-Time Centered Space".
        EXPLANATION:  https://life123.science/diffusion

        :param time_step:       Delta time over which to carry out this single diffusion step;
                                    if too large, an Exception will be raised.
        :param species_index:   ID (in the form of an integer index) of the chemical species under consideration
        :param delta_x:         Distance between consecutive bins

        :return:                A 1-D Numpy array with the CHANGE in concentration for the given species across all bins
        """
        assert self.system is not None, "BioSim1D.diffuse_step_single_species(): Must first initialize the system"
        assert not self.chem_data.missing_diffusion_rate(), "BioSim1D.diffuse_step_single_species(): Must first set the diffusion rates"
        assert self.sealed == True, "BioSim1D.diffuse_step_single_species(): For now, there's no provision for exchange with the outside"

        increment_vector = np.zeros(self.n_bins, dtype=float)       # One element per bin

        if self.n_bins == 1:
            return increment_vector                                 # There's nothing to do in the case of just 1 bin!

        diff = self.chem_data.get_diffusion_rate(species_index=species_index)     # The diffusion rate of the specified single species

        assert not self.is_excessive(time_step, diff, delta_x), \
            f"diffuse_step_single_species(): Excessive large time_step ({time_step}). Should be < {self.max_time_step(diff, delta_x)}"


        # Carry out a 1-D convolution operation, with a tile of size 3 (or 2 if only 2 bins)
        #print(f"Diffusing species # {species_index}")

        max_bin_number = self.n_bins - 1    # Bin numbers range from 0 to max_bin_number, inclusive

        effective_diff = diff * time_step   # We're calling this quantity "Effective Diffusion" (NOT a standard term)
        if delta_x != 1:
            effective_diff /= (delta_x**2)


        for i in range(self.n_bins):    # Bin number, ranging from 0 to max_bin_number, inclusive
            #print(f"Processing bin number {i}")
            current_conc = self.system[species_index , i]   # Concentration in the center of the convolution tile

            if i == 0 :                     # Special case for the first bin (no left neighbor)
                increment_vector[i] = effective_diff * (self.system[species_index , 1] - current_conc)
            elif i == max_bin_number :      # Special case for the last bin (no right neighbor)
                increment_vector[i] = effective_diff * (self.system[species_index , i - 1] - current_conc)
            else:
                increment_vector[i] = effective_diff * \
                                        (self.system[species_index , i + 1] - current_conc
                                         + self.system[species_index , i - 1] - current_conc)

        return increment_vector



    def diffuse_step_single_species_5_1_stencils(self, time_step: float, species_index=0, delta_x=1) -> np.array:
        """
        Similar to diffuse_step_single_species(), but using a "5+1 stencil";
        i.e. spatial derivatives are turned into finite elements using 5 adjacent bins instead of 3.

        For more info, see diffuse_step_single_species()

        IMPORTANT: the actual system concentrations are NOT changed.

        :param time_step:       Delta time over which to carry out this single diffusion step;
                                    if too large, an Exception will be raised.
        :param species_index:   ID (in the form of an integer index) of the chemical species under consideration
        :param delta_x:         Distance between consecutive bins

        :return:                A Numpy array with the CHANGE in concentration for the given species across all bins
        """
        assert self.system is not None, "Must first initialize the system"
        assert self.n_bins > 0, "Must first set the number of bins"
        assert self.n_bins >= 3, "For very small number of bins, use another method"
        assert not self.chem_data.missing_diffusion_rate(), "Must first set the diffusion rates"
        assert self.sealed == True, "For now, there's no provision for exchange with the outside"


        increment_vector = np.zeros(self.n_bins, dtype=float)   # One element per bin

        if self.n_bins == 1:
            return increment_vector                             # There's nothing to do in the case of just 1 bin!

        diff = self.chem_data.get_diffusion_rate(species_index=species_index)     # The diffusion rate of the specified single species

        # TODO: this Upper Bound is based on a *different* method, and should be made more specific to this method
        #assert not self.is_excessive(time_step, diff, delta_x), \
            #f"Excessive large time_fraction. Should be < {self.max_time_step(diff, delta_x)}"


        # Carry out a 1-D convolution operation, with a tile of size 5

        max_bin_number = self.n_bins - 1     # Bin numbers range from 0 to max_bin_number, inclusive

        effective_diff = diff * time_step   # We're calling this quantity "Effective Diffusion" (NOT a standard term)
        if delta_x != 1:
            effective_diff /= (delta_x**2)

        #print("effective_diff: ", effective_diff)

        # The coefficients for the "Central Differences" for the spatial 2nd partial derivative,
        #   to "accuracy 4" (using 5 term: 2 left neighbors and 2 right neighbors)
        C2 = -1/12
        C1 = 4/3
        C0 = - 2.5

        leftmost = self.system[species_index , 0]
        rightmost = self.system[species_index , max_bin_number]

        for i in range(self.n_bins):    # Bin number, ranging from 0 to max_bin_number, inclusive
            #print(f"Processing bin number {i}")
            C_i = self.system[species_index , i]

            # The boundary conditions, at left and right edges of the system,
            # state that the flux is zero across the boundaries
            if i == 0:                     # Special cases for the first 2 bins
                C_i_minus_2 = leftmost
                C_i_minus_1 = leftmost
            elif i == 1:
                C_i_minus_2 = leftmost
                C_i_minus_1 = self.system[species_index , i-1]
            else:
                C_i_minus_2 = self.system[species_index , i-2]
                C_i_minus_1 = self.system[species_index , i-1]

            if i == max_bin_number:      # Special cases for the last 2 bins
                C_i_plus_1 = rightmost
                C_i_plus_2 = rightmost
            elif i == max_bin_number - 1:
                C_i_plus_1 = self.system[species_index , i+1]
                C_i_plus_2 = rightmost
            else:
                C_i_plus_1 = self.system[species_index , i+1]
                C_i_plus_2 = self.system[species_index , i+2]

            #print("The 5 bins under consideration: ", C_i_minus_2, C_i_minus_1, C_i, C_i_plus_1, C_i_plus_2)
            # Compute the "Central Differences" for the 2nd partial derivative, to "accuracy 4"
            increment_vector[i] = effective_diff * \
                                      (  C2 * C_i_minus_2
                                       + C1 * C_i_minus_1
                                       + C0 * C_i
                                       + C1 * C_i_plus_1
                                       + C2 * C_i_plus_2)

        return increment_vector




    def is_excessive(self, time_step, diff_rate, delta_x) -> bool:
        """
        Use a loose heuristic to determine if the requested time step is too long,
        given the diffusion rate and delta_x.
        This is also based on the "Von Neumann stability analysis"
        (an explanation can be found at: https://www.youtube.com/watch?v=QUiUGNwNNmo)

        :param time_step:
        :param diff_rate:
        :param delta_x:
        :return:
        """
        if time_step > self.max_time_step(diff_rate, delta_x):
            return True
        else:
            return False



    def max_time_step(self, diff_rate, delta_x) -> float:
        """
        Determine a reasonable upper bound on the time step, for the given diffusion rate and delta_x
        This is also based on the "Von Neumann stability analysis"
        (an explanation can be found at: https://www.youtube.com/watch?v=QUiUGNwNNmo)

        :param diff_rate:
        :param delta_x:
        :return:
        """
        return delta_x**2 * self.time_step_threshold/diff_rate



    def react(self, total_duration=None, time_step=None, n_steps=None, snapshots=None, silent=False) -> None:
        """
        Update the system concentrations as a result of all the reactions in all bins - taking
        the presence of membranes into account, if applicable.

        CAUTION : NO diffusion is performed. Use this function
                  only if you intend to do reactions without diffusion!

        The duration and granularity of the reactions is specified with 2 out of the 3 parameters:
            total_duration, time_step, n_steps

        For each bin, or each membrane-separated side of bin, (or combined group of bins - not currently implemented),
        process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        Optionally, save some data from the individual reaction steps

        TODO: in case of any Exception, the state of the system is still valid, as of the time before this call

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param snapshots:       OPTIONAL dict that may contain any the following keys:
                                        -"frequency"
                                        -"sample_bin" (Required integer; if not present, no snapshots are taken)
                                        -"species" (NOT YET IMPLEMENTED)
                                        -"initial_caption" (default blank. NOT YET IMPLEMENTED)
                                        -"final_caption" (default blank. NOT YET IMPLEMENTED)
                                    If provided, take a system snapshot after running a multiple
                                    of "frequency" run steps (default 1, i.e. at every step.)
                                    EXAMPLE: snapshots={"frequency": 2, "sample_bin": 0}
        :return:                None
        """
        time_step, n_steps = self.reaction_dynamics.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        # TODO: this is an old system being phased out
        if snapshots:
            frequency = snapshots.get("frequency", 1)         # If not present, it will be 1
            sample_bin = snapshots.get("sample_bin", None)    # If not present, it will be None
        else:
            frequency = None
            sample_bin = None


        for i in range(n_steps):
            self.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            # Update the state of the system
            self.system += self.delta_reactions           # Matrix operation to update all the concentrations
            if self.uses_membranes():
                self.system_B += self.delta_reactions_B   # Matrix operation to update all the concentrations

            self.system_time += time_step

            self.capture_snapshot(step_count=i+1)       # Save historical values (if enabled)
                                                        # It's i+1 because we save the conc. values at the END of the step

            # Preserve some of the data, as requested  TODO: this is an old system being phased out
            if snapshots and ((i+1)%frequency == 0) and (sample_bin is not None):
                self.add_snapshot(self.bin_snapshot(bin_address = snapshots["sample_bin"]))


        if not silent:
            # Print out a summary, at the termination of the run
            print(f"System Time is now: {self.system_time:,.5g}")



    def reaction_step(self, delta_time: float) -> None:
        """
        Compute and store the incremental concentration changes in all bins,
        from all reactions,
        for a single time step of duration delta_time.

        The incremental concentration changes are stored in the class variable
        "delta_reactions", which contains a Numpy array that gets cleared and set.

        IMPORTANT: the actual system concentrations are NOT changed.

        For each bin, process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        TODO: parallelize the computation over the separate bins
        TODO: explore looping over reactions first, and then over bins

        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                            concentration won't vary significantly during this span
        :return:            None (note: the class variable "delta_reactions" gets updated)
        """
        assert self.reaction_dynamics is not None, \
            "reaction_step(): must first set the Reactions object"

        self.delta_reactions = np.zeros((self.n_species, self.n_bins), dtype=float)
        if self.uses_membranes():
            self.delta_reactions_B = np.zeros((self.n_species, self.n_bins), dtype=float)

        # For each bin
        for bin_n in range(self.n_bins):     # Bin number, ranging from 0 to max_bin_number, inclusive
            if self.debug:
                print(f"BioSim1D.reaction_step(): processing the all the reactions in bin number {bin_n}")

            # Obtain the Delta-concentration for each species, for this bin
            conc_array = self.bin_snapshot_array(bin_address=bin_n)
            #print(f"\conc_array in bin {bin_n}: ", conc_array)

            # Obtain the Delta-conc for each species, for bin number bin_n (a NumPy array)
            increment_vector, _, _ = self.reaction_dynamics.reaction_step_common(delta_time=delta_time, conc_array=conc_array,
                                                                                 variable_steps=False)  # Using fixed time steps

            # Replace the "bin_n" column of the self.delta_reactions matrix with the contents of the vector increment_vector
            self.delta_reactions[:, bin_n] = np.array([increment_vector])

            #print(self.delta_reactions)


        # Now process the "other side" of membrane-containing bins, if applicable
        # TODO: the implementation of this part is expected to change soon!
        if self.uses_membranes():
            for bin_n in self.bins_with_membranes():
                # Obtain the Delta-concentration for each species, for this bin
                conc_array = self.bin_snapshot_array(bin_address=bin_n, trans_membrane=True)
                #print(f"\n Post-membrane side conc_dict in bin {bin_n}: ", conc_dict)

                # Obtain the Delta-conc for each species, for bin number bin_n (a NumPy array)
                increment_vector, _, _ = self.reaction_dynamics.reaction_step_common(delta_time=delta_time, conc_array=conc_array,
                                                                                     variable_steps=False)  # Using fixed time steps

                # Replace the "bin_n" column of the self.delta_reactions_B matrix with the contents of the vector increment_vector
                self.delta_reactions_B[:, bin_n] = np.array([increment_vector])



    def react_diffuse(self, total_duration=None, time_step=None, n_steps=None, delta_x = 1) -> None:
        """
        It expects 2 out of the following 3 arguments:  total_duration, time_step, n_steps
        Perform a series of reaction and diffusion constant time steps.

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each constant time step
        :param n_steps:         The desired number of constant steps
        :param delta_x:         Distance between consecutive bins
        :return:                None
        """
        time_step, n_steps = self.reaction_dynamics.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        for i in range(n_steps):
            # TODO: split off the diffusion step and the reaction steps to different computing cores
            self.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            self.diffuse_step(time_step, delta_x=delta_x)
            # Merge into the concentrations of the various bins/chemical species pairs,
            # the increments concentrations computed separately by the reaction and the diffusion steps
            self.system += self.delta_reactions     # Matrix operation to update all the concentrations
                                                    #   from the reactions
            self.system += self.delta_diffusion     # Matrix operation to update all the concentrations
                                                    #   from the diffusion
            self.system_time += time_step

            self.capture_snapshot(step_count=i)     # Save historical values (if enabled)






    #####################################################################################################

    '''                                    ~   VISUALIZATION   ~                                           '''

    def ________VISUALIZATION________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def visualize_system(self, title_prefix=None, colors=None) -> pgo.Figure:
        """
        Visualize the current state of the system of all the chemicals as a combined line plot,
        using plotly

        :param title_prefix:[OPTIONAL] A string to prefix to the auto-generated title
        :param colors:      [OPTIONAL] If None, then use the registered colors (if specified),
                                or the hardwired defaults as a last resort
        :return:            A Plotly "Figure" object
        """
        title = f"System snapshot at time t={self.system_time:.8g}"
        if title_prefix:
            title = title_prefix + ".  " + title

        chem_labels = self.chem_data.get_all_labels()
        n_chem = len(chem_labels)

        if colors is None:
            # Attempt to make use of the previously-registered colors, if available
            colors = self.chem_data.get_registered_colors(chem_labels)
            if colors is None:
                # Fall back to default colors
                colors = Colors.assign_default_colors(n_chem)


        if n_chem == 1:
            chem_name = chem_labels[0]      # The only chemical

            fig = px.line(y=self.lookup_species(species_name=chem_name),
                          title= title,
                          color_discrete_sequence = colors,
                          labels={"y": f"[{chem_name}]", "x":"Bin number"},)
        else:
            fig = px.line(data_frame=self.system_snapshot(), y=chem_labels,
                          title= title,
                          color_discrete_sequence = colors,
                          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})                        

        return fig



    def system_heatmap(self, chem_labels=None, title_prefix ="", row_height=150,
                           colors=None) -> pgo.Figure:
        """
        Produce a heatmap, and return it as a plotly Figure object

        :param chem_labels: [OPTIONAL] NOT YET USED.  For now, ALL chemicals get shown
        :param title_prefix:[OPTIONAL] A string to prefix to the auto-generated title
        :param row_height:  [OPTIONAL] Height of each of the heatmap rows;
                                if too small to fit the heatmap, it will automatically be made larger
        :param colors:      [OPTIONAL] If None, then use the registered colors (if specified);
                                in case of absence of registered colors, use black if just 1 chemical,
                                or the hardwired default colors otherwise
        :return:            A Plotly "Figure" object with a Heatmap
        """
        # TODO: generalize to sub-selection of chemicals
        title = f"System snapshot at time t={self.system_time:.8g}"
        if title_prefix:
            title = title_prefix + ".  " + title

        '''
        # ANOTHER APPROACH TO CREATE A PLOTLY HEATMAP, using imshow()
        fig = px.imshow(self.system_snapshot().T,
                title= title,
                labels=dict(x="Bin number", y="Chem. species", color="Concentration"),
                text_auto='.2f', color_continuous_scale="gray_r")         # text_auto=’auto’

        # Insert a little spacing between adjacent cells horizontally, and substantially more vertical
        fig.update_traces(xgap=3)       # Alt way to specify it: fig.data[0].xgap=3
        fig.update_traces(ygap=8)       # Vertical separation is larger because individual rows belong to different chemicals
        '''

        if chem_labels is None:
            chem_labels = self.chem_data.get_all_labels()

        conc_matrix = self.system_snapshot().T
        hm = pgo.Heatmap(z=conc_matrix,
                            y=chem_labels,
                            colorscale='gray_r', colorbar={'title': 'Concentration'},
                            xgap=3, ygap=8, texttemplate = '%{z:.2f}',
                            hovertemplate= 'Conc.: %{z}<br>Bin number: %{x}<br>Chem: %{y}<extra>%{y}</extra>')

        # Create the Figure object
        fig = pgo.Figure(data=hm)


        if colors is None:
            # Attempt to make use of the previously-registered colors, if available
            colors = self.chem_data.get_registered_colors(chem_labels)
            if (colors is None) and len(chem_labels) > 1:
                # Fall back to default colors (but don't bother if just 1 chemical)
                colors = Colors.assign_default_colors(len(chem_labels))

        min_height = 235        # Too small leads to an obscure JS error message
        if len(chem_labels) == 1:
            row_height = min_height

        if colors:
            # "delta" affect the vertical extension of the color-coded edges; fine-tune for visual esthetics
            if row_height <= 200:
                delta = 0.46
            elif row_height <= 350:
                delta = 0.48
            else:
                delta = 0.49

            # Add color-coded edges using shapes
            for i, col in enumerate(colors):
                for j in range(self.n_bins):
                    fig.add_shape(
                        type="rect",
                        x0=j - 0.5,   x1=j + 0.5,         # Make x0  more neg to extend to left; make x1 more pos to extend to right
                        y0=i - delta, y1=i + delta,       # Make y0  more neg to extend to bottom; make y1 more pos to extend to top
                        line=dict(color=col, width=3),    # Edge color and width
                        xref="x", yref="y"
                    )

            # Add annotations for row color-coding
            for i, col in enumerate(colors):
                fig.add_shape(
                    type="rect",
                    x0 = -0.5,     x1 = -0.4,       # Slightly outside the heatmap
                    y0 = i - delta, y1 = i + delta,
                    xref="x", yref="y",
                    line=dict(width=0),
                    fillcolor=col
                )


        height = max(min_height, row_height * len(chem_labels))

        fig.update_layout(title=title,
                          xaxis={'title': 'Bin number'}, yaxis={'title': 'Chem. species'},
                          yaxis_autorange="reversed",
                          height=height)

        return fig



    def single_species_heatmap(self, species_index: int, heatmap_pars: dict, graphic_component, header=None) -> None:
        """
        Send to the HTML log, a heatmap representation of the concentrations of
        the single requested species.  Note: if using in Jupyterlab, this image will NOT be displayed there

        IMPORTANT: must first call GraphicLog.config(), or an Exception will be raised

        :param species_index:       Index identifying the species of interest
        :param heatmap_pars:        A dictionary of parameters for the heatmap
        :param graphic_component:   A string with the name of the graphic module to use.  EXAMPLE: "vue_heatmap_11"
        :param header:              Optional string to display just above the heatmap
        :return:                    None
        """
        if not GraphicLog.is_initialized():
            raise Exception("Prior to calling single_species_heatmap(), "
                            "need to initialize the graphics module with a call to GraphicLog.config()")

        if header:
            log.write(f"{header}", style=log.h1, newline=False, also_print=False)


        #
        # Prepare the data for the heatmap
        #

        # List of concentrations in the various bins, for the requested species
        species_concentrations = list(self.lookup_species(species_index))
        #print(species_concentrations)

        all_data = {
            "y_labels": [f"Mol {species_index}"],

            # Data for the heatmap, by rows, starting with the bottom one
            "heatmap_data": [species_concentrations],

            # Set the range of values in the heatmap bins
            "range_min": heatmap_pars["range"][0],
            "range_max": heatmap_pars["range"][1],

            # Set the dimensions and margins of the heatmap
            "outer_width": heatmap_pars["outer_width"],
            "outer_height": heatmap_pars["outer_height"],
            "margins": heatmap_pars["margins"]
        }

        # Send the plot to the HTML log file.
        # The version of the heatmap Vue component specified
        # in the call to GraphicLog.config() will be used
        GraphicLog.export_plot(all_data, graphic_component, unpack=True)



    def single_species_line_plot(self, species_index: int, plot_pars: dict, graphic_component, header=None) -> None:
        """
        Send to the HTML log, a line plot representation of the concentrations of
        the single requested species.  To plot more than 1 species, use line_plot() instead

        IMPORTANT: must first call GraphicLog.config(), or an Exception will be raised

        :param species_index:       Index identifying the species of interest
        :param plot_pars:           A dictionary of parameters for the plot
        :param graphic_component:   A string with the name of the graphic module to use.  EXAMPLE: "vue_curves_3"
        :param header:              Optional string to display just above the plot
        :return:                    None
        """
        if not GraphicLog.is_initialized():
            raise Exception("Prior to calling single_species_line_plot(), "
                            "need to initialize the graphics module with a call to GraphicLog.config()")

        if header:
            log.write(f"{header}", style=log.h1, newline=False)


        #
        # Prepare the data for the plot
        #

        # List of concentrations in the various bins, for the requested species
        species_concentrations = list(self.lookup_species(species_index))
        #print(species_concentrations)

        all_data = {
            "y_labels": [f"Mol {species_index}"],

            # Concentration data for the plots (for now just 1 chemical species), in index order
            "data": species_concentrations,

            # Set the range of values for the y-scale of the plot
            "range_min": plot_pars["range"][0],
            "range_max": plot_pars["range"][1],

            # Set the dimensions and margins of the plot
            "outer_width": plot_pars["outer_width"],
            "outer_height": plot_pars["outer_height"],
            "margins": plot_pars["margins"]
        }

        # Send the plot to the HTML log file.
        # The version of the heatmap Vue component specified
        # in the call to GraphicLog.config() will be used
        GraphicLog.export_plot(all_data, graphic_component, unpack=True)



    def line_plot(self, plot_pars: dict, graphic_component, header=None, color_mapping=None) -> None:
        """
        Send to the HTML log, a line plot representation of the concentrations of
        all the chemical species species.

        IMPORTANT: must first call GraphicLog.config(), or an Exception will be raised

        :param plot_pars:           A dictionary of parameters (such as "outer_width") for the plot
        :param graphic_component:   A string with the name of the graphic module to use.  EXAMPLE: "vue_curves_4"
        :param header:              [OPTIONAL] String to display just above the plot
        :param color_mapping:       [OPTIONAL] Dict mapping index numbers to color names or RBG hex values.
                                        If not provided, the colors associated to the chemicals are used, if any
        :return:                    None
        """
        #TODO: offer an option to limit which chemical species to display

        if not GraphicLog.is_initialized():
            raise Exception("Prior to calling line_plot(), "
                            "need to initialize the graphics module with a call to GraphicLog.config()")

        if header:
            # Display the requested header just above the plot, in the log file
            log.write(f"{header}", style=log.h1, newline=False)


        #
        # Prepare the data for the plot
        #

        # Concentration data for the plots
        #       outer level : order of chemical-species index,
        #       inner level : in bin index order from left to right
        species_concentrations = [list(self.lookup_species(i)) for i in range(self.n_species)]
        #print(species_concentrations)

        all_data = {
            "curve_labels": self.chem_data.get_all_labels(),

            # Concentration data for the plots
            "plot_data": species_concentrations,

            # Set the range of values for the y-scale of the plot
            "range_min": plot_pars["range"][0],
            "range_max": plot_pars["range"][1],

            # Set the dimensions and margins of the plot
            "outer_width": plot_pars["outer_width"],
            "outer_height": plot_pars["outer_height"],
            "margins": plot_pars["margins"]
        }

        # If a color mapping was provided, add it to the data
        if color_mapping:
            all_data["color_mapping"] = color_mapping
        elif col_map := self.chem_data.get_color_mapping_by_index():
            all_data["color_mapping"] = col_map


        # Send the plot to the HTML log file.
        # The version of the heatmap Vue component specified
        # in the call to GraphicLog.config() will be used
        GraphicLog.export_plot(all_data, graphic_component, unpack=True)



    def plot_history_single_bin(self, bin_address :int, colors=None, title=None, smoothed=False) -> pgo.Figure:
        """
        Using plotly, draw the plots of chemical concentration values over time at the specified bin,
        based on the historical data that was saved when running simulations.

        Note: if this plot is to be later combined with others, use PlotlyHelper.combine_plots()

        :param bin_address: A single bin address (an integer)
        :param colors:      [OPTIONAL] List of CSS color names for each of the heatmaps.
                                If provided, its length must match that of the data;
                                if None, then use the registered colors (if specified),
                                or the hardwired defaults as a last resort
        :param title:       [OPTIONAL] Label for the top of the plot.  If not passed, a default is used
        :param smoothed:    [OPTIONAL] If True, a spline is used to smooth the lines;
                                otherwise (default), line segments are used
        :return:            A plotly "Figure" object; an Exception is raised if no historical data is found
        """
        # TODO: add more options

        self.assert_valid_bin(bin_address)

        if title is None:
            if self.chem_data.number_of_chemicals() == 1:
                chem_title = f"chemical `{self.chem_data.get_label(0)}`"    # The label of the only chemical in the system
            else:
                chem_title = "all chemicals"

            title = f"Concentration changes with time of {chem_title} at bin {bin_address}"

        # Retrieve the historical data
        df = self.conc_history.bin_history(bin_address = bin_address)
        if type(df) == str:         # No data was found
            raise Exception(df)

        #chem_labels = list(df.columns).remove("SYSTEM TIME")  # All the column names, except the independent var
        chem_labels = self.conc_history.restrict_chemicals  # The chemicals for which history was kept

        if colors is None:  # Attempt to make use of the previously-registered colors, if available
            colors = self.chem_data.get_registered_colors(chem_labels)

        return PlotlyHelper.plot_pandas(df, x_var="SYSTEM TIME", y_label="Concentration",
                                        colors=colors, legend_header="Chemical", title=title,
                                        smoothed=smoothed)






    #########################################################################
    #                                                                       #
    #                       FOURIER ANALYSIS                                #
    #                                                                       #
    #########################################################################

    def frequency_analysis(self, species_name: str, threshold = 0.001, n_largest = None) -> pd.DataFrame:
        """
        Return the individual frequencies, and their relative amplitudes,
        in the concentration values of the specified chemical species.
        A Discrete Fourier Transform is used for the computation.

        :param species_name:    The name of the chemical whose concentration we want to analyze
        :param threshold:       Minimum amplitudes of the frequency components to be considered non-zero
                                    (NOTE: these are the raw values returned by the DFT - not the normalized ones.)
        :param n_largest:       If specified, only the rows with the given number of largest amplitudes gets returned
                                    (if there are fewer rows to start with, they all get returned)

        :return:                A Pandas dataframe with 2 columns, "Frequency" and "Relative Amplitude";
                                    amplitudes are relative the the smallest nonzero frequency (which is taken to be 1.0)
                                    EXAMPLE:
                                                   Frequency  Relative Amplitude
                                                0        0.0                 3.0
                                                1        2.0                 1.0
                                                2        4.0                 0.5
                                                3        8.0                 0.2
        """

        conc_samples = self.lookup_species(species_name=species_name)

        #import plotly.express as px
        #fig = px.line(y=conc_samples)
        #fig.show()

        # Perform a DFT, to extract the frequency components
        # The size of the computed arrays xf and yf is  (self.n_bins/2 + 1) if self.n_bins is even
        # or (self.n_bins + 1) /2 if self.n_bins is odd
        xf = rfftfreq(self.n_bins, 1 / self.n_bins)
        yf = rfft(conc_samples)
        magnitude_yf = np.abs(yf)   # The magnitudes of the complex numbers in yf

        #print(xf)
        #print(magnitude_yf)

        #print("Size of returned array: ", len(xf))

        #fig = px.line(x=xf, y=magnitude_yf)
        #fig.show()

        above_threshold = magnitude_yf > threshold  # Numpy array of booleans,
        # indicating the positions in magnitude_yf array
        # where the value exceeds the given threshold
        #print(above_threshold)

        amps = magnitude_yf[above_threshold]        # Above-threshold amplitudes
        freqs = xf[above_threshold]                 # Frequencies corresponding to above-threshold amplitudes

        #print(freqs)
        #print(amps)

        if np.allclose(freqs[0], 0.):
            # If there's a "DC bias" (zero-frequency component)
            amps[0] /= 2.   # Amplitudes of the DC components show up differently;
            # see https://www.sjsu.edu/people/burford.furman/docs/me120/FFT_tutorial_NI.pdf
            baseline_amp = amps[1]  # Use the first non-zero frequency as the baseline value
        else:
            baseline_amp = amps[0]  # Use the first non-zero frequency as the baseline value

        scaled_amps = amps / baseline_amp   # Amplitudes scaled by the amplitude
        #   of the first non-zero frequency

        # Assemble and return as a Pandas dataframe the frequencies with amplitudes above the threshold,
        #   together with their relative amplitudes
        df = pd.DataFrame()
        df['Frequency'] = freqs
        df['Relative Amplitude'] = scaled_amps

        if (n_largest is not None) and (n_largest < len(df)):
            return df.nlargest(n=n_largest, columns="Relative Amplitude")

        return df
