import numpy as np
import pandas as pd
import math
from scipy.fft import rfft, rfftfreq    # Fast Fourier Transforms to extract frequency components
from scipy.stats import norm
from typing import Union, List, Tuple
from life123.movies import CollectionTabular
from life123.uniform_compartment import UniformCompartment
import plotly.express as px
from life123.html_log import HtmlLog as log
from life123.visualization.graphic_log import GraphicLog



class BioSim1D:
    """
    1D simulations of diffusion and reactions, 
    with an early partial implementation of membranes
    """


    def __init__(self, n_bins=None, chem_data=None, reactions=None):
        """

        :param n_bins:      The number of compartments (bins) to use in the simulation
        :param chem_data:   (OPTIONAL) Object of class "ReactionData";
                                if not specified, it will get extracted from the "UniformCompartment" class
        :param reactions:   (OPTIONAL) Object of class "UniformCompartment";
                                if not specified, it'll get instantiated here   TODO: maybe no longer necessary
        """
        self.debug = False

        self.n_bins = 0         # Number of spatial compartments (bins) used in the simulation

        self.n_species = 1      # The number of (non-water) chemical species   TODO: phase out?

        self.chem_data = None   # Object of type "ReactionData", with info on the individual chemicals and their reactions

        self.system_length = None   # System extension, from the middle of the leftmost bin to the middle of the rightmost one.
                                    #  The de-facto default value, though not used, is (n_bins-1)

        self.global_Dx = 1      # Used in cases when not using ad-hoc local changes in x-scale resolution

        self.system = None      # Concentration data in the System we're simulating, for all the chemicals
                                #   NumPy array of floats, of dimension: (n_species x n_bins).
                                #   Each row represents a species

        self.system_earlier = None  # NOT IN CURRENT USE.  Envisioned for simulations where the past 2 time states are used
                                    # to compute the state at the next time step

        self.membranes = None   # NumPy boolean array of dimension n_bins;
                                #   each boolean value marks the presence of a membrane in that bin.
                                #   Any given bin is expected to either contain, or not contain, a single membrane passing thru it
                                #   TODO - explore possible alternative: list of bin addresses the contains a membrane,
                                #          or other sparse-matrix representation

        self.A_fraction = None  # NumPy array of floats, of dimension n_bins;
                                #   each value records the fraction (between 0. and 1.) of the LEFT one
                                #   of the 2 parts into which the membrane splits the bin.
                                #   0. if N/A

        self.system_B = None    # Just like "system", but only for the "A" fractions of the bin, where applicable;
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

        self.reaction_dynamics = None       # Object of class "UniformCompartment"

        self.history = CollectionTabular()       # To store user-selected snapshots of (parts of) the system,
                                            #   whenever requested by the user.
                                            #   Note that we're using the "tabular" format - friendly to Pandas

        self.system_time = None             # Global time of the system, from initialization on

        if (n_bins is not None) or (chem_data is not None) or (reactions is not None):
            self.initialize_system(n_bins=n_bins, chem_data=chem_data, reactions=reactions)




    #########################################################################
    #                                                                       #
    #                           SYSTEM-WIDE                                 #
    #                                                                       #
    #########################################################################

    def initialize_system(self, n_bins: int, chem_data=None, reactions=None) -> None:
        """
        Initialize all concentrations to zero.
        Membranes, if present, need to be set later.

        TODO?: maybe allow optionally passing n_species in lieu of chem_data,
              and let it create and return the "Chemicals" object in that case

        :param n_bins:      The number of compartments (bins) to use in the simulation
        :param chem_data:   (OPTIONAL) Object of class "ReactionData";
                                if not specified, it will get extracted from the "UniformCompartment" class
        :param reactions:   (OPTIONAL) Object of class "UniformCompartment";
                                if not specified, it'll get instantiated here
        :return:            None
        """
        assert n_bins >= 1, "The number of bins must be at least 1"

        assert chem_data is not None or reactions is not None, \
            "BioSim1D: at least one of the arguments `chem_data` and `reactions` must be set"
            # TODO: maybe drop this requirement?  And then set the system matrix later on?

        if chem_data:
            self.chem_data = chem_data
        else:
            self.chem_data = reactions.chem_data

        if reactions:
            self.reaction_dynamics = reactions
        else:
            self.reaction_dynamics = UniformCompartment(chem_data=chem_data)

        self.n_bins = n_bins

        self.n_species = chem_data.number_of_chemicals()

        assert self.n_species >= 1, \
            "At least 1 chemical species must be declared prior to calling initialize_system()"

        # Initialize all concentrations to zero
        self.system = np.zeros((self.n_species, n_bins), dtype=float)

        self.system_time = 0             # "Start the clock"



    def system_size(self) -> int:
        """
        Note: the bin numbers will range between 0 and system_size - 1
        :return:    The number of bins in the system
        """
        return self.n_bins



    def reset_system(self) -> None:
        """
        WARNING - THIS IS VERY PARTIAL.  TODO: expand, or drop (not sure if really needed anymore)
        :return:
        """
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
        TODO: deal with membranes, and anything else needed to later restore the complete system state

        :return:    A dict of (for now a part of) the elements needed to later restore the complete system state
        """
        return {"system": self.system.copy(), "system_time": self.system_time}



    
    def restore_system(self, new_state: dict) -> None:
        """
        Replace (some, for now, of) the various parts the System's internal state.
        For details of the data structure, see the class variable "system"
        TODO: membranes aren't yet managed. System length and global_Dx are currently not modified

        :param new_state:   Numpy array containing the desired new System's internal state
        :return:
        """
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



    
    def system_snapshot(self) -> pd.DataFrame:
        """
        Return a snapshot of all the concentrations of all the species, across all bins
        as a Pandas dataframe
        TODO: make allowance for membranes

        :return:    A Pandas dataframe: each row is a bin,
                        and each column a chemical species
        """
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



    
    def set_all_uniform_concentrations(self, conc_list: Union[list, tuple]) -> None:
        """
        Set the concentrations of all species at once, uniformly across all bins
        :param conc_list:   List or tuple of concentration values for each of the chemical species
        :return:            None
        """
        assert len(conc_list) == self.chem_data.number_of_chemicals(), \
            f"set_all_uniform_concentrations(): the argument must be a list or tuple of size {self.chem_data.number_of_chemicals()}"

        for i, conc in enumerate(conc_list):
            self.set_uniform_concentration(species_index=i, conc=conc)



    
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
        if species_name is not None:
            species_index = self.chem_data.get_index(species_name)
        else:
            self.chem_data.assert_valid_species_index(species_index)

        self.assert_valid_bin(bin_address)


        assert conc >= 0., \
            f"set_bin_conc(): the concentration must be a positive number or zero (the requested value was {conc})"

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


    
    def inject_conc_to_bin(self, bin_address: int, species_index: int, delta_conc: float, zero_clip = False) -> None:
        """
        Add the requested concentration to the cell with the given address, for the specified chem species

        :param bin_address:     The zero-based bin number of the desired cell
        :param species_index:   Zero-based index to identify a specific chemical species
        :param delta_conc:      The concentration to add to the specified location
        :param zero_clip:       If True, any requested increment causing a concentration dip below zero, will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        self.assert_valid_bin(bin_address)

        if (self.system[species_index, bin_address] + delta_conc) < 0. :
            # Take special action if the requested change would make the bin concentration negative
            if zero_clip:
                self.system[species_index, bin_address] = 0
                return
            else:
                raise Exception("inject_conc_to_bin(): The requested concentration change would result in a negative final value")

        # Normal scenario, not leading to negative values for the final concentration
        self.system[species_index, bin_address] += delta_conc




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
            self.inject_conc_to_bin(bin_address = x, species_index = species_index,
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



    ########  DIMENSION-RELATED  ################
    
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




    #########################################################################
    #                                                                       #
    #                              TO VIEW                                  #
    #                                                                       #
    #########################################################################

    def assert_valid_bin(self, bin_address: int) -> None:
        """
        Raise an Exception if the given bin number isn't valid

        :param bin_address:  An integer that ought to be between 0 and (self.n_bins-1), inclusive
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


    
    def bin_concentration(self, bin_address: int, species_index=None, species_name=None, trans_membrane=False) -> float:
        """
        Return the concentration at the requested bin of the specified species

        :param bin_address:     The bin number
        :param species_index:   The index order of the chemical species of interest
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :param trans_membrane:  If True, consider the "other side" of the bin, i.e. the portion across the membrane
        :return:                A concentration value at the indicated bin, for the requested species
        """
        if species_name is not None:
            species_index = self.chem_data.get_index(species_name)

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



    def show_system_snapshot(self) -> None:
        """

        :return:    None
        """
        print(f"SYSTEM SNAPSHOT at time {self.system_time}:")
        print(self.system_snapshot())



    
    def describe_state(self, concise=False) -> None:
        """
        A simple printout of the state of the system, for now useful only for small systems

        TODO: The goal for ASCII-based printouts involving membranes is something like.  Maybe use Pandas?

           _____________________
        A: |20|18|12|8( )  2| 6|    Diff rate: 0.2 (M: 0.01)
        B: |30| 2|56|4( )3.5|12|    Diff rate: 0.8 (M: 0.3)
           ---------------------
             0  1  2       3  4

        :param concise: If True, only produce a minimalist printout with just the concentration values
        :return:        None
        """
        print(f"SYSTEM STATE at Time t = {self.system_time:,.8g}:")

        if concise:             # A minimalist printout...
            print(self.system)   # ...only showing the concentration data (a Numpy array)
            if self.uses_membranes():
                print("Membranes:")
                print(self.system_B)
            return

        # If we get thus far, it's a FULL printout

        print(f"{self.n_bins} bins and {self.n_species} species:")

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



    def show_membranes(self, n_decimals=1) -> str:
        """
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




    #########################################################################
    #                                                                       #
    #                        CHANGE RESOLUTIONS                             #
    #                                                                       #
    #########################################################################
    
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




    #########################################################################
    #                                                                       #
    #                               DIFFUSION                               #
    #                                                                       #
    #########################################################################

    
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




    #########################################################################
    #                                                                       #
    #                               REACTIONS                               #
    #                                                                       #
    #########################################################################

    def react(self, total_duration=None, time_step=None, n_steps=None, snapshots=None) -> None:
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

        # TODO: validation; also, implement "species" option for snapshots
        first_snapshot = True
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
            # Preserve some of the data, as requested
            if snapshots and ((i+1)%frequency == 0) and (sample_bin is not None):
                self.add_snapshot(self.bin_snapshot(bin_address = snapshots["sample_bin"]))



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

        #number_reactions = self.reaction_dynamics.number_of_reactions()

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




    #########################################################################
    #                                                                       #
    #                         REACTION-DIFFUSION                            #
    #                                                                       #
    #########################################################################

    def react_diffuse(self, total_duration=None, time_step=None, n_steps=None, delta_x = 1) -> None:
        """
        It expects 2 of the arguments:  total_duration, time_step, n_steps
        Perform a series of reaction and diffusion time steps.

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param delta_x:         Distance between consecutive bins
        :return:                None
        """
        time_step, n_steps = self.reaction_dynamics.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        for i in range(n_steps):
            # TODO: split off the reaction step and the diffusion step to 2 different computing cores
            self.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            self.diffuse_step(time_step, delta_x=delta_x)
            # Merge into the concentrations of the various bins/chemical species pairs,
            # the increments concentrations computed separately by the reaction and the diffusion steps
            self.system += self.delta_reactions     # Matrix operation to update all the concentrations
                                                    #   from the reactions
            self.system += self.delta_diffusion     # Matrix operation to update all the concentrations
                                                    #   from the diffusion
            self.system_time += time_step




    #####################################################################################################

    '''                                    ~   VISUALIZATION   ~                                           '''

    def ________VISUALIZATION________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def visualize_system(self, caption=None, colors=None) -> None:
        """
        Visualize the current state of the system as a line plot,
        using plotly

        :param caption: Optional caption to prefix to the default one
        :param colors:  [Not yet used]
        :return:        None
        """
        title = f"System snapshot at time t={self.system_time}"
        if caption:
            title = caption + ".  " + title

        if self.n_species == 1:
            chem_name = self.chem_data.get_label(species_index=0)      # The only chemical in the system

            fig = px.line(y=self.lookup_species(species_name=chem_name),
                          title= title,
                          labels={"y": f"[{chem_name}]", "x":"Bin number"},)
            fig.show()
        else:
            print("NOT YET IMPLEMENTED")
            '''
            #TODO: generalize
            fig = px.line(data_frame=self.system_snapshot(), y=["A"],
                          title= f"Diffusion. System snapshot at time t={self.system_time}",
                          color_discrete_sequence = ['red'],
                          labels={"value":"concentration", "variable":"Chemical", "index":"Bin number"})                        
            fig.show()
            '''



    def single_species_heatmap(self, species_index: int, heatmap_pars: dict, graphic_component, header=None) -> None:
        """
        Send to the HTML log, a heatmap representation of the concentrations of
        the single requested species.

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
        TODO: offer an option to limit which chemical species to display

        IMPORTANT: must first call GraphicLog.config(), or an Exception will be raised

        :param plot_pars:           A dictionary of parameters (such as "outer_width") for the plot
        :param graphic_component:   A string with the name of the graphic module to use.  EXAMPLE: "vue_curves_4"
        :param header:              OPTIONAL string to display just above the plot
        :param color_mapping:       OPTIONAL dict mapping index numbers to color names or RBG hex values
        :return:                    None
        """
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


        # Send the plot to the HTML log file.
        # The version of the heatmap Vue component specified
        # in the call to GraphicLog.config() will be used
        GraphicLog.export_plot(all_data, graphic_component, unpack=True)




    #########################################################################
    #                                                                       #
    #                              HISTORY                                  #
    #                                                                       #
    #########################################################################
    
    def add_snapshot(self, data_snapshot: dict, caption ="") -> None:
        """
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
        Retrieve and return a Pandas dataframe with the system history that had been saved
        using add_snapshot()

        :return:        a Pandas dataframe
        """
        return self.history.get_dataframe()




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
