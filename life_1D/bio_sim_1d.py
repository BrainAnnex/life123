import numpy as np
import pandas as pd
import math
from scipy.fft import rfft, rfftfreq    # Fast Fourier Transforms to extract frequency components
from typing import Union, List, Tuple
from modules.movies.movies import Movie
from modules.reactions.reactions import Reactions
from modules.html_log.html_log import HtmlLog as log
from modules.visualization.graphic_log import GraphicLog


class BioSim1D:
    """
    Note: for at least the time being, this class does NOT get instantiated
    """

    #####################
    #  Class variables  #
    #####################

    verbose = False

    n_bins = 0          # Number of spacial compartments (bins) used in the simulation

    n_species = 1       # The number of (non-water) chemical species   TODO: phase out?

    chem_data = None    # Object of type "Chemicals", with info on the individual chemicals,
                        #   incl. their names and diffusion rates

    system_length = None    # System extension, from the middle of the leftmost bin to the middle of the rightmost one.
                            #  The de-facto default value, though not used, is (n_bins-1)

    global_Dx = 1           # Used in cases when not using ad-hoc local changes in x-scale resolution

    system = None           # Concentration data in the System we're simulating, for all the chemicals
                            #   NumPy array of floats, of dimension: (n_species x n_bins).
                            #   Each row represents a species

    system_earlier = None   # NOT IN CURRENT USE.  Envisioned for simulations where the past 2 time states are used
                            # to compute the state at the next time step

    membranes = None        # NumPy boolean array of dimension n_bins;
                            #   each boolean value marks the presence of a membrane in that bin.
                            #   Any given bin is expected to either contain, or not contain, a single membrane passing thru it
                            #   TODO - explore possible alternative: list of bin addresses the contains a membrane,
                            #          or other sparse-matrix representation

    A_fraction = None       # NumPy array of floats, of dimension n_bins;
                            #   each value records the fraction (between 0. and 1.) of the LEFT one
                            #   of the 2 parts into which the membrane splits the bin.
                            #   0. if N/A

    system_B = None         # Just like "system", but only for the "A" fractions of the bin, where applicable;
                            #   0. if N/A

    delta_diffusion = None  # Buffer for the concentration changes from diffusion step (n_species x n_bins)
    delta_reactions = None  # Buffer for the concentration changes from reactions step (n_species x n_bins)

    delta_reactions_B = None    # Same as above, but for the "other sides" of bins with membranes
                                # TODO: explore sparse-matrix representations

    sealed = True           # If True, no exchange with the outside; if False, immersed in a "bath"

    # Only applicable if "sealed" is False:
    bath_concentrations = None      # A NumPy array for each species
    container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container

    time_step_threshold = 0.33333   # This is used to set an Upper Bound on the single time steps
                                    #   in the diffusion process.
                                    #   See explanation in file overly_large_single_timesteps.py

    all_reactions = None            # Object of class "Reactions"

    history = Movie(tabular=True)   # To store user-selected snapshots of (parts of) the system,
                                    #   whenever requested by the user

    system_time = None              # Global time of the system, from initialization on

    debug = False



    #########################################################################
    #                                                                       #
    #                           SYSTEM-WIDE                                 #
    #                                                                       #
    #########################################################################

    @classmethod
    def initialize_system(cls, n_bins: int, chem_data=None, reactions=None) -> None:
        """
        Initialize all concentrations to zero.
        Membranes, if present, need to be set later.

        TODO?: maybe allow optionally passing n_species in lieu of chem_data,
              and let it create and return the "Chemicals" object in that case

        :param n_bins:      The number of compartments (bins) to use in the simulation
        :param chem_data:   (OPTIONAL) Object of class "Chemicals";
                                if not specified, it will get extracted from the "Reactions" class
        :param reactions:   (OPTIONAL) Object of class "Reactions";
                                if not specified, it'll get instantiated here
        :return:            None
        """
        assert n_bins >= 1, "The number of bins must be at least 1"

        assert chem_data is not None or reactions is not None, \
            "BioSim1D: at least one of the arguments `chem_data` and `reactions` must be set"
            # TODO: maybe drop this requirement?  And then set the system matrix later on?

        if chem_data:
            cls.chem_data = chem_data
        else:
            cls.chem_data = reactions.chem_data

        if reactions:
            cls.all_reactions = reactions
        else:
            cls.all_reactions = Reactions(chem_data=chem_data)

        cls.n_bins = n_bins

        cls.n_species = chem_data.n_species

        assert cls.n_species >= 1, \
            "At least 1 chemical species must be declared prior to calling initialize_system()"

        # Initialize all concentrations to zero
        cls.system = np.zeros((cls.n_species, n_bins), dtype=float)

        cls.system_time = 0             # "Start the clock"



    @classmethod
    def reset_system(cls) -> None:
        """
        WARNING - THIS IS VERY PARTIAL.  TODO: expand, or switch to instantiated class
        :return:
        """
        cls.system = None
        cls.system_earlier = None
        cls.system_B = None
        cls.membranes = None
        cls.A_fraction = None
        cls.chem_data = None
        cls.all_reactions = None
        cls.history = None

        cls.n_bins = 0
        cls.system_length = None
        cls.global_Dx = 1
        cls.system_time = 0

        cls.delta_reactions = None
        cls.delta_reactions_B = None



    @classmethod
    def save_system(cls) -> dict:
        """
        For now, just return a copy of cls.system, with a "frozen" snapshot of the current system state
        TODO: deal with membranes, and anything else needed to later restore the complete system state

        :return:    A dict of (for now a part of) the elements needed to later restore the complete system state
        """
        return {"system": cls.system.copy(), "system_time": cls.system_time}



    @classmethod
    def restore_system(cls, new_state: dict) -> None:
        """
        Replace (some, for now, of) the various parts the System's internal state.
        For details of the data structure, see the class variable "system"
        TODO: membranes aren't yet managed. System length and global_Dx are currently not modified

        :param new_state:   Numpy array containing the desired new System's internal state
        :return:
        """
        cls.system = new_state["system"]
        cls.n_species, cls.n_bins = cls.system.shape     # Extract from the new state
        assert cls.n_species == cls.chem_data.n_species, \
            f"restore_system(): inconsistency in the number of chemical species in the specified new state ({cls.n_species}) " \
            f"vs. what's stored in the `Chemicals` object ({cls.chem_data.n_species})"

        cls.system_time = new_state["system_time"]



    @classmethod
    def replace_system(cls, new_state: np.array) -> None:
        """
        Replace the System's internal state.
        For details of the data structure, see the class variable "system"
        IMPORTANT: membranes aren't handled. System length and global_Dx are currently not modified

        :param new_state:   Numpy array containing the desired new System's internal state
        :return:
        """
        cls.system = new_state
        cls.n_species, cls.n_bins = new_state.shape     # Extract from the new state
        assert cls.n_species == cls.chem_data.n_species, \
            "replace_system(): inconsistency in the number of chemical species vs. what's stored in the `Chemicals` object"



    @classmethod
    def system_snapshot(cls) -> pd.DataFrame:
        """
        Return a snapshot of all the concentrations of all the species,
        as a Pandas dataframe
        TODO: make allowance for membranes

        :return:    A Pandas dataframe: each row is a bin,
                        and each column a chemical species
        """
        all_chem_names = cls.chem_data.get_all_names()
        if cls.system is None:
            return pd.DataFrame(columns = all_chem_names)   # Empty dataframe

        matrix = cls.system.T

        df = pd.DataFrame(matrix, columns = all_chem_names)

        return df



    @classmethod
    def compare_states(cls, state1: np.array, state2: np.array, verbose=False):
        """

        :param state1:
        :param state2:
        :param verbose:
        :return:
        """
        diff = state1 - state2
        abs_diff = abs(diff)
        print("Max of unsigned absolute differences: ", np.max(abs_diff))

        rel_diff = diff / state1

        if verbose:
            print("Relative differences: ", rel_diff)

        print("Max of unsigned relative differences: ", np.max(abs(rel_diff)))

        print("Mean of relative differences: ", np.mean(rel_diff))
        print("Median of relative differences: ", np.median(rel_diff))
        print("Standard deviation of relative differences: ", np.std(rel_diff))

        print("np.allclose with lax tolerance?  (rtol=1e-01, atol=1e-01) : ",
              np.allclose(state1 , state2, rtol=1e-01, atol=1e-01))
        print("np.allclose with mid tolerance?  (rtol=1e-02, atol=1e-03) : ",
              np.allclose(state1 , state2, rtol=1e-02, atol=1e-03))
        print("np.allclose with tight tolerance?  (rtol=1e-03, atol=1e-05) : ",
              np.allclose(state1 , state2, rtol=1e-03, atol=1e-05))
        print("np.allclose with extra-tight tolerance?  (rtol=1e-05, atol=1e-08) : ",
              np.allclose(state1 , state2, rtol=1e-05, atol=1e-08))





    #########################################################################
    #                                                                       #
    #                SET/MODIFY CONCENTRATIONS (or membranes)               #
    #                                                                       #
    #########################################################################

    @classmethod
    def set_uniform_concentration(cls, conc: float, species_index=None, species_name=None) -> None:
        """
        Assign the given concentration to all the bins of the specified species (identified by its index or name.)
        Any previous values get over-written

        :param conc:            The desired value of chemical concentration for the above species
        :param species_index:   Zero-based index to identify a specific chemical species
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :return:                None
        """
        if species_name is not None:
            assert cls.chem_data, f"BioSim1D.set_uniform_concentration(): must first call initialize_system()"
            species_index = cls.chem_data.get_index(species_name)
        else:
            cls.chem_data.assert_valid_index(species_index)

        assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"

        cls.system[species_index] = np.full(cls.n_bins, conc, dtype=float)

        if cls.uses_membranes():
            # If membranes are present, also set "side B" of the bins with membranes, to the same concentration
            # TODO: see if it can be done as vector operation
            for bin_address in cls.bins_with_membranes():
                if cls.membranes[bin_address]:
                    cls.system_B[species_index, bin_address] = conc



    @classmethod
    def set_all_uniform_concentrations(cls, conc_list: Union[list, tuple]) -> None:
        """
        Set the concentrations of all species at once, uniformly across all bins
        :param conc_list:   List or tuple of concentration values for each of the chemical species
        :return:            None
        """
        assert len(conc_list) == cls.chem_data.n_species, \
            f"set_all_uniform_concentrations(): the argument must be a list or tuple of size {cls.chem_data.n_species}"

        for i, conc in enumerate(conc_list):
            cls.set_uniform_concentration(species_index=i, conc=conc)



    @classmethod
    def set_bin_conc(cls, bin_address: int, conc: float, species_index=None, species_name=None,
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
            assert cls.chem_data, f"BioSim1D.set_bin_conc(): must first call initialize_system()"
            species_index = cls.chem_data.get_index(species_name)
        else:
            cls.chem_data.assert_valid_index(species_index)

        cls.assert_valid_bin(bin_address)


        assert conc >= 0., \
            f"set_bin_conc(): the concentration must be a positive number or zero (the requested value was {conc})"

        if across_membrane or both_sides:
            assert cls.system_B is not None, \
                "set_bin_conc(): the `other_side` option cannot be used unless membranes are first set"
            cls.system_B[species_index, bin_address] = conc

        if both_sides or (not across_membrane):
            cls.system[species_index, bin_address] = conc



    @classmethod
    def set_species_conc(cls, conc_list: Union[list, tuple, np.ndarray], species_index=None, species_name=None) -> None:
        """
        Assign the requested list of concentration values to all the bins, in order, for the specified species.

        :param conc_list:       A list, tuple or Numpy array with the desired concentration value to assign to the specified location
        :param species_index:   Zero-based index to identify a specific chemical species
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :return:                None
        """
        if species_name is not None:
            assert cls.chem_data, f"BioSim1D.set_species_conc(): must first call initialize_system()"
            species_index = cls.chem_data.get_index(species_name)
        else:
            cls.chem_data.assert_valid_index(species_index)

        assert (type(conc_list) == list) or (type(conc_list) == tuple) or (type(conc_list) == np.ndarray), \
            f"set_species_conc(): the argument `conc_list` must be a list, tuple or Numpy array; the passed value was {type(conc_list)})"

        assert len(conc_list) == cls.n_bins, \
            "set_species_conc(): the argument `conc_list` must be a list of concentration values for ALL the various bins (wrong length)"

        cls.system[species_index] = conc_list



    @classmethod
    def inject_conc_to_bin(cls, bin_address: int, species_index: int, delta_conc: float, zero_clip = False) -> None:
        """
        Add the requested concentration to the cell with the given address, for the specified chem species

        :param bin_address:     The zero-based bin number of the desired cell
        :param species_index:   Zero-based index to identify a specific chemical species
        :param delta_conc:      The concentration to add to the specified location
        :param zero_clip:       If True, any requested increment causing a concentration dip below zero, will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        cls.assert_valid_bin(bin_address)

        if (cls.system[species_index, bin_address] + delta_conc) < 0. :
            if zero_clip:
                cls.system[species_index, bin_address] = 0
                return
            else:
                raise Exception("The requested concentration change would result in a negative final value")

        # Normal scenario, not leading to negative values for the final concentration
        cls.system[species_index, bin_address] += delta_conc



    @classmethod
    def inject_sine_conc(cls, species_name, frequency, amplitude, bias=0, phase=0, zero_clip = False) -> None:
        """
        Add to the concentrations of the specified chem species a sinusoidal signal across all bins

        Note:   A sine wave of the form   A sin(B x - C)
                has an amplitude of A, a period of 2Pi/B and a right phase shift of C (in radians)

        In Mathematica:  Plot[Sin[B x - C] /. {B -> 2 Pi, C -> 0} , {x, 0, 1}, GridLines -> Automatic]

        :param species_name:    The name of the chemical species whose concentration we're modifying
        :param frequency:       Number of waves along the length of the system
        :param amplitude:       Amplitude of the Sine wave.  Note that peak-to-peak values are double the amplitude
        :param bias:            Amount to be added to all values (akin to "DC bias" in electrical circuits)
        :param phase:           In degrees: phase shift to the RIGHT.  EXAMPLE: 180 to flip the Sine curve
        :param zero_clip:       If True, any requested increment causing a concentration dip below zero,
                                will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        assert cls.chem_data, f"BioSim1D.inject_sine_conc(): must first call initialize_system()"
        species_index = cls.chem_data.get_index(species_name)

        period = cls.n_bins / frequency
        #print("period: ", period)

        phase_radians = phase * math.pi / 180

        B = 2*math.pi / period
        #print("B: ", B)

        for x in range(cls.n_bins):
            conc = amplitude * math.sin(B*x - phase_radians) + bias
            #print(conc)
            cls.inject_conc_to_bin(bin_address = x, species_index = species_index,
                                   delta_conc = conc, zero_clip = zero_clip)



    @classmethod
    def frequency_analysis(cls, species_name: str, threshold = 0.001, n_largest = None) -> pd.DataFrame:
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

        conc_samples = cls.lookup_species(species_name=species_name)

        #import plotly.express as px
        #fig = px.line(y=conc_samples)
        #fig.show()

        # Perform a DFT, to extract the frequency components
        # The size of the computed arrays xf and yf is  (cls.n_bins/2 + 1) if cls.n_bins is even
        # or (cls.n_bins + 1) /2 if cls.n_bins is odd
        xf = rfftfreq(cls.n_bins, 1 / cls.n_bins)
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




    ########  DIMENSION-RELATED  ################

    @classmethod
    def set_dimensions(cls, length) -> None:
        """
        Set the overall length of the system.
        Doing so, will permit to convert bin numbers to positional values

        :param length:
        :return:
        """
        assert (type(length) == float) or (type(length) == int), "set_dimensions(): length must be a number"
        assert length > 0, "set_dimensions(): length must be positive"

        cls.system_length = length
        cls.global_Dx = length / (cls.n_bins-1)


    @classmethod
    def x_coord(cls, bin_address):
        """"
        Return the x coordinate of the middle of the specified bin.
        By convention, for the leftmost bin, it's zero,
        and for the rightmost, it's the overall length of the system
        """
        assert cls.system_length, "x_coord(): must first call set_dimensions()"

        return bin_address * cls.global_Dx




    ########  MEMBRANE-RELATED  ################

    @classmethod
    def uses_membranes(cls) -> bool:
        """
        Return True if membranes are part of the system

        :return:
        """
        return cls.membranes is not None


    @classmethod
    def bins_with_membranes(cls) -> [int]:
        """

        :return:
        """
        # Create and return a list of bin numbers whose entry in cls.membranes is True
        return [bin_number for bin_number in range(cls.n_bins)
                if cls.membranes[bin_number]]



    @classmethod
    def set_membranes(cls, membrane_pos: Union[List, Tuple]) -> None:
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
        cls.membranes = np.zeros(cls.n_bins, dtype=bool)
        cls.A_fraction = np.zeros(cls.n_bins, dtype=float)
        if cls.system_B is None:
            cls.system_B = np.zeros((cls.n_species, cls.n_bins), dtype=float)


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

            cls.assert_valid_bin(bin_number)
            cls.membranes[bin_number] = True
            cls.A_fraction[bin_number] = fraction

            # Equalize all the concentrations across the compartments on both sides of the membrane
            if cls.chem_data:
                for chem_index in range(cls.chem_data.number_of_chemicals()):
                    # TODO: do as vector operation
                    conc = cls.bin_concentration(bin_address=bin_number, species_index=chem_index)
                    cls.set_bin_conc(bin_address=bin_number, species_index=chem_index, conc=conc, across_membrane=True)




    #########################################################################
    #                                                                       #
    #                              TO VIEW                                  #
    #                                                                       #
    #########################################################################

    @classmethod
    def assert_valid_bin(cls, bin_address: int) -> None:
        """
        Raise an Exception if the given bin number isn't valid

        :param bin_address:  An integer that ought to be between 0 and (cls.n_bins-1), inclusive
        :return:            None
        """
        assert type(bin_address) == int, \
            f"BioSim1D: the requested bin address ({bin_address}) is not an integer; its type is {type(bin_address)}"

        if bin_address < 0 or bin_address >= cls.n_bins:
            raise Exception(f"BioSim1D: the requested bin address ({bin_address}) is out of bounds for the system size; "
                            f"allowed range is [0-{cls.n_bins-1}], inclusive")




    @classmethod
    def lookup_species(cls, species_index=None, species_name=None, trans_membrane=False) -> np.array:
        """
        Return the NumPy array of concentration values across the all bins (from left to right)
        for the single specified chemical species

        :param species_index:   The index order of the chemical species of interest
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :param trans_membrane:  It True, consider only the "other side" of the bins, i.e. the portion across the membrane
                                    (it will be zero for bins without membrane)
        :return:                A NumPy array of concentration values across the bins (from left to right);
                                    the size of the array is the number of bins
        """
        if species_name is not None:
            assert cls.chem_data, f"BioSim1D.lookup_species(): must first call initialize_system()"
            species_index = cls.chem_data.get_index(species_name)

        cls.chem_data.assert_valid_index(species_index)

        if trans_membrane:
            return cls.system_B[species_index]
        else:
            return cls.system[species_index]



    @classmethod
    def bin_concentration(cls, bin_address: int, species_index=None, species_name=None, trans_membrane=False) -> float:
        """
        Return the concentration at the requested bin of the specified species

        :param bin_address:     The bin number
        :param species_index:   The index order of the chemical species of interest
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :param trans_membrane:  If True, consider the "other side" of the bin, i.e. the portion across the membrane
        :return:                A concentration value at the indicated bin, for the requested species
        """
        if species_name is not None:
            assert cls.chem_data, f"BioSim1D.bin_concentration(): must first call initialize_system()"
            species_index = cls.chem_data.get_index(species_name)

        cls.chem_data.assert_valid_index(species_index)

        if trans_membrane:
            return cls.system_B[species_index, bin_address]
        else:
            return cls.system[species_index, bin_address]



    @classmethod
    def bin_snapshot(cls, bin_address: int) -> dict:
        """
        Extract the concentrations of all the chemical species at the specified bin,
        as a dict whose keys are the names of the species
        EXAMPLE:  {'A': 10.0, 'B': 50.0}

        :param bin_address: An integer with the bin number
        :return:            A dict of concentration values; the keys are the names of the species
        """
        assert type(bin_address) == int, "bin_snapshot(): the argument must be an integer"

        d = {}
        for species_index in range(cls.n_species):
            name = cls.chem_data.get_name(species_index)
            conc = cls.bin_concentration(bin_address, species_index)
            d[name] = conc

        return d



    @classmethod
    def show_system_snapshot(cls) -> None:
        """"
        """
        print(f"SYSTEM SNAPSHOT at time {cls.system_time}:")
        print(cls.system_snapshot())



    @classmethod
    def describe_state(cls, concise=False) -> None:
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
        print(f"SYSTEM STATE at Time t = {cls.system_time}:")

        if concise:             # A minimalist printout...
            print(cls.system)   # ...only showing the concentration data (a Numpy array)
            if cls.uses_membranes():
                print("Membranes:")
                print(cls.system_B)
            return

        # If we get thus far, it's a FULL printout

        print(f"{cls.n_bins} bins and {cls.n_species} species:")

        # Show a line of line of data for each chemical species in turn
        for species_index in range(cls.n_species):
            name = cls.chem_data.get_name(species_index)
            if name:    # If a name was provided, show it
                name = f" ({name})"
            else:
                name = ""

            if cls.membranes is None:
                all_conc = cls.system[species_index]
            else:
                all_conc = "|"
                for bin_no in range(cls.n_bins):
                    all_conc += str(cls.bin_concentration(bin_no, species_index=species_index))
                    if cls.membranes[bin_no]:
                        # Add a symbol for the membrane, and the additional membrane concentration data (on the "other side")
                        all_conc += "()"    # To indicate a membrane
                        all_conc += str(cls.bin_concentration(bin_no, species_index=species_index, trans_membrane=True))
                    all_conc += "|"

            if cls.chem_data.diffusion_rates is None:
                print(f"  Species {species_index}{name}. Diff rate: NOT SET. Conc: {all_conc}")
            else:
                print(f"  Species {species_index}{name}. Diff rate: {cls.chem_data.diffusion_rates[species_index]}. Conc: {all_conc}")




    @classmethod
    def show_membranes(cls, n_decimals=1) -> str:
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
        if cls.membranes is None:
            print("No membranes present.  Call set_membranes() to set them")
            return ""

        # Prepare the middle line
        box_contents = "|"
        for bin_no, val in enumerate(cls.membranes):
            if val:
                fraction = round(cls.A_fraction[bin_no], n_decimals)
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

    @classmethod
    def increase_spacial_resolution(cls, factor:int) -> None:
        """
        Increase the spacial resolution of the system by cloning and repeating
        each bin, by the specified number of times.
        Replace the System's internal state.

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
        new_state = np.repeat(cls.system, factor, axis=1)
        cls.replace_system(new_state)



    @classmethod
    def decrease_spacial_resolution(cls, factor:int) -> None:
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
        assert cls.n_bins % factor == 0, f"The number of bins (currently {cls.n_bins}) must be a multiple of the requested scaling factor"

        reduced_n_bins = int(cls.n_bins / factor)
        # The result matrix will have the same number of chemical species, but fewer bins
        new_state = np.zeros((cls.n_species, reduced_n_bins), dtype=float)

        for i in range(reduced_n_bins):
            start_col = factor * i      # The start column will initially be 0, and will get incremented by factor
            col_group = cls.system[ : , start_col:start_col+factor] # Extract a submatrix containing the number of columns specified by "factor",
                                                                    # starting with the column specified by start_col
            compressed_col_group = np.sum(col_group, axis=1, keepdims=True) / factor    # Create a single column that is the average of the columns in the group
            new_state[:, [i]] = compressed_col_group                                       # Store the newly-computed column of averages in the appropriate place
                                                                                        # in the result matrix
        cls.replace_system(new_state)



    @classmethod
    def smooth_spacial_resolution(cls) -> None:
        """
        EXAMPLE: if the system is
                        [[10., 20., 30.]
                         [ 2., 8.,   4.]]
                then the result will be
                        [[10., 15., 20., 25., 30.]
                         [ 2.,  5.,  8.,  6.,  4.]]

        :return:
        """
        n_bins = cls.n_bins
        new_n_bins = n_bins * 2 - 1     # The final number of bins
        new_state = np.zeros((cls.n_species, new_n_bins), dtype=float)

        for start_col in range(n_bins-1):                           # The start column will be between 0 and (cls.n_bins-2), inclusive
            col_group = cls.system[ : , start_col:start_col+2]      # Extract a submatrix containing 2 columns,
                                                                    # starting with the one in position start_col
            avg_col = np.sum(col_group, axis=1, keepdims=True) / 2. # Create a single column that is the average of the columns in the group

            new_state[:, [2*start_col]] = cls.system[ : , start_col:start_col+1]   # Set one column, from the first column in the group
            new_state[:, [2*start_col+1]] = avg_col                                # Set the next column with the newly-computed column of averages


        new_state[ : , -1:] = cls.system[ : , -1:]                 # Set the last column of result to the last column of cls.system

        cls.replace_system(new_state)





    #########################################################################
    #                                                                       #
    #                               DIFFUSION                               #
    #                                                                       #
    #########################################################################

    @classmethod
    def is_excessive(cls, time_step, diff_rate) -> bool:
        """
        Use a loose heuristic to determine if the requested time step is too long,
        given the diffusion rate

        :param time_step:
        :param diff_rate:
        :return:
        """
        value = time_step * diff_rate

        if value > cls.time_step_threshold:
            return True
        else:
            return False


    @classmethod
    def max_time_step(cls, diff_rate) -> float:
        """
        Determine a reasonable upper bound on the time step, for the given diffusion rate
        :param diff_rate:
        :return:
        """
        return cls.time_step_threshold/diff_rate



    @classmethod
    def diffuse(cls, total_duration=None, time_step=None, n_steps=None, delta_x = 1) -> dict:
        """
        Uniform-step diffusion, with 2 out of 3 criteria specified:
            1) until reaching, or just exceeding, the desired time duration
            2) using the given time step
            3) carrying out the specified number of steps

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param delta_x:         Distance between consecutive bins
        :return:                A dictionary with data about the status of the operation
        """
        time_step, n_steps = cls.all_reactions.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        for i in range(n_steps):
            if cls.debug:
                if (i < 2) or (i >= n_steps-2):
                    print(f"    Performing diffusion step {i}...")
                elif i == 2:
                    print("    ...")

            cls.diffuse_step(time_step, delta_x=delta_x)
            cls.system += cls.delta_diffusion     # Matrix operation to update all the concentrations
            cls.system_time += time_step

        if cls.debug:
            print(f"\nSystem after Delta time {total_duration}, at end of {n_steps} steps of size {time_step}:")
            cls.describe_state(concise=True)
            print()

        status = {"steps": n_steps}
        return status



    @classmethod
    def diffuse_step(cls, time_step, delta_x = 1) -> None:
        """
        Diffuse all the species by the given time step:
        clear the delta_diffusion array, and then re-compute it from all the species.

        IMPORTANT: the actual system concentrations are NOT changed.

        :param time_step:   Time step over which to carry out the diffusion
                            If too large - as determined by the method is_excessive() - an Exception will be raised
        :param delta_x:     Distance between consecutive bins
        :return:            None
        """
        # TODO: parallelize the independent computations

        cls.delta_diffusion = np.zeros((cls.n_species, cls.n_bins), dtype=float)

        for species_index in range(cls.n_species):

            increment_vector = cls.diffuse_step_single_species(time_step, species_index=species_index, delta_x=delta_x)
            #print("Increment vector is: ", increment_vector)

            # For each bin, update the concentrations from the buffered increments
            cls.delta_diffusion[species_index] = increment_vector      # Vector operation to a row of the matrix delta_diffusion



    @classmethod
    def diffuse_step_single_species(cls, time_step: float, species_index=0, delta_x = 1) -> np.array:
        """
        Diffuse the specified single species, for the specified time step, across all bins,
        and return an array of the changes in concentration ("Delta concentration")
        for the given species across all bins.

        IMPORTANT: the actual system concentrations are NOT changed.

        We're assuming an isolated environment, with nothing diffusing thru the "walls"

        EXPLANATION:  https://life123.science/diffusion

        :param time_step:       Time step over which to carry out the diffusion.
                                If too large, an Exception will be raised.
        :param species_index:   ID (in the form of an integer index) of the chemical species under consideration
        :param delta_x:         Distance between consecutive bins

        :return:                A Numpy array with the change in concentration for the given species across all bins
        """
        assert cls.system is not None, "Must first initialize the system"
        assert cls.n_bins > 0, "Must first set the number of bins"
        assert cls.chem_data.diffusion_rates is not None, "Must first set the diffusion rates"
        assert cls.sealed == True, "For now, there's no provision for exchange with the outside"

        increment_vector = np.zeros(cls.n_bins, dtype=float)   # One element per bin

        if cls.n_bins == 1:
            return increment_vector                 # There's nothing to do in the case of just 1 cell!

        diff = cls.chem_data.diffusion_rates[species_index]   # The diffusion rate of the specified single species

        assert not cls.is_excessive(time_step, diff), f"Excessive large time_fraction. Should be < {cls.max_time_step(diff)}"


        # Carry out a convolution operation, with a tile of size 3 (or 2 if only 2 bins)
        #print(f"Diffusing species # {species_index}")

        max_bin_number = cls.n_bins - 1     # Bin numbers range from 0 to max_bin_number, inclusive

        effective_diff = diff * time_step
        if delta_x != 1:
            effective_diff /= (delta_x**2)


        for i in range(cls.n_bins):    # Bin number, ranging from 0 to max_bin_number, inclusive
            #print(f"Processing bin number {i}")
            current_conc = cls.system[species_index , i]

            if i == 0 :                     # Special case for the first bin (no left neighbor)
                increment_vector[i] = effective_diff * (cls.system[species_index , 1] - current_conc)
            elif i == max_bin_number :      # Special case for the last bin (no right neighbor)
                increment_vector[i] = effective_diff * (cls.system[species_index , i - 1] - current_conc)
            else:
                increment_vector[i] = effective_diff * \
                                        (cls.system[species_index , i + 1] - current_conc
                                         + cls.system[species_index , i - 1] - current_conc)

        return increment_vector




    #########################################################################
    #                                                                       #
    #                               REACTIONS                               #
    #                                                                       #
    #########################################################################

    @classmethod
    def set_reactions(cls, reactions) -> None:
        """
        Register the given set of reactions.
        Anything previously registered, will be over-written.

        :param reactions:   Object of class "Reactions"
        :return:            None
        """
        cls.all_reactions = reactions



    @classmethod
    def react(cls, total_duration=None, time_step=None, n_steps=None, snapshots=None) -> None:
        """
        Update the system concentrations as a result of all the reactions in all bins - taking
        the presence of membranes into account, if applicable.
        CAUTION : NO diffusion is performed.

        For each bin, or each membrane-separated side of bin, (or combined group of bins - not currently implemented),
        process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        TODO: in case of any Exception, the state of the system is still valid, as of the time before this call

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param snapshots:       OPTIONAL dict that may contain any the following keys:
                                    "frequency", "sample_bin", "sample_species"
                                    If provided, take a system snapshot after running a multiple of "frequency" runs.
                                    EXAMPLE: snapshots={"frequency": 2, "sample_bin": 0}
        :return:                None
        """
        time_step, n_steps = cls.all_reactions.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        # TODO: validation; also, implement "sample_species" option for snapshots
        if snapshots is None:
            frequency = None
        else:
            frequency = snapshots.get("frequency", 1)

        for i in range(n_steps):
            cls.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            # Update the state of the system
            cls.system += cls.delta_reactions           # Matrix operation to update all the concentrations
            if cls.uses_membranes():
                cls.system_B += cls.delta_reactions_B   # Matrix operation to update all the concentrations

            cls.system_time += time_step
            # Preserve some of the data, as requested
            if (frequency is not None) and ((i+1)%frequency == 0):
                cls.save_snapshot(cls.bin_snapshot(bin_address = snapshots["sample_bin"]))



    @classmethod
    def reaction_step(cls, delta_time: float) -> None:
        """
        Clear and compute the delta_reactions array (a class variable),
        based on all the reactions in all bins.
        IMPORTANT: the actual system concentrations are NOT changed.

        For each bin, process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        TODO: parallelize the computation over the separate bins
        TODO: explore looping over reactions first, and then over bins

        :param delta_time:
        :return:            None
        """
        assert cls.all_reactions is not None, \
            "reaction_step(): must first set the Reactions object"

        #number_reactions = cls.all_reactions.number_of_reactions()

        cls.delta_reactions = np.zeros((cls.n_species, cls.n_bins), dtype=float)
        if cls.uses_membranes():
            cls.delta_reactions_B = np.zeros((cls.n_species, cls.n_bins), dtype=float)

        # For each bin
        for bin_n in range(cls.n_bins):     # Bin number, ranging from 0 to max_bin_number, inclusive
            if cls.verbose:
                print(f"reaction_step(): processing the all the reactions in bin number {bin_n}")

            # Obtain the Delta-concentration for each species, for this bin
            conc_dict = {species_index: cls.system[species_index , bin_n]
                         for species_index in range(cls.n_species)}
            #print(f"\nconc_dict in bin {bin_n}: ", conc_dict)

            # Obtain the Delta-conc for each species, for bin number bin_n (a NumPy array)
            increment_vector = cls.all_reactions.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=delta_time)

            # Replace the "bin_n" column of the cls.delta_reactions matrix with the contents of the vector increment_vector
            cls.delta_reactions[:, bin_n] = np.array([increment_vector])

            #print(cls.delta_reactions)


        # Now process the "other side" of membrane-containing bins, if applicable
        if cls.uses_membranes():
            for bin_n in cls.bins_with_membranes():
                # Obtain the Delta-concentration for each species, for this bin
                conc_dict = {species_index: cls.bin_concentration(bin_address=bin_n, species_index=species_index, trans_membrane=True)
                             for species_index in range(cls.n_species)}
                #print(f"\n Post-membrane side conc_dict in bin {bin_n}: ", conc_dict)

                # Obtain the Delta-conc for each species, for bin number bin_n (a NumPy array)
                increment_vector = cls.all_reactions.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=delta_time)

                # Replace the "bin_n" column of the cls.delta_reactions_B matrix with the contents of the vector increment_vector
                cls.delta_reactions_B[:, bin_n] = np.array([increment_vector])




    #########################################################################
    #                                                                       #
    #                         REACTION-DIFFUSION                            #
    #                                                                       #
    #########################################################################

    @classmethod
    def react_diffuse(cls, total_duration=None, time_step=None, n_steps=None, delta_x = 1) -> None:
        """
        It expects 2 of the arguments:  total_duration, time_step, n_steps
        Perform a series of reaction and diffusion time steps.

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param delta_x:         Distance between consecutive bins
        :return:                None
        """
        time_step, n_steps = cls.all_reactions.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        for i in range(n_steps):
            # TODO: split off the reaction step and the diffusion step to 2 different computing cores
            cls.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            cls.diffuse_step(time_step, delta_x=delta_x)
            # Merge into the concentrations of the various bins/chemical species pairs,
            # the increments concentrations computed separately by the reaction and the diffusion steps
            cls.system += cls.delta_reactions   # Matrix operation to update all the concentrations
                                                #   from the reactions
            cls.system += cls.delta_diffusion   # Matrix operation to update all the concentrations
                                                #   from the diffusion
            cls.system_time += time_step



    #########################################################################
    #                                                                       #
    #                           VISUALIZATION                               #
    #                                                                       #
    #########################################################################

    @classmethod
    def single_species_heatmap(cls, species_index: int, heatmap_pars: dict, graphic_component, header=None) -> None:
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
        species_concentrations = list(cls.lookup_species(species_index))
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
        # The version of the heatmap Vue component specified in the call to GraphicLog.config() will be used
        GraphicLog.export_plot(all_data, graphic_component)



    @classmethod
    def single_species_line_plot(cls, species_index: int, plot_pars: dict, graphic_component, header=None) -> None:
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
        species_concentrations = list(cls.lookup_species(species_index))
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
        # The version of the heatmap Vue component specified in the call to GraphicLog.config() will be used
        GraphicLog.export_plot(all_data, graphic_component)



    @classmethod
    def line_plot(cls, plot_pars: dict, graphic_component, header=None, color_mapping=None) -> None:
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
        species_concentrations = [list(cls.lookup_species(i)) for i in range(cls.n_species)]
        #print(species_concentrations)

        all_data = {
            "curve_labels": cls.chem_data.get_all_names(),

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
        # The version of the heatmap Vue component specified in the call to GraphicLog.config() will be used
        GraphicLog.export_plot(all_data, graphic_component)




    #########################################################################
    #                                                                       #
    #                                HISTORY                                #
    #                                                                       #
    #########################################################################


    @classmethod
    def save_snapshot(cls, data_snapshot, caption = "") -> None:
        """

        :param data_snapshot:
        :param caption:
        :return:
        """
        cls.history.append(pars=cls.system_time,
                           data_snapshot = data_snapshot, caption=caption)


    @classmethod
    def get_history(cls):
        """

        :return:
        """
        return cls.history.get()
