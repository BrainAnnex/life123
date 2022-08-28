import numpy as np
import math
import pandas as pd
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

    system = None       # Concentration data in the System we're simulating, for all the chemicals
                        #   NumPy array of reals, of dimension: (n_species x n_bins).
                        #   Each row represents a species

    system_earlier = None   # NOT IN CURRENT USE.  Envisioned for simulations where the past 2 time states are used
                            # to compute the state at the next time step

    membranes = None        # NumPy boolean array of dimension n_bins;
                            #   each boolean value marks the presence of a membrane in that bin.
                            #   Any given bin is expected to either contain, or not contain, a single membrane passing thru it

    A_fraction = None       # NumPy array of reals, of dimension n_bins;
                            #   each value records the fraction (between 0. and 1.) of the LEFT one
                            #   of the 2 parts into which the membrane splits the bin.
                            #   0. if N/A

    system_B = None         # Just like "system", but only for the "A" fractions of the bin, where applicable;
                            #   0. if N/A

    delta_diffusion = None  # Buffer for the concentration changes from diffusion step (n_species x n_bins)
    delta_reactions = None  # Buffer for the concentration changes from reactions step (n_species x n_bins)

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
    #                    SET/MODIFY CONCENTRATIONS (or membranes)           #
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

        # Initialize all concentrations to zero
        cls.system = np.zeros((cls.n_species, n_bins), dtype=float)

        cls.system_time = 0             # "Start the clock"




    @classmethod
    def replace_system(cls, new_state: np.array) -> None:
        """
        Replace the System's internal state.
        For details of the data structure, see the class variable "system"

        :param new_state:   Numpy array containing the desired new System's internal state
        :return:
        """
        cls.system = new_state
        cls.n_species, cls.n_bins = new_state.shape     # Extract from the new state
        assert cls.n_species == cls.chem_data.n_species, \
            "replace_system(): inconsistency in the number of chemical species vs. the Chemicals object"



    @classmethod
    def set_uniform_concentration(cls, conc: float, species_index=None, species_name=None) -> None:
        """
        Assign the given concentration to all the bin of the specified species (identified by its index).
        Any previous values get over-written

        :param conc:            The desired value of chemical concentration for the above species
        :param species_index:   Zero-based index to identify a specific chemical species
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :return:                None
        """
        if species_name is not None:
            assert cls.chem_data, f"BioSim1D.set_uniform_concentration(): must first call initialize_system()"
            species_index = cls.chem_data.get_index(species_name)

        cls.chem_data.assert_valid_index(species_index)

        assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"

        cls.system[species_index] = np.full(cls.n_bins, conc, dtype=float)


    @classmethod
    def set_all_uniform_concentrations(cls, conc_list: [float]) -> None:
        """
        Set the concentrations of all species at once, uniformly across all bins
        :param conc_list:   List of concentration values for each of the chemical species
        :return:            None
        """
        assert len(conc_list) == cls.chem_data.n_species, \
            f"set_all_uniform_concentrations(): the argument must be a list of size {cls.chem_data.n_species}"

        for i, conc in enumerate(conc_list):
            cls.set_uniform_concentration(species_index=i, conc=conc)



    @classmethod
    def set_bin_conc(cls, bin: int, species_index: int, conc: float) -> None:
        """
        Assign the requested concentration value to the bin with the given index, for the specified species

        :param bin:             The zero-based bin number of the desired compartment
        :param species_index:   Zero-based index to identify a specific chemical species
        :param conc:            The desired concentration value to assign to the specified location
        :return:                None
        """
        assert bin < cls.n_bins, f"set_bin_conc(): the requested cell index ({bin}) must be in the range [0 - {cls.n_bins - 1}], inclusive"
        cls.chem_data.assert_valid_index(species_index)

        assert conc >= 0., f"set_bin_conc(): the concentration must be a positive number or zero (the requested value was {conc})"

        cls.system[species_index, bin] = conc


    @classmethod
    def set_species_conc(cls, species_index: int, conc_list: Union[list, tuple, np.ndarray]) -> None:
        """
        Assign the requested list of concentration values to all the bins, in order, for the specified species

        :param species_index:   Zero-based index to identify a specific chemical species
        :param conc_list:       A list, tuple or Numpy array with the desired concentration value to assign to the specified location
        :return:                None
        """
        cls.chem_data.assert_valid_index(species_index)
        assert (type(conc_list) == list) or (type(conc_list) == tuple) or (type(conc_list) == np.ndarray), \
            f"set_species_conc(): the argument `conc_list` must be a list, tuple or Numpy array; the passed value was {type(conc_list)})"

        assert len(conc_list) == cls.n_bins, \
            "set_species_conc(): the argument `conc_list` must be a list of concentration values for ALL the various bins (wrong length)"

        cls.system[species_index] = conc_list



    @classmethod
    def inject_conc_to_bin(cls, bin: int, species_index: int, delta_conc: float, zero_clip = True) -> None:
        """
        Add the requested concentration to the cell with the given index, for the specified species

        :param bin:             The zero-based bin number of the desired cell
        :param species_index:   Zero-based index to identify a specific chemical species
        :param delta_conc:      The concentration to add to the specified location
        :param zero_clip:       If True, any requested increment causing a concentration dip below zero, will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        assert bin < cls.n_bins, f"The requested cell index ({bin}) must be in the range [0 - {cls.n_bins - 1}]"

        if (cls.system[species_index, bin] + delta_conc) < 0. :
            if zero_clip:
                cls.system[species_index, bin] = 0
            else:
                raise Exception("The requested concentration change would result in a negative final value")

        # Normal scenario, not leading to negative values for the final concentration
        cls.system[species_index, bin] += delta_conc



    ########  EXPERIMENTAL!  MEMBRANE-RELATED  ################

    @classmethod
    def set_membranes(cls, membrane_pos: Union[List[int], Tuple[int]]) -> None:
        """
        Set the class variable "membranes", a NumPy boolean array of dimension n_cells;
        each element marks the presence of a membrane in that bin.

        IMPORTANT: any previously-set membrane information is lost.

        :param membrane_pos:    A list or tuple of indexes of bins that contain membranes
        :return:                None
        """
        cls.membranes = np.zeros(cls.n_bins, dtype=bool)

        # Loop over all the values passed as a list or tuple
        for bin_number in membrane_pos:
            cls.assert_valid_bin_number(bin_number)
            cls.membranes[bin_number] = True




    #########################################################################
    #                                                                       #
    #                              TO VIEW                                  #
    #                                                                       #
    #########################################################################

    @classmethod
    def assert_valid_bin_number(cls, bin_number: int) -> None:
        """
        Raise an Exception if the given bin number isn't valid
        :param bin_number:  An integer that ought to be between 0 and (cls.n_bins-1), inclusive
        :return:            None
        """
        if bin_number < 0 or bin_number >= cls.n_bins:
            raise Exception(f"BioSim1D: the requested bin number ({bin_number}) is out of bounds for the system size; "
                            f"allowed range is [0-{cls.n_bins-1}], inclusive")




    @classmethod
    def lookup_species(cls, index: int) -> np.array:
        """
        Return the NumPy array of concentration values across the bins (from left to right)
        for the single specified species

        :param index:   The index order of the chemical species of interest
        :return:        A NumPy array of concentration values across the bins (from left to right)
        """
        cls.chem_data.assert_valid_index(index)
        return cls.system[index]


    @classmethod
    def bin_concentration(cls, bin_address: int, species_index: int) -> float:
        """
        Return the concentration at the requested bin of the specified species

        :param bin_address:     The bin index
        :param species_index:   The index order of the chemical species of interest
        :return:                A concentration value at the indicated bin, for the requested species
        """
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
    def system_snapshot(cls) -> pd.DataFrame:
        """
        Return a snapshot of all the concentrations of all the species,
        as a Pandas dataframe

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
    def show_system_snapshot(cls) -> None:
        """"
        """
        print(f"SYSTEM SNAPSHOT at time {cls.system_time}:")
        print(cls.system_snapshot())



    @classmethod
    def describe_state(cls, concise=False) -> None:
        """
        A printout of the state of the system, for now useful only for small systems

        :param concise: If True, only produce a minimalist printout with just the concentration values
        :return:        None
        """
        print(f"SYSTEM STATE at Time t = {cls.system_time}:")

        if concise:             # A minimalist printout...
            print(cls.system)   # ...only showing the concentration data (a Numpy array)
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

            if cls.chem_data.diffusion_rates is None:
                print(f"  Species {species_index}{name}. Diff rate: NOT SET. Conc: ", cls.system[species_index])
            else:
                print(f"  Species {species_index}{name}. Diff rate: {cls.chem_data.diffusion_rates[species_index]}. Conc: ",
                      cls.system[species_index])




    @classmethod
    def show_membranes(cls) -> str:
        """
        A simple-minded early method to visualize where the membranes are.

        EXAMPLE (with 2 membrane on the right part of a 5-bin system):
                    ___________
                    | | |*| |*|
                    -----------
        :return:
        """
        box_width = 2 * cls.n_bins + 1

        box = "\n"
        box += "_" * box_width + "\n"   # The top of the box

        # Prepare the middle line
        box_contents = "|"
        for val in cls.membranes:
            if val:
                box_contents += "*|"
            else:
                box_contents += " |"

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
    def increase_spacial_resolution(cls, factor:int) -> np.array:
        """
        Increase the spacial resolution of the system by cloning and repeating
        each bin, by the specified number of times

        EXAMPLE: if the system is
                        [[11. 12. 13.]
                         [ 5. 15. 25.]]
                and factor=2, then the result will be
                        [[11. 11. 12. 12. 13. 13.]
                         [ 5.  5. 15. 15. 25. 25.]]

        :param factor:  Number of bins into which to split each bin (replicating their concentration values)
        :return:        The derived system state
        """
        assert type(factor) == int, "The argument `factor` must be an integer"
        return np.repeat(cls.system, factor, axis=1)


    @classmethod
    def decrease_spacial_resolution(cls, factor:int) -> np.array:
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
        result = np.zeros((cls.n_species, reduced_n_bins), dtype=float)

        for i in range(reduced_n_bins):
            start_col = factor * i      # The start column will initially be 0, and will get incremented by factor
            col_group = cls.system[ : , start_col:start_col+factor] # Extract a submatrix containing the number of columns specified by "factor",
                                                                    # starting with the column specified by start_col
            compressed_col_group = np.sum(col_group, axis=1, keepdims=True) / factor    # Create a single column that is the average of the columns in the group
            result[:, [i]] = compressed_col_group                                       # Store the newly-computed column of averages in the appropriate place
                                                                                        # in the result matrix

        return result



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
    def diffuse(cls, total_duration=None, time_step=None, n_steps=None) -> dict:
        """
        Uniform-step diffusion, with 2 out of 3 criteria specified:
            1) until reaching, or just exceeding, the desired time duration
            2) using the given time step
            3) carrying out the specified number of steps

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
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

            cls.diffuse_step(time_step)
            cls.system += cls.delta_diffusion     # Matrix operation to update all the concentrations
            cls.system_time += time_step

        if cls.debug:
            print(f"\nSystem after Delta time {total_duration}, at end of {n_steps} steps of size {time_step}:")
            cls.describe_state(concise=True)
            print()

        status = {"steps": n_steps}
        return status



    @classmethod
    def diffuse_step(cls, time_step) -> None:
        """
        Diffuse all the species by the given time step:
        clear the delta_diffusion array, and then re-compute it from all the species.

        IMPORTANT: the actual system concentrations are NOT changed.

        :param time_step:   Time step over which to carry out the diffusion.
                            If too large, an Exception will be raised.
        :return:            None
        """
        # TODO: parallelize the independent computations

        cls.delta_diffusion = np.zeros((cls.n_species, cls.n_bins), dtype=float)

        for species_index in range(cls.n_species):

            increment_vector = cls.diffuse_step_single_species(time_step, species_index=species_index)
            #print("Increment vector is: ", increment_vector)

            # For each bin, update the concentrations from the buffered increments
            cls.delta_diffusion[species_index] = increment_vector      # Vector operation to a row of the matrix delta_diffusion



    @classmethod
    def diffuse_step_single_species(cls, time_step: float, species_index=0) -> np.array:
        """
        Diffuse the specified single species, for the specified time step, across all bins,
        and return an array of the changes in concentration ("Delta concentration")
        for the given species across all bins.

        IMPORTANT: the actual system concentrations are NOT changed.

        We're assuming an isolated environment, with nothing diffusing thru the "walls"

        EXPLANATION:  https://life123.science/diffusion

        TODO: allow to specify the Delta_x (for now, taken to always be 1)

        :param time_step:       Time step over which to carry out the diffusion.
                                If too large, an Exception will be raised.
        :param species_index:   ID (in the form of an integer index) of the chemical species under consideration
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
        Update the system concentrations as a result of all the reactions in all bins.
        CAUTION : NO diffusion is taken into account.

        For each bin, process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        TODO: in case of any Exception, the state of the system is still valid, as of the time before this call

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param snapshots:       OPTIONAL dict with the keys: "frequency", "sample_bin", "sample_species"
                                    If provided, take a system snapshot after running a multiple of "frequency" runs.
                                    EXAMPLE: snapshots={"frequency": 2, "sample_bin": 0}
        :return:                None
        """
        time_step, n_steps = cls.all_reactions.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        # TODO: validation; implement sample_species
        if snapshots is None:
            frequency = None
        else:
            frequency = snapshots.get("frequency", 1)

        for i in range(n_steps):
            cls.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            cls.system += cls.delta_reactions   # Matrix operation to update all the concentrations
            cls.system_time += time_step
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

        # For each bin
        for bin_n in range(cls.n_bins):     # Bin number, ranging from 0 to max_bin_number, inclusive
            if cls.verbose:
                print(f"reaction_step(): processing the all the reactions in bin number {bin_n}")

            # Obtain the Delta-concentration for each species, for this bin
            conc_dict = {species_index: cls.system[species_index , bin_n]
                         for species_index in range(cls.n_species)}
            #print(f"\nconc_dict in bin {bin_n}: ", conc_dict)


            # Obtain the Delta-conc for each species, for bin number bin_n
            # OLD APPROACH
            #increment_vector = cls.single_bin_reaction_step(bin_n, delta_time, number_reactions)

            # NEW APPROACH
            increment_vector = cls.all_reactions.single_compartment_reaction_step(conc_dict=conc_dict, delta_time=delta_time)


            # Replace the "bin_n" column of the cls.delta_reactions matrix with the contents of the vector increment_vector
            #cls.delta_reactions[:, bin_n] = increment_vector.transpose()
            cls.delta_reactions[:, bin_n] = np.array([increment_vector])

            #print(cls.delta_reactions)




    #########################################################################
    #                                                                       #
    #                         REACTION-DIFFUSION                            #
    #                                                                       #
    #########################################################################

    @classmethod
    def react_diffuse(cls, total_duration=None, time_step=None, n_steps=None) -> None:
        """
        It expects 2 of the arguments:  total_duration, time_step, n_steps
        Perform a series of reaction and diffusion time steps.

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :return:                None
        """
        time_step, n_steps = cls.all_reactions.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        for i in range(n_steps):
            # TODO: split off the reaction step and the diffusion step to 2 different computing cores
            cls.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            cls.diffuse_step(time_step)
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
