import numpy as np
import pandas as pd
from typing import Union, List, Tuple
from src.modules.movies.movies import MovieTabular
from src.modules.reactions.reaction_dynamics import ReactionDynamics


class BioSim2D:
    """
    2D simulations of diffusion and reactions
    """

    def __init__(self, n_bins=None, chem_data=None, reactions=None):
        """

        :param n_bins:      A pair with the bin size in the x- and y- coordinates
        :param chem_data:
        :param reactions:
        """
        self.debug = False

        self.n_bins_x = 0       # Number of x-direction spatial compartments (bins) used in the simulation
        self.n_bins_y = 0       # Number of y-direction spatial compartments (bins) used in the simulation

        self.n_species = 1      # The number of (non-water) chemical species    TODO: phase out?

        self.chem_data = None   # Object of type "ReactionData", with info on the individual chemicals and their reactions

        self.system = None      # Concentration data in the System we're simulating, for all the chemicals
                                #   NumPy array of dimension (n_species x n_bins_x x n_bins_y)
                                #   Each plane represents a species

        # The following buffers are (n_species x n_bins_x x n_bins_y)
        self.delta_diffusion = None  # Buffer for the concentration changes from diffusion step
        self.delta_reactions = None  # Buffer for the concentration changes from reactions step

        self.sealed = True           # If True, no exchange with the outside; if False, immersed in a "bath"

        # Only applicable if "sealed" is False:
        self.bath_concentrations = None      # A NumPy array for each species
        self.container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container


        self.reaction_dynamics = None        # Object of class "ReactionDynamics"

        self.system_time = None              # Global time of the system, from initialization on

        if (n_bins is not None) or (chem_data is not None) or (reactions is not None):
            self.initialize_system(n_bins=n_bins, chem_data=chem_data, reactions=reactions)




    #########################################################################
    #                                                                       #
    #                           SYSTEM-WIDE                                 #
    #                                                                       #
    #########################################################################
    
    def initialize_system(self, n_bins: (int, int), chem_data=None, reactions=None) -> None:
        """
        Initialize all concentrations to zero.

        :param n_bins:     The number of compartments (bins) to use in the simulation,
                                in the x- and y- dimensions, as a pair of integers
        :param chem_data:   (OPTIONAL) Object of class "Chemicals";
                                if not specified, it will get extracted from the "Reactions" class
        :param reactions:   (OPTIONAL) Object of class "Reactions";
                                if not specified, it'll get instantiated here

        :return:            None
        """
        assert type(n_bins) == tuple, "BioSim2D: the argument `n_bins` must be a pair of integers"
        (n_cells_x, n_cells_y) = n_bins
        assert n_cells_x >= 1, "The number of bins must be at least 1 in each dimension"
        assert n_cells_y >= 1, "The number of bins must be at least 1 in each dimension"

        assert chem_data is not None or reactions is not None, \
            "BioSim2D: at least one of the arguments `chem_data` and `reactions` must be set"

        if chem_data:
            self.chem_data = chem_data
        else:
            self.chem_data = reactions.chem_data

        if reactions:
            self.reaction_dynamics = reactions
        else:
            self.reaction_dynamics = ReactionDynamics(reaction_data=chem_data)

        self.n_bins_x = n_cells_x
        self.n_bins_y = n_cells_y

        self.n_species = chem_data.n_species

        # Initialize all concentrations to zero
        self.system = np.zeros((self.n_species, n_cells_x, n_cells_y), dtype=float)

        self.system_time = 0             # "Start the clock"



    def system_size(self) -> (int, int):
        """
        Return a pair of integers with the system size in the x- and y-dimensions
        Note: the bin numbers will range between 0 and system_size - 1

        :return:    The pair (x-dimension, y-dimension)
        """
        return (self.n_bins_x, self.n_bins_y)



    def system_snapshot(self, species_index=0) -> pd.DataFrame:
        """
        Return a snapshot of all the concentrations of the given species, across all bins,
        as a Pandas dataframe

        :return:    A Pandas dataframe with the concentration data for the single specified chemical;
                    rows and columns correspond to the system's rows and columns
        """
        # TODO: validation for species_index

        matrix = self.system[species_index]

        df = pd.DataFrame(matrix)

        return df




    #########################################################################
    #                                                                       #
    #                     SET/MODIFY CONCENTRATIONS                         #
    #                                                                       #
    #########################################################################

    def set_bin_conc(self, bin_x: int, bin_y: int, species_index: int, conc: float) -> None:
        """
        Assign the requested concentration value to the cell with the given index, for the specified species

        :param bin_x:           The zero-based bin number of the desired cell
        :param bin_y:           The zero-based bin number of the desired cell
        :param species_index:   Zero-based index to identify a specific chemical species
        :param conc:            The desired concentration value to assign to the specified location
        :return:                None
        """
        #assert bin < self.n_bins, f"The requested cell index ({bin}) must be in the range [0 - {self.n_bins - 1}]"
        assert species_index < self.n_species, \
            f"The requested species index ({bin}) must be in the range [0 - {self.n_species - 1}]"

        assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"

        self.system[species_index, bin_x, bin_y] = conc



    
    def set_bin_conc_all_species(self, bin_x: int, bin_y: int, conc_list: [float]) -> None:
        """
        Assign the requested concentration values to the cell with the given index,
        for all the species in their index order

        :param bin_x:       The zero-based bin number of the desired cell
        :param bin_y:       The zero-based bin number of the desired cell
        :param conc_list:   A list with the desired concentration values to assign to the specified location
        :return:            None
        """
        #assert bin < self.n_bins, f"The requested cell index ({bin}) must be in the range [0 - {self.n_bins - 1}]"

        for i, conc in enumerate(conc_list):
            assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"
            self.system[i, bin_x, bin_y] = conc



    def set_species_conc(self, conc_data: Union[list, tuple, np.ndarray], species_index=None, species_name=None) -> None:
        """
        For the single specified species, assign the requested list of concentration values to all the bins,
        in row order first (top to bottom) and then in column order.

        EXAMPLE:  set_species_conc([[1, 2, 3], [4, 5, 6]], species_index=0)
                  will set the system state, for the specified chemical, to:
                            [1, 2, 3]
                            [4, 5, 6]

        :param conc_data:       A list, tuple or Numpy array with the desired concentration values
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
            raise Exception("BioSim2D.set_species_conc(): must provide a `species_name` or `species_index`")
        else:
            self.chem_data.assert_valid_species_index(species_index)

        assert (type(conc_data) == list) or (type(conc_data) == tuple) or (type(conc_data) == np.ndarray), \
                    f"BioSim2D.set_species_conc(): the argument `conc_list` must be a list, tuple or Numpy array; " \
                    f"the passed value was of type {type(conc_data)})"

        if type(conc_data) == np.ndarray:
            assert conc_data.shape == (self.n_bins_x, self.n_bins_y), \
                f"set_species_conc(): the numpy array `conc_list` must have dimensions {(self.n_bins_x, self.n_bins_y)}; instead, it was {conc_data.shape}"
        else:
            # Verify that the list or tuple corresponds to a matrix of the correct size
            assert all((type(ele) == list or type(ele) == tuple) for ele in conc_data), \
                f"BioSim2D.set_species_conc(): the argument `conc_list` must represent a matrix: all its elements must be lists or tuples"
            assert len(conc_data) == self.n_bins_x, \
                f"BioSim2D.set_species_conc(): the argument `conc_list` must represent a matrix with {self.n_bins_x} rows (found {len(conc_data)})"
            assert all(len(ele) == self.n_bins_y for ele in conc_data), \
                f"BioSim2D.set_species_conc(): the argument `conc_list` must represent a matrix with {self.n_bins_y} columns"

            conc_data = np.array(conc_data, dtype=float)

        # Verify that none of the concentrations are negative
        assert conc_data.min() >= 0, \
            f"BioSim2D.set_species_conc(): concentrations cannot be negative (values like {conc_data.min()} aren't permissible)"

        # Update the system state
        self.system[species_index] = conc_data



    def inject_conc_to_bin(self, bin_address: (int, int), species_index: int, delta_conc: float, zero_clip = False) -> None:
        """
        Add the requested concentration to the cell with the given address, for the specified chem species

        :param bin_address:     A pair with the zero-based bin numbers of the desired cell, in the x- and y-coordinates
        :param species_index:   Zero-based index to identify a specific chemical species
        :param delta_conc:      The concentration to add to the specified location
        :param zero_clip:       If True, any requested increment causing a concentration dip below zero, will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        #self.assert_valid_bin(bin_address)     #TODO: define

        bin_x, bin_y = bin_address  # Unpack the bin address

        if (self.system[species_index, bin_x, bin_y] + delta_conc) < 0. :
            # Take special action if the requested change would make the bin concentration negative
            if zero_clip:
                self.system[species_index, bin_x, bin_y] = 0
                return
            else:
                raise Exception("inject_conc_to_bin(): The requested concentration change would result in a negative final value")

        # Normal scenario, not leading to negative values for the final concentration
        self.system[species_index, bin_x, bin_y] += delta_conc




    #########################################################################
    #                                                                       #
    #                              TO VIEW                                  #
    #                                                                       #
    #########################################################################
    
    def describe_state(self, concise=False) -> None:
        """
        For each chemical species, show its name (or index, if name is missing),
        followed by the matrix of concentrations values for that chemical

        :param concise:     Not yet used
        :return:            None
        """
        #np.set_printoptions(linewidth=125)
        print(f"SYSTEM STATE at Time t = {self.system_time}:")
        for species_index in range(self.n_species):
            chem_name = self.chem_data.get_name(species_index)
            if chem_name is None:
                print(f"Species {species_index}:")      # Use the index, if the name isn't available
            else:
                print(f"Species `{chem_name}`:")

            #print(self.system[species_index])
            print(self.system_snapshot(species_index))



    def lookup_species(self, species_index=None, species_name=None, copy=False) -> np.array:
        """
        Return the NumPy array of concentration values across the all bins
        (from top to bottom, and then left to right),
        for the single specified chemical species.
        NOTE: what is being returned NOT a copy, unless specifically requested

        :param species_index:   The index order of the chemical species of interest
        :param species_name:    (OPTIONAL) If provided, it over-rides the value for species_index
        :param copy:            If True, an independent numpy array will be returned: a *copy* rather than a view
        :return:                A NumPy 2-D array of concentration values across the bins
                                    (from top to bottom, and then left to right);
                                    the size of the array is (n_bins_y x n_bins_x)
        """
        if species_name is not None:
            species_index = self.chem_data.get_index(species_name)
        else:
            self.chem_data.assert_valid_species_index(species_index)

        species_conc = self.system[species_index]

        if copy:
            return species_conc.copy()
        else:
            return species_conc



    def bin_snapshot_array(self, bin_address: (int, int)) -> np.array:
        """
        Extract the concentrations of all the chemical species at the specified bin,
        as a Numpy array in the index order of the species
        EXAMPLE: np.array([10., 50.)]

        :param bin_address:     A pair with the zero-based bin numbers of the desired cell, in the x- and y-coordinates
        :return:                A Numpy array  of concentration values, in the index order of the species
        """
        bin_x, bin_y = bin_address      # Unpack the bin address


        return self.system[:, bin_x, bin_y]




    #########################################################################
    #                                                                       #
    #                               DIFFUSION                               #
    #                                                                       #
    #########################################################################


    def diffuse(self, total_duration=None, time_step=None, n_steps=None, h=1, algorithm="5_point") -> dict:
        """
        Uniform-step diffusion, with 2 out of 3 criteria specified:
            1) until reaching, or just exceeding, the desired time duration
            2) using the given time step
            3) carrying out the specified number of steps

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param h:               Distance between consecutive bins in both the x- and y-directions
                                    (For now, they must be equal)
        :param algorithm:       (OPTIONAL) String with a name specifying the method to use to solve the diffusion equation.
                                    Currently available options: "5_point"
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

            self.diffuse_step(time_step, h=h, algorithm=algorithm)
            self.system += self.delta_diffusion     # Array operation to update all the concentrations
            self.system_time += time_step

        if self.debug:
            print(f"\nSystem after Delta time {total_duration}, at end of {n_steps} steps of size {time_step}:")
            self.describe_state(concise=True)
            print()

        status = {"steps": n_steps}
        return status



    def diffuse_step(self, time_step, h=1, algorithm="5_point") -> None:
        """
        Diffuse all the species for the given time step, across all bins;
        clear the delta_diffusion array, and then re-compute it from all the species.

        IMPORTANT: the actual system concentrations are NOT changed.

        :param time_step:   Time step over which to carry out the diffusion
                            If too large - as determined by the method is_excessive() - an Exception will be raised
        :param h:           Distance between consecutive bins in both the x- and y-directions
                                (For now, they must be equal)
        :param algorithm:   String with a name specifying the method to use to solve the diffusion equation.
                                Currently available options: "5_point"
        :return:            None (the array in the class variable "delta_diffusion" gets set)
        """
        # TODO: parallelize the independent computations

        # 3-D array of incremental changes at every bin, for each chemical species
        self.delta_diffusion = np.zeros((self.n_species, self.n_bins_x, self.n_bins_y), dtype=float)

        for species_index in range(self.n_species):
            increment_matrix = self.diffuse_step_single_species(time_step=time_step, h=h, species_index=species_index, algorithm=algorithm)

            #print("Increment vector is: ", increment_vector)

            # For each bin, update the concentrations from the buffered increments
            self.delta_diffusion[species_index] = increment_matrix      # Matrix operation to a plane of the 3-D array delta_diffusion



    def diffuse_step_single_species(self, time_step: float, h: float, species_index=0, algorithm="5_point") -> np.array:
        """
        Diffuse the specified single chemical species, for the given time step, across all bins,
        and return a 2-D array of the changes in concentration ("Delta concentration")
        for the given species across all bins.

        IMPORTANT: the actual system concentrations are NOT changed.

        We're assuming an isolated environment, with nothing diffusing thru the "walls"

        EXPLANATION of the methodology:  https://life123.science/diffusion

        TODO: also test on tiny systems smaller than 3x3

        :param time_step:       Delta time over which to carry out this single diffusion step;
                                    TODO: add - if too large, an Exception will be raised.
        :param species_index:   ID (in the form of an integer index) of the chemical species under consideration
        :param h:               Distance between consecutive bins, ASSUMED the same in both directions
        :param algorithm:       String with a name specifying the method to use to solve the diffusion equation.
                                    Currently available options: "5_point"

        :return:                A 2-D Numpy array with the CHANGE in concentration for the given species across all bins
        """
        assert self.system is not None, "BioSim2D.diffuse_step_single_species(): Must first initialize the system"
        assert not self.chem_data.missing_diffusion_rate(), "BioSim2D.diffuse_step_single_species(): Must first set the diffusion rates"
        assert self.sealed == True, "BioSim2D.diffuse_step_single_species(): For now, there's no provision for exchange with the outside"

        increment_matrix = np.zeros((self.n_bins_x, self.n_bins_y), dtype=float)   # One element per bin

        if self.n_bins_x and self.n_bins_y == 1:
            return increment_matrix                                 # There's nothing to do in the case of just 1 bin!

        diff = self.chem_data.get_diffusion_rate(species_index)     # The diffusion rate of the specified single species

        #assert not self.is_excessive(time_step, diff, delta_x), \  # TODO: implement
            #f"Excessive large time_step ({time_step}). Should be < {self.max_time_step(diff, delta_x)}"

        # We're calling the following quantity "Effective Diffusion" (NOT a standard term)
        effective_diff = diff * time_step / (h ** 2)

        if algorithm == "5_point":
            self.convolution_5_point_stencil(increment_matrix, species_index, effective_diff)
        else:
            raise Exception(f"diffuse_step_single_species(): Unknown algorithm: `{algorithm}`")

        return increment_matrix



    def convolution_5_point_stencil(self, increment_matrix, species_index, effective_diff) -> None:
        """
        Carry out a 2-D convolution operation on increment_matrix,
        with a tile of size 3 that implements a 5-point stencil

        TODO: maybe pass self.system[species_index] as argument
        TODO: move the x effective_diff to the calling function

        :param increment_matrix:
        :param species_index:
        :param effective_diff:
        :return:                None (increment_matrix gets modified)
        """
        max_bin_x = self.n_bins_x - 1    # Bin numbers range from 0 to max_bin_x, inclusive (in x-direction)
        max_bin_y = self.n_bins_y - 1    # Bin numbers range from 0 to max_bin_y, inclusive (in y-direction)

        for i in range(self.n_bins_x):          # From 0 to max_bin_x, inclusive
            for j in range(self.n_bins_y):      # From 0 to max_bin_y, inclusive
                C_ij = self.system[species_index, i, j]     # Concentration in the center of the convolution tile

                # The boundary conditions state that the flux is zero across boundaries
                if i == 0:
                    C_left = self.system[species_index, 0, j]
                else:
                    C_left = self.system[species_index, i-1, j]

                if i == max_bin_x:
                    C_right = self.system[species_index, max_bin_x, j]
                else:
                    C_right = self.system[species_index, i+1, j]

                if j == 0:
                    C_above = self.system[species_index, i, 0]
                else:
                    C_above = self.system[species_index, i, j-1]

                if j == max_bin_y:
                    C_below = self.system[species_index, i, max_bin_y]
                else:
                    C_below = self.system[species_index, i, j+1]


                # Convolution with a 5-point stencil
                increment_matrix[i, j] = effective_diff * \
                                         (C_above + C_below + C_left + C_right - 4 * C_ij)





    #########################################################################
    #                                                                       #
    #                               REACTIONS                               #
    #                                                                       #
    #########################################################################

    
    def react(self, total_duration=None, time_step=None, n_steps=None, snapshots=None) -> None:
        """
        Update the system concentrations as a result of all the reactions in all bins.
        CAUTION : NO diffusion is taken into account.

        The duration and granularity of the reactions is specified with 2 out of the 3 parameters:
            total_duration, time_step, n_steps

        For each bin, process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        TODO: in case of any Exception, the state of the system is still valid, as of the time before this call

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :param snapshots:       NOT YET USED
                                OPTIONAL dict with the keys: "frequency", "sample_bin", "sample_species"
                                    If provided, take a system snapshot after running a multiple of "frequency" runs
        :return:                None
        """

        time_step, n_steps = self.reaction_dynamics.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)
        #if snapshots is None:
            #frequency = None
        #else:
            #frequency = snapshots.get("frequency", 1)

        for i in range(n_steps):
            self.reaction_step(time_step)         # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            self.system += self.delta_reactions   # Matrix operation to update all the concentrations
            self.system_time += time_step
            #if (frequency is not None) and ((i+1)%frequency == 0):
                #self.save_snapshot(self.bin_snapshot(bin_address = snapshots["sample_bin"]))



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

        :param delta_time:  The time duration of the reaction step - assumed to be small enough that the
                            concentration won't vary significantly during this span
        :return:            None (note: the class variable "delta_reactions" gets updated)
        """
        assert self.reaction_dynamics is not None, \
            "BioSim2D.reaction_step(): must first set the Reactions object"

        self.delta_reactions = np.zeros((self.n_species, self.n_bins_x, self.n_bins_y), dtype=float)

        # For each bin
        for bin_n_x in range(self.n_bins_x):
            for bin_n_y in range(self.n_bins_y):
                if self.debug:
                    print(f"BioSim2D.reaction_step(): processing the all the reactions "
                          f"in bin number ({bin_n_x}, {bin_n_y})")

                # Obtain the Delta-concentration for each species, for this bin
                conc_array = self.bin_snapshot_array(bin_address=(bin_n_x, bin_n_y))
                if self.debug:
                    print(f"\nconc_dict in bin ({bin_n_x}, {bin_n_y}): ", conc_array)


                # Obtain the Delta-conc for each species, for the current bin
                increment_vector, _, _ = self.reaction_dynamics.reaction_step_orchestrator(delta_time=delta_time, conc_array=conc_array)
                                                                                                    #delta_time=delta_time)

                # Replace the appropriate column of the self.delta_reactions matrix
                # with the contents of the vector increment_vector
                self.delta_reactions[:, bin_n_x, bin_n_y] = np.array([increment_vector])

                if self.debug:
                    print(self.delta_reactions)
