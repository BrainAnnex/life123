import numpy as np
import pandas as pd
from typing import Union
import plotly.express as px
import plotly.graph_objects as pgo
from life123.uniform_compartment import UniformCompartment
from life123.history import HistoryBinConcentration
from life123.visualization.plotly_helper import PlotlyHelper
from life123.visualization.colors import Colors



class BioSim2D:
    """
    2D simulations of diffusion and reactions
    """

    def __init__(self, x_bins :int, y_bins :int, chem_data=None, reaction_handler=None):
        """
        :param x_bins:          The bin size in the x-coordinates.  Notice that this is the number of COLUMNS in the data matrix
        :param y_bins:          The bin size in the y-coordinates.  Notice that this is the number of ROWS in the data matrix

        [At least one of the 2 following arguments must be provided]
        :param chem_data:       [OPTIONAL] Object of class "ChemData";
                                    if not specified, it will get extracted
                                    from the "UniformCompartment" class (if passed to the next argument)
        :param reaction_handler:[OPTIONAL] Object of class "UniformCompartment";
                                    if not specified, it'll get instantiated here
        """
        self.debug = False

        self.n_bins_x = 0       # Number of x-direction spatial compartments (bins) used in the simulation
        self.n_bins_y = 0       # Number of y-direction spatial compartments (bins) used in the simulation

        self.n_species = 1      # The number of (non-water) chemical species    TODO: phase out?

        self.chem_data = None   # Object of type "ChemData", with info on the individual chemicals and their reactions

        self.reactions = None   # Object of type "Reactions", with info on all the reactions

        self.reaction_dynamics = None   # Object of class "UniformCompartment"

        self.system = None      # Concentration data in the System we're simulating, for all the chemicals:
                                #   NumPy array of dimension (n_species) x (n_bins_x) x (n_bins_y)
                                #   Each plane (sub-array along the 0-th axis) represents a chemical species;
                                #   e.g, self.system[0] is the matrix for the 0-th chemical.
                                #   IMPORTANT: the x-dimension and y-dimension values are stored as, respectively,
                                #              rows and columns in that matrix.  That means that the matrix is the TRANSPOSE
                                #              of the xy grid layout.  Visualization methods will need to make the transformation!
                                #              See: https://github.com/BrainAnnex/life123/discussions/75

        # The following buffers are of size (n_species x n_bins_x x n_bins_y)
        self.delta_diffusion = None  # Buffer for the concentration changes from diffusion step
        self.delta_reactions = None  # Buffer for the concentration changes from reactions step

        self.sealed = True           # If True, no exchange with the outside; if False, immersed in a "bath"

        # Only applicable if "sealed" is False:
        self.bath_concentrations = None      # A NumPy array for each species
        self.container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container

        self.system_time = None              # Global time of the system, from initialization on

        self._initialize_system(x_bins=x_bins, y_bins=y_bins, chem_data=chem_data, reaction_handler=reaction_handler)

        self.conc_history = HistoryBinConcentration(active=False)




    #########################################################################
    #                                                                       #
    #                           SYSTEM-WIDE                                 #
    #                                                                       #
    #########################################################################

    def _initialize_system(self, x_bins :int, y_bins: int, chem_data=None, reaction_handler=None) -> None:
        """
        Initialize all concentrations to zero.

        :param x_bins:      The bin size in the x-coordinates.  Notice that this is the number of COLUMNS in the data matrix
        :param y_bins:      The bin size in the y-coordinates.  Notice that this is the number of ROWS in the data matrix
        :param chem_data:   (OPTIONAL) Object of class "ChemData";
                                if not specified, it will get extracted from the "Reactions" class
        :param reaction_handler:   (OPTIONAL) Object of class "Reactions";
                                if not specified, it'll get instantiated here

        :return:            None
        """
        assert type(x_bins) == int, "BioSim2D() instantiation: the argument `x_bins` must be an integer"
        assert type(y_bins) == int, "BioSim2D() instantiation: the argument `y_bins` must be an integer"

        assert x_bins >= 1 and y_bins >= 1, \
            "BioSim2D() instantiation: The number of bins must be at least 1 in each dimension"

        assert chem_data is not None or reaction_handler is not None, \
            "BioSim2D() instantiation: at least one of the arguments `chem_data` or `reaction_handler` must be set"

        if chem_data:
            self.chem_data = chem_data
        else:
            self.chem_data = reaction_handler.chem_data

        if reaction_handler:
            self.reaction_dynamics = reaction_handler
        else:
            self.reaction_dynamics = UniformCompartment(chem_data=self.chem_data)

        self.reactions = self.reaction_dynamics.get_reactions()     # TODO: Maybe use self.get_reactions()

        self.n_bins_x = x_bins
        self.n_bins_y = y_bins

        self.n_species = self.chem_data.number_of_chemicals()

        assert self.n_species >= 1, \
            "BioSim2D() instantiation: At least 1 chemical species must be declared prior to instantiating class"

        # Initialize all bin concentrations to zero
        self.system = np.zeros((self.n_species, x_bins, y_bins), dtype=float)

        self.system_time = 0             # "Start the clock"




    #####################################################################################################

    '''                                    ~   VIEW/UPDATE SYSTEM   ~                                           '''

    def ________VIEW_SYSTEM________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

    def system_size(self) -> (int, int):
        """
        Return a pair of integers with the system sizes in the x- and y-dimensions
        Note: the bin numbers start with 0

        :return:    The pair (x-dimension, y-dimension)
        """
        return (self.n_bins_x, self.n_bins_y)



    def get_system_time(self) -> float:
        """
        Return the current system time
        :return:
        """
        return self.system_time



    def get_chem_data(self):
        """

        :return:    An Object of type "ChemData"
        """
        return self.chem_data



    def get_reactions(self):
        """
        Return all the associated reactions

        :return:    Object ot type "Reactions" (with data about all the reactions)
        """
        return self.reactions



    def system_snapshot_arr_xy(self, chem_label=None, chem_index=None) -> np.ndarray:
        """
        Return a snapshot of all the concentrations of the given chemical species,
        across ALL BINS, as a Numpy array in XY coordinates.
        IMPORTANT: the rows of the matrix store the x-coordinates, and the columns store the y-coordinates;
                   so, this matrix is the TRANSPOSE of the X-Y cartesian representation
        If a Pandas dataframe is desired, for a single chemical, use system_snapshot()

        :param chem_label:  String with the label to identify the chemical of interest
        :param chem_index:  Integer to identify the chemical of interest.  Cannot specify both chem_label and chem_index
        :return:            A 2-D Numpy array of concentration values in XY coordinates
        """
        assert (chem_label is not None) or (chem_index is not None), \
            "system_snapshot_xy(): at least one of the args `chem_label` or `chem_index` must be provided"

        if chem_label is not None:
            assert chem_index is None, "system_snapshot_xy(): cannot pass both arguments `chem_label` and `chem_index`"
            chem_index = self.chem_data.get_index(chem_label)
        else:
            assert chem_index is not None, "system_snapshot_xy(): must pass one of the arguments `chem_label` or `chem_index`"
            self.chem_data.assert_valid_species_index(chem_index)

        matrix = self.system[chem_index].T    # A 2-D Numpy array with the chemical data in XY dimensions (notice the transpose)

        return matrix



    def system_snapshot(self, chem_label=None, chem_index=None, cartesian=True) -> pd.DataFrame:
        """
        Return a snapshot of all the concentrations of the given chemical species,
        across ALL BINS, in a grid XY coordinate system, as a Pandas dataframe.
        The columns of the dataframe are the x-coordinates, increasing to the right.
        The rows of the the dataframe are the y-coordinates, by default increasing in the up direction (earlier rows)

        :param chem_label:  String with the label to identify the chemical of interest
        :param chem_index:  Integer to identify the chemical of interest.  Cannot specify both chem_label and chem_index
        :param cartesian:   If True (default) a Cartesian grid coordinate is used, with y-bin numbers increasing up
        :return:            A Pandas dataframe with the concentration data for the single specified chemical;
                            rows and columns correspond to the system's rows and columns
        """
        # TODO: offer the option to pass multiple labels
        # TODO: fix inconsistency vs. 1D case (single chem vs. all)

        matrix = self.system_snapshot_arr_xy(chem_label=chem_label, chem_index=chem_index)  # In XY coordinates

        df = pd.DataFrame(matrix)

        if cartesian:
            return df[::-1]     # For other ways to do this: https://www.geeksforgeeks.org/how-to-reverse-row-in-pandas-dataframe/

        return df



    def bin_concentration(self, bin_address: (int, int), species_index=None, species_label=None) -> float:
        """
        Return the concentration at the requested bin of the specified chemical species

        :param bin_address:     A pair of integers identifying the bin of interest
        :param species_index:   The index order of the chemical species of interest
        :param species_label:   [OPTIONAL] If provided, it over-rides the value for species_index
        :return:                A concentration value at the indicated bin, for the requested species
        """
        if species_label is not None:
            species_index = self.chem_data.get_index(species_label)

        self.chem_data.assert_valid_species_index(species_index)

        xbin, ybin = bin_address

        return self.system[species_index, xbin, ybin]



    def describe_state(self, cartesian=True) -> None:
        """
        A minimalist view of all the chemical concentrations across the 2D system.
        For each chemical species, show its name (or index, if name is missing),
        followed by the matrix of its concentrations values

        :param cartesian:   If True (default) a Cartesian grid coordinate is used, with y-bin numbers increasing up
        :return:            None
        """
        # Set the maximum width of each row, when printing a Pandas dataframe, to fit more content
        pd.set_option('display.width', 150)

        print(f"SYSTEM STATE at Time t = {self.system_time:.8g}:")
        for species_index in range(self.n_species):
            chem_name = self.chem_data.get_label(species_index)
            if chem_name is None:
                print(f"Species {species_index}:")      # Use the index, if the name isn't available
            else:
                print(f"Species `{chem_name}`:")

            print(self.system_snapshot(chem_index=species_index, cartesian=cartesian))



    def lookup_species(self, species_index=None, species_name=None, copy=False) -> np.array:
        """
        Return the NumPy array of concentration values across ALL THE BINS
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



    def selected_concentrations(self, bins, chem_labels=None) -> dict:
        """
        Extract and return the concentration values of one or more (use None for all) chemicals,
        in one or more bins.
        The value is returned as a dictionary where the keys are bin addresses, and the values are dicts of
        concentration values for the various chemicals (identified by their labels)
            EXAMPLE:
                    {   (0,0): {"A": 1.3, "B": 3.9},
                        (3,5): {"A": 4.6, "B": 2.7}
                    }

            TODO: alternate returned structure - a Pandas dataframe such as
                    BIN ADDRESS      A       B
                        (0,0)       1.3     3.9
                        (3,5)       4.6     2.7

        :param bins:    Bin address (pair of integers), or list of bin addresses. Use None to indicate all
        :param chem_labels:  Chemical label, or list of labels. Use None to indicate all
        :return:        A dict indexed by bin address
        """
        if bins is None:
            bins = [(x,y) for x in range(self.n_bins_x) for y in range(self.n_bins_y)]  # All bins
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





    #####################################################################################################

    '''                                    ~   UPDATE SYSTEM   ~                                           '''

    def ________UPDATE_SYSTEM________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def assert_valid_bin(self, bin_address: (int,int)) -> None:
        """
        Raise an Exception if the given bin address isn't valid

        :param bin_address: A pair of integers identifying a bin in the 2D system
        :return:            None
        """
        assert type(bin_address) == tuple, \
            f"BioSim2D: the requested bin address `{bin_address}` is not a pair of integers; its type is {type(bin_address)}"

        assert len(bin_address) == 2, \
            f"BioSim2D: the requested bin address {bin_address} is not a pair of integers; it has length {len(bin_address)}"

        bin_x, bin_y =  bin_address

        assert type(bin_x) == int, \
            f"BioSim2D: the requested bin address {bin_address} isn't composed of integers; " \
            f"the types in the pair are ({type(bin_x)}, {type(bin_y)})"

        assert type(bin_y) == int, \
            f"BioSim2D: the requested bin address {bin_address} isn't composed of integers; " \
            f"the types in the pair are ({type(bin_x)}, {type(bin_y)})"

        assert (bin_x >= 0) and (bin_x < self.n_bins_x),\
            f"BioSim2D: the requested bin address {bin_address} is out of bounds for the system size; " \
            f"allowed range for the x-bin is [0-{self.n_bins_x-1}], inclusive"

        assert (bin_y >= 0) and (bin_y < self.n_bins_y),\
            f"BioSim2D: the requested bin address {bin_address} is out of bounds for the system size; " \
            f"allowed range for the y-bin is [0-{self.n_bins_y-1}], inclusive"



    def set_bin_conc(self, bin_address: (int, int), chem_label :str, conc: float) -> None:
        """
        Assign the requested concentration value to the bin with the given address,
        for the specified chemical species

        :param bin_address: A pair with the zero-based bin numbers of the desired cell, in the x- and y-coordinates
        :param chem_label:  String with the label to identify the chemical of interest
        :param conc:        The desired concentration value to assign to the specified location
        :return:            None
        """
        self.assert_valid_bin(bin_address)

        bin_x, bin_y = bin_address  # Unpack the bin address

        assert conc >= 0., \
            f"set_bin_conc(): The concentration must be a positive number or zero (the provided value was {conc})"

        species_index = self.chem_data.get_index(chem_label)
        self.system[species_index, bin_x, bin_y] = conc



    def inject_conc_to_bin(self, bin_address: (int, int), chem_index: int, delta_conc: float, zero_clip = False) -> None:
        """
        Add the requested concentration to the cell with the given address, for the specified chemical species

        :param bin_address: A pair with the zero-based bin numbers of the desired cell, in the x- and y-coordinates
        :param chem_index:  Zero-based index to identify a specific chemical species
        :param delta_conc:  The concentration to add to the specified location
        :param zero_clip:   If True, any requested increment causing a concentration dip below zero, will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        #TODO: also allow chem label, as done in 1D
        self.assert_valid_bin(bin_address)

        bin_x, bin_y = bin_address  # Unpack the bin address

        if (self.system[chem_index, bin_x, bin_y] + delta_conc) < 0. :
            # Take special action if the requested change would make the bin concentration negative
            if zero_clip:
                self.system[chem_index, bin_x, bin_y] = 0
                return
            else:
                raise Exception("inject_conc_to_bin(): The requested concentration change would result in a negative final value")

        # Normal scenario, not leading to negative values for the final concentration
        self.system[chem_index, bin_x, bin_y] += delta_conc



    def set_bin_conc_all_species(self, bin_address: (int, int), conc_list: [float]) -> None:
        """
        Assign the requested concentration values to the cell with the given index,
        for all the chemical species in their index order

        :param bin_address: A pair with the zero-based bin numbers of the desired cell, in the x- and y-coordinates
        :param conc_list:   A list with the desired concentration values to assign to the specified location
        :return:            None
        """
        self.assert_valid_bin(bin_address)

        bin_x, bin_y = bin_address  # Unpack the bin address

        for i, conc in enumerate(conc_list):
            assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"
            self.system[i, bin_x, bin_y] = conc



    def set_species_conc(self, conc_data: Union[list, tuple, np.ndarray], species_index=None, species_name=None) -> None:
        """
        For the single specified chemical species, assign the requested list of concentration values to all the bins,
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
        arr = self.system_snapshot_arr_xy(chem_label=chem_label, chem_index=chem_index)
        total = np.sum(arr)
        return np.allclose(expected, total)





    #####################################################################################################

    '''                                      ~   HISTORY   ~                                          '''

    def ________HISTORY________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def enable_history(self, bins=None, frequency=1, chem_labels=None, take_snapshot=False) -> None:
        """
        Request history capture, with the specified parameters.
        If history was already enabled, this function can be used to alter its capture parameters.

        :param bins:            Bin address (pair of integers), or list of bin addresses. Use None to indicate all
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



    def get_bin_history(self, bin_address :(int,int)) -> pd.DataFrame:
        """
        Get the concentration history at the given bin(s) of all the chemicals
        whose history was requested by a call of enable_history()

        :param bin_address: A single bin address (a pair of integers)
        :return:            A Pandas data frame
        """
        return self.conc_history.bin_history(bin_address = bin_address)



    def capture_snapshot(self, step_count=None) -> None:
        """
        Add to the system history (if enabled) a snapshot of (part of) the current data

        :param step_count:
        :return:            None
        """
        if not self.conc_history.to_capture(step_count):
            return

        data_snapshot = self.selected_concentrations(bins=self.conc_history.restrict_bins,
                                                     chem_labels=self.conc_history.restrict_chemicals)
        '''
           EXAMPLE of data_snapshot:
                { (0,0): {"A": 1.3, "B": 3.9},
                  (3,5): {"A": 4.6, "B": 2.7}
                }        
        '''
        self.conc_history.save_snapshot(step_count=step_count, system_time=self.system_time,
                                        data_snapshot=data_snapshot)





    #####################################################################################################

    '''                                    ~   SIMULATIONS   ~                                           '''

    def ________SIMULATIONS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################

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
            self.capture_snapshot(step_count=i)     # Save historical values (if enabled)

        if self.debug:
            print(f"\nSystem after Delta time {total_duration}, at end of {n_steps} steps of size {time_step}:")
            self.describe_state()
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
                                (for now, they must be equal)
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
        assert self.system is not None, \
                "BioSim2D.diffuse_step_single_species(): Must first initialize the system"
        assert not self.chem_data.missing_diffusion_rate(), \
                "BioSim2D.diffuse_step_single_species(): Must first set the diffusion rates"
        assert self.sealed == True, \
                "BioSim2D.diffuse_step_single_species(): For now, there's no provision for exchange with the outside"

        increment_matrix = np.zeros((self.n_bins_x, self.n_bins_y), dtype=float)   # One element per bin

        if self.n_bins_x and self.n_bins_y == 1:
            return increment_matrix                                 # There's nothing to do in the case of just 1 bin!

        diff = self.chem_data.get_diffusion_rate(species_index=species_index)     # The diffusion rate of the specified single species

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

        :param increment_matrix:A Numpy matrix of the correct size
        :param species_index:   Integer to identify the chemical of interest
        :param effective_diff:
        :return:                None (increment_matrix gets modified)
        """
        #TODO: maybe pass self.system[species_index] as argument
        #TODO: move the multiplication by effective_diff to the calling function

        max_bin_x = self.n_bins_x - 1    # Bin numbers range from 0 to max_bin_x, inclusive (in x-direction)
        max_bin_y = self.n_bins_y - 1    # Bin numbers range from 0 to max_bin_y, inclusive (in y-direction)

        for i in range(self.n_bins_x):          # From 0 to max_bin_x, inclusive
            for j in range(self.n_bins_y):      # From 0 to max_bin_y, inclusive
                C_ij = self.system[species_index, i, j]     # Concentration in the center of the convolution tile

                # The boundary conditions state that the flux is zero across boundaries;
                # we can attain that by giving identical concentrations to the hypothetical neighboring bins
                # across boundaries, both up/down and left/right (diagonals values are not used in this algorithm)
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



    def react(self, total_duration=None, time_step=None, n_steps=None) -> None:
        """
        Update the system concentrations as a result of all the reactions in all bins.
        CAUTION : NO diffusion is taken into account.

        The duration and granularity of the reactions is specified with 2 out of the 3 parameters:
            total_duration, time_step, n_steps

        For each bin, process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :return:                None
        """
        # TODO: in case of any Exception, the state of the system is still valid, as of the time before this call

        time_step, n_steps = self.reaction_dynamics.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        for i in range(n_steps):
            self.reaction_step(time_step)         # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            self.system += self.delta_reactions   # Matrix operation to update all the concentrations
            self.system_time += time_step

            self.capture_snapshot(step_count=i)     # Save historical values (if enabled)



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
                increment_vector, _, _ = self.reaction_dynamics.reaction_step_common(delta_time=delta_time, conc_array=conc_array,
                                                                                     variable_steps=False)

                # Replace the appropriate column of the self.delta_reactions matrix
                # with the contents of the vector increment_vector
                self.delta_reactions[:, bin_n_x, bin_n_y] = np.array([increment_vector])

                if self.debug:
                    print(self.delta_reactions)



    def react_diffuse(self, total_duration=None, time_step=None, n_steps=None, h = 1) -> None:
        """
        It expects 2 out of the following 3 arguments:  total_duration, time_step, n_steps
        Perform a series of reaction and diffusion constant time steps.

        :param total_duration:  The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each constant time step
        :param n_steps:         The desired number of constant steps
        :param h:               Distance between consecutive bins in both the x- and y-directions
                                    (for now, they must be equal)
        :return:                None
        """
        time_step, n_steps = self.reaction_dynamics.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)

        for i in range(n_steps):
            # TODO: split off the diffusion step and the reaction steps to different computing cores
            self.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            self.diffuse_step(time_step, h=h)
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

    def heatmap_single_chem_greyscale(self, chem_label :str, title_prefix = "",
                                      width=None, height=550, cartesian=True) -> pgo.Figure:
        """
        Create and return a greyscale heatmap (a Plotly Figure object) of the 2D concentrations
        of the specified single chemical, using the current system data.

        :param chem_label:  Label to identify the chemical of interest
        :param title_prefix:[OPTIONAL] A string to prefix to the auto-generated Heatmap title
        :param width:       [OPTIONAL] If not specified, all the available space width is used
        :param height:      [OPTIONAL] Height of the heatmap graphics
        :param cartesian:   If True (default) a Cartesian grid coordinate is used, with y-bin numbers increasing up
        :return:            A Plotly "Figure" object
        """
        title = f"System state at time t={self.system_time:.5g} for `{chem_label}`"
        if title_prefix:
            title = f"{title_prefix}.  {title}"

        # Get the concentration data for the requested chemical
        df = self.system_snapshot(chem_label=chem_label, cartesian=False)   # Note: no need to reverse the rows in the dataframe;
                                                                            #       the heatmap will take care of that if requested

        origin = "lower" if cartesian else "upper"      # Position of the [0, 0] bin, in the upper left or lower left corner
        # Create the plotly heatmap
        fig = px.imshow(img=df,
                        title=title,
                        labels={"x": "x bin", "y": "y bin", "color": "Conc."},
                        text_auto=".2f", color_continuous_scale="gray_r",
                        width=width, height=height, origin=origin)

        # Insert a little spacing between adjacent cells
        fig.update_traces(xgap=2)       # Alt way to specify it: fig.data[0].xgap=2
        fig.update_traces(ygap=2)

        return fig



    def heatmap_single_chem(self, chem_label :str, title_prefix = "", width=None, height=550,
                            color=None, cartesian=True) -> pgo.Figure:
        """
        Create and return a heatmap (a Plotly Figure object) of the 2D concentration
        of the specified single chemical, using the current system data.
        A gradient of the specified color is used (grayscale if not color specified.)

        Note: if requesting a greyscale, this function is almost identical to heatmap_single_chem_greyscale(),
              but shows a little more info when the mouse hovers over the heatmap bins

        :param chem_label:  Label to identify the chemical of interest
        :param title_prefix:[OPTIONAL] A string to prefix to the auto-generated Heatmap title
        :param height:      [OPTIONAL] Height of the heatmap graphics
        :param width:       [OPTIONAL] If not specified, all the available space width is used
        :param color:       [OPTIONAL] The standard (CSS) name of the color to use for the gradient;
                                if not specified, a grayscale is used
        :param cartesian:   If True (default) a Cartesian grid coordinate is used, with y-bin numbers increasing up
        :return:            A Plotly "Figure" object
        """
        title = f"System state at time t={self.system_time:.5g} for `{chem_label}`"
        if title_prefix:
            title = f"{title_prefix}.  {title}"

        if color is None:
            color_scale = "gray_r"
        else:
            lighter_color = Colors.lighten_color(color, factor=.96)
            color_scale = [
                [0.0, lighter_color],   # Light tint
                [1.0, color],           # Full color
            ]

        # Get the concentration data for the requested chemical
        df = self.system_snapshot(chem_label=chem_label, cartesian=False)   # Note: no need to reverse the rows in the dataframe;
                                                                            #       the heatmap will take care of that if requested

        # Create the Heatmap plotly object
        hm = pgo.Heatmap(z=df,
                         colorscale=color_scale,
                         xgap=2, ygap=2,
                         hovertemplate='Conc.: %{z}<br>x bin: %{x}<br>y bin: %{y}<extra>' + chem_label + '</extra>',
                         texttemplate = '%{z:.2f}',
                         colorbar_title="Conc."
                         )

        # Create the Figure object
        fig = pgo.Figure(data=hm)

        # Update layout
        fig.update_layout(
            title=title,
            height=height,
            width=width,
            xaxis_title='x bin',
            yaxis_title='y bin'
        )

        if not cartesian:
            # Invert the y-axis for the heatmap (the default is the y-axis increasing down)
            fig.update_layout(
                yaxis_autorange="reversed"
            )

        return fig



    def system_heatmaps(self, chem_labels=None, title_prefix = "", height=None, colors=None, cartesian=True) -> pgo.Figure:
        """
        Prepare and return a Plotly Figure object containing a grid of heatmaps (up to a max of 12)

        :param chem_labels: [OPTIONAL] List of Labels to identify the chemicals of interest;
                                if not specified, it means ALL chemicals.
                                The max number of heatmaps that can be shown together is 12
        :param title_prefix:[OPTIONAL] A string to prefix to the auto-generated title for the grid of heatmaps
        :param height:      [OPTIONAL] Height of the overall grid of heatmaps
        :param colors:      [OPTIONAL] List of CSS color names for each of the heatmaps.
                                If provided, its length must match that of the data;
                                if None, then use the registered colors (if specified),
                                or the hardwired defaults as a last resort
        :param cartesian:   If True (default) a Cartesian grid coordinate is used, with y-bin numbers increasing up
        :return:            A Plotly "Figure" object
        """
        if chem_labels is None:
            chem_labels = self.chem_data.get_all_labels()

        # Get the concentration data for all the requested chemicals
        data = [self.system_snapshot(chem_label=chem, cartesian=False)
                    for chem in chem_labels]
        # Note: no need to reverse the rows in the dataframes; the heatmaps will take care of that, if requested

        if (n_chem := len(chem_labels)) == 1:
            title = f"System state at time t={self.system_time:.5g} for `{chem_labels[0]}`"
        else:
            title = f"System state at time t={self.system_time:.5g} for {n_chem} chemicals:"

        if title_prefix:
            title = f"{title_prefix}.  {title}"

        if colors is None:  # Attempt to make use of the previously-registered colors, if available
            colors = self.chem_data.get_registered_colors(chem_labels)

        return PlotlyHelper.heatmap_grid(array_list=data, labels=chem_labels, title=title,
                                         height=height, colors=colors, z_name="Conc.", max_n_cols=4,
                                         cartesian=cartesian)



    def plot_history_single_bin(self, bin_address :(int,int), colors=None, title=None, smoothed=False) -> pgo.Figure:
        """
        Using plotly, draw the plots of chemical concentration values over time at the specified bin,
        based on historical data that was saved when running simulations.
        The plotting will involve the chemicals for which history-keeping was requested for this given bin.

        Note: if this plot is to be later combined with others, use PlotlyHelper.combine_plots()
              EXAMPLE:
                    from life123 import PlotlyHelper
                    p1 = plot_history(various args, show=False)
                    p2 = plot_history(various args, show=False)
                    PlotlyHelper.combine_plots([p1, p2], other optional args)

        :param bin_address: A single bin address (a pair of integers)
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
                chem_label = f"chemical `{self.chem_data.get_label(0)}`"    # The label of the only chemical in the system
            else:
                chem_label = "all chemicals"

            title = f"Concentration changes with time of {chem_label} at bin (x={bin_address[0]}, y={bin_address[1]})"

        df = self.conc_history.bin_history(bin_address = bin_address)
        if type(df) == str:         # No data was found
            raise Exception(df)

        if colors is None:  # Attempt to make sure of the previously-registered colors, if available
            colors = self.chem_data.get_registered_colors(self.conc_history.restrict_chemicals)

        return PlotlyHelper.plot_pandas(df, x_var="SYSTEM TIME", y_label="Concentration",
                                        colors=colors, legend_header="Chemical", title=title,
                                        smoothed=smoothed)
