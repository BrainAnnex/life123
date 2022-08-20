import numpy as np
from modules.movies.movies import Movie
from modules.reactions.reactions import Reactions


class BioSim2D:
    """
    Note: for at least the time being, this class doesn't get instantiated
    """

    #####################
    #  Class variables  #
    #####################

    n_bins_x = 0       # Number of x-direction spacial compartments (bins) used in the simulation
    n_bins_y = 0       # Number of y-direction spacial compartments (bins) used in the simulation

    n_species = 1       # The number of (non-water) chemical species    TODO: phase out?

    chem_data = None    # Object of type "Chemicals", with info on the individual chemicals, incl. their names

    system = None       # Concentration data in the System we're simulating, for all the chemicals
                        # NumPy array of dimension (n_species x n_bins_x x n_bins_y)
                        # Each plane represents a species

    # The following buffers are (n_species x n_bins_x x n_bins_y)
    delta_diffusion = None  # Buffer for the concentration changes from diffusion step
    delta_reactions = None  # Buffer for the concentration changes from reactions step

    sealed = True       # If True, no exchange with the outside; if False, immersed in a "bath"

    # Only applicable if "sealed" is False:
    bath_concentrations = None      # A NumPy array for each species
    container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container


    all_reactions = None            # Object of class "Reactions"

    system_time = None              # Global time of the system, from initialization on

    debug = False



    #########################################################################
    #                                                                       #
    #                     SET/MODIFY CONCENTRATIONS                         #
    #                                                                       #
    #########################################################################

    @classmethod
    def initialize_system(cls, n_bins: (int, int), chem_data=None, reactions=None) -> None:
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
            cls.chem_data = chem_data
        else:
            cls.chem_data = reactions.chem_data

        if reactions:
            cls.all_reactions = reactions
        else:
            cls.all_reactions = Reactions(chem_data=chem_data)

        cls.n_bins_x = n_cells_x
        cls.n_bins_y = n_cells_y

        cls.n_species = chem_data.n_species

        # Initialize all concentrations to zero
        cls.system = np.zeros((cls.n_species, n_cells_x, n_cells_y), dtype=float)

        cls.system_time = 0             # "Start the clock"



    @classmethod
    def set_bin_conc(cls, bin_x: int, bin_y: int, species_index: int, conc: float) -> None:
        """
        Assign the requested concentration value to the cell with the given index, for the specified species

        :param bin_x:           The zero-based bin number of the desired cell
        :param bin_y:           The zero-based bin number of the desired cell
        :param species_index:   Zero-based index to identify a specific chemical species
        :param conc:            The desired concentration value to assign to the specified location
        :return:                None
        """
        #assert bin < cls.n_bins, f"The requested cell index ({bin}) must be in the range [0 - {cls.n_bins - 1}]"
        assert species_index < cls.n_species, f"The requested species index ({bin}) must be in the range [0 - {cls.n_species - 1}]"

        assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"

        cls.system[species_index, bin_x, bin_y] = conc



    @classmethod
    def set_bin_conc_all_species(cls, bin_x: int, bin_y: int, conc_list: [float]) -> None:
        """
        Assign the requested concentration values to the cell with the given index,
        for all the species in their index order

        :param bin_x:       The zero-based bin number of the desired cell
        :param bin_y:       The zero-based bin number of the desired cell
        :param conc_list:   A list with the desired concentration values to assign to the specified location
        :return:            None
        """
        #assert bin < cls.n_bins, f"The requested cell index ({bin}) must be in the range [0 - {cls.n_bins - 1}]"

        for i, conc in enumerate(conc_list):
            assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"
            cls.system[i, bin_x, bin_y] = conc




    #########################################################################
    #                                                                       #
    #                              TO VIEW                                  #
    #                                                                       #
    #########################################################################

    @classmethod
    def describe_state(cls, concise=False) -> None:
        """

        :param concise:
        :return:
        """
        print(f"SYSTEM STATE at Time t = {cls.system_time}:")
        for species_index in range(cls.n_species):
            print(f"Species `{cls.chem_data.get_name(species_index)}`:")
            print(cls.system[species_index])




    #########################################################################
    #                                                                       #
    #                               DIFFUSION                               #
    #                                                                       #
    #########################################################################

    # TODO: TBA




    #########################################################################
    #                                                                       #
    #                               REACTIONS                               #
    #                                                                       #
    #########################################################################

    @classmethod
    def react(cls, total_duration=None, time_step=None, n_steps=None, snapshots=None) -> None:
        """
        Update the system concentrations as a result of all the reactions in all bins.
        CAUTION : NO diffusion is taken into account.

        For each bin, process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        TODO: in case of any Exception, the state of the system is still valid, as of the time before this call

        :param total_duration:
        :param time_step:
        :param n_steps:
        :param snapshots:       OPTIONAL dict with the keys: "frequency", "sample_bin", "sample_species"
                                    If provided, take a system snapshot after running a multiple of "frequency" runs
        :return:                None
        """

        time_step, n_steps = cls.all_reactions.specify_steps(total_duration=total_duration,
                                                             time_step=time_step,
                                                             n_steps=n_steps)
        #if snapshots is None:
            #frequency = None
        #else:
            #frequency = snapshots.get("frequency", 1)

        for i in range(n_steps):
            cls.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            cls.system += cls.delta_reactions   # Matrix operation to update all the concentrations
            cls.system_time += time_step
            #if (frequency is not None) and ((i+1)%frequency == 0):
                #cls.save_snapshot(cls.bin_snapshot(bin_address = snapshots["sample_bin"]))



    @classmethod
    def reaction_step(cls, delta_time: float) -> None:
        """
        Clear and compute the delta_reactions array (a class variable),
        based on all the reactions in all bins.
        IMPORTANT: the actual system concentrations are NOT changed.

        For each bin, process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        :param delta_time:
        :return:            None
        """
        assert cls.all_reactions is not None, \
            "reaction_step(): must first set the Reactions object"

        cls.delta_reactions = np.zeros((cls.n_species, cls.n_bins_x, cls.n_bins_y), dtype=float)

        # For each bin
        for bin_n_x in range(cls.n_bins_x):
            for bin_n_y in range(cls.n_bins_y):
                if cls.debug:
                    print(f"reaction_step(): processing the all the reactions in bin number ({bin_n_x}, {bin_n_y})")

                # Obtain the Delta-concentration for each species, for this bin
                conc_dict = {species_index: cls.system[species_index , bin_n_x, bin_n_y]
                             for species_index in range(cls.n_species)}
                if cls.debug:
                    print(f"\nconc_dict in bin ({bin_n_x}, {bin_n_y}): ", conc_dict)


                # Obtain the Delta-conc for each species, for the current bin
                increment_vector = cls.all_reactions.single_compartment_reaction_step(conc_dict=conc_dict,
                                                                                      delta_time=delta_time)

                # Replace the appropriate column of the cls.delta_reactions matrix
                # with the contents of the vector increment_vector
                cls.delta_reactions[:, bin_n_x, bin_n_y] = np.array([increment_vector])

                if cls.debug:
                    print(cls.delta_reactions)
