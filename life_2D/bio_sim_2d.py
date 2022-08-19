import numpy as np


class BioSim2D:
    """
    Note: for at least the time being, this class doesn't get instantiated
    """

    #####################
    #  Class variables  #
    #####################

    n_bins_x = 0       # Number of x-direction spacial compartments (bins) used in the simulation
    n_bins_y = 0       # Number of y-direction spacial compartments (bins) used in the simulation

    n_species = 1       # The number of (non-water) chemical species

    system = None       # Concentration data in the System we're simulating, for all the chemicals
                        # NumPy array of dimension (n_species x n_cells_x x n_cells_y)
                        # Each plane represents a species

    diffusion_rates = None  # NumPy array of diffusion rates for the various species

    sealed = True       # If True, no exchange with the outside; if False, immersed in a "bath"

    # Only applicable if "sealed" is False:
    bath_concentrations = None      # A NumPy array for each species
    container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container


    all_reactions = None            # Object of class "Reactions"

    system_time = None              # Global time of the system, from initialization on




    #########################################################################
    #                                                                       #
    #                     SET/MODIFY CONCENTRATIONS                         #
    #                                                                       #
    #########################################################################

    @classmethod
    def initialize_system(cls, n_bins: (int, int), chem_data, reactions=None) -> None:
        """
        Initialize all concentrations to zero.

        :param n_bins:     The number of compartments (bins) to use in the simulation,
                                in the x- and y- dimensions, as a pair of integers
        :param chem_data:   An object of class "Chemicals"
        :param reactions:   (OPTIONAL) Object of class "Reactions".  It may also be set later

        :return:            None
        """
        (n_cells_x, n_cells_y) = n_bins
        assert n_cells_x >= 1, "The number of bins must be at least 1 in either dimension"
        assert n_cells_y >= 1, "The number of bins must be at least 1 in either dimension"

        cls.n_bins_x = n_cells_x
        cls.n_bins_y = n_cells_y
        cls.n_species = chem_data.n_species

        # Initialize all concentrations to zero
        cls.system = np.zeros((cls.n_species, n_cells_y, n_cells_x), dtype=float)

        cls.diffusion_rates = None
        cls.names = None
        cls.chem_data = chem_data

        if reactions:
            cls.all_reactions = reactions

        cls.system_time = 0             # "Start the clock"



    @classmethod
    def set_diffusion_rates(cls, diff_list: list) -> None:
        """
        Set the diffusion rates of all the chemical species at once

        :param diff_list:   List of diffusion rates, in index order
        :return:            None
        """
        assert cls.n_species > 0, "Must first call initialize_system()"
        assert len(diff_list) == cls.n_species, \
            "The number of items in the diffusion list must equal the registered number of species"
        cls.diffusion_rates = np.array(diff_list, dtype=float)
