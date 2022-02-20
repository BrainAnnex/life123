import numpy as np


class BioSim3D:
    """
    Note: for at least the time being, this class doesn't get instantiated
    """

    #####################
    #  Class variables  #
    #####################

    n_cells_x = 0       # Number of x-direction spacial compartments (bins) used in the simulation
    n_cells_y = 0       # Number of y-direction spacial compartments (bins) used in the simulation
    n_cells_z = 0       # Number of z-direction spacial compartments (bins) used in the simulation

    n_species = 1       # The number of (non-water) chemical species

    univ = None         # NumPy array of dimension (n_species x n_cells_x x n_cells_y x n_cells_z)
    # Each block represents a species

    diffusion_rates = None  # NumPy array of diffusion rates for the various species

    sealed = True       # If True, no exchange with the outside; if False, immersed in a "bath"

    # Only applicable if "sealed" is False:
    bath_concentrations = None      # A NumPy array for each species
    container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container



    @classmethod
    def initialize_universe(cls, n_cells: (int, int, int), n_species: int) -> None:
        """

        :param n_cells:     The number of compartments (bins) to use in the simulation,
                                in the x-, y- and z- dimensions, as a triplet of integers
        :param n_species:   The number of (non-water) chemical species.  It must be at least 1
        :return:            None
        """
        pass