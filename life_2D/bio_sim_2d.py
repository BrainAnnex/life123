import numpy as np


class BioSim2D:
    """
    Note: for at least the time being, this class doesn't get instantiated
    """

    # Class variables

    n_cells_x = 0       # Number of x-direction spacial compartments (bins) used in the simulation
    n_cells_y = 0       # Number of y-direction spacial compartments (bins) used in the simulation

    n_species = 1       # The number of (non-water) chemical species

    univ = None         # NumPy array of dimension (n_species x n_cells_x x n_cells_y)
                        # Each plane represents a species

    diffusion_rates = None  # NumPy array of diffusion rates for the various species



    @classmethod
    def initialize_universe(cls, n_cells: (int, int), n_species: int) -> None:
        """

        :param n_cells:     The number of compartments (bins) to use in the simulation,
                                in the x- and y- dimensions, as a pair of integers
        :param n_species:   The number of (non-water) chemical species.  It must be at least 1
        :return:            None
        """
        (n_cells_x, n_cells_y) = n_cells
        assert n_cells_x >= 1, "The number of cells must be at least 1 in either dimension"
        assert n_cells_y >= 1, "The number of cells must be at least 1 in either dimension"
        assert n_species >= 1, "The number of (non-water) chemical species must be at least 1"

        cls.n_cells_x = n_cells_x
        cls.n_cells_y = n_cells_y
        cls.n_species = n_species

        cls.univ = np.zeros((n_cells_y, n_cells_x, n_species), dtype=float)



    @classmethod
    def set_diffusion_rates(cls, diff_list: list) -> None:
        """
        Set the diffusion rates of all the chemical species at once

        :param diff_list:   List of diffusion rates, in index order
        :return:            None
        """
        assert cls.n_species > 0, "Must first call initialize_universe()"
        assert len(diff_list) == cls.n_species, \
            "The number of items in the diffusion list must equal the registered number of species"
        cls.diffusion_rates = np.array(diff_list, dtype=float)
