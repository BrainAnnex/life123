import numpy as np


class BioSim:

    # Class variables
    n_cells = 0         # Number of spacial compartments used in the simulation
                        #   For now fixed; in later iterations, it will be adjustable

    n_species = 1       # The number of (non-water) chemical species


    @classmethod
    def initialize_universe(cls, n_cells: int, n_species: int) -> None:
        """

        :param n_cells:     The umber of spacial compartments to use in the simulation
        :param n_species:   The number of (non-water) chemical species.  It must be at least 1
        :return:            None
        """
        cls.n_cells = n_cells

        assert n_species >= 1, "The number of (non-water) chemical species must be at least 1"
        cls.n_species = n_species
