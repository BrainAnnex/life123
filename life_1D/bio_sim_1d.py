import numpy as np


class BioSim1D:
    """
    Note: for at least the time being, this class doesn't get instantiated
    """

    #####################
    #  Class variables  #
    #####################

    n_cells = 0         # Number of spacial compartments (bins) used in the simulation

    n_species = 1       # The number of (non-water) chemical species

    univ = None         # "Universe"/Container/System
                        # NumPy array of dimension (n_species x n_cells).
                        # Each row represents a species

    diffusion_rates = None  # NumPy array of diffusion rates for the various species

    sealed = True       # If True, no exchange with the outside; if False, immersed in a "bath"

    # Only applicable if "sealed" is False:
    bath_concentrations = None      # A NumPy array for each species
    container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container



    @classmethod
    def initialize_universe(cls, n_cells: int, n_species: int) -> None:
        """

        :param n_cells:     The number of compartments (bins) to use in the simulation
        :param n_species:   The number of (non-water) chemical species.  It must be at least 1
        :return:            None
        """
        assert n_cells >= 1, "The number of cells must be at least 1"
        assert n_species >= 1, "The number of (non-water) chemical species must be at least 1"


        cls.n_cells = n_cells
        cls.n_species = n_species

        cls.univ = np.zeros((n_species, n_cells), dtype=float)



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



    @classmethod
    def set_uniform_concentration(cls, species_index: int, conc: float) -> None:
        """
        Assign the given concentration to all the cells of the specified species (identified by its index).
        Any previous values get over-written

        :param species_index:   Zero-based index to identify a specific chemical species
        :param conc:            The desired value of chemical concentration for the above species
        :return:                None
        """
        assert (species_index >= 0) and (species_index < cls.n_species), \
                    f"The species_index must be an integer between 0 and {cls.n_species - 1}"

        assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"

        cls.univ[species_index] = np.full(cls.n_cells, conc, dtype=float)



    @classmethod
    def set_conc_of_cell(cls, cell_index: int, species_index: int, conc: float) -> None:
        """
        Assign the requested concentration value to the cell with the given index, for the specified species

        :param cell_index:      The zero-based bin number of the desired cell
        :param species_index:   Zero-based index to identify a specific chemical species
        :param conc:            The desired concentration value to assign to the specified location
        :return:                None
        """
        assert cell_index < cls.n_cells, f"The requested cell index ({cell_index}) must be in the range [0 - {cls.n_cells - 1}]"

        assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"

        cls.univ[species_index, cell_index] = conc



    @classmethod
    def inject_conc_to_cell(cls, cell_index: int, species_index: int, delta_conc: float, zero_clip = True) -> None:
        """
        Add the requested concentration to the cell with the given index, for the specified species

        :param cell_index:      The zero-based bin number of the desired cell
        :param species_index:   Zero-based index to identify a specific chemical species
        :param delta_conc:      The concentration to add to the specified location
        :param zero_clip:       If True, any requested increment causing a concentration dip below zero, will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        assert cell_index < cls.n_cells, f"The requested cell index ({cell_index}) must be in the range [0 - {cls.n_cells - 1}]"

        if (cls.univ[species_index, cell_index] + delta_conc) < 0. :
            if zero_clip:
                cls.univ[species_index, cell_index] = 0
            else:
                raise Exception("The requested concentration change would result in a negative final value")

        # Normal scenario, not leading to negative values for the final concentration
        cls.univ[species_index, cell_index] += delta_conc




    ####    PERFORM DIFFUSION    ####

    @classmethod
    def diffuse_step_single_species(cls, species_index=0, time_fraction=1.0):
        """
        Diffuse the specified single species, for 1 time step

        :param species_index:   ID (in the form of an integer index) of the chemical species under consideration
        :param time_fraction:   Fractional time step
        :return:                A NumPy array of floats
        """
        assert species_index == 0, "Index must be zero: for now, only 1 chemical species is implemented"
        assert cls.diffusion_rates is not None, "The diffusion rate not set yet"
        assert cls.n_cells > 0, "Must first set the number of cells"

        if cls.n_cells == 1:
            return cls.univ     # There's nothing to do in the case of just 1 cell!

        diff = cls.diffusion_rates[species_index]   # The diffusion rate of the specified single species

        effective_diff = diff * time_fraction

        assert effective_diff <= 0.3333, f"Excessive large time_fraction. Should be <= {0.3333/diff}"
                                            # A value of 0.5 would be such that 2 adjacent cells - in the
                                            # absence of other neighbors - would equilibrate in concentration.
                                            # A value of 1/3 would be such that 3 adjacent cells - in the
                                            # absence of other neighbors - would equilibrate in concentration.

        initial_conc = cls.univ
        new_conc = np.zeros(cls.n_cells, dtype=float)
        max_bin_number = cls.n_cells - 1     # Bin numbers range from 0 to max_bin_number


        # Carry out a convolution operation with a tile of size 3.
        # TODO: try a convolution in place, to avoid creating a new NumPy array
        for i in range(cls.n_cells):
            #print(f"Processing i={i}")

            if i == 0 :                 # Special case for the first bin (no left neighbor)
                new_conc[i] += initial_conc[i] + effective_diff * (initial_conc[i + 1] - initial_conc[i])
            elif i == max_bin_number :  # Special case for the last bin (no right neighbor)
                new_conc[i] += initial_conc[i] + effective_diff * (initial_conc[max_bin_number - 1] - initial_conc[max_bin_number])
            else:
                new_conc[i] += initial_conc[i] \
                               + effective_diff * (initial_conc[i + 1] - initial_conc[i]) \
                               + effective_diff * (initial_conc[i - 1] - initial_conc[i])

        cls.univ = new_conc

        return new_conc
