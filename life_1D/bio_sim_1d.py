import numpy as np
import math
from modules.html_log.html_log import HtmlLog as log



class BioSim1D:
    """
    Note: for at least the time being, this class doesn't get instantiated
    """

    #####################
    #  Class variables  #
    #####################

    n_bins = 0          # Number of spacial compartments (bins) used in the simulation

    n_species = 1       # The number of (non-water) chemical species

    univ = None         # "Universe"/Container/System
                        # NumPy array of dimension: (n_species x n_cells).
                        # Each row represents a species

    diffusion_rates = None  # NumPy array of diffusion rates for the various species

    sealed = True       # If True, no exchange with the outside; if False, immersed in a "bath"

    # Only applicable if "sealed" is False:
    bath_concentrations = None      # A NumPy array for each species
    container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container

    time_step_threshold = 0.33333   # This is used to set an Upper Bound on the single time steps
                                    # in the diffusion process.
                                    # See explanation in file overly_large_single_timesteps.py



    @classmethod
    def initialize_universe(cls, n_bins: int, n_species: int) -> None:
        """
        Initialize the entire system, EXCEPT the diffusion rates.  TODO: offer the option to do here
        All concentrations set to zero.

        :param n_bins:      The number of compartments (bins) to use in the simulation
        :param n_species:   The number of (non-water) chemical species.  It must be at least 1
        :return:            None
        """
        assert n_bins >= 1, "The number of bins must be at least 1"
        assert n_species >= 1, "The number of (non-water) chemical species must be at least 1"

        cls.n_bins = n_bins
        cls.n_species = n_species

        cls.univ = np.zeros((n_species, n_bins), dtype=float)
        cls.diffusion_rates = None



    @classmethod
    def lookup_species(cls, index: int) -> np.array:
        """
        Return the array of concentration values for the single specified species

        :param index:
        :return:
        """
        assert 0 <= index < cls.n_species, f"The species index must be in the range [0-{cls.n_species - 1}]"
        return cls.univ[index]



    @classmethod
    def describe_state(cls, show_diffusion_rates=False, concise=False) -> None:
        """
        A simple printout of the state of the system, for now useful only for small systems

        :param show_diffusion_rates:
        :param concise: Only applicable if show_diff is False.
                        If True, don't include a header with the number of bins and number of species
        :return:
        """
        if show_diffusion_rates:
            if cls.diffusion_rates is None:
                print("Diffusion rates not yet set!")
                cls.describe_state(show_diffusion_rates=False)
            else  :
                for species in range(cls.n_species):
                    print(f"Species {species}. Diff rate: {cls.diffusion_rates[species]}. Conc: ", cls.univ[species])

        else:
            if concise:
                print(cls.univ)
            else:
                if cls.n_species == 1:
                    print(f"{cls.n_bins} bins and 1 species: ", cls.univ)
                else:
                    print(f"{cls.n_bins} bins and {cls.n_species} species:\n", cls.univ)



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

        cls.univ[species_index] = np.full(cls.n_bins, conc, dtype=float)



    @classmethod
    def set_bin_conc(cls, bin: int, species_index: int, conc: float) -> None:
        """
        Assign the requested concentration value to the cell with the given index, for the specified species

        :param bin:             The zero-based bin number of the desired cell
        :param species_index:   Zero-based index to identify a specific chemical species
        :param conc:            The desired concentration value to assign to the specified location
        :return:                None
        """
        assert bin < cls.n_bins, f"The requested cell index ({bin}) must be in the range [0 - {cls.n_bins - 1}]"
        assert species_index < cls.n_species, f"The requested species index ({bin}) must be in the range [0 - {cls.n_species - 1}]"

        assert conc >= 0., f"The concentration must be a positive number or zero (the requested value was {conc})"

        cls.univ[species_index, bin] = conc


    @classmethod
    def set_species_conc(cls, species_index: int, conc_list: list) -> None:
        """
        Assign the requested concentration value to the cell with the given index, for the specified species

        :param species_index:   Zero-based index to identify a specific chemical species
        :param conc_list:            The desired concentration value to assign to the specified location
        :return:                None
        """
        assert species_index < cls.n_species, f"The requested species index ({bin}) must be in the range [0 - {cls.n_species - 1}]"

        cls.univ[species_index] = conc_list



    @classmethod
    def inject_conc_to_cell(cls, bin: int, species_index: int, delta_conc: float, zero_clip = True) -> None:
        """
        Add the requested concentration to the cell with the given index, for the specified species

        :param bin:      The zero-based bin number of the desired cell
        :param species_index:   Zero-based index to identify a specific chemical species
        :param delta_conc:      The concentration to add to the specified location
        :param zero_clip:       If True, any requested increment causing a concentration dip below zero, will make the concentration zero;
                                otherwise, an Exception will be raised
        :return:                None
        """
        assert bin < cls.n_bins, f"The requested cell index ({bin}) must be in the range [0 - {cls.n_bins - 1}]"

        if (cls.univ[species_index, bin] + delta_conc) < 0. :
            if zero_clip:
                cls.univ[species_index, bin] = 0
            else:
                raise Exception("The requested concentration change would result in a negative final value")

        # Normal scenario, not leading to negative values for the final concentration
        cls.univ[species_index, bin] += delta_conc



    #########################################################################
    #                                                                       #
    #                        PERFORM DIFFUSION                              #
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

        :param diff_rate:
        :return:
        """
        return cls.time_step_threshold/diff_rate



    @classmethod
    def diffuse(cls, time_duration=None, time_step=None, n_steps=None, verbose=False) -> dict:
        """
        Uniform-step diffusion, until reach, or just exceeding, the desired time duration.

        :param time_duration:
        :param time_step:
        :param n_steps:
        :param verbose:
        :return:                A dictionary with data about the status of the operation
        """
        assert (not time_duration or not time_step or not n_steps), \
                        "Cannot specify all 3 arguments: time_duration, time_step, n_steps"

        assert (time_duration and time_step) or (time_duration and n_steps) or (time_step and n_steps), \
                        "Must provide exactly 2 arguments from:  time_duration, time_step, n_steps"

        if not time_step:
            time_step = time_duration / n_steps

        if not n_steps:
            n_steps = math.ceil(time_duration / time_step)

        for i in range(n_steps):
            if verbose:
                if (i < 2) or (i >= n_steps-2):
                    print(f"    Performing diffusion step {i}...")
                elif i == 2:
                    print("    ...")

            cls.diffuse_step(time_step)

        if verbose:
            print(f"\nSystem after Delta time {time_duration}, at end of {n_steps} steps of size {time_step}:")
            cls.describe_state(concise=True)
            print()

        status = {"steps": n_steps}
        return status



    @classmethod
    def diffuse_step(cls, time_step) -> None:
        """
        Diffuse all the species by the given time step

        :param time_step:   Time step over which to carry out the diffusion.
                            If too large, an Exception will be raised.
        :return:            None
        """
        # TODO: parallelize the independent computations
        for species_index in range(cls.n_species):
            cls.diffuse_step_single_species(time_step, species_index=species_index)



    @classmethod
    def diffuse_step_single_species(cls, time_step, species_index=0) -> None:
        """
        Diffuse the specified single species, for the specified time step.
        We're assuming an isolated environment, with nothing diffusing thru the "walls"

        :param time_step:       Time step over which to carry out the diffusion.
                                If too large, an Exception will be raised.
        :param species_index:   ID (in the form of an integer index) of the chemical species under consideration
        :return:                None
        """
        assert cls.univ is not None, "Must first initialize the system"
        assert cls.n_bins > 0, "Must first set the number of bins"
        assert cls.diffusion_rates is not None, "Must first set the diffusion rates"
        assert cls.sealed == True, "For now, there's no provision for exchange with the outside"

        if cls.n_bins == 1:
            return                  # There's nothing to do in the case of just 1 cell!

        diff = cls.diffusion_rates[species_index]   # The diffusion rate of the specified single species

        assert not cls.is_excessive(time_step, diff), f"Excessive large time_fraction. Should be < {cls.max_time_step(diff)}"


        max_bin_number = cls.n_bins - 1     # Bin numbers range from 0 to max_bin_number, inclusive


        # Carry out a convolution operation in place, with a tile of size 3 (or 2 if only 2 bins)

        effective_diff = diff * time_step

        #print(f"Diffusing species # {species_index}")
        prev_conc = 0   # Not necessary; just to stop complaints from the code analyzer

        for i in range(cls.n_bins):    # Bin number, ranging from 0 to max_bin_number, inclusive
            #print(f"Processing bin number {i}")
            current_conc = cls.univ[species_index , i]

            if i == 0 :                     # Special case for the first bin (no left neighbor)
                cls.univ[species_index , i] +=  \
                                effective_diff * (cls.univ[species_index , i + 1] - current_conc)
            elif i == max_bin_number :      # Special case for the last bin (no right neighbor)
                cls.univ[species_index , i] += \
                                effective_diff * (prev_conc - current_conc)
            else:
                cls.univ[species_index , i] +=  \
                                 effective_diff * (cls.univ[species_index , i + 1] - current_conc) \
                               + effective_diff * (prev_conc - current_conc)

            prev_conc = current_conc




    #########################################################################
    #                                                                       #
    #                           VISUALIZATION                               #
    #                                                                       #
    #########################################################################

    @classmethod
    def single_species_heatmap(cls, species_index: int, heatmap_pars: dict, header=None) -> None:
        """
        Send to the log a heatmap representation of the concentrations of
        the specified species

        :param species_index:   Index identifying the species of interest
        :param heatmap_pars:    A dictionary of parameters for the heatmap
        :param header:          Optional string to display just above the heatmap
        :return:
        """
        if header:
            log.write(f"{header}", style=log.h1, newline=False)

        my_groups = [str(i) for i in range(cls.n_bins)]
        #print()
        #print(my_groups)

        my_data = [{"group": str(i), "variable": "Mol 0", "value": str(cls.univ[species_index, i])}
                   for i in range(cls.n_bins)]
        #print(my_data)

        all_data = {
            "my_groups": my_groups,
            "my_vars": [f"Mol {species_index}"],
            "my_data": my_data,
            "range_min": heatmap_pars["range"][0],
            "range_max": heatmap_pars["range"][1],
            "outer_width": heatmap_pars["outer_width"],
            "outer_height": heatmap_pars["outer_height"],
            "margins": heatmap_pars["margins"]
        }

        log.export_plot_Vue(data=all_data,
                            component_name="vue-heatmap-9",
                            component_file="../../../modules/visualization/vue_components/heatmap9.js")
