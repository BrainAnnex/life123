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

    chem_data = None    # Object with info on the individual chemicals, incl. their names

    system = None       # Concentration data in the System we're simulating, for all the chemicals
                        # NumPy array of dimension: (n_species x n_cells).
                        # Each row represents a species

    system_earlier = None   # NOT IN CURRENT USE.  Envisioned for simulations where the past 2 time states are used
                            # to compute the state at the next time step

    delta_diffusion = None  # Buffer for the concentration changes from diffusion step (n_species x n_cells)
    delta_reactions = None  # Buffer for the concentration changes from reactions step (n_species x n_cells)

    sealed = True       # If True, no exchange with the outside; if False, immersed in a "bath"

    # Only applicable if "sealed" is False:
    bath_concentrations = None      # A NumPy array for each species
    container_diffusion = None      # A NumPy array for each species: diffusion rate in/out of the container

    time_step_threshold = 0.33333   # This is used to set an Upper Bound on the single time steps
                                    # in the diffusion process.
                                    # See explanation in file overly_large_single_timesteps.py

    all_reactions = None             # Object of class "Reactions"      # TODO: add a setter method

    verbose = False



    #########################################################################
    #                                                                       #
    #                              TO VIEW                                  #
    #                                                                       #
    #########################################################################


    @classmethod
    def lookup_species(cls, index: int) -> np.array:
        """
        Return the NumPy array of concentration values across the bins (from left to right)
        for the single specified species

        :param index:
        :return:
        """
        assert 0 <= index < cls.n_species, f"The species index must be in the range [0-{cls.n_species - 1}]"
        return cls.system[index]


    @classmethod
    def bin_concentration(cls, bin_address: int, species_index: int):
        """
        Return the concentration at the requested bin of the specified species

        :param bin_address:
        :param species_index:
        :return:
        """
        return cls.system[species_index, bin_address]



    @classmethod
    def describe_state(cls, show_diffusion_rates=False, concise=False) -> None:
        """
        A simple printout of the state of the system, for now useful only for small systems

        :param show_diffusion_rates:    NO LONGER USED.  TODO: phase out
        :param concise: Only applicable if show_diffusion_rates is False.
                        If True, don't include a header with the number of bins and number of species
        :return:        None
        """
        if concise:     # A minimalist printout
            print(cls.system)
            return

        if cls.n_species == 1:
            print(f"{cls.n_bins} bins and 1 species: ")
        else:
            print(f"{cls.n_bins} bins and {cls.n_species} species:\n")


        for species_index in range(cls.n_species):
            name = cls.chem_data.get_name(species_index)
            if name:
                name = f" ({name})"

            if cls.chem_data.diffusion_rates is None:
                print(f"  Species {species_index}{name}. Diff rate: NOT SET. Conc: ", cls.system[species_index])
            else:
                print(f"  Species {species_index}{name}. Diff rate: {cls.chem_data.diffusion_rates[species_index]}. Conc: ",
                      cls.system[species_index])





    #########################################################################
    #                                                                       #
    #                     SET/MODIFY CONCENTRATIONS                         #
    #                                                                       #
    #########################################################################

    @classmethod
    def initialize_system(cls, n_bins: int, chem_data, reactions=None) -> None:
        """
        Initialize all concentrations to zero.

        TODO: maybe allow optionally passing n_species in lieu of chem_data,
              and let it create and return the "Chemicals" object in that case

        :param n_bins:      The number of compartments (bins) to use in the simulation
        :param chem_data:   An object of class "Chemicals"
        :param reactions:   (OPTIONAL) Object of class "Reactions".  It may also be set later

        :return:            None
        """
        assert n_bins >= 1, "The number of bins must be at least 1"
        assert chem_data.n_species >= 1, "The number chemical species set in the `chem_data` object must be at least 1"

        cls.n_bins = n_bins
        cls.n_species = chem_data.n_species

        cls.system = np.zeros((cls.n_species, n_bins), dtype=float)

        cls.diffusion_rates = None
        cls.names = None
        cls.chem_data = chem_data

        if reactions:
            cls.all_reactions = reactions


    @classmethod
    def replace_system(cls, new_state: np.array) -> None:
        """
        Replace the System's internal state

        :param new_state:
        :return:
        """
        cls.system = new_state
        cls.n_species, cls.n_bins = new_state.shape     # Extract from the new state



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

        cls.system[species_index] = np.full(cls.n_bins, conc, dtype=float)


    @classmethod
    def set_all_uniform_concentrations(cls, conc_list: [float]) -> None:
        """
        Set the concentrations of all species at once, uniformly across all bins
        :param conc_list:
        :return:
        """
        if cls.chem_data.n_species:
            assert len(conc_list) == cls.chem_data.n_species, \
                f"The argument to 'set_all_uniform_concentrations()' must be a list of size {cls.chem_data.n_species}"

        for i, conc in enumerate(conc_list):
            cls.set_uniform_concentration(species_index=i, conc=conc)



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

        cls.system[species_index, bin] = conc


    @classmethod
    def set_species_conc(cls, species_index: int, conc_list: list) -> None:
        """
        Assign the requested list of concentration values to all the bins, in order, for the specified species

        :param species_index:   Zero-based index to identify a specific chemical species
        :param conc_list:            The desired concentration value to assign to the specified location
        :return:                None
        """
        assert species_index < cls.n_species, f"The requested species index ({bin}) must be in the range [0 - {cls.n_species - 1}]"

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
    def diffuse(cls, time_duration=None, time_step=None, n_steps=None, verbose=False) -> dict:
        """
        Uniform-step diffusion, with 2 out of 3 criteria specified:
            1) until reaching, or just exceeding, the desired time duration
            2) using the given time step
            3) carrying out the specified number of steps

        :param time_duration:
        :param time_step:
        :param n_steps:
        :param verbose:
        :return:                A dictionary with data about the status of the operation
        """
        # TODO: factor out this part, in common to diffuse() and react()
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
            cls.system += cls.delta_diffusion     # Matrix operation

        if verbose:
            print(f"\nSystem after Delta time {time_duration}, at end of {n_steps} steps of size {time_step}:")
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
        :return:
        """
        cls.all_reactions = reactions



    @classmethod
    def react(cls, time_duration=None, time_step=None, n_steps=None) -> None:
        """
        Update the system concentrations as a result of all the reactions in all bins.

        For each bin, process all the reactions in it - based on
        the INITIAL concentrations (prior to this reaction step),
        which are used as the basis for all the reactions.

        TODO: in case of any Exception, the state of the system is still valid, as of the time before this call

        :param time_duration:
        :param time_step:
        :param n_steps:
        :return:                None
        """
        # TODO: factor out this part, in common to diffuse() and react()
        assert (not time_duration or not time_step or not n_steps), \
            "Cannot specify all 3 arguments: time_duration, time_step, n_steps"

        assert (time_duration and time_step) or (time_duration and n_steps) or (time_step and n_steps), \
            "Must provide exactly 2 arguments from:  time_duration, time_step, n_steps"

        if not time_step:
            time_step = time_duration / n_steps

        if not n_steps:
            n_steps = math.ceil(time_duration / time_step)

        for i in range(n_steps):
            cls.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            cls.system += cls.delta_reactions   # Matrix operation



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
        number_reactions = cls.all_reactions.number_of_reactions()

        cls.delta_reactions = np.zeros((cls.n_species, cls.n_bins), dtype=float)

        # For each bin
        for bin_n in range(cls.n_bins):     # Bin number, ranging from 0 to max_bin_number, inclusive
            if cls.verbose:
                print(f"Processing the reaction in bin number {bin_n}")

            # Obtain the Delta-conc for each species, for bin number bin_n
            increment_vector = cls.single_bin_reaction_step(bin_n, delta_time, number_reactions)

            # Replace the "bin_n" column of the cls.delta_reactions matrix with the contents of the vector increment_vector
            cls.delta_reactions[:, bin_n] = increment_vector.transpose()

        #print(cls.delta_reactions)



    @classmethod
    def single_bin_reaction_step(cls, bin_n: int, delta_time: float, number_reactions: int) -> np.array:
        """
        For the given bin, do a single reaction time step for ALL the reactions in it -
        based on the INITIAL concentrations in the bin (prior to this reaction step),
        which are used as the basis for all the reactions.

        IMPORTANT: the actual system concentrations are NOT changed.

        :param bin_n:
        :param delta_time:
        :param number_reactions:
        :return:                    The increment vector for all the chemical species concentrations for this bin
        """

        # Compute the forward and back conversions of all the reactions
        delta_fwd_list, delta_back_list = cls.compute_rates(bin_n, delta_time, number_reactions)
        if cls.verbose:
            print(f"    delta_fwd_list: {delta_fwd_list} | delta_back_list: {delta_back_list}")


        increment_vector = np.zeros((cls.n_species, 1), dtype=float)       # One element per species

        # For each reaction, adjust the concentrations of the reactants and products,
        # based on the forward and back rates of the reaction
        for i in range(number_reactions):
            if cls.verbose:
                print(f"    adjusting the species concentrations based on reaction number {i}")

            # TODO: turn into a more efficient single step, as as:
            #(reactants, products) = cls.all_reactions.unpack_terms(i)
            reactants = cls.all_reactions.get_reactants(i)
            products = cls.all_reactions.get_products(i)

            # Adjust the concentrations

            #   The reactants decrease based on the forward reaction,
            #             and increase based on the reverse reaction
            for r in reactants:
                stoichiometry, species_index, order = r
                increment_vector[species_index] += stoichiometry * (- delta_fwd_list[i] + delta_back_list[i])

                if (cls.system[species_index , bin_n] + increment_vector[species_index]) < 0:
                    raise Exception(f"The given time interval ({delta_time}) leads to negative concentrations in reactions: make it smaller!")

                #cls.univ[species_index , bin_n] += increment_vector[species_index]


            #   The products increase based on the forward reaction,
            #             and decrease based on the reverse reaction
            for p in products:
                stoichiometry, species_index, order = p
                increment_vector[species_index] += stoichiometry * (delta_fwd_list[i] - delta_back_list[i])

                if (cls.system[species_index , bin_n] + increment_vector[species_index]) < 0:
                    raise Exception(f"The given time interval ({delta_time}) leads to negative concentrations in reactions: make it smaller!")

                #cls.univ[species_index , bin_n] += increment_vector[species_index]

        # END for

        return increment_vector



    @classmethod
    def compute_rates(cls, bin_n: int, delta_time: float, number_reactions: int) -> (list, list):
        """
        For each of the reactions, compute its forward and back "contributions" (rates, multiplied by delta_time)

        :param bin_n:
        :param delta_time:
        :param number_reactions:
        :return:                A pair of lists (List of forward conversions, List of reverse conversions);
                                    each list has 1 entry per reaction, in the index order of the reactions
                                TODO: simply return the list of their differences!
                                TODO: also, make a note of large relative increments (to guide future time-step choices)
        """
        delta_fwd_list = []             # It will have 1 entry per reaction
        delta_back_list = []            # It will have 1 entry per reaction
        for i in range(number_reactions):
            if cls.verbose:
                print(f"    evaluating the rates for reaction number {i}")

            # TODO: turn into a more efficient single step, as as:
            #(reactants, products, fwd_rate_coeff, back_rate_coeff) = cls.all_reactions.unpack_reaction(i)
            reactants = cls.all_reactions.get_reactants(i)
            products = cls.all_reactions.get_products(i)
            fwd_rate_coeff = cls.all_reactions.get_forward_rate(i)
            back_rate_coeff = cls.all_reactions.get_reverse_rate(i)

            delta_fwd = delta_time * fwd_rate_coeff         # TODO: save, to avoid re-computing at each bin
            for r in reactants:
                stoichiometry, species_index, order = r
                conc = cls.system[species_index , bin_n]      # TODO: make more general for 2D and 3D
                delta_fwd *= conc ** order      # Raise to power

            delta_back = delta_time * back_rate_coeff       # TODO: save, to avoid re-computing at each bin
            for p in products:
                stoichiometry, species_index, order = p
                conc = cls.system[species_index , bin_n]      # TODO: make more general for 2D and 3D
                delta_back *= conc ** order     # Raise to power

            #print(f"    delta_fwd: {delta_fwd} | delta_back: {delta_back}")
            delta_fwd_list.append(delta_fwd)
            delta_back_list.append(delta_back)

        return( (delta_fwd_list, delta_back_list) )



    #########################################################################
    #                                                                       #
    #                         REACTION-DIFFUSION                            #
    #                                                                       #
    #########################################################################

    @classmethod
    def react_diffuse(cls, time_duration=None, time_step=None, n_steps=None) -> None:
        """
        It expects 2 of the arguments:  time_duration, time_step, n_steps
        Perform a series of reaction and diffusion time steps.

        :param time_duration:   The overall time advance (i.e. time_step * n_steps)
        :param time_step:       The size of each time step
        :param n_steps:         The desired number of steps
        :return:                None
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
            # TODO: split off the reaction step and the diffusion step to 2 different computing cores
            cls.reaction_step(time_step)        # TODO: catch Exceptions in this step; in case of failure, repeat with a smaller time_step
            cls.diffuse_step(time_step)
            # Merge into the concentrations of the various bins/chemical species pairs,
            # the increments concentrations computed separately by the reaction and the diffusion steps
            cls.system += cls.delta_reactions     # Matrix operation
            cls.system += cls.delta_diffusion     # Matrix operation




    #########################################################################
    #                                                                       #
    #                           VISUALIZATION                               #
    #                                                                       #
    #########################################################################

    @classmethod
    def single_species_heatmap(cls, species_index: int, heatmap_pars: dict, header=None) -> None:
        """
        Send to the HTML log, a heatmap representation of the concentrations of
        the single requested species.
        This version uses the Vue component heatmap11.js

        :param species_index:   Index identifying the species of interest
        :param heatmap_pars:    A dictionary of parameters for the heatmap
        :param header:          Optional string to display just above the heatmap
        :return:
        """
        if header:
            log.write(f"{header}", style=log.h1, newline=False)

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

        log.export_plot_Vue(data=all_data,
                            component_name="vue-heatmap-11",
                            component_file="../../../modules/visualization/vue_components/heatmap11.js")
