###  NOT IN CURRENT USE!

class BioSim1D:

    @classmethod
    def reaction_step_NOT_YET_USED(cls, delta_time: float) -> None:
        """
        For each reaction, carry it out in all bins - based on the initial concentrations,
        which are used as the basis for all the reactions.
        And finally add up all the individual contributions and update the concentrations.

        :param delta_time:
        :return:
        """
        number_reactions = cls.all_reactions.number_of_reactions()

        cumulative_increments = np.zeros((cls.n_species, cls.n_bins), dtype=float)  # A clone of the system state shape, with all zeros

        # For each reaction
        for reaction_id in range(number_reactions):     # Ranging from 0 to number_reactions-1
            if cls.verbose:
                print(f"Processing the {reaction_id}-th reaction")

            cls.single_reaction_NOT_YET_USED(reaction_id, delta_time)


    @classmethod
    def single_reaction_NOT_YET_USED(cls, reaction_id, delta_time):
        # For each bin_number
        for bin_n in range(cls.n_bins):     # Bin number, ranging from 0 to max_bin_number, inclusive
            if cls.verbose:
                print(f"Processing the reaction in bin number {bin_n}")



    @classmethod
    def reaction_step_old(cls, time_step: float) -> None:
        """
        NOTE: THIS WILL NOT WORK CORRECTLY FOR MULTIPLE REACTIONS, but should be correctable if adding an array of DELTA_conc
        :param time_step:
        :return:
        """

        number_reactions = cls.all_reactions.number_of_reactions()

        # TODO: loop over the bins FIRST!  -> Done in new reaction_step()
        for i in range(number_reactions):
            print(f"Evaluating reaction number {i}")
            reactants = cls.all_reactions.get_reactants(i)
            products = cls.all_reactions.get_products(i)
            fwd_rate = cls.all_reactions.get_forward_rate(i)
            back_rate = cls.all_reactions.get_reverse_rate(i)

            for bin_n in range(cls.n_bins):    # Bin number, ranging from 0 to max_bin_number, inclusive
                #print(f"    processing the reaction in bin_number number {bin_n}")
                delta_fwd = time_step * fwd_rate
                for r in reactants:
                    stoichiometry, species_index, order = r
                    conc = cls.univ[species_index , bin_n]
                    delta_fwd *= conc ** order      # Raise to power

                delta_back = time_step * back_rate
                for p in products:
                    stoichiometry, species_index, order = p
                    conc = cls.univ[species_index , bin_n]
                    delta_back *= conc ** order     # Raise to power

                print(f"    delta_fwd: {delta_fwd} | delta_back: {delta_back}")

                # Adjust the concentrations based on the forward reaction;
                #   the reactants decrease and the products increase
                for r in reactants:
                    stoichiometry, species_index, order = r
                    cls.univ[species_index , bin_n] -= delta_fwd * stoichiometry

                for p in products:
                    stoichiometry, species_index, order = p
                    cls.univ[species_index , bin_n] += delta_fwd * stoichiometry

                # Adjust the concentrations based on the back reaction;
                #   the reactants increase and the products decrease
                for r in reactants:
                    stoichiometry, species_index, order = r
                    cls.univ[species_index , bin_n] += delta_back * stoichiometry

                for p in products:
                    stoichiometry, species_index, order = p
                    cls.univ[species_index , bin_n] -= delta_back * stoichiometry
