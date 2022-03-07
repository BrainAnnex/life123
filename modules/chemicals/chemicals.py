import numpy as np


class Chemicals:
    """
    Object with info on the individual chemicals, incl. their Names and Diffusion rates
    """

    def __init__(self, n_species=0, diffusion_rates=None, names=None):
        """

        :param n_species:
        :param diffusion_rates:
        :param names:
        """
        self.diffusion_rates = None     # NumPy array of diffusion rates for the various species
        self.names = None               # List of the names of the various species
        self.name_dict = {}             # To map assigned names to their index

        self.n_species = n_species      # OPTIONAL.  This can be automatically set up by various calls,
                                        #            such as set_diffusion_rates() or set_names()
                                        # Over-written by diffusion_rates, if present

        if diffusion_rates:
            self.set_diffusion_rates(diffusion_rates)

        if names:
            self.set_names(names)



    def set_diffusion_rates(self, diff_list: list) -> None:
        """
        Set the diffusion rates of all the chemical species at once

        :param diff_list:   List of diffusion rates, in index order
        :return:            None
        """
        if self.n_species:
            assert len(diff_list) == self.n_species, \
                f"The number of items in the diffusion list must equal the previously declared number of species ({self.n_species})"
        else:
            self.n_species = len(diff_list)

        self.diffusion_rates = np.array(diff_list, dtype=float)



    def set_names(self, name_list)  -> None:
        """
        Set the names of all the chemical species, given in index order

        :param name_list:   List of the names of the chemical species, in index order
        :return:            None
        """
        if self.n_species:
            assert len(name_list) == self.n_species, \
                "fThe number of items in the list of name must equal the previously declared number of species ({self.n_species})"
        else:
            self.n_species = len(name_list)

        self.names = name_list

        for i, name in enumerate(name_list):
            self.name_dict[name] = i



    def get_name(self, species_index: int) -> str:
        """
        Return the name of the species with the given index.
        If no name was assigned, return an empty string.

        :param species_index:
        :return:                The name of the species with the given index if present,
                                    or an empty string if not
        """
        try:
            return self.names[species_index]
        except Exception as ex:
            return ""


    def get_index(self, species_name: str) -> int:
        """
        Return the index of the species with the given name.

        :param species_name:
        :return:                The index of the species with the given name
        """
        return self.name_dict[species_name]
