from typing import Union
import numpy as np


class Chemicals:
    """
    Object with info on the individual chemicals, incl. their Names and Diffusion rates.

    The chemicals are assigned an index position (starting from zero)
    based on order with which they were first added.
    """

    def __init__(self, n_species=0, diffusion_rates=None, names=None):
        """

        :param n_species:       The number of chemicals - exclusive of water
                                # NOTE: the diffusion rates and names, if both provided, must be in the same order
        :param diffusion_rates: A list with the diffusion rates of the chemicals
        :param names:           A list with the names of the chemicals
        """
        self.diffusion_rates = None     # NumPy array of diffusion rates for the various species
        self.names = None               # List of the names of the various species
        self.name_dict = {}             # To map assigned names to their positional index (in the list of chemicals)
                                        #           This is automatically set and maintained

        self.n_species = n_species      # OPTIONAL.  This can be automatically set up by various calls,
                                        #            such as set_diffusion_rates() or set_names()
                                        # Over-written by diffusion_rates, if present

        if diffusion_rates:
            self.set_diffusion_rates(diffusion_rates)

        if names:
            self.set_names(names)



    def set_diffusion_rates(self, diff_list: list) -> None:
        """
        Set the diffusion rates of all the chemical species, given in index order

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

        # Create a dictionary to map the assigned names to their positional index
        for i, name in enumerate(name_list):
            self.name_dict[name] = i



    def get_name(self, species_index: int) -> Union[str, None]:
        """
        Return the name of the species with the given index.
        If no name was assigned, return None.

        :param species_index:   An integer (starting with zero) corresponding to the
                                    original order with which the chemical species were first added
        :return:                The name of the species with the given index if present,
                                    or None if not
        """
        assert type(species_index) == int, \
            f"Chemicals.get_name(): the argument `species_index` must be an integer (value passed was: {species_index})"

        assert species_index >= 0, \
            f"Chemicals.get_name(): the argument `species_index` must be a non-negative integer (value passed was: {species_index})"

        try:
            return self.names[species_index]
        except Exception as ex:
            return None


    def get_index(self, species_name: str) -> Union[int, None]:
        """
        Return the index of the species with the given name,
        or None if not found

        :param species_name:
        :return:                The index of the species with the given name
        """
        return self.name_dict.get(species_name, None)
