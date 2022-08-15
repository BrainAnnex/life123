from typing import Union
import numpy as np


class Chemicals:
    """
    Object with info on the individual chemicals, incl. their Names and Diffusion rates.

    The chemicals are assigned an index position (starting from zero)
    based on order with which they were first added.
    """

    def __init__(self, diffusion_rates=None, names=None):
        """
        If diffusion rates and names are both provided, they must have the same count,
        and the passed elements must appear in the same order.
        It's ok not to pass any data yet, and later add it.

        :param diffusion_rates: A list or tuple with the diffusion rates of the chemicals
        :param names:           A list with the names of the chemicals
        """
        self.diffusion_rates = None     # NumPy array of diffusion rates for the various species
        self.names = None               # List of the names of the various species
        self.name_dict = {}             # To map assigned names to their positional index (in the list of chemicals)
                                        #       This is automatically set and maintained

        self.n_species = 0              # The number of chemicals - exclusive of water

        if diffusion_rates and names:
            assert len(diffusion_rates) == len(names), \
                "Chemicals instantiation: the supplied number of diffusion_rates and names don't match"

        if diffusion_rates:
            self.set_diffusion_rates(diffusion_rates)

        if names:
            self.set_names(names)



    def number_of_chemicals(self) -> int:
        # Return the number of registered chemicals - exclusive of water
        return self.n_species



    def set_diffusion_rates(self, diff_list: list) -> None:
        """
        Set the diffusion rates of all the chemical species, given in index order.
        If chemicals names have already been registered, then the number of diffusion rates
        must match that of the names.

        :param diff_list:   List or tuple of diffusion rates, in index order
        :return:            None
        """
        # Validate
        assert self.diffusion_rates is None, \
            f"Chemicals.set_diffusion_rates(): can be invoked only if no diffusion rates were previously set"

        arg_type = type(diff_list)
        assert arg_type == list or arg_type == tuple,   \
            f"Chemicals.set_diffusion_rates(): the diffusion rates must be a list or tuple.  What was passed was of type {arg_type}"

        if self.names is not None:
            assert len(diff_list) == len(self.names), \
                f"Chemicals.set_diffusion_rates(): the number of the passed diffusion rates ({len(diff_list)}) " \
                f"doesn't match that of the registered chemicals ({len(self.names)})"


        self.diffusion_rates = np.array(diff_list, dtype=float)

        self.n_species = len(diff_list)



    def set_names(self, name_list)  -> None:
        """
        Set the names of all the chemical species, given in index order.
        If diffusion rates have already been registered, then the number of names
        must match that of the diffusion rates.

        :param name_list:   List or tuple of the names of the chemical species, in index order
        :return:            None
        """
        # Validate
        assert self.names is None, \
            f"Chemicals.set_names(): can be invoked only if no names for the chemical species were previously set"

        arg_type = type(name_list)
        assert arg_type == list or arg_type == tuple, \
            f"Chemicals.set_names(): the names must be a list or tuple.  What was passed was of type {arg_type}"

        if self.diffusion_rates is not None:
            assert len(name_list) == len(self.diffusion_rates), \
                f"Chemicals.set_names(): the number of the passed names ({len(name_list)}) " \
                f"doesn't match that of the registered diffusion rates ({len(self.diffusion_rates)})"


        self.names = name_list

        self.n_species = len(name_list)

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


    def get_all_names(self) -> [str]:
        """
        Return a list with the names of all the chemical species, in their index order

        :return:    A list of strings
        """
        return self.names



    def get_index(self, name: str) -> Union[int, None]:
        """
        Return the index of the chemical species with the given name,
        or None if not found

        :param name:    Name of the chemical species of interest
        :return:        The index of the species with the given name
        """
        return self.name_dict.get(name, None)



    def get_diffusion_rate(self, species_index: int) -> Union[str, None]:
        """
        Return the diffusion rate of the species with the given index.
        If no name was assigned, return None.

        :param species_index:   An integer (starting with zero) corresponding to the
                                    original order with which the chemical species were first added
        :return:                The value of the diffusion rate for the species with the given index if present,
                                    or None if not
        """
        assert type(species_index) == int, \
            f"Chemicals.get_diffusion_rate(): the argument `species_index` must be an integer (value passed was: {species_index})"

        assert species_index >= 0, \
            f"Chemicals.get_diffusion_rate(): the argument `species_index` must be a non-negative integer (value passed was: {species_index})"

        try:
            return self.diffusion_rates[species_index]
        except Exception as ex:
            return None



    def add_chemical(self, name: str, diffusion_rate: float) -> None:
        """
        Register one more chemical species, with a name and a diffusion rate.
        This can only be done if an EQUAL number of names and diffusion rates were set;
        if not the case, use set_names() or set_diffusion_rates() instead

        :param name:            Name of the chemical species to add to this object
        :param diffusion_rate:  Float value for the diffusion rate (in water) of this chemical
        :return:                None
        """
        if self.names is None or self.diffusion_rates is None or \
            len(self.names) != np.size(self.diffusion_rates):
            raise Exception("Chemicals.add_chemical() may only be used when an equal number of names and diffusion rates were set")

        assert type(diffusion_rate) == int or type(diffusion_rate) == float, \
            "Chemicals.add_chemical(): the diffusion rate argument, if provided, must be a number"
        assert diffusion_rate >= 0., \
            "Chemicals.add_chemical(): the diffusion rate argument cannot be negative"
        assert type(name) == str, \
            f"Chemicals.add_chemical(): a name must be provided, as a string value.  Value passed was {name}"

        self.diffusion_rates = np.append(self.diffusion_rates, diffusion_rate)   # Append to the NumPy array

        self.name_dict[name] = len(self.names)      # EXAMPLE: append to dictionary the entry 'some name': 123
        self.names.append(name)

        self.n_species += 1