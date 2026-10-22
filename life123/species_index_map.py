from life123.species_registry import SpeciesRegistry
import numpy as np
import numpy.typing as npt



class SpeciesIndexMap:
    """
    Manage Species Index for high-level modules
    that keep a system state that utilizes a subset of the overall species in a registry.

    A pair of indexes is used to reconcile the species id's to their index position in the system state array

    ALTERNATE NAME IDEAS:  * IndexSpecies   * SpeciesIndex
    """


    def __init__(self, species_ids :set|list|tuple =None, species_registry :SpeciesRegistry =None):
        """

        :param species_ids:         [OPTIONAL] Set, list or tuple of the ID's of species whose indexing needs to be managed;
                                        more can be added later
        :param species_registry:    [OPTIONAL] Object of type "SpeciesRegistry";
                                        if passed, species id's are enforced to exist in it
        """
        # TODO: maybe make `species_registry` required?

        # Pair of indexes to reconcile the species id's to their index position in the system state array
        self.index_to_species: list[str] = []           # EXAMPLE: ["Species A", "Species X"]
        self.species_to_index: dict[str, int] = {}      # EXAMPLE: {"Species A": 0, "Species X": 1}

        self.species_registry = species_registry

        if (species_ids is None) and (species_registry is not None):
            species_ids = species_registry.get_all_species_ids(sort=True)
            # TODO: if both supplied, verify compatibility

        if species_ids is not None:
            self.add_species(species_ids)



    def __str__(self):
        return f"""Object with the following 2 mappings:
            species_to_index: {self.species_to_index}
            index_to_species: {self.index_to_species}
               """



    def number_of_mapped_species(self) -> int:
        """
        Number of species being mapped
        (generally matching, or possibly larger, than those kept in the system state)

        :return:
        """
        return len(self.index_to_species)



    def index_of(self, species_id :str, enforce=True) -> int|None:
        """
        Locate and return the array index (meant for system-state arrays)
        associated to the given species ID.

        :param species_id:  A species ID
        :param enforce:
        :return:
        """
        species_index = self.species_to_index.get(species_id, None)

        if enforce:
            if species_index is None:
                raise IndexError(f'IndexSpecies.index_of(): no species with id "{species_id}" is currently registered ' \
                                 f'in the system-state array. \nDid you include it?')

        return species_index


    def species_at(self, species_index :int) -> str:
        """
        Locate and return the species ID associated to the given index
        (typically, index used for storage in system-state arrays)

        :param species_index:
        :return:                A species ID
        """
        assert species_index >= 0, \
            f"IndexSpecies.species_at(): the species index must be a non-negative number"

        try:
            species_id = self.index_to_species[species_index]   # Note: negative indices are allowed in python
        except IndexError:
            raise IndexError(f"IndexSpecies.species_at(): there is no species linked "
                            f"to index value {species_index} in the system-state array.  \n"
                            f"Did you add all the species to the index?")

        return species_id



    def get_all_species_ids(self, sort=False) -> list[str]:
        """
        Return a list with the id's of all the species,
        NOT in any particular order

        :param sort:    [OPTIONAL]
        :return:        A list of strings with the species id's,
                            NOT in any particular order (unless `sort` is True)
        """
        if sort:
            return sorted(list(self.index_to_species))

        return list(self.index_to_species)     # The dictionary keys as a list



    def add_species(self, species_ids : str | set[str] | list[str] | tuple[str, ...]) -> int:
        """
        Specify that the given group of species id's should be included in the index management
        Note: the calling module will typically need to expand the system state array accordingly

        :param species_ids: String, set, list or tuple of the ID's of species whose indexing needs to be managed;
                                no harm in also including species that are already managed (they will be ignored).
                                If a set is passed, it will be sorted in alphabetical order;
                                if you want a different specific order, pass a list or tuple in the order you want it
                                (the order has only UX significance, for displays, etc)
        :return:            The number of newly-managed species
        """
        if type(species_ids) is str:
            species_id_list = [species_ids]
        elif type(species_ids) is set:
            species_id_list = sorted(list(species_ids))     # The sorting is just for UX reasons
        else:
            species_id_list = species_ids

        number_added = 0
        for i, sp_id in enumerate(species_id_list):
            if self.species_to_index.get(sp_id) is not None:
                continue        # Already indexed this species

            if self.species_registry:
                assert self.species_registry.species_exists(sp_id), \
                        f"add_species(): attempting to index a species (`{sp_id}`) that doesn't exist in the SpeciesRegistry; \n" \
                        f"make sure to first add it there, with add_species()"

            new_index = len(self.index_to_species)
            self.index_to_species.append(sp_id)
            self.species_to_index[sp_id] = new_index

            number_added += 1

        return number_added



    def clear_index(self) -> None:
        """
        A fresh start on all indexes

        :return:    None
        """
        self.index_to_species = []
        self.species_to_index = {}



    def concentrations_to_array(self, conc: dict[str, float], dtype: npt.DTypeLike = float) -> np.ndarray:
        """
        Convert a dictionary of species concentrations into a Numpy array,
        following the species ordering in `self.index_to_species`

        :param conc:    Dictionary mapping species IDs to concentrations
        :param dtype:   [OPTIONAL] Numpy dtype for the resulting array (default: float)
        :return:        Concentration array ordered according to `self.index_to_species`
        """
        # TODO: pytest
        return np.array(
            [conc[species_id] for species_id in self.index_to_species],
            dtype=dtype
        )
