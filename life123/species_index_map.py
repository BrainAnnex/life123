
class SpeciesIndexMap:
    """
    Manage Species Index for high-level modules
    that keep a system state that utilizes a subset of the overall species in a registry.

    A pair of indexes is used to reconcile the species id's to their index position in the system state array

    ALTERNATE NAME IDEAS:  * IndexSpecies   * SpeciesIndex
    """


    def __init__(self, species_ids=None):
        """

        :param species_ids: [OPTIONAL] Set, list or tuple of the ID's of species whose indexing needs to be managed;
                                more can be added later
        """
        # Pair of indexes to reconcile the species id's to their index position in the system state array
        self.index_to_species: list[str] = []           # EXAMPLE: ["Species A", "Species X"]
        self.species_to_index: dict[str, int] = {}      # EXAMPLE: {"Species A": 0, "Species X": 1}

        if species_ids is not None:
            self.add_species(species_ids)




    def __str__(self):
        return f"""Object with the following 2 mappings:
            species_to_index: {self.species_to_index}
            index_to_species: {self.index_to_species}
               """



    def number_of_system_species(self) -> int:
        """
        Number of species being simulated (and kept in the system state)

        :return:
        """
        return len(self.index_to_species)



    def index_of(self, species_id :str, enforce=True) -> int:
        """
        Locate and return the array index (meant for system-state arrays)
        associated to the given species ID.

        :param species_id:  A species ID
        :return:
        """
        species_index = self.species_to_index.get(species_id, None)

        if enforce:
            assert species_index is not None, \
                f'IndexSpecies.index_of(): no species with id "{species_id}" is currently registered ' \
                f'in the system-state array. \nDid you include it?'

        return species_index


    def species_at(self, species_index :int) -> str:
        """
        Locate and return the species ID associated to the given index
        (typically, index used for storage in system-state arrays)

        :param species_index:
        :return:                A species ID
        """
        try:
            species_id = self.index_to_species[species_index]
        except IndexError:
            raise IndexError(f"IndexSpecies.species_at(): there is no species linked "
                            f"to index value {species_index} in the system-state array.  \n"
                            f"Did you add all the species?")

        return species_id



    def add_species(self, species_ids : set[str] | list[str] | tuple[str, ...]) -> int:
        """
        Specify that the given group of species id's should be included in the index management
        Note: the calling module will typically need to expand the system state array accordingly

        :param species_ids: Set, list or tuple of the ID's of species whose indexing needs to be managed;
                                no harm in also including species that are already managed (they will be ignored)
        :return:            The number of newly-managed species
        """
        species_id_list = sorted(list(species_ids))     # The sorting is just for UX reasons
                                                        # TODO: maybe only sort sets, and follow user's order in lists/tuples

        number_added = 0
        for i, sp_id in enumerate(species_id_list):
            if self.species_to_index.get(sp_id) is not None:
                continue        # Already indexed this species

            new_index = len(self.index_to_species)
            self.index_to_species.append(sp_id)
            self.species_to_index[sp_id] = new_index

            number_added += 1

        return number_added



    def clear_index(self):
        self.index_to_species = []
        self.species_to_index = {}
