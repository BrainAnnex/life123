
class IndexSpecies:
    """
    Manage Species Index for high-level modules
    that keep a system state that utilizes a subset of the overall species in a registry.

    A pair of indexes is used to reconcile the species id's to their index position in the system state array

    ALTERNATE NAME IDEAS:  * SpeciesIndexMap   * SpeciesIndex
    """


    def __init__(self):

        # Pair of indexes to reconcile the species id's to their index position in the system state array
        self.index_to_species: list[str] = []           # EXAMPLE: ["Species A", "Species X"]
        self.species_to_index: dict[str, int] = {}      # EXAMPLE: {"Species A": 0, "Species X": 1}



    def number_of_system_species(self) -> int:
        """
        Number of species being simulated (and kept in the system state)
        :return:
        """
        return len(self.index_to_species)



    def locate_species_index(self, species_id :str) -> int:
        """

        :param species_id:
        :return:
        """
        species_index = self.species_to_index.get(species_id, None)

        assert species_index is not None, \
            f'UniformCompartment.locate_species_index(): no species with id "{species_id}" is currently registered ' \
            f'in the system-state array. \nDid you add reaction including it?'

        return species_index


    def locate_species_id(self, species_index :int) -> str:
        """

        :param species_index:
        :return:
        """
        try:
            species_id = self.index_to_species[species_index]
        except IndexError:
            raise Exception(f"locate_species_id(): there is no species linked "
                            f"to index value {species_index} in the system-state array.  \n"
                            f"Maybe you didn't add all the reactions?")

        return species_id



    def add_species(self, species_id_set : set[str] | list[str]) -> int:
        """
        Specify a set to species id's to add to the index management

        :param species_id_set:  Set of ID's of species whose indexing needs to be managed;
                                    no action taken for any species that is already managed
        :return:                The number of newly-managed species
        """
        species_id_list = sorted(list(species_id_set))      # The sorting is just for UX reasons

        number_added = 0
        for i, sp_id in enumerate(species_id_list):
            if self.species_to_index.get(sp_id) is not None:
                continue        # Already indexed this species

            new_index = len(self.index_to_species)
            self.index_to_species.append(sp_id)
            self.species_to_index[sp_id] = new_index

            number_added += 1

            """
            # Expand the system state array for concentrations. TODO: do it for all the newly-added species at once
            if self.system is None:
                self.system = np.array([0], dtype='d')      # float64      TODO: allow users to specify the type
            else:
                self.system = np.pad(self.system, (0, 1))
            """

        return number_added



    def clear_index(self):
        self.index_to_species = []
        self.species_to_index = {}
