import numpy as np
from life123.species_registry import SpeciesRegistry
from life123.species_index_map import SpeciesIndexMap




class UniformCompartment:
    """
    Used to simulate the dynamics of reactions (in a single compartment.)
    This might be thought of as a "zero-dimensional system".

    This class is currently covering multiple roles
    that are expected to get separated in the future:

        1. reaction-system/model holder (utilizing class "ReactionRegistry")
        2. holder of state space
        3. scratch workspace for one reaction calculation -> BEING RELOCATED TO class "ReactionSimulator"
        4. time integrator (simulator for the kinetics of the reactions) -> BEING RELOCATED TO class "ReactionSimulator"

    Note that UniformCompartment is treated as a "reaction handler" by spatial modules
    that offer reaction-diffusion (BioSim1D, BioSim2D, BioSim3D)
    """

    def __init__(self, species_index_map=None, species_registry=None,
                 macromolecules=None, sync_func=None, species_list=None, names=None):
        """
        Note: AT MOST 1 of the following 2 arguments can be passed

        :param species_registry:[OPTIONAL 1] Object of type "SpeciesRegistry" (with data about the species).
                                If passed, cannot pass either of the args `reactions` (an object that contains `species_data`)
                                nor `names` (names are contained in the `species_data` object)

        :param names:           [OPTIONAL 2] A single name, or list or tuple of names, of the species.
                                (note: providing a list allows one to make the species appear in a particular desired order.)
                                If passed, cannot pass `species_data` (both those object contain the chemical names)

        :param macromolecules:      [OPTIONAL] Object of class "Macromolecule"
        :param species_index_map:   [OPTIONAL] Object of type "IndexSpecies"
        """
        self.sync_func = sync_func                  # Optional external function to invoke for "lazy synchronization"

        self.species_registry = species_registry    # Object of type "SpeciesRegistry" (with data about all the species)

        self.species_index_map = species_index_map  # Indexes to reconcile the species id's
                                                    # to their index position in the system state array

        assert (species_index_map is not None) or (species_registry), \
            "UniformCompartment() instantiation: at least one " \
            "of the arguments `species_index_map` or `species_registry` must be passed"

        if self.species_registry is None:
            # Use the species registry in the object species_index_map
            self.species_registry = self.species_index_map.species_registry
        else:
            # If both `species_registry` and `species_index_map` were passed, check for compatibility
            if self.species_index_map and self.species_index_map.species_registry:
                assert species_registry == self.species_index_map.species_registry,\
                    "incompatibility between objects passed by arguments `species_registry` and `species_index_map`"

        if self.species_index_map is None:
            self.species_index_map = SpeciesIndexMap(species_registry=self.species_registry)
            # Index all the species found in the species registry
            all_species = self.species_registry.get_all_species_ids()
            self.species_index_map.add_species(all_species)

        if self.species_index_map is None:

            system_state_size = 0
        else:
            system_state_size = self.species_index_map.number_of_mapped_species()

        self.system_state :np.ndarray = np.zeros(system_state_size, dtype='d')   # float64    CONCENTRATION SYSTEM STATE
                            # TODO: allow users to specify the type
                            # Concentration data in the single compartment we're simulating, for all the chemicals
                            # A 1-d Numpy array of the concentrations (floats), in their index order;
                            # the array size is the total number of chemical species.
                            # Each entry is the concentration of a species whose id is indexed by the IndexSpecies module
                            # Note that this is the counterpart - with 1 less dimension - of the array by the same name
                            #       in the class BioSim1D


        self.macromolecules = macromolecules    # Object of type "MacroMolecules"

        self.macro_system = {}      # A dict mapping macromolecule names to their counts
                                    # EXAMPLE:  {"M1": 1, "M2": 3, "M3": 1}

        self.macro_system_state = {}  # The counterpart of the system data for macro-molecules, if present
                                    # Binding fractions of the applicable transcription factors, for all the macro-molecules,
                                    # over the previous time step, indexed by macromolecule and by binding site number
                                    # EXAMPLE:   {"M1": {1: ("A", 0.2), 2: ("F", 0.93), ("P", 0.)},
                                    #             "M2": {1: ("C", 0.5), 2: ("F", 0.1)},
                                    #             "M3": {} }

                                    # For background, see: https://www.annualreviews.org/doi/10.1146/annurev-cellbio-100617-062719


        if macromolecules:
            ligands = macromolecules.get_ligands()     # Set of species id's
            self.add_species(ligands)





    #####################################################################################################

    '''                              ~   TO READ SYSTEM DATA   ~                                      '''

    def ________TO_READ_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def get_species_registry(self) -> SpeciesRegistry:
        """
        Return the "SpeciesRegistry" object being used

        :return:    Object of type "SpeciesRegistry"
        """
        return self.species_registry



    def number_of_system_species(self) -> int:
        """
        Number of species being simulated (and kept in the system state)
        :return:
        """
        number = len(self.system_state)
        assert number <= self.species_index_map.number_of_mapped_species(), \
            f"number_of_system_species(): the system size ({number}) " \
            f"exceeds the number of species mapped in the species index ({self.species_index_map.number_of_mapped_species()})"

        return number



    def get_system_conc(self) -> np.array:
        """
        Retrieve the concentrations of ALL the chemicals as a Numpy array

        :return:        A Numpy array with the concentrations of ALL the chemicals,
                        in their index order
                            EXAMPLE:  array([12.3, 4.56, 0.12])
                                      The 0-th chemical has concentration 12.3, and so on...
        """
        return self.system_state



    def get_species_conc(self, species_id :str) -> float:
        """
        Return the current system concentration of the given species, specified by its id.
        If no species with that id exists, an Exception is raised

        Synchronization, if registered, will take place

        :param species_id:  The id of a species
        :return:            The current system concentration of the above species
        """
        try:
            species_index = self.species_index_map.index_of(species_id, enforce=True)
        except IndexError:
            if self.sync_func:
                self.sync_func()

            # Try again, post-synchronization
            species_index = self.species_index_map.index_of(species_id, enforce=True)

        return self.system_state[species_index]



    def get_conc_dict(self, species_ids=None) -> dict:
        """
        Retrieve the concentrations of the requested species (by default all),
        as a dictionary indexed by the species id

        Synchronization, if registered, will take place

        :param species_ids: [OPTIONAL] List, tuple or set of the id's of the species of interest;
                                by default, return all

        :return:            A dictionary, possibly empty, indexed by the species id's, of the concentration values;
                                EXAMPLE: {"A": 1.2, "D": 4.67}
        """
        if species_ids is None:
            return {self.species_index_map.species_at(index): self.system_state[index]
                    for index, conc in enumerate(self.system_state)}
        else:
            assert isinstance(species_ids, (list, tuple, set)), \
                f"UniformCompartment.get_conc_dict(): the argument `species` must be a list, tuple or set" \
                f" (it was of type {type(species_ids)})"

            conc_dict = {}
            for sp_id in species_ids:
                # This will also take care of synchronization, if needed
                conc_dict[sp_id] = self.get_species_conc(species_id=sp_id)


            return conc_dict







    #####################################################################################################

    '''                             ~   TO SET SYSTEM DATA   ~                                   '''

    def ________TO_SET_DATA________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def add_species(self, species_ids : str | set[str] | list[str] | tuple[str, ...]) -> int:
        """
        Add the requested species, or group of species, to the system state,
        with values initialized to zero.
        The existing part of the system state won't be affected.

        :param species_ids: String, set, list or tuple of the ID's of species whose indexing needs to be managed;
                                no harm in also including species that are already managed (they will be ignored).
                                If a set is passed, it will be sorted in alphabetical order;
                                if you want a different specific order, pass a list or tuple in the order you want it
                                (the order has only UX significance, for displays, etc)

        :return:                The number of newly-managed species
        """
        number_added = self.species_index_map.add_species(species_ids)
        #print("number_added: ", number_added)

        if number_added > 0:
            # Expand the system state array for concentrations
            self.system_state = np.pad(self.system_state, (0, number_added))

        return number_added



    def set_system_state(self, conc : np.ndarray | list | tuple) -> None:
        """
        Completely set the entire system state, from the given concentration array or list/tuple.
        The array indexes must be present in `self.species_index_map` (i.e. have an assigned
        species correspondence).
        IMPORTANT: make sure that the array is in the right order, according to the current species index map.
                   A safer function to use is set_conc(), where species id's are used, instead of their array indexes

        :param conc:  Numpy array or list/tuple of concentration values
        :return:            None
        """
        assert len(conc) <= self.species_index_map.number_of_mapped_species(), \
            f"set_system_state(): the number of the passed concentration values ({len(conc)}) " \
            f"exceeds the size of the species index ({self.species_index_map.number_of_mapped_species()});\n" \
            f"try using set_conc() instead"

        if self.sync_func:
            self.sync_func()

        assert min(conc) >= 0, \
            f"UniformCompartment.set_system_state(): values meant to be species concentrations cannot be negative " \
            f"(such as the passed value {min(conc)})"

        if isinstance(conc, (list, tuple)):
            self.system_state = np.array(conc)
        else:

            self.system_state = conc.copy()    # Clone



    def set_single_conc(self, conc, species_id=None, species_index=None,
                        allow_growth=True, snapshot=None) -> None:
        """
        Set the concentrations of 1 species.

        Note: if both `species_id` and `species_index` and  are provided,an Exception is raised

        :param conc:            A non-negative number with the desired concentration value
                                    of the species specified below.  (Any previous value will get over-written)
        :param species_id:      [OPTIONAL] The id of the species of interest.
        :param species_index:   [OPTIONAL] An integer that indexes the species of interest (numbering starts at 0)
                                    If both `species_name` and `species_index` are provided, an Exception is raised.
                                    At least one of `species_name` or `species_index` must be specified

        :param allow_growth:    [OPTIONAL] Only applicable if `species_index` is specified.
                                    If True (default), allow the system state to expand if species_index exceeds the system size
        :param snapshot:        [OPTIONAL] A function that accepts an argument called "caption"
        :return:                None
        """
        # Validate the arguments
        assert conc >= 0, \
            f"UniformCompartment.set_single_conc(): species concentrations cannot be negative (value passed: {conc})"

        if species_index is not None:
            assert species_id is None, \
                "set_single_conc(): cannot pass both arguments `species_id` and `species_index`"

            assert species_index < self.species_index_map.number_of_mapped_species(), \
                f"set_single_conc(): the requested species index ({species_index}) " \
                f"exceeds the range of the species index (0-{self.species_index_map.number_of_mapped_species()-1});\n" \
                f"try using the argument `species_id` instead"


        if self.sync_func:
            self.sync_func()


        if species_id is not None:
            species_index = self.species_index_map.index_of(species_id)

        elif species_index is not None:
            # Validate the value
            self.species_index_map.species_at(species_index)    # This will raise an IndexError if invalid

        else:
            raise Exception("UniformCompartment.set_single_conc(): at least one "
                            "of the arguments `species_name` or `species_index` must be provided")


        # If no array position already exists
        # for the species whose concentration we're setting
        if species_index >= len(self.system_state):
            assert allow_growth, \
                f"UniformCompartment.set_single_conc(): the given index ({species_index}) cannot fit in the current " \
                f"system size ({len(self.system_state)}); use the argument `allow_growth=True` to automatically grow the system-state array"

            number_of_entries_to_add = species_index - len(self.system_state) + 1
            self.system_state = np.pad(self.system_state, (0, number_of_entries_to_add))


        self.system_state[species_index] = conc

        if snapshot is not None:
            # Invoke a passed external function
            snapshot(caption=f"Set concentration of `{self.species_index_map.species_at(species_index)}`")



    def set_conc(self, conc :dict, snapshot_funct=None) -> None:
        """
        Set the concentrations of a group of species

        :param conc:            A dict indexed by the species id, for the species of interest.
                                    Anything not specified will be left undisturbed
                                    EXAMPLE:  {"A": 12.4, "B": 0.23, "E": 2.6}

                                Note: any previous values will get over-written

        :param snapshot_funct:  [OPTIONAL] A function that accepts an argument called "caption"

        :return:                None
        """
        assert type(conc) is dict, \
            "set_conc(): argument `conc` must be a dict"

        for sp_id, conc_value in conc.items():
            self.set_single_conc(conc=conc_value, species_id=sp_id, snapshot=None)  # This function will also take care
                                                                                    # of optional synchronization

        if snapshot_funct is not None:
            # Invoke a passed external function
            snapshot_funct(caption="Set concentration")






    #####################################################################################################

    '''                                 ~  MACROMOLECULE STRUCTURE ~                                   '''

    def ________MACROMOLECULE_DYNAMICS________(DIVIDER):
        pass        # Used to get a better structure view in IDEs such asPycharm
    #####################################################################################################


    def set_macromolecules(self, data=None) -> None:
        """
        Specify the macromolecules, and their counts, to be included in the system.
        The fractional occupancy is set to 0, for all binding sites of all the specified macromolecules.
        Any previous data gets over-written.

        Note: to set a single fractional occupancy value, use set_occupancy()

        :param data:    A dict mapping macromolecule names to their counts
                            EXAMPLE:  {"M1": 1, "M2": 3, "M3": 1}
                        If any of the requested macromolecules isn't registered, an Exception will be raised
                        If data=None, then the set of registered macromolecules is used,
                            and all their counts are set to 1
        :return:        None.
                        The object variables self.macro_system and self.macro_system_state get set
        """
        if data is None:
            # Use the registered macromolecules, and set all counts to 1
            data = {}
            for mm in self.macromolecules.get_macromolecules():
                data[mm] = 1        # EXAMPLE, after this operation: data = {"M1": 1}

        self.macro_system = data


        self.macro_system_state = {}    # Reset
        for mm in data.keys():          # For each macromolecule in our system
            binding_sites_and_ligands = self.macromolecules.get_binding_sites_and_ligands(mm)     # EXAMPLE: {1: "A", 2: "C"}
            d = {}
            for (site_number, ligand) in binding_sites_and_ligands.items():
                d[site_number] = (ligand, 0.)           # All "binding occupancy fractions" are set to 0.

            self.macro_system_state[mm] = d



    def get_occupancy(self, macromolecule :str, site_number :int) -> float:
        """
        Get the fractional occupancy at the given binding site of the specified macromolecule.

        :param macromolecule:           Id of a species that represents a large molecule
        :param site_number:             Integer to identify a binding site on the macromolecule
        :return:                        A number between 0. and 1., representing the fractional occupancy
        """
        assert self.macro_system_state != {}, \
            "get_occupancy(): The system state for macromolecules has not been set yet;  " \
            "use set_macromolecules() or set_occupancy()"

        assert macromolecule in self.macro_system_state, \
            f"get_occupancy(): No occupancy data yet set for macromolecule `{macromolecule}`"

        assert site_number in self.macro_system_state[macromolecule], \
            f"get_occupancy(): No occupancy data yet set for site number {site_number} " \
            f"of macromolecule `{macromolecule}`"

        (ligand, fractional_occupancy) = self.macro_system_state[macromolecule][site_number]
        return fractional_occupancy



    def set_occupancy(self, macromolecule :str, site_number :int, fractional_occupancy :float) -> None:
        """
        Set the fractional occupancy at the given binding site of the specified macromolecule,
        using the requested value.
        If the specified macromolecule hasn't yet been added to the dynamical system state,
        automatically add it with count 1

        :param macromolecule:           Id of a species that represents a large molecule
        :param site_number:             Integer to identify a binding site on the macromolecule
        :param fractional_occupancy:    A number between 0. and 1., inclusive
        :return:                        None
        """
        assert 0. <= fractional_occupancy <= 1., \
            f"set_occupancy(): the value for the fractional occupancy must be a number between 0. and 1. " \
            f"Value given: {fractional_occupancy}"

        ligand = self.macromolecules.get_ligand_name(macromolecule=macromolecule, site_number=site_number)

        # If the specified macromolecule hasn't yet been added to the dynamic system, automatically add it
        # with count 1
        if self.macro_system_state == {}:
            self.set_macromolecules({macromolecule: 1})

        self.macro_system_state[macromolecule][site_number] = (ligand, fractional_occupancy)
