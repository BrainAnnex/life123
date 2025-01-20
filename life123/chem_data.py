from typing import Union, List, NamedTuple, Set

import numpy as np
import pandas as pd



class ChemCore:
    """
    Core data about the chemical species, such as their names, labels
    and indexes (position in their listing)

    Note: End users will typically utilize the class ChemData, which extends this one
    """

    def __init__(self):

        self.chemical_data = [] # Basic data for all chemicals, *except* water and macro-molecules.
                                # Each list entry represents 1 chemical,
                                # and is a dict required to contain the keys "label" and "name";
                                # it may also contain other arbitrary other keys ("notes" is a commonly-used one)
                                # EXAMPLE: [{"label": "A",   "name": "A"},
                                #           {"label": "NAD", "name": "Nicotinamide adenine dinucleotide", "note": "some note"},
                                #           {"label": "B",   "name": "B"]
                                # The position in this list is referred to as the "INDEX" of that chemical;
                                #          in out example above, "A" has index 0, while "NAD" has index 1
                                # Labels must be unique; likewise, names must be unique.
                                # The name of any chemical cannot be the same as the label of a different one.
                                # TODO: maybe use a Pandas dataframe


        self.label_dict = {}    # To map assigned chemical labels to their positional index
                                # (in the ordered list of chemicals, self.chemical_data);
                                # this is automatically set and maintained
                                # EXAMPLE: {"A": 0, "B": 1, "C": 2}

        self.color_dict = {}    # To map assigned chemical labels
                                # to optionally specified color values to use in visualizations
                                # EXAMPLE: {"A": "red", "B": "d07a19"}




    def number_of_chemicals(self) -> int:
        """
        Return the number of registered chemicals - exclusive of water and of macro-molecules

        :return:    The number of registered chemicals - exclusive of water and of macro-molecules
        """
        return len(self.chemical_data)



    def assert_valid_species_index(self, species_index: int) -> None:
        """
        Raise an Exception if the specified species_index (meant to identify a chemical) isn't valid

        :param species_index:   An integer that indexes the chemical of interest (numbering starts at 0)
        :return:                None
        """
        n_species = self.number_of_chemicals()
        assert type(species_index) == int,  f"The specified species index ({species_index}) must be an integer: " \
                                            f"instead, it has type {type(species_index)}"

        assert 0 <= species_index < n_species, \
            f"The specified species index ({species_index}) is not in the expected range [0 - {n_species - 1}], inclusive.  " \
            f"I.e., there is no chemical species assigned to this index"



    def get_label_mapping(self) -> dict:
        """
        Return a dict with all the mappings of the chemical names to the registered index

        :return:
        """
        return self.label_dict



    def get_index(self, label :str) -> int:
        """
        Return the index of the chemical species with the given label.
        Indexes are the integers assigned, in autoincrement order,
        at the time each chemical is first registered.
        If not found, an Exception is raised

        :param label:   Label of the chemical species of interest
        :return:        The index of the species with the given name
        """
        index = self.label_dict.get(label)
        assert index is not None, \
            f"ChemData.get_index(): No chemical species named `{label}` was found"

        return index



    def label_exists(self, label :str) -> bool:
        """
        Return True if the chemical with the given name exists (i.e. was registered),
        or False otherwise

        :param label:   The label of a chemical
        :return:        True if it exists, or False otherwise
        """
        return label in self.label_dict.keys()



    def get_label(self, species_index: int) -> str:
        """
        Return the label of the species with the given index.

        :param species_index:   An integer (starting with zero) corresponding to the
                                    original order with which the chemical species were first added
        :return:                The name of the species with the given index.
                                    If missing or blank, an Exception is raised
        """
        self.assert_valid_species_index(species_index)

        name = self.chemical_data[species_index].get("label")    # If "label" is not present, None will be returned

        assert name, \
            f"get_label(): A chemical species with the requested index ({species_index}) is present, but it lacks a name"

        return name



    def get_all_labels(self) -> [str]:
        """
        Return a list with the labels of all the chemical species,
        in their index order of registration.
        If any label is missing or blank, an Exception is raised

        :return:    A list of strings with the chemical names,
                        in their registered index order
        """
        all_labels = []
        for i, c in enumerate(self.chemical_data):
            label = c.get("label", None)
            assert label is not None, \
                f"get_all_labels(): missing or blank chemical name in index position {i}"

            all_labels.append(label)

        return all_labels



    def all_chemicals(self) -> pd.DataFrame:
        """
        Returns a Pandas dataframe with all the known information
        about all the registered chemicals (not counting macro-molecules),
        in their registration index order

        :return:    A Pandas dataframe
        """
        df = pd.DataFrame(self.chemical_data)

        # Add a column for plot colors, if any were registered
        if self.color_dict != {}:
            color_df = pd.DataFrame(list(self.color_dict.items()), columns=["label", "plot_color"])
            df = pd.merge(df, color_df, on="label", how="left") # All rows from df and matching rows from color_df
            df["plot_color"] = df["plot_color"].fillna("")      # Turn any NaN in the 'plot_color' column into a blank

        return df



    def add_chemical(self, name :str, label=None, note=None, plot_color=None, **kwargs) -> int:
        """
        Register a new chemical species, with a name
        and (optionally) :
            - a label (typically, a short version of the name, or a stand-in for it)
            - a note (will be stored in a "note" column)
            - a color value used in visualizations
            - any other named argument(s) that the user wishes to store (i.e. arbitrary-named arguments)

        EXAMPLE:  add_chemical("Nicotinamide adenine dinucleotide",
                               label = "NAD", note = "my note about this substrate", CAS_number = "CAS 53-84-9")

        Note: if also wanting to set the diffusion rate in a single function call,
              use ChemData.add_chemical_with_diffusion() instead

        :param name:        Name of the new chemical species to register; names must be unique -
                               an Exception will be raised if the name was already registered as a name or as a label
        :param label:       [OPTIONAL] Typically, a short version of the name, or a stand-in for it.
                                If provided, it must be unique, and cannot be identical to the name of another chemical;
                                if not provided, the name will be used as a label
        :param note:        [OPTIONAL] String with note to attach to the chemical
        :param plot_color:  [OPTIONAL] String with color value to attach to the chemical for visualizations
        :param kwargs:      [OPTIONAL] Dictionary of named arguments (with any desired names)
        :return:            The integer index assigned to the newly-added chemical
        """
        # Validate
        assert type(name) == str, \
            f"add_chemical(): a chemical's name must be provided, as a string value.  " \
            f"What was passed ({name}) was of type {type(name)}"

        assert name != "", \
            "add_chemical(): the chemical's name cannot be a blank string"

        assert name not in self.label_dict.keys(), \
            f"add_chemical(): the requested name (`{name}`) is ALREADY registered"

        for chem in self.chemical_data:
            assert name != chem["name"], \
                f"add_chemical(): the requested name (`{name}`) is ALREADY registered"


        if not label:
            label = name    # Use the (full) name as a label, if no full label is provided
        else:
            assert type(label) == str, \
                f"add_chemical(): a chemical's label, if provided, must be a string value.  " \
                f"What was passed ({label}) was of type {type(label)}"
            assert label not in self.label_dict.keys(), \
                f"add_chemical(): the requested name (`{label}`) is ALREADY registered"
            for chem in self.chemical_data:
                assert label != chem["name"], \
                    f"add_chemical(): the requested label (`{label}`) is ALREADY registered as the name of another chemical"

        index = len(self.chemical_data)     # The next available index number (list position)
        self.label_dict[label] = index

        # Prepare a dict with all the available data
        d = {"name": name, "label": label}

        if plot_color is not None:
            self.color_dict[label] = plot_color

        if note is not None:
            d["note"] = note

        d.update(kwargs)        # Merge dictionary kwargs into d

        self.chemical_data.append(d)

        return index



    def get_plot_color(self, label :str) -> Union[str, None]:
        """
        Return the name of the plot color previously associated to the given chemical

        :param label:   A string to identify a particular chemical
        :return:        The name of the associated plot color, or None if not found
        """
        return self.color_dict.get(label, None)



    def get_registered_colors(self, chem_labels) -> [str]:
        """
        Return a list of the colors registered for the specified chemicals

        :param chem_labels: List of chemical labels
        :return:            List of color names, with as many entries as the chemicals of interest;
                                any of the entries might be None
        """
        # Attempt to use the colors registered for individual chemicals, if present
        registered_colors = []
        for label in chem_labels:
            stored_color = self.get_plot_color(label)     # Will be None if no color was registered for this chemical
            registered_colors.append(stored_color)

        return registered_colors        # List of color names, with as many entries as the chemicals of interest;
                                        # any of the entries might be None



    def get_color_mapping_by_label(self) -> dict:
        """
        EXAMPLE: {"A": "red", "B": "orange", "C": "green"}

        :return:    A dict of plot colors, indexed by chemical labels
        """
        return self.color_dict


    def get_color_mapping_by_index(self) -> dict:
        """
        EXAMPLE: {0: 'red', 1: 'orange', 2: 'green'}

        :return:    A dict of plot colors, indexed by chemical index
        """
        return   {self.label_dict[k]: v for k, v in self.color_dict.items()}


    def get_all_colors(self) -> list:
        """

        :return:    A list of all colors, in the index position of their respective chemicals
        """
        return list(self.color_dict.values())






###############################################################################################################
###############################################################################################################


class Diffusion(ChemCore):
    """
    Extends its parent class, to manage diffusion-related data

    End users will typically utilize the class ChemData, which extends this one
    """

    def __init__(self):

        super().__init__()          # Invoke the constructor of its parent class

        self.diffusion_rates = {}   # Values for the diffusion rates, indexed by chemical "label".
                                    # All values must be non-negative numbers.
                                    # Only chemicals with an assigned diffusion rate will be present here.
                                    # EXAMPLE: {"A": 6.4, "B": 12.0}



    def add_chemical_with_diffusion(self, name :str, diff_rate :Union[float, int], label=None, note=None, **kwargs) -> int:
        """
        Register a new chemical species, with a name, a diffusion rate (in water),
        and (optionally) :
            - a note
            - any other named argument(s) that the user wishes to store (i.e. arbitrary named arguments)

        EXAMPLE:  add_chemical(label = "P1", diff_rate = 0.1,
                               note = "my note about P1", name = "protein P1")

        Note: if no diffusion is to be set, can also simply use ChemData.add_chemical()

        :param name:       Label of the new chemical species to register;
                                an Exception will be raised if the name was already registered
        :param diff_rate:   Floating-point number with the diffusion rate coefficient (in water) of this chemical
        :param label:       [OPTIONAL] Typically, a short version of the name, or a stand-in for it;
                                if not provided, the name will be used as a label
        :param note:        [OPTIONAL] String with note to attach to the chemical
        :param kwargs:      [OPTIONAL] Dictionary of named arguments (with any desired names)
        :return:            The integer index assigned to the newly-added chemical
        """
        # Validate the diffusion_rate, if present; other arguments get validated by add_chemical()
        self.assert_valid_diffusion(diff_rate)

        index = self.add_chemical(name=name, label=label, note=note, **kwargs)

        if label is None:
            label = name

        self.set_diffusion_rate(label=label, diff_rate=diff_rate)

        return index



    def set_diffusion_rate(self, label :str, diff_rate :Union[float, int]) -> None:
        """
        Set the diffusion rate of the given chemical species (identified by its name)

        :param label:       Label of a chemical species
        :param diff_rate:   Diffusion rate (in water) for the above chemical
        :return:            None
        """
        self.assert_valid_diffusion(diff_rate)
        self.diffusion_rates[label] = diff_rate



    def assert_valid_diffusion(self, diff) -> None:
        """
        Raise an Exception if the specified diffusion value isn't valid.
        Valid values are non-negative numbers (integer, float or numpy integers/floats)

        :param diff:    Diffusion rate
        :return:        None
        """
        assert isinstance(diff, (float, int, np.integer, np.floating)), \
            f"assert_valid_diffusion(): The value for the diffusion rate ({diff}) is not a valid number; it is of type {type(diff)}"

        assert diff >= 0., \
            f"assert_valid_diffusion(): the diffusion rate ({diff}) cannot be negative"



    def get_diffusion_rate(self, species_index=None, name=None) -> Union[float, int, None]:
        """
        Return the diffusion rate of the specified chemical species.
        If no value was assigned (but the chemical exists), return None.

        :param name:            Name of the chemical of interest
        :param species_index:   Alternate way to specify the chemical, using its zero-based index (order
                                    in which it was registered);
                                    `name` and `species_index` cannot be both specified, or an Exception will be raised
        :return:                The value of the diffusion rate for the species with the given index if present,
                                    or None if not
        """
        assert (name is None) or (species_index is None), \
            "get_diffusion_rate(): cannot specify BOTH `name` and `species_index`"

        if name is None:
            name = self.get_label(species_index)
        else:
            assert self.label_exists(name), \
                f"get_diffusion_rate(): No chemical named `{name}` exists"

        return self.diffusion_rates.get(name)      # If not present, None is returned



    def get_all_diffusion_rates(self) -> list:
        """
        Return a list of the diffusion rates of all the chemicals,
        in the order of their indexes.

        If any value is missing, None is used for it

        :return:    A list of numbers (or None values) with the diffusion rates
        """
        return [self.diffusion_rates.get(name) for name in self.get_all_labels()]
        # If any value is not present, None is used for it



    def missing_diffusion_rate(self) -> bool:
        """
        Determine whether any of the registered chemicals has a missing diffusion rates

        :return:    True if any of the diffusion rates (for all the registered chemicals) is missing;
                        False otherwise
        """
        if len(self.diffusion_rates) < self.number_of_chemicals():
            return True

        for name, diff in self.diffusion_rates.items():
            if diff is None:
                return True     # TODO: this should never occur

        return False



###############################################################################################################
###############################################################################################################


class ChemicalAffinity(NamedTuple):
    """
    Used for binding of ligands to macromolecules (e.g. Transcription Factors to DNA)
    """
    chemical: str   # Name of ligand
    Kd: float       # Dissociation constant; inversely related to binding affinity
                    # Note: dissociation constants are for now assumed to be constant,
                    #       regardless of what other (nearby) sites are occupied by ligands



class Macromolecules(Diffusion):
    """
    Extends its parent class to manage modeling of large molecules (such as DNA)
    with multiple binding sites (for example, for Transcription Factors)

    End users will typically utilize the class ChemData, which extends this one
    """

    def __init__(self):

        super().__init__()          # Invoke the constructor of its parent class


        self.macromolecules = []    # List of names.  EXAMPLE: ["M1", "M2"]
                                    # The position in the list is referred to as the
                                    #   "index" of that macro-molecule
                                    # Names are enforced to be unique

        self.binding_sites = {}     # A dict whose keys are macromolecule names.
                                    # The values are in turn dicts, indexed by binding-site number.
        # EXAMPLE:
        #       {"M1": {1: ChemicalAffinity("A", 2.4), 2: ChemicalAffinity("C", 5.1)},
        #        "M2": {1: ChemicalAffinity("C", 9.1), 2: ChemicalAffinity("B", 0.3),
        #               3: ChemicalAffinity("A", 1.8), 4: ChemicalAffinity("C", 2.3)}
        #        }
        #       where "M1", "M2" are macro-molecules,
        #           and "A", "B", "C" are bulk chemicals (such as transcription factors),
        #           all previously-declared;
        #           the various ChemicalAffinity's are NamedTuples (objects)
        #           storing a ligand name and its dissociation constant at that site.

        # TODO: maybe make a new class for a SINGLE macromolecule (akin to what done for reactions)

        # Info on Binding Site Affinities : https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6787930/




    def add_macromolecules(self, names: Union[str, List[str]]) -> None:
        """
        Register one or more macromolecule species, specified by their name(s)
        Note: this is a register of names, NOT of dynamical information
              about counts of macromolecules in the system (which is the domain of the class UniformCompartment)

        :param names:   A string, or list of strings, with the name(s) of the macromolecule(s)
        :return:        None.  The object attribute self.macro_molecules will get modified
        """
        if type(names) == str:  # If a single name was given, turn it into a list
            names = [names]

        for m in names:
            if m in self.macromolecules:
                # Warn of redundant attempt to re-register an existing macromolecule name
                print(f"WARNING: Macromolecule `{m}` was already registered.  Skipped...")
            else:
                self.macromolecules.append(m)  # Grow the list of registered macromolecules



    def get_macromolecules(self) -> [str]:
        """
        Return a list of the names of all the registered macromolecules

        :return:    A (possibly empty) list of the names of all the registered macromolecules
        """
        return self.macromolecules



    def set_binding_site_affinity(self, macromolecule: str, site_number: int, ligand: str, Kd) -> None:
        """
        Set the values of the binding affinity of the given macromolecule, at the indicated site on it,
        for the specified chemical species.

        Any previously-set value (for that macromolecule/site_number/chemical) will get over-written.

        IMPORTANT: Only 1 chemical (ligand type) can be associated to a particular site of a given macromolecule; attempting
        to associate another one will result in error.  In case multiple ligands can bind to the same physical site on
        the macromolecule, simply assign multiple site numbers for each of them

        NOTE: at present, no allowance is made for the fact that if the macromolecule is already bound
              to chemical "X" then its affinity of "Y" might be different than in the absence of "X"

        :param macromolecule:   Name of a macromolecule; if not previously-declared,
                                    it will get added to the list of registered macromolecules
        :param site_number:     Unique integer to identify a binding site on the macromolecule
        :param ligand:          Name of a previously-declared (bulk) chemical;
                                    if not found, an Exception will be raised       TODO: inconsistent with "macromolecule" arg
        :param Kd:              Dissociation constant, in units of concentration (typically microMolar).
                                    Note that the dissociation constant is inversely proportional to the binding affinity
        :return:                None
        """
        assert ligand in self.get_all_labels(), \
            f"set_binding_site_affinity(): no chemical named `{ligand}` found; use add_chemical() first"

        assert type(site_number) == int, \
            f"set_binding_site_affinity(): the argument `site_number` must be an integer"

        if macromolecule not in self.macromolecules:
            self.add_macromolecules([macromolecule])

        if self.binding_sites.get(macromolecule) is None:
            self.binding_sites[macromolecule] = {}          # Initialize the dict of binding sites for this macromolecule

        binding_data = self.binding_sites[macromolecule]    # This will be a dict whose key is the site_number

        if site_number in binding_data:
            existing_affinity_data = binding_data[site_number]
            if existing_affinity_data.chemical != ligand:
                raise Exception(f"set_binding_site_affinity(): "
                                f"site number {site_number} of macromolecule `{macromolecule}` was previously associated to chemical `{existing_affinity_data.chemical}` "
                                f"(attempting to set an affinity value for chemical `{ligand}`)")

        binding_data[site_number] = ChemicalAffinity(chemical=ligand, Kd=Kd)



    def get_binding_site_affinity(self, macromolecule: str, site_number: int) -> ChemicalAffinity:
        """
        Return the value of the binding affinity of the given macromolecule.
        If no value was previously set, an Exception is raised

        :param macromolecule:   Name of a macromolecule; if not found, an Exception will get raised
        :param site_number:     Integer to identify a binding site on the macromolecule
        :return:                The NamedTuple (ligand name, dissociation constant)
                                    if no value was previously set, an Exception is raised
        """
        assert macromolecule in self.macromolecules, \
            f"get_binding_site_affinity(): no macromolecule named `{macromolecule}` found"

        assert type(site_number) == int, \
            f"get_binding_site_affinity(): the argument `site_number` must be an integer"

        binding_data = self.binding_sites.get(macromolecule)    # This will be a dict whose key is the site_number.  (None if not found)
        assert binding_data is not None, "not found"

        chem_affinity = binding_data.get(site_number)           # Will be None if not found
        assert chem_affinity is not None, "not found"

        return chem_affinity



    def get_binding_sites(self, macromolecule) -> [int]:
        """
        Get a list of all the binding-site numbers of the given macromolecule.
        If the requested macromolecule isn't registered, an Exception will be raised

        :param macromolecule:   The name of a macromolecule
        :return:                A (possibly empty) list of integers, representing the binding-site numbers.
                                    EXAMPLE: [1, 2]
        """
        binding_site_info = self.binding_sites.get(macromolecule)   # EXAMPLE: {1: ChemicalAffinity("A", 2.4), 2: ChemicalAffinity("C", 5.1)}
        #   will be None if no binding sites are found
        if binding_site_info is None:
            assert macromolecule in self.get_macromolecules(), \
                f"get_binding_sites(): the requested macromolecule ({macromolecule}) isn't registered; use add_macromolecules()"
            return []
        return list(binding_site_info)                              # EXAMPLE: [1, 2]



    def get_binding_sites_and_ligands(self, macromolecule) -> dict:
        """
        Return a mapping (python dict) from binding-site number to ligand species, for the given macromolecule
        If the requested macromolecule isn't registered, an Exception will be raised

        :param macromolecule:   The name of a macromolecule
        :return:                A dict whose keys are binding-site numbers and values are their respective ligands
                                    EXAMPLE: {1: "A", 2: "C"}
        """
        binding_site_info = self.binding_sites.get(macromolecule)   # EXAMPLE: {1: ChemicalAffinity("A", 2.4), 2: ChemicalAffinity("C", 5.1)}
                                                                    #   It will be None if no binding sites are found
        if binding_site_info is None:
            assert macromolecule in self.get_macromolecules(), \
                f"get_binding_sites(): the requested macromolecule ({macromolecule}) isn't registered; use add_macromolecules()"
            return {}

        d = {}
        for (site_number, affinity_obj) in binding_site_info.items():
            ligand = affinity_obj.chemical
            d[site_number] = ligand

        return d        # EXAMPLE: {1: "A", 2: "C"}



    def get_ligand_name(self, macromolecule: str, site_number: int) -> str:
        """
        Return the name of the ligand associated to the specified site
        on the given macromolecule.
        If not found, an Exception is raised

        :param macromolecule:   The name of a macromolecule
        :param site_number:     Integer to identify a binding site on the macromolecule
        :return:                The name of the ligand (chemical species)
        """
        binding_site_info = self.binding_sites.get(macromolecule)   # EXAMPLE: {1: ChemicalAffinity("A", 2.4), 2: ChemicalAffinity("C", 5.1)}
                                                                    #   It will be None if no binding sites are found
        if binding_site_info is None:
            assert macromolecule in self.get_macromolecules(), \
                f"get_ligand_name(): the requested macromolecule (`{macromolecule}`) isn't registered; use add_macromolecules()"

            raise Exception(f"get_ligand_name(): no binding sites are defined on macromolecule {macromolecule}")

        ligand_data = binding_site_info.get(site_number)        # EXAMPLE:  ChemicalAffinity("A", 2.4)
                                                                # It will be None if this binding site isn't present
        assert ligand_data is not None, \
            f"get_ligand_name(): binding site {site_number} is not found on macromolecule {macromolecule}"

        return  ligand_data.chemical



    def show_binding_affinities(self) -> None:
        """
        Print out the Dissociation Constant for each Binding Site in each Macromolecule

        :return:    None
        """
        for mm in self.get_macromolecules():
            #binding_sites_and_ligands = self.get_binding_sites_and_ligands(mm)
            print(mm, " :")
            for site_number in self.get_binding_sites(mm):
                aff = self.get_binding_site_affinity(macromolecule=mm, site_number=site_number)
                print(f"   Site {site_number} - Kd (dissociation const) for {aff.chemical} : {aff.Kd}")



    def reset_macromolecule(self, macromolecule) -> None:
        """
        Erase all data for the specified macromolecule

        :param macromolecule:   The name of a macromolecule
        :return:                None
        """
        if self.binding_sites.get(macromolecule) is not None:
            del self.binding_sites[macromolecule]



    def clear_macromolecules(self) -> None:
        """
        Reset all macromolecules to their original state

        :return:    None
        """
        self.macromolecules = []
        self.binding_sites = {}




###############################################################################################################
###############################################################################################################


class ChemData(Macromolecules):
    """
    Data about all the chemicals and (if applicable) reactions,
    including:
        - names
        - diffusion rates
        - macro-molecules Binding Site Affinities (for Transcription Factors)


    Notes:  - for now, the temperature is assumed constant everywhere, and unvarying (or very slowly varying)

            - we're using a "daisy chain" of classes extending the previous one, starting from ChemCore
              and ending in this user-facing class:
                    ChemCore <- Diffusion <- Macromolecules <- ChemData
    """
    def __init__(self, names=None, labels=None, diffusion_rates=None, plot_colors=None):
        """
        Any non-None arguments MUST all have the same length (and appear in the same order), or all be scalars.
        It's ok to skip passing any data at instantiation, and later add it, with calls to add_chemical().
        Reactions, if applicable, need to be added later by means of calls to add_reaction()
        Macro-molecules, if applicable, need to be added later.

        If no names nor labels are provided, but diffusion rate or plot colors are given,
        the strings "Chemical 1", "Chemical 2", ..., are used

        :param names:           [OPTIONAL] A single name, or list or tuple of names, of the chemicals.
                                    If not provided, the names are made equal to the labels.
                                    If neither names nor labels is provided, "Chemical 1", "Chemical 2", ..., are used
                                    (as many as the diffusion rates or plot colors)

        :param labels:          [OPTIONAL] A single label, or list or tuple of labels, of the chemicals,
                                    in the same order as the names (if provided).
                                    If not provided, the labels are made equal to the names.
                                    If neither names nor labels is provided, "Chemical 1", "Chemical 2", ..., are used
                                    (as many as the diffusion rates or plot colors)

        :param diffusion_rates: [OPTIONAL] A non-negative number, or a list/tuple/Numpy array with the diffusion rates of the chemicals,
                                    in the same order as the names/labels (if provided).

        :param plot_colors:     [OPTIONAL] A single name, or list or tuple of names, of the plotting colors for the chemicals,
                                    in the same order as the names/labels (if provided).
        """
        # TODO: allow a way to optionally pass macromolecules as well

        super().__init__()       # Invoke the constructor of its parent class


        non_none_args = [arg for arg in (names, labels, diffusion_rates, plot_colors) if arg is not None]

        if len(non_none_args) == 0:
            return      # All args are None; nothing else to do


        if names is not None:
            if (dt := type(names)) == str:
                names = [names]
            else:
                assert dt == list or dt == tuple, \
                    f"ChemData instantiation: the `names` argument must be a list or tuple, or a string for a single chemical.  " \
                    f"What was passed ({names}) was of type {dt}"

        if labels is not None:
            if (dt := type(labels)) == str:
                labels = [labels]
            else:
                assert dt == list or dt == tuple, \
                    f"ChemData instantiation: the `labels` argument must be a list or tuple, or a string for a single chemical.  " \
                    f"What was passed ({labels}) was of type {dt}"

        if plot_colors is not None:
            if (dt := type(plot_colors)) == str:
                plot_colors = [plot_colors]
            else:
                assert dt == list or dt == tuple, \
                    f"ChemData instantiation: the `plot_colors` argument must be a list or tuple, or a string for a single chemical.  " \
                    f"What was passed ({plot_colors}) was of type {dt}"

        if diffusion_rates is not None:
            if (dt := type(diffusion_rates)) == int or dt == float:
                diffusion_rates = [diffusion_rates]
            else:
                assert dt == list or dt == tuple  or dt == np.ndarray, \
                    f"ChemData instantiation(): the `diffusion_rates` argument must be a list or tuple or Numpy array, or a number for a single chemical.  " \
                    f"What was passed ({diffusion_rates}) was of type {dt}"


        # We need to re-compute non_none_args because some of the original arguments may have been modified by now
        non_none_args = [arg for arg in (names, labels, diffusion_rates, plot_colors) if arg is not None]

        # Assert that all lengths are the same
        if len(non_none_args) > 1:  # No need to check if only one or no argument is present
            lengths = []
            for arg in non_none_args:
                    lengths.append(len(arg))

            assert all(l == lengths[0] for l in lengths), \
                "ChemData instantiation: all the non-None arguments must have the same length"


        # Populate the data structure
        if names is None and labels is None:
            if diffusion_rates is not None:
                n_species = len(diffusion_rates)
            else:
                n_species = len(plot_colors)        # plot_colors cannot be None, otherwise all args would be (already excluded)

            names = [f"Chemical {i+1}" for i in range(n_species)]   # The strings "Chemical 1", "Chemical 2", ..., will be used


        if labels is None:
            labels = names

        if names is None:
            names = labels

        for i, chem_name in enumerate(names):
            assert type(chem_name) == str, \
                f"ChemData instantiation: all the names must be strings.  The passed value ({chem_name}) is of type {type(chem_name)}"

            l = labels[i]
            assert type(l) == str, \
                f"ChemData instantiation: all the labels must be strings.  The passed value ({l}) is of type {type(l)}"

            if diffusion_rates is None:
                self.add_chemical(name=chem_name, label=l)
            else:
                diff = diffusion_rates[i]
                self.assert_valid_diffusion(diff)
                self.add_chemical(name=chem_name, label=l)
                self.set_diffusion_rate(label=l, diff_rate=diff)

            if plot_colors is not None:
                color = plot_colors[i]
                assert type(color) == str, \
                    f"ChemData instantiation: all the colors must be strings.  The passed value ({color}) is of type {type(color)}"
                self.color_dict[l] = color
