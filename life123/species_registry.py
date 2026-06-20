from dataclasses import dataclass, field, asdict
from itertools import islice
from life123.visualization.colors import Colors
from typing import NamedTuple
import string
import numpy as np
import pandas as pd



@dataclass(slots=True)    # (slots=True) has the effect of prohibiting non-listed fields, and of making the class more efficient
class Species:
    """
    EXAMPLE:
        Species(
            id="metK",
            categories=["protein", "enzyme"],
            molecular_weight=43000,
            ec_number="2.5.1.6",
            metadata={"compartment": "cytosol"}
        )
    """
    id: str
    name: str  = ""                                 # EXAMPLE: "methionine adenosyltransferase"
    label: str  = ""                                # Primarily for use in visualization;
                                                    #   typically, a short version of the name, or a stand-in for it
                                                    #   EXAMPLE: "metK"
    sort_order: float = 0                           # By default, zero-based order of registration
    diffusion_rate: float|None = None
    charge: int = 0
    formula: str = ""
    ec_number: str = ""                             # EXAMPLE: "2.5.1.6"
    cas_number: str = ""                            # EXAMPLE: "53-84-9"
                                                    #   See: https://commonchemistry.cas.org
    chebi: str = ""
    locus_tag: str = ""                             # EXAMPLE: "0432"
    molecular_weight: float|None = None
    categories: list = field(default_factory=list)  # EXAMPLE: ["protein", "enzyme"]
    metadata: dict = field(default_factory=dict)    # A dict of any extra fields. EXAMPLE: {"compartment": "cytosol"}
    annotation: str = ""                            # Freeform notes
    plot_color: str = ""                            # For available names, see colors.py
                                                    #   EXAMPLES: "darkorange" , "d07a19"
    """
    EXAMPLES of categories:   
        small_molecule
        ion
        enzyme
        protein
        protein_complex
        dna
        rna
        lipid
        carbohydrate
        metabolite
        drug
        macromolecular_assembly
        nucleotide
        purine
    """

    def __setattr__(self, name :str, value):
        """
        Note: here we do validation, and set a few defaults,
              both during instantiation as well as when modifying an existing Species

        :param name:
        :param value:
        :return:
        """
        if (name == "label") or (name == "name"):
            # Special treatment for the properties (fields) "label" and "name"
            if (value is None):
                value = self.id             # Use the `id` value if missing
            else:
                assert type(value)==str, \
                    f"The field '{name}', if specified, must be a string; the value given ({value}) is of type {type(value)}"
                if value.strip() == "":     # Detect if blank
                    value = self.id         # Use the `id` value if missing

        else:
            self._validate(field_name=name, value=value)


        # Carry out the actual setting of the property
        object.__setattr__(self, name, value)



    def _validate(self, field_name :str, value) -> None:
        """
        Note: "label" and "name" are handled separately; not here

        :param field_name:
        :param value:
        :return:            None
        """
        # TODO: no longer needed by SpeciesRegistry class; can now be moved into Species class
        # TODO: add more validation

        match field_name:
            case "diffusion_rate":
                # Acceptable values: None, or a non-negative number (possibly Numpy)
                if (value is not None):
                    assert isinstance(value, (float, int, np.integer, np.floating)), \
                        f"The value for `{field_name}` ({value}) is not a valid number; the value given is of type {type(value)}"

                    assert value >= 0, \
                        f"`{field_name}` cannot be negative; the value given: {value}"


            case "categories":
                # Acceptable values: lists
                assert type(value) == list, \
                    f"`{field_name}` must be a list; the value given ({value}) is of type {type(value)}"





############################################################################################

class SpeciesRegistry:
    """

    """
    def __init__(self, id=None, species=None, n_species=None,
                 diffusion_rates=None, plot_colors=None, **kwargs):
        """

        Note: all arguments are optional.  AT MOST, 1 of the following 3 may be passed

        :param id:          [OPTIONAL] A string, or list or tuple of them.  The same value(s) will be used as id, name and label

        :param species:     [OPTIONAL] An object of type "Species", or a list or tuple of them.
                                Their `sort_order` attributes will be modified, to be a zero-base auto-increment: 0, 1, 2, ...

        :param n_species:   [OPTIONAL] The desired number of species.  Their id's, names and labels are automatically
                                assigned as: "A", "B", ..., "Z", "Z2", "Z3", ...

        :param diffusion_rates: [OPTIONAL] OBSOLETE
        :param plot_colors:     [OPTIONAL] OBSOLETE

        :param **kwargs:        [OPTIONAL] Other named arguments.
                                    Their names must match the property names of the class "Species".
                                    If present, one of the 3 core arguments (`id` or `species` or `n_species`) must also be present
                                    EXAMPLE: diffusion_rate=[5, 2], plot_color=["blue", "red"]
        """
        # Raise Exception if any of the obsolete arguments are being used
        assert diffusion_rates is None, \
            "SpeciesRegistry() instantiation: argument `diffusion_rates` is obsolete; use the singular `diffusion_rate` instead"

        assert plot_colors is None, \
            "SpeciesRegistry() instantiation: argument `plot_colors` is obsolete; use the singular `plot_color` instead"


        self.count = 0                      # Species registration index

        self.by_id: dict[str, Species] = {} # Species classes indexed by `id`.  dict[str, Species]   EXAMPLE:
        """
        {
            "A":    Species(id="A", name="my chem", label="A", "notes="example", ...
                           ),
            "metK": Species(id="metK", name="methionine adenosyltransferase", label="metK", "note="a lot of metadata!",
                            ec_number="2.5.1.6", locus_tag="0432", molecular_weight: 43000,
                            categories=["protein", "enzyme"], ...                               
                            )
        }
        """

        # TODO: maybe the 2 indexes below really belong to the simulator
        #self.id_to_index: dict[str, int] = {}       # EXAMPLE: {"Species A": 0, "Species B": 1}
        #self.index_to_species: list[Species] = []   # EXAMPLE: ["Species A", "Species B"]


        # Determine which arguments, if any, aren't None
        passed_arg_values = [arg
                             for arg in (id, species, n_species) if arg is not None]

        #if len(passed_arg_values) == 0:
        #    return      # All the arguments are None; nothing else to do

        assert len(passed_arg_values) <= 1, \
            f"SpeciesRegistry() instantiation: cannot pass more than 1 of the arguments `id`, `species`, `n_species`.  " \
            f"{len(passed_arg_values)} arguments were passed"


        # `id` argument
        if id is not None:
            if type(id) == str:
                id = [id]
            else:
                assert (type(id) == list) or (type(id) == tuple), \
                    f"SpeciesRegistry() instantiation: the `id` argument, if provided, must be a string, list or tuple"

            for id in id:
                self.add_species(id=id)


        # `species` argument
        if species is not None:
            if type(species) == Species:
                species = [species]
            else:
                assert type(species) == list or type(species) == tuple, \
                    f"SpeciesRegistry() instantiation: the `species` argument, if provided, must be a Species object, or list/tuple of them"

            for species in species:
                assert type(species) == Species, \
                    "SpeciesRegistry() instantiation: the `species` argument must be a list or tuple of `Species` objects"

                species.sort_order = self.count
                self.by_id[species.id] = species
                self.count += 1


        # `n_species` argument
        if n_species is not None:
            assert (type(n_species) == int) and (n_species > 0), \
                f"SpeciesRegistry() instantiation: the `n_species` argument, if provided, must be a positive integer"

            id = self._generate_generic_names(n_species)     # Generates the strings "A", "B", ..., "Z", "Z2", "Z3", ...

            for id in id:
                self.add_species(id=id)


        # Process any other named arguments that might be present
        #print('kwargs: ', kwargs)           # EXAMPLE: {'diffusion_rate': [123], 'plot_color': ['red']}

        number_species = self.number_of_species()
        if number_species == 0:
            assert len(kwargs) == 0, \
                f"SpeciesRegistry(): you cannot pass special arguments {kwargs} , unless you also pass one of " \
                f"the 3 core arguments (`id` or `species` or `n_species`)"

        for k, v in kwargs.items():
            #print(f"key: {k} -> value: {v}")    # EXAMPLE: key: plot_color -> value: ('blue', 'turquoise')

            if type(v) != list and type(v) != tuple:
                v = [v]         # Turn single values into lists

            assert len(v) == number_species, \
                f"SpeciesRegistry() instantiation: the number of elements in the special argument `{k}` ({len(v)}) " \
                f"doesn't match the requested number of species ({number_species})"

            self.set_all_values(field_name=k, values=v)



    def __repr__(self):
        msg = """This is on object of type life123.species_registry.SpeciesRegistry
                 To display is, use the function as_recordset() or as_dataframe()
                 EXAMPLE:   my_species_registry.as_dataframe()
              """
        return  msg






    #####################################################################################################

    '''                                 ~   GET VALUES  ~                                             '''

    def ________GET_VALUES________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def get_species(self, id :str) -> Species:
        """
        Retrieve, and return, the species with the given id

        :param id:  String with a unique value to identify a species
        :return:    An object of type "Species"
        """
        s = self.by_id.get(id)

        if s is None:   # If not found
            # Check whether the problem was a bad `id` format, for better error reporting
            assert type(id) == str, \
                f"add_species(): the argument `id` must be a string; it is of type {type(id)}"

            raise Exception(f"get_species(): no species with id `{id}` found")

        return s



    def get_species_id(self, species_index :int) -> str:
        """
        Return the id of the species with the given registration index.

        :param species_index:   An integer (starting with zero) corresponding to the
                                    original order with which the species were first registered
        :return:                The id of the chemical species with the given index;
                                    if missing, an Exception is raised
        """
        # TODO: this will become replacement for the old get_label()
        # TODO: it'd be better if the calling module kept a lookup table instead :  index -> species_id

        if len(self.by_id) == 0:
            raise Exception("get_species_id(): no lookup is possible, because there are no species currently registered")

        assert type(species_index) == int, \
            f"get_species_id(): argument `species_index` must be an integer; it was of type {type(species_index)}"

        assert (species_index >= 0) and (species_index < len(self.by_id)), \
            f"get_species_id(): argument `species_index` must be between 0 and {len(self.by_id) - 1} (inclusive).  " \
            f"The value passed was {species_index}"

        # Extract the species_index-th key efficiently
        return next(islice(iter(self.by_id), species_index, species_index + 1))



    def get_all_species_ids(self) -> list[str]:
        """
        Return a list with the id's of all the species,
        in their index order of registration.

        :return:    A list of strings with the species id's,
                        in their registered index order
        """
        # TODO: this is to replace the old get_all_labels()
        return list(self.by_id)     # The dictionary keys as a list



    def get_species_index(self, id :str) -> int:
        """
        Return the index of the species with the given id.
        Indexes are the integers assigned, in autoincrement order,
        at the time each species is first registered.
        If not found, an Exception is raised

        :param id:  String with a unique value to identify a species
        :return:    The integer index of the species with the given id
                        (the order in which it was added to the registry)
        """
        # TODO: replacement for the old get_index()
        # TODO: consider maintaining a lookup - or better yet doing without this feature!
        index = next((i for i, k in enumerate(self.by_id) if k == id), None)

        assert index is not None, \
            f"get_species_index(): No species with id `{id}` was found"

        return index



    def assert_valid_species_index(self, index :int) -> None:
        """
        Raise an Exception if the specified index (meant to identify a  species) isn't valid

        :param index:  An integer that indexes the species of interest (numbering starts at 0)
        :return:       None
        """
        #TODO: this replaces assert_valid_chem_index()
        #TODO: pytest
        n_species = self.number_of_species()
        assert type(index) == int, f"assert_valid_species_index(): The specified index ({index}) must be an integer: " \
                                   f"instead, it has type {type(index)}"

        assert 0 <= index < n_species, \
            f"assert_valid_species_index(): The specified index ({index}) is not in the expected range [0 - {n_species - 1}], inclusive.  " \
            f"I.e., there is no species assigned to this index"



    def get_value(self, species_id :str, field_name :str):
        """
        If the species isn't found, or if it lacks the requested field (property) name,
        an Exception is raised

        :param species_id:  String with a unique value to identify a species
        :param field_name:  The name of the property of interest.  EXAMPLE: "formula"
        :return:            The value of the above property, for the specified species.
                                (Note that None may be returned if that happens to be the property value)
        """
        s = self.get_species(species_id)
        try:
            value = getattr(s, field_name)
        except:
            raise Exception(f"get_value(): species `{species_id}` has no attribute '{field_name}'")

        return value



    def get_all_values(self, field_name :str):
        """
        Retrieve all the values of the specified field (property),
        across all species in order of their registration index

        :param field_name:  The name of the property of interest.  EXAMPLE: "formula"
        :return:            A list of all values.  Some may be None
        """
        all_values = []
        for key, val in self.by_id.items():
            value = self.get_value(species_id=key, field_name=field_name)
            all_values.append(value)

        return all_values



    def max_value(self, field_name :str):
        """

        :param field_name:  The name of the property of interest.  EXAMPLE: "diffusion_rate"
        :return:            The max of all values.
                                If the species registry is empty, or if any value is None,
                                an Exception will be raised
        """
        assert self.number_of_species() > 0, \
            "get_max_value(): no values found, because no species are present"

        try:
            max_value = max(self.get_all_values(field_name=field_name))
        except Exception as ex:
            if self.has_missing_values(field_name=field_name):      # For better error message
                raise Exception("get_max_value(): cannot determine a max, because some values are missing")
            else:
                raise Exception(f"get_max_value(): {ex}")       # Just pass thru the Exception message


        return max_value



    def has_missing_values(self, field_name :str) -> bool:
        """
        Determine whether any of the registered species has any missing value for the specified field (property)

        :return:            True if any of the value is missing;
                                False otherwise
        :param field_name:  The name of the property of interest.  EXAMPLE: "diffusion_rate"
        :return:
        """
        for s in self.by_id.values():
            value = getattr(s, field_name, None)
            if value is None:
                return True

        return False



    def number_of_species(self) -> int:
        """
        Return the number of registered species

        :return:    The number of registered species
        """
        return len(self.by_id)



    def as_recordset(self) -> list[dict]:
        """
        Return all the registry data in the form of a list of dictionaries,
        where each dictionary contains the data for a species

        :return: EXAMPLE:
                    [
                        {'id': 'A', 'name': 'A', 'label': 'A', .... }
                    ]
        """
        recordset = []
        # Loop over values in the dictionary self.by_id
        for v in self.by_id.values():
            d = asdict(v)
            recordset.append(d)

        return recordset



    def as_dataframe(self, sort=False) -> pd.DataFrame:
        """
        Returns a Pandas dataframe with all the known information
        about all the registered species,
        including id's, labels, names, diffusion rates, plot colors, etc.

        Missing values with appear as NaN in the dataframe

        :param sort:    [OPTIONAL] If True, the dataframe is sorted by its standard column 'sort_order'.
                            Default, False
        :return:        A Pandas dataframe with info about all the registered species
        """
        df = pd.DataFrame(self.as_recordset())   # Get labels, names and optional extra fields
        # Note: missing diffusion values are returned as None, and will be turned into NaN in the dataframe

        if sort:
            # Sort by 'sort_order' in place (ascending order)
            df.sort_values(by='sort_order', inplace=True, ignore_index=True)

        return df





    #####################################################################################################

    '''                                 ~   SET VALUES  ~                                             '''

    def ________SET_VALUES________(DIVIDER):
        pass        # Used to get a better structure view in IDEs
    #####################################################################################################


    def add_species(self, id :str, skip_duplicates=False, **kwargs) -> Species|None:
        """
        Register a new species, with an id and (optionally) any of the parameters listed for the class Species,
        such as name, label, annotation, diffusion_rate, ec_number, etc.

            - a color value used for this chemical species in visualizations
            - any other named argument(s) that the user wishes to store (i.e. arbitrary-named arguments)

        EXAMPLE:  add_species(id = "NAD", name = "Nicotinamide adenine dinucleotide",
                              annotation = "my note about this substrate",
                              cas_number = "53-84-9", plot_color="paleturquoise")

        :param id:              String with a unique value to identify this species being added.
                                    If `name` or `label` aren't passed as arguments, they will be
                                    made equal to the id
        :param skip_duplicates: [OPTIONAL] Normally, duplicate labels or duplicate names raise Exception;
                                    however, if True, any such duplicate is silently skipped
        :param kwargs:          [OPTIONAL] Any other argument to pass;
                                    for list and explanations, see under the class `Species`
        :return:                An object of class `Species`, representing the newly-registered species;
                                    however, if it already exists (and skip_duplicates=True), None is returned
        """
        assert type(id) == str, \
            f"add_species(): the argument `id` must be a string; it is of type {type(id)}"

        if id in self.by_id:
            # Duplicate is detected
            if skip_duplicates:
                return None
            else:
                raise Exception(f"add_species(): the requested species id (`{id}`) ALREADY exists")


        s = Species(id=id, sort_order=self.count, **kwargs)
        self.by_id[id] = s
        self.count += 1

        return s



    def set_value(self, species_id :str, field_name=None, value=None, **kwargs) -> None:
        """

        :param species_id:  String with a unique value to identify a species
        :param field_name:  The name of the property of interest.  EXAMPLE: "formula"
        :param value:       The value that we wish to set for the above property;
                                None is an acceptable value
        :param kwargs:      [OPTIONAL] Other named arguments.
                                Their names must match the property names of the class "Species".
                                If present, cannot have arguments `field_name` nor `value`
                                EXAMPLE: diffusion_rate=123.
        :return:            None
        """
        s = self.get_species(species_id)
        # TODO: validate for unique name and label, if applicable

        if field_name is not None:      # Setting by `field_name` and `value`
            # Note: value could be None
            assert kwargs == {}, \
                f"set_value(): Cannot specify both `field_name` and `{list(kwargs)[0]}`"
            setattr(s, field_name, value)       # In-place modification
            return

        # If we get here, `field_name` is None
        assert len(kwargs) > 0, "set_value(): Must also pass the argument `value`"

        # We're setting by special named arguments
        assert value is None, \
            f"set_value(): Cannot specify both `value` and `{list(kwargs)[0]}`"

        for k, v in kwargs.items():
            #print(f"key: {k} -> value: {v}")    # EXAMPLE: key: plot_color -> value: ('blue', 'turquoise')
            setattr(s, k, v)       # In-place modification



    def set_all_values(self, field_name :str, values :list) -> None:
        """
        Set values across ALL species, in order of their registration index,
        for the specified field (property)

        :param field_name:  The name of the property of interest.  EXAMPLE: "formula"
        :param values:      A list of values compatible with the above field
        :return:            None
        """
        assert type(values) == list or type(values) == tuple, \
            "set_values(): the argument `values` must be a list or tuple"

        assert len(values) == self.number_of_species(), \
            f"set_values(): the number of entries in argument `values` ({len(values)}) " \
            f"doesn't match the number of registered species ({self.number_of_species()})"

        index = 0
        for key, val in self.by_id.items():
            self.set_value(species_id=key, field_name=field_name, value=values[index])
            index += 1



    def assign_colors(self, species_ids=None) -> list[str]:
        """
        Return a list of the colors registered for the requested chemical species
        or, if not specified, for all chemicals.
        If any species lacks a registered color, assign new ones from a palette of default colors,
        and register the new assignments for future use.

        :param species_ids: List of species id's; use None to mean ALL the species
        :return:            List of color names, with as many entries as the species of interest,
                                and in their same order
                                EXAMPLE:   ["blue", "yellow", "cyan"]
        """
        if species_ids is None:
            species_ids = self.get_all_species_ids()
        elif type(species_ids) == str:
            species_ids = [species_ids]
        else:
            assert type(species_ids) == list, \
                f"assign_colors(): argument `species_ids` must be a string or a list of strings; " \
                f"instead, it was {type(species_ids)}"


        # Attempt to use the colors registered for the individual species;
        # if not already present, assign new ones
        registered_colors = []
        fallback_colors = None      # A list that will be created if needed
        new_color_index = 0
        for label in species_ids:
            stored_color = self.get_value(species_id=label, field_name="plot_color")     # Will be an empty string if no color was registered for this species
            if stored_color:
                registered_colors.append(stored_color)
            else:   # In case no previously-registered color for this chemical
                if fallback_colors is None:
                    fallback_colors = Colors.assign_default_colors(len(species_ids))

                new_color = fallback_colors[new_color_index]
                new_color_index += 1
                registered_colors.append(new_color)
                self.set_value(species_id=label, field_name="plot_color", value=new_color)


        return registered_colors        # List of color names, with as many entries as the species of interest



    def get_color_mapping_by_index(self) -> dict:
        """
        EXAMPLE: {0: 'red', 1: 'orange', 2: 'green'}

        :return:    A dict of plot colors, indexed by species index
        """
        return  {i: s.plot_color
                        for i, s in enumerate(self.by_id.values())
                }





    #############  PRIVATE  #############

    def _generate_generic_names(self, n :int) -> [str]:
        """
        Generate a list of n elements using the strings "A", "B", ..., "Z", "Z2", "Z3", ...
        (i.e., if there aren't enough letters in the English alphabet, proceed with "Z2", "Z3", ...)

        :param n:   Number of desired names
        :return:    List of n names of the form ["A", "B", ..., "Z", "Z2", "Z3", ...]
        """
        letters = list(string.ascii_uppercase)      # ['A', 'B', ..., 'Z']

        if n <= 26:
            return letters[:n]

        alphanumeric = [f"Z{i-25}" for i in range(27, n+1)]     # Note that 27 gets mapped to ["Z2"],
                                                                #      28 to ["Z2", "Z3"], etc.

        return letters + alphanumeric





############################################################################################


class ChemicalAffinity(NamedTuple):
    """
    Used for binding of ligands to macromolecules (e.g. Transcription Factors to DNA)
    """
    chemical: str   # Name of ligand
    Kd: float       # Dissociation constant; inversely related to binding affinity
                    # Note: dissociation constants are for now assumed to be constant,
                    #       regardless of what other (nearby) sites are occupied by ligands




class Macromol():
    """
    Manage the modeling of large molecules (such as DNA)
    with multiple binding sites (for example, for Transcription Factors).

    CAUTION: this class is currently in flux, and likely to change.
    """

    def __init__(self):
        # TODO: probably pass the SpeciesRegistry as (optional?) argument

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




    def add_macromolecules(self, names: str|list[str]) -> None:
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
        :param ligand:          Id of a previously-declared (bulk) chemical;
                                    if not found, an Exception will be raised       TODO: inconsistent with "macromolecule" arg
        :param Kd:              Dissociation constant, in units of concentration (typically microMolar).
                                    Note that the dissociation constant is inversely proportional to the binding affinity
        :return:                None
        """
        # TODO: restore the validation below
        #assert ligand in self.get_all_labels(), \
        #    f"set_binding_site_affinity(): no chemical named `{ligand}` found; use add_species() first"

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
