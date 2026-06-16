from dataclasses import dataclass, field, asdict
from itertools import islice
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
    def __post_init__(self):
        # Validate some arguments

        if (self.diffusion_rate is not None):
            assert isinstance(self.diffusion_rate, (float, int, np.integer, np.floating)), \
                f"The value for the diffusion_rate ({self.diffusion_rate}) is not a valid number; it is of type {type(self.diffusion_rate)}"

            assert self.diffusion_rate >= 0, "diffusion_rate, if provided, cannot be negative"

        assert type(self.categories) == list, \
            "The `categories` argument, if provided, must be a list"


        # TODO: add more validation


        # Assign default values to some arguments, based on available arguments
        if self.name == "":
            self.name = self.id

        if self.label == "":
            self.label = self.id




############################################################################################

class SpeciesRegistry:
    """

    """
    def __init__(self, ids=None, species=None, n_species=None):
        """

        [AT MOST, 1 of the following arguments may be passed:]
        :param ids:         [OPTIONAL] A string, or list or tuple of them.  The same value(s) will be used as id, name and label

        :param species:     [OPTIONAL] An object of type "Species", or a list or tuple of them.
                                Their `sort_order` attributes will be modified, to be a zero-base auto-increment: 0, 1, 2, ...

        :param n_species:   [OPTIONAL] The desired number of species.  Their id's, names and labels are automatically
                                assigned as: "A", "B", ..., "Z", "Z2", "Z3", ...
        """

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
                            for arg in (ids, species, n_species) if arg is not None]

        if len(passed_arg_values) == 0:
            return      # All the arguments are None; nothing else to do

        assert len(passed_arg_values) == 1, \
            f"SpeciesRegistry(): cannot pass more than 1 of the arguments `ids`, `species`, `n_species`.  " \
            f"{len(passed_arg_values)} arguments were passed"


        # `ids` argument
        if ids is not None:
            if type(ids) == str:
                ids = [ids]
            else:
                assert (type(ids) == list) or (type(ids) == tuple), \
                    f"SpeciesRegistry(): the `ids` argument, if provided, must be a string, list or tuple"

            for id in ids:
                self.add_species(id=id)


        # `species` argument
        if species is not None:
            if type(species) == Species:
                species = [species]
            else:
                assert type(species) == list or type(species) == tuple, \
                    f"SpeciesRegistry(): the `species` argument, if provided, must be a Species object, or list/tuple of them"

            for species in species:
                assert type(species) == Species, \
                    "SpeciesRegistry(): the `species` argument must be a list or tuple of `Species` objects"

                species.sort_order = self.count
                self.by_id[species.id] = species
                self.count += 1


        # `n_species` argument
        if n_species is not None:
            assert (type(n_species) == int) and (n_species > 0), \
                f"SpeciesRegistry(): the `n_species` argument, if provided, must be a positive integer"

            ids = self._generate_generic_names(n_species)     # Generates the strings "A", "B", ..., "Z", "Z2", "Z3", ...

            for id in ids:
                self.add_species(id=id)



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
                "add_species(): the argument `id` must be a string"

            raise Exception(f"get_species(): no species with id `{id}` found")

        return s



    def get_species_id(self, species_index :int) -> str:
        """
        Return the id of the species with the given index.

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



    def number_of_species(self) -> int:
        """
        Return the number of registered species

        :return:    The number of registered species
        """
        return len(self.by_id)



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
            "add_species(): the argument `id` must be a string"

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
