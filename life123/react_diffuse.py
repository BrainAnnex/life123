from life123.index_species import IndexSpecies
from life123.uniform_compartment import UniformCompartment
from life123 import BioSim1D, BioSim2D, BioSim3D
from life123 import SpeciesRegistry


class ReactDiffuse:
    """
    WARNING : this is an experimental class, only partially implemented.

    Top-level user entry point for:
        1. Diffusion (without reaction), in 1-, 2- or 3-D
        2. Reaction in a perfectly-mixed uniform compartment (no diffusion)
        3. Reaction-Diffusion, in 1-, 2- or 3-D
    """

    def __init__(self, type :str, species_registry=None, species_ids=None, n_bins=None):
        """

        :param type:                One of "uniform", "1d", "2d" or "3d"
        :param species_registry:    Object of type "SpeciesRegistry";
                                        if not specified, it will be created
        """
        ALLOWED_TYPES = ["uniform", "1d", "2d" or "3d"]
        assert type in ALLOWED_TYPES, f"Unknown type: '{type}'"

        self.module = None

        # Object of type "SpeciesRegistry", with info on the individual species
        if species_registry is None:
            self.species_registry = SpeciesRegistry()
        else:
            self.species_registry = species_registry


        self.index_species = IndexSpecies()
        if species_ids is not None:
            self.index_species.add_species(species_ids)


        match type:
            case "uniform":
                self.module = UniformCompartment(index_species=self.index_species)


            case "1d":
                assert n_bins is not None, "Must pass a value for argument `n_bins`"
                if not species_ids:     # None, or an empty list
                    all_registry_species = self.species_registry.get_all_species_ids()
                    assert len(all_registry_species) > 1, \
                        "No species were specified, either thru argument `species_ids` or" \
                        "thru the argument `species_registry`"
                    self.index_species.add_species(all_registry_species)

                self.module = BioSim1D(n_bins=n_bins, species_data=species_registry, index_species=self.index_species)
