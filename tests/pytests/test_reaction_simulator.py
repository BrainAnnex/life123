import math

import numpy as np

from life123.reaction_simulator import ReactionSimulator, VariableTimeSteps
from life123.species_registry import SpeciesRegistry
from life123.reaction_registry import ReactionRegistry
from life123.species_index_map import SpeciesIndexMap
from life123.history import HistoryUniformConcentration, HistoryReactionRate
from life123.uniform_compartment import UniformCompartment




def update_concentrations(conc, delta_conc) -> None:
    """
    UTILITY HELPER FUNCTION.  TODO: eventually move to one of the libraries

    Update the values of the dict `conc` based on the increments in `delta_conc` for the corresponding keys

    :param conc:
    :param delta_conc:
    :return:            None
    """
    for k in conc:
        conc[k] += delta_conc.get(k, 0)     # Missing values default to zero





#####################################################################################################

'''                            ~  Class ReactionSimulator  ~                                     '''

def ________ReactionSimulator________(DIVIDER):
    pass         # Used to get a better structure view in IDEs such asPycharm
#####################################################################################################


def test_CONSTRUCTOR():

    # Run 1
    ind = SpeciesIndexMap(species_ids=["A", "B", "C"])
    sim = ReactionSimulator(species_index_map=ind)

    assert type(sim.species_registry) is SpeciesRegistry
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B', 'C']

    assert type(sim.reaction_registry) is ReactionRegistry
    assert sim.reaction_registry.number_of_reactions() == 0

    assert type(sim.species_index_map) is SpeciesIndexMap
    assert sim.species_index_map == ind

    assert sim.method == "forward_euler"
    assert sim.system_time == 0
    assert sim.previous_system is None
    assert sim.exact is False
    assert math.isclose(sim.temp, 298.15)
    assert sim.diagnostics_enabled is False
    assert sim.diagnostics is None
    assert sim.diagnostic_data_snapshot == {}
    assert sim.reaction_step_diagnostics is None
    assert type(sim.conc_history) is HistoryUniformConcentration
    assert sim.system_rxn_rates == {}
    assert type(sim.adaptive_steps) is VariableTimeSteps

    assert sim._reaction_registry_version == 0

    assert type(sim.uniform_compartment) is UniformCompartment
    assert len(sim.uniform_compartment.system) == 3
    assert np.allclose(sim.uniform_compartment.system, np.array([0, 0, 0]))


    # Run 2
    sr = SpeciesRegistry(n_species=4)
    sim = ReactionSimulator(species_registry=sr)

    assert type(sim.species_registry) is SpeciesRegistry
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'D']

    assert type(sim.reaction_registry) is ReactionRegistry
    assert sim.reaction_registry.number_of_reactions() == 0

    assert type(sim.species_index_map) is SpeciesIndexMap
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C', 'D']

    assert np.allclose(sim.uniform_compartment.system, np.array([0, 0, 0, 0]))


    # Run 3
    rxns = ReactionRegistry()
    rxns.add_reaction(reactants="R", products="P", reaction_model="mass action")
    assert rxns.number_of_reactions() == 1
    sim = ReactionSimulator(reaction_registry=rxns)

    assert type(sim.species_registry) is SpeciesRegistry
    assert sim.species_registry.get_all_species_ids(sort=True) == ['P', 'R']

    assert type(sim.reaction_registry) is ReactionRegistry
    assert sim.reaction_registry.number_of_reactions() == 1

    assert type(sim.species_index_map) is SpeciesIndexMap
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['P', 'R']

    assert np.allclose(sim.uniform_compartment.system, np.array([0, 0]))


    # Run 4
    rxns = ReactionRegistry()
    assert rxns.number_of_reactions() == 0
    assert rxns.version == 0
    sim = ReactionSimulator(reaction_registry=rxns)
    assert np.allclose(sim.uniform_compartment.system, np.array([]))

    rxns.add_reaction(reactants="R", products="P", reaction_model="mass action")
    assert rxns.number_of_reactions() == 1
    assert rxns.version == 1    # A mutation in the ReactionRegistry object just occurred

    assert type(sim.species_registry) is SpeciesRegistry
    assert sim.species_registry.get_all_species_ids(sort=True) == ['P', 'R']    # Automatically added when the reaction was added

    assert type(sim.reaction_registry) is ReactionRegistry
    assert sim.reaction_registry.number_of_reactions() == 1

    assert type(sim.species_index_map) is SpeciesIndexMap
    # The SpeciesIndexMap object has not yet registered the new species
    assert sim._reaction_registry_version == 0  # The stale version
    assert sim.species_index_map.number_of_mapped_species() == 0
    assert sim.uniform_compartment.number_of_system_species() == 0

    sim.set_system_conc([10, 20])
    assert sim._reaction_registry_version == 1  # The above update triggered a re-synchronization

    assert sim.species_index_map.get_all_species_ids(sort=True) == ['P', 'R']

    assert np.allclose(sim.uniform_compartment.system, np.array([10, 20]))

    # TODO next: add the reaction late on, thru the sim object
