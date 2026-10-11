import pytest
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

    assert type(sim.system) is UniformCompartment
    assert len(sim.system.system_state) == 3
    assert np.allclose(sim.system.system_state, np.array([0, 0, 0]))


    # Run 2
    sr = SpeciesRegistry(n_species=4)
    sim = ReactionSimulator(species_registry=sr)

    assert type(sim.species_registry) is SpeciesRegistry
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'D']

    assert type(sim.reaction_registry) is ReactionRegistry
    assert sim.reaction_registry.number_of_reactions() == 0

    assert type(sim.species_index_map) is SpeciesIndexMap
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C', 'D']

    assert np.allclose(sim.system.system_state, np.array([0, 0, 0, 0]))


    # Run 3
    rxns = ReactionRegistry()
    assert rxns.version == 0

    rxns.add_reaction(reactants=["A", "B"], products="C", reaction_model="mass action")
    assert rxns.number_of_reactions() == 1
    assert rxns.version == 1    # A mutation in the ReactionRegistry object just occurred

    sim = ReactionSimulator(reaction_registry=rxns)

    assert type(sim.species_registry) is SpeciesRegistry
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B', 'C']

    assert type(sim.reaction_registry) is ReactionRegistry
    assert sim.reaction_registry.number_of_reactions() == 1

    assert type(sim.species_index_map) is SpeciesIndexMap
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B', 'C']

    assert sim.system.number_of_system_species() == 3
    assert np.allclose(sim.system.system_state, np.array([0, 0, 0]))

    # 2nd reaction added
    rxns.add_reaction(reactants="R", products="P", reaction_model="mass action")
    assert rxns.number_of_reactions() == 2
    assert rxns.version == 2    # A mutation in the ReactionRegistry object just occurred
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'P', 'R']  # Auto-updated when rxn added

    assert sim._reaction_registry_version == 1                                      # Stale
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B', 'C']  # Stale
    assert np.allclose(sim.system.system_state, np.array([0, 0, 0]))                # Unchanged
    assert sim.system.number_of_system_species() == 3                               # Unchanged

    with pytest.raises(Exception):
        sim.system.set_system_state([10, 20, 30, 40, 50])   # It's still unknown
                                                            # what species the last 2 array elements correspond to!

    sim.system.set_system_state([10])    # This action will trigger a "lazy" re-synchronization
    assert sim._reaction_registry_version == 2                                                  # Now up-to-date
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'P', 'R']    # Now up-to-date
    assert np.allclose(sim.system.system_state, np.array([10]))                                 # A smaller system state
    assert sim.system.number_of_system_species() == 1                                           # A smaller system state

    sim.system.set_system_state([10, 20, 30, 40, 50])       # Now that the species index got updated, this is a legit operation
    assert np.allclose(sim.system.system_state, np.array([10, 20, 30, 40, 50]))                 # A full system state
    assert sim.system.number_of_system_species() == 5                                           # A full system state

    # 3rd reaction added
    rxns.add_reaction(reactants="P", products=["A", "F"], reaction_model="mass action")
    assert rxns.number_of_reactions() == 3
    assert rxns.version == 3    # A mutation in the ReactionRegistry object just occurred
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'F', 'P', 'R']  # Auto-updated when rxn added

    assert sim._reaction_registry_version == 2                                                  # Stale
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'P', 'R']    # Stale
    #print(sim.species_index_map.species_to_index)  # {'B': 0, 'A': 1, 'C': 2, 'P': 3, 'R': 4}  # Order might vary
    assert np.allclose(sim.system.system_state, np.array([10, 20, 30, 40, 50]))                 # Unchanged
    assert sim.system.number_of_system_species() == 5                                           # Unchanged

    with pytest.raises(Exception):
        sim.system.set_single_conc(conc=123, species_index=5)   # It's still unknown
                                                                # what species corresponds to index 5

    sim.system.set_single_conc(conc=123, species_id="F")    # This action will trigger a "lazy" re-synchronization
    assert sim._reaction_registry_version == 3                                                      # Now up-to-date
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'F', 'P', 'R']   # Now up-to-date
    assert np.allclose(sim.system.system_state, np.array([10, 20, 30, 40, 50, 123]))                # Expanded
    assert sim.system.number_of_system_species() == 6                                               # Expanded

    # 4th reaction added
    rxns.add_reaction(reactants="F", products=["Z"], reaction_model="mass action")
    assert rxns.number_of_reactions() == 4
    assert rxns.version == 4    # A mutation in the ReactionRegistry object just occurred
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'F', 'P', 'R', 'Z']  # Auto-updated when rxn added

    assert sim._reaction_registry_version == 3                                                      # Stale
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'F', 'P', 'R']   # Stale
    assert np.allclose(sim.system.system_state, np.array([10, 20, 30, 40, 50, 123]))                # Unchanged
    assert sim.system.number_of_system_species() == 6                                               # Unchanged

    conc_value = sim.system.get_species_conc(species_id="Z")    # This action will trigger a "lazy" re-synchronization
    assert conc_value == 0      # System state gets expanded with zeroes
    assert sim._reaction_registry_version == 4                                                          # Now up-to-date
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'F', 'P', 'R', 'Z']  # Now up-to-date
    assert np.allclose(sim.system.system_state, np.array([10, 20, 30, 40, 50, 123, 0]))                 # Expanded
    assert sim.system.number_of_system_species() == 7                                                   # Expanded

    # TODO: add a reaction directly thru sim
