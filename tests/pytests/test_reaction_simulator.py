import pytest
import math
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal

from life123.reaction_simulator import ReactionSimulator, VariableTimeSteps, ReactionStep
from life123.species_registry import SpeciesRegistry, MacroMolecules
from life123.reaction_registry import ReactionRegistry
from life123.species_index_map import SpeciesIndexMap
from life123.history import HistoryUniformConcentration, HistoryReactionRate
from life123.uniform_compartment import UniformCompartment
from life123.reactions import ReactionDefinition
from life123.reaction_simulator import ExcessiveTimeStepHard
from life123.kinetic_models import Custom_Model
from life123.reaction_kinetics import ReactionKinetics
from life123.collections import CollectionTabular
from life123.diagnostics import Diagnostics
from tests.utilities.comparisons import *




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


def test_CONSTRUCTOR_1():
    # This series of tests focuses on instantiation using arg `species_index_map`
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


def test_CONSTRUCTOR_2():
    # This series of tests focuses on instantiation using arg `species_registry`
    sr = SpeciesRegistry(n_species=4)
    sim = ReactionSimulator(species_registry=sr)

    assert type(sim.species_registry) is SpeciesRegistry
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'D']

    assert type(sim.reaction_registry) is ReactionRegistry
    assert sim.reaction_registry.number_of_reactions() == 0

    assert type(sim.species_index_map) is SpeciesIndexMap
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C', 'D']

    assert np.allclose(sim.system.system_state, np.array([0, 0, 0, 0]))


def test_CONSTRUCTOR_3():
    # This series of tests focuses on instantiation using arg `reaction_registry`
    rxns = ReactionRegistry()
    assert rxns.version == 0

    # 1st reaction, prior to instantiating ReactionSimulator
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
    rxns.add_reaction(reactants="F", products="Z", reaction_model="mass action")
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

    # 5th reaction added - this time thru the ReactionSimulation object
    sim.add_reaction(reactants="J", products="K", reaction_model="mass action")
    assert rxns.number_of_reactions() == 5
    assert rxns.version == 5    # A mutation in the ReactionRegistry object just occurred
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'F', 'J', 'K', 'P', 'R', 'Z']  # Auto-updated when rxn added

    # Nothing is stale this time around!
    assert sim._reaction_registry_version == 5
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'F', 'J', 'K', 'P', 'R', 'Z']    # Up-to-date
    assert np.allclose(sim.system.system_state, np.array([10, 20, 30, 40, 50, 123, 0, 0, 0]))           # Expanded for the 2 new species
    assert sim.system.number_of_system_species() == 9                                                   # Expanded for the 2 new species


    sp_extra = SpeciesRegistry()
    with pytest.raises(Exception):
        # Argument incompatibility
        ReactionSimulator(reaction_registry=rxns, species_registry=sp_extra)



def test_CONSTRUCTOR_4():
    # This series of tests focuses on instantiation using arg `species_ids`
    with pytest.raises(Exception):
        ReactionSimulator(species_ids=123)

    sim = ReactionSimulator(species_ids=["A", "B"])

    assert type(sim.species_registry) is SpeciesRegistry
    assert sim.species_registry.get_all_species_ids(sort=True) == ['A', 'B']

    assert type(sim.reaction_registry) is ReactionRegistry
    assert sim.reaction_registry.number_of_reactions() == 0

    assert type(sim.species_index_map) is SpeciesIndexMap
    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B']

    assert sim.system.number_of_system_species() == 2
    assert np.allclose(sim.system.system_state, np.array([0, 0]))


def test_CONSTRUCTOR_5():
    # This series of tests focuses on instantiation with no arguments

    sim = ReactionSimulator()

    #TODO: continue



def test_specify_steps():
    sim = ReactionSimulator()

    with pytest.raises(Exception):
        # Too few arguments
        sim.specify_steps()
        sim.specify_steps(duration=15.1)
        sim.specify_steps(initial_step=0.2)
        sim.specify_steps(n_steps=30)
        # Too many arguments
        sim.specify_steps(duration=15.1, initial_step=0.2, n_steps=30)

    assert sim.specify_steps(initial_step=0.5, n_steps=24) == (0.5, 24)
    assert sim.specify_steps(duration=12.0, initial_step=0.5) == (0.5, 24)
    assert sim.specify_steps(duration=12.0, n_steps=24) == (0.5, 24)



def test_single_compartment_react():

    # TODO: test `self.number_neg_concs` and `self.number_soft_aborts`

    # Test based on experiment "cycles_1"  (with 3 reactions)
    species_registry = SpeciesRegistry(ids=["A", "B", "C", "E_high", "E_low"])
    rxns = ReactionRegistry(species_registry=species_registry)

    # Unimolecular reaction A <-> B, mostly in forward direction (favored energetically)
    rxns.add_reaction(reactants="A", products="B",
                      reaction_model="mass action", kinetic_parameters={"kF": 9., "kR": 3.})

    # Unimolecular reaction B <-> C, also favored energetically
    rxns.add_reaction(reactants="B", products="C",
                      reaction_model="mass action", kinetic_parameters={"kF": 8., "kR": 4.})

    # Reaction C + E_High <-> A + E_Low, also favored energetically, but kinetically slow.
    # HYPOTHETICALLY treated as a mass-action reaction
    rxns.add_reaction(reactants=["C" , "E_high"], products=["A", "E_low"],
                      reaction_model="mass action", kinetic_parameters={"kF": 1., "kR": 0.2})


    sim = ReactionSimulator(reaction_registry=rxns, enable_diagnostics=True)

    assert sim.species_index_map.get_all_species_ids(sort=True) == ['A', 'B', 'C', 'E_high', 'E_low']
    assert sim.reaction_registry.number_of_reactions() == 3
    assert sim.system_time == 0


    initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}
    sim.set_conc(conc=initial_conc, snapshot=True)

    assert np.allclose(sim.system.system_state, [100., 0., 0., 1000. , 0.])
    assert np.allclose(sim.system.system_state, [100., 0., 0., 1000. , 0.])


    sim.single_compartment_react(initial_step=0.0005, target_end_time=0.0035, variable_steps=False)
    assert np.allclose(sim.system.system_state, sim.system.system_state)

    run1 = sim.system.get_system_conc()

    assert np.allclose(sim.system_time, 0.0035)
    assert np.allclose(run1, [9.69252541e+01, 3.05696280e+00, 1.77831454e-02, 9.99980686e+02, 1.93144884e-02])  # Line 7, below

    df = sim.diagnostics.conc_data.get_diagnostic_conc_data()
    row_0 = df.iloc[0]      # Line 0, below
    expected_row = pd.Series([0.0, 100.0, 0.0, 0.0, 1000.0, 0.0, "Set concentration"],
                             index =['TIME', 'A', 'B', 'C', 'E_high', 'E_low', 'caption'])
    assert row_0.equals(other=expected_row)
    """
         TIME           A         B  ...       E_high     E_low            caption
    0  0.0000  100.000000  0.000000  ...  1000.000000  0.000000  Set concentration
    1  0.0005   99.550000  0.450000  ...  1000.000000  0.000000                   
    2  0.0010   99.102700  0.895500  ...  1000.000000  0.000000                   
    3  0.0015   98.658981  1.336540  ...   999.999100  0.000900                   
    4  0.0020   98.219251  1.773164  ...   999.996870  0.003130                   
    5  0.0025   97.783686  2.205413  ...   999.993108  0.006892                   
    6  0.0030   97.352350  2.633332  ...   999.987725  0.012275                   
    7  0.0035   96.925254  3.056963  ...   999.980686  0.019314
    """
    assert sim.diagnostics.conc_data.explain_time_advance(return_times=True, silent=True) == \
               ([0.0, 0.0035], [0.0005])



    # Now repeat the process, step-by-step
    #print("\n\nNEW RUN:")
    sim2 = ReactionSimulator(reaction_registry=rxns, enable_diagnostics=True)
    assert sim2.system_time == 0

    sim2.set_conc(conc=initial_conc, snapshot=True)

    assert np.allclose(sim2.system.system_state, [100., 0., 0., 1000. , 0.])

    for i in range(7):
        sim2.single_compartment_react(initial_step=0.0005, n_steps=1, variable_steps=False)
        assert math.isclose(sim2.system_time, 0.0005*(i+1))
        #print("   SYSTEM TIME NOW: ", sim2.system_time, "\n")


    run2 = sim2.system.get_system_conc()
    assert np.allclose(run2, run1)                  # Same result as before
    assert np.allclose(sim2.system_time, 0.0035)    # Same as before

    df2 = sim2.diagnostics.conc_data.get_diagnostic_conc_data()
    row_0 = df2.iloc[0]      # The initial row
    assert row_0.equals(other=expected_row)     # expected_row is the same as for the first run

    #print(sim2.diagnostics.explain_time_advance(return_times=True))
    assert sim2.diagnostics.conc_data.explain_time_advance(return_times=True, silent=True) == \
               ([0.0, 0.0035], [0.0005])



def test_single_compartment_react__catch_neg_conc():
    # Based on "Run 3" of experiment "negative_concentrations_1

    species_registry = SpeciesRegistry(ids=["U", "X", "S"])
    rxns = ReactionRegistry(species_registry=species_registry)

    # 1st reaction: 2 S <-> U , HYPOTHETICALLY with 1st-order kinetics for all species (mostly forward)
    r = ReactionDefinition(reactants=[(2, "S")], products="U", species_registry=species_registry,
                           reaction_model="custom",
                           kinetic_parameters={"kF": 8., "kR": 2., "rate_function": Custom_Model.kinetic_rate_first_order})
    rxns.register_reaction(r)

    sim = ReactionSimulator(reaction_registry=rxns, method="forward_euler",
                            enable_diagnostics=True)

    # 2nd reaction: elementary reaction S <-> X (mostly forward)
    rxns.add_reaction(reactants="S", products="X",
                      reaction_model="mass action", kinetic_parameters={"kF": 6., "kR": 3.})

    sim.system.set_conc(conc={"U": 50., "X": 100., "S": 0.})
    assert sim.species_index_map.index_to_species == ['S', 'U', 'X']
    assert sim.species_index_map.species_to_index == {'S': 0, 'U': 1, 'X': 2}


    # Perform a first reaction step
    sim.single_compartment_react(initial_step=0.25, n_steps=1, variable_steps=False)

    assert np.allclose(sim.system_time, 0.25)
    assert np.allclose(sim.system.system_state, [ 125., 25.,  25. ])   # "S", "U", "X"
    assert sim.system.get_conc_as_dict() == {'S': 125, 'U': 25, 'X': 25}


    # Attempt a second reaction step
    with pytest.raises(ExcessiveTimeStepHard) as ex:
        # This step would make [S] negative
        sim.single_compartment_react(initial_step=0.25, n_steps=1, variable_steps=False)
    msg = ex.value.details
    #details = ex.value.details
    #assert details["function"] == "reaction_step_common_fixed_step"
    #assert details["delta_time"] == 0.25
    #assert details["caption"] == "aborted: neg. conc. in `S` from rxn # 0"
    #assert details["system_time"] == 0.25
    #assert details["rxn_index"] == 0
    #assert details["rate"] == 950   # 8. * 125 - 2. * 25
    #assert details["species_id"] == "S"
    #assert details["previous_function"] == "_validate_increment"
    assert msg == """reaction_step_common_fixed_step(): unable to complete the reaction step.  System time is now: 0.25 
Try REDUCING the time step, or switching to variable time steps. 
DETAILS received from the low-level function `_validate_increment()`: 
      The attempted time step (0.25) would lead to a NEGATIVE concentration in the species `S` from the reaction `2 S <-> U` (rxn # 0)
      Baseline concentration value of `S` : 125 at system time 0.25; requested change (NOT carried out): -475, due to a reaction rate of 950"""

    # Nothing has changed, since that last step wasn't actually taken
    assert np.allclose(sim.system_time, 0.25)
    assert np.allclose(sim.system.system_state, [ 125., 25.,  25. ])

    # A smaller step saves the day!
    sim.single_compartment_react(initial_step=0.03, n_steps=1, variable_steps=False)
    assert np.allclose(sim.system_time, 0.28)    # 0.25 + 0.03
    assert np.allclose(sim.system.system_state, [47.75, 53.5,  45.25])



def test_single_compartment_react__variable_steps():
    # Based on experiment "variable_steps_1"

    species_registry = SpeciesRegistry(ids=["U", "X", "S"])
    rxns = ReactionRegistry(species_registry=species_registry)

    # 1st reaction:  2 S <-> U , HYPOTHETICALLY with 1st-order kinetics for all species (mostly forward)
    r = ReactionDefinition(reactants=[(2, "S")], products="U", species_registry=species_registry,
                           reaction_model="custom",
                           kinetic_parameters={"kF": 8., "kR": 2., "rate_function": Custom_Model.kinetic_rate_first_order})
    rxns.register_reaction(r)

    #uc = UniformCompartment(species_data=species_registry, preset=None, reactions=rxns)
    sim = ReactionSimulator(reaction_registry=rxns, method="forward_euler",
                            enable_diagnostics=True, preset=None)

    # 2nd reaction: elementary reaction S <-> X (mostly forward)
    rxns.add_reaction(reactants="S", products="X",
                      reaction_model="mass action", kinetic_parameters={"kF": 6., "kR": 3.})

    sim.system.set_conc(conc={"U": 50., "X": 100., "S": 0.}, snapshot_funct=sim.capture_conc_snapshot)
    assert sim.species_index_map.index_to_species == ['S', 'U', 'X']
    assert sim.species_index_map.species_to_index == {'S': 0, 'U': 1, 'X': 2}


    sim.adaptive_steps.set_thresholds(norm="norm_A", low=0.25, high=0.64, abort=1.44)
    sim.adaptive_steps.set_step_factors(abort=0.5, downshift=0.5, upshift=2.0)

    sim.single_compartment_react(initial_step=0.01, target_end_time=0.2,
                                 variable_steps=True)

    df = sim.get_history()
    print(df)
    assert len(df) == 23

    assert np.allclose(df.iloc[0][['SYSTEM TIME', 'S', 'U', 'X']].to_numpy(dtype='float16'),
                       [0.0000,  0.000000,  50.000000,  100.000000])
    assert df.iloc[0]["caption"] == "Set concentration"

    assert np.allclose(df.iloc[1][['SYSTEM TIME', 'S', 'U', 'X']].to_numpy(dtype='float16'),
                       [0.0050,  2.500000,  49.500000,  98.500000], rtol=1e-03)     # Notice the halved step size
    assert df.iloc[1]["caption"] == "1st reaction step"
    assert df.iloc[1]["step"] == "1"

    assert np.allclose(df.iloc[2][['SYSTEM TIME', 'S', 'U', 'X']].to_numpy(dtype='float16'),
                       [0.0075,  3.596250,  49.302500,  97.798750], rtol=1e-03)
    assert df.iloc[2]["step"] == "2"

    assert np.allclose(df.iloc[22][['SYSTEM TIME', 'S', 'U', 'X']].to_numpy(dtype='float16'),
                       [0.2050,  19.671183,  55.600598,   69.127620], rtol=1e-03)
    assert df.iloc[22]["caption"] == "last reaction step"
    assert df.iloc[22]["step"] == "22"



def test__single_compartment_react_main_loop_1():
    # FIXED steps

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})
    assert ind.index_to_species == ["A", "B"]

    sim = ReactionSimulator(system=np.array([10., 50.]), species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler")
    sim.system_time = 8

    new_count, recommended_next_step = sim._single_compartment_react_main_loop(time_step=0.1, variable_steps=False,
                                                                               step_count=0, n_steps=1000)
    assert new_count == 1
    assert math.isclose(recommended_next_step, 0.1)
    assert np.allclose(sim.system.system_state, [17., 43.])
    assert np.allclose(sim.previous_system, [10., 50.])
    assert math.isclose(sim.system_time, 8.1)

    # Check the concentration history
    df = sim.conc_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 8.1, "A": 17.0, "B": 43.0, "step": "1", "caption": ""}   # Note: "step" is a string!
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)

    # Check the rate history
    expected_rates = {0: -70.0}
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, expected_rates)
    df = sim.rate_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 8, "rxn0_rate": -70.0, "step": "0"}
    # Note: "step" is a string!  System time and step refer to the START of the simulation step
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)



    # RESET to initial concentrations, and simulate a step just shy of excessive
    sim.system.set_system_state([10., 50.])
    sim.system_time = 0
    sim.rate_history = HistoryReactionRate(active=True)         # We're not instantiating a new "ReactionSimulator"; so, we must reset some values
    sim.conc_history = HistoryUniformConcentration(active=True)
    # Pick a time step just below 50/70 (where 70 is the initial reaction rate),
    # to bring the second species concentration to almost zero
    new_count, recommended_next_step = sim._single_compartment_react_main_loop(time_step=0.7142857142, variable_steps=False,
                                                                               step_count=5, n_steps=1000)
    assert new_count == 6
    assert math.isclose(recommended_next_step, 0.7142857142)
    assert np.allclose(sim.system.system_state, [60., 0.])       # This time step was so large that it converted all B to A !
    assert math.isclose(sim.system_time, 0.7142857142)
    assert np.allclose(sim.previous_system, [10., 50.])

    # Check the concentration history
    df = sim.conc_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 0.7142857142, "A": 60.0, "B": 0.0, "step": "6", "caption": ""}   # Note: "step" is a string!
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)

    # Check the rate history
    expected_rates = {0: -70.0}
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, expected_rates)
    df = sim.rate_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 0, "rxn0_rate": -70.0, "step": "5"}
    # Note: "step" is a string!  System time and step refer to the START of the simulation step
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)



    # RESET to initial concentrations, and attempt an excessive step
    sim.system.set_system_state([10., 50.])
    sim.system_time = 20
    # Pick a time step just above the previous time, to tip [B] into negative values
    with pytest.raises(ExcessiveTimeStepHard) as ex:
        sim._single_compartment_react_main_loop(time_step=0.714285715, variable_steps=False,
                                                step_count=1, n_steps=1000)
    msg = ex.value.details
    #details = ex.value.details
    #assert details["function"] == "reaction_step_common_fixed_step"
    #assert math.isclose(details["delta_time"], 0.714285715)
    #assert details["caption"] == 'aborted: neg. conc. in `B` from rxn # 0'
    #assert details["system_time"] == 20
    #assert details["rate"] == -70.0
    #assert details["rxn_index"] == 0
    #assert details["previous_function"] == "_validate_increment"
    assert msg == """reaction_step_common_fixed_step(): unable to complete the reaction step.  System time is now: 20 
Try REDUCING the time step, or switching to variable time steps. 
DETAILS received from the low-level function `_validate_increment()`: 
      The attempted time step (0.714286) would lead to a NEGATIVE concentration in the species `B` from the reaction `A <-> B` (rxn # 0)
      Baseline concentration value of `B` : 50 at system time 20; requested change (NOT carried out): -50, due to a reaction rate of -70"""

    # The following values are all unchanged
    assert sim.system_time == 20
    assert np.allclose(sim.system.system_state, [10., 50.])
    assert np.allclose(sim.previous_system, [10., 50.])



def test__single_compartment_react_main_loop_2_a():
    # VARIABLE steps

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})
    assert ind.index_to_species == ["A", "B"]

    ### First simulation
    sim = ReactionSimulator(system=np.array([10., 50.]), species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler", preset="fast")  # Note the preset
    sim.system_time = 8

    new_count, step_recommended = sim._single_compartment_react_main_loop(time_step=0.02, variable_steps=True,
                                                                          step_count=0, n_steps=1000)
    assert new_count == 1
    assert math.isclose(step_recommended, 0.02)     # "Stay on course" (i.e. keep same step size)
    assert np.allclose(sim.system.system_state, [11.4, 48.6])
    assert np.allclose(sim.previous_system, [10., 50.])
    assert math.isclose(sim.system_time, 8.02)

    # Inspect the `sim.reaction_step_diagnostics` data object
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.02)
    assert sim.reaction_step_diagnostics.decision_data == {'action': 'stay', 'operation': 'stay',
                                                           'step_factor': 1,
                                                           'applicable_norms': 'ALL'}
    assert compare_dicts(sim.reaction_step_diagnostics.norms, {'norm_A': 0.98, 'norm_B': 0.14})
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})

    # Check the concentration history
    df = sim.conc_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 8.02, "A": 11.4, "B": 48.6, "step": "1", "caption": ""}   # Note: "step" is a string!
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)

    # Check the rate history
    expected_rates = {0: -70.0}
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, expected_rates)
    df = sim.rate_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 8, "rxn0_rate": -70.0, "step": "0"}
    # Notes: "step" is a string!  System time and step refer to the START of the simulation step
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)



    ### New simulation, with re-instantiated "ReactionSimulator" object, and somewhat large time step (but not so large as to cause hard errors)
    initial_system = np.array([10., 50.])
    sim = ReactionSimulator(system=initial_system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler", preset="fast")  # Note the preset
    sim.system_time = 17

    new_count, step_recommended = sim._single_compartment_react_main_loop(time_step=0.04, variable_steps=True,
                                                                          step_count=5, n_steps=1000)
    assert new_count == 6
    assert math.isclose(step_recommended, 0.0192)   # "Go smaller" in next round
    assert np.allclose(sim.system.system_state, [11.68, 48.32])
    assert np.allclose(sim.previous_system, [10., 50.])
    assert math.isclose(sim.system_time, 17.024)    # A step of just 0.024 was actually taken

    # Inspect the `sim.reaction_step_diagnostics` data object
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.024)
    assert sim.reaction_step_diagnostics.decision_data == {'action': 'high', 'operation': 'downshift',
                                                           'step_factor': 0.8,
                                                           'applicable_norms': ['norm_A']}
    assert compare_dicts(sim.reaction_step_diagnostics.norms, {'norm_A': 1.4112, 'norm_B': 0.168})
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})

    # Check the concentration history
    df = sim.conc_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 17.024, "A": 11.68, "B": 48.32, "step": "6", "caption": ""}   # Note: "step" is a string!
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)

    # Check the rate history
    expected_rates = {0: -70.0}
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, expected_rates)
    df = sim.rate_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 17, "rxn0_rate": -70.0, "step": "5"}
    # Notes: "step" is a string!  System time and step refer to the START of the simulation step
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)



    ### New simulation, with re-instantiated "ReactionSimulator" object, and excessive large time step (so large as to cause a backtracking)
    initial_system = np.array([10., 50.])
    sim = ReactionSimulator(system=initial_system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler", preset="fast")  # Note the preset
    sim.system_time = 666

    new_count, step_recommended = sim._single_compartment_react_main_loop(time_step=0.8, variable_steps=True,
                                                                          step_count=13, n_steps=1000)  # , explain_variable_steps=(-1, 1000)

    assert new_count == 14
    assert math.isclose(step_recommended, 0.0186624)   # "Go much smaller" in next round
    assert np.allclose(sim.system.system_state, [11.306368, 48.693632])
    assert np.allclose(sim.previous_system, [10., 50.])
    assert math.isclose(sim.system_time, 666.0186624)    # A step of just 0.0186624 was actually taken

    # Inspect the `sim.reaction_step_diagnostics` data object
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.0186624)    # 0.8 * 0.5 * (0.6)^6, from multiple re-tries
    assert sim.reaction_step_diagnostics.decision_data == {'action': 'stay', 'operation': 'stay',
                                                           'step_factor': 1,
                                                           'applicable_norms': 'ALL'}
    assert compare_dicts(sim.reaction_step_diagnostics.norms, {'norm_A': 0.853298675712, 'norm_B': 0.1306368})
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})

    # Check the concentration history
    df = sim.conc_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 666.0186624, "A": 11.306368, "B": 48.693632, "step": "14", "caption": ""}   # Note: "step" is a string!
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)

    # Check the rate history
    expected_rates = {0: -70.0}
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, expected_rates)
    df = sim.rate_history.get_history().get_dataframe()
    row_expected = {"SYSTEM TIME": 666, "rxn0_rate": -70.0, "step": "13"}
    # Notes: "step" is a string!  System time and step refer to the START of the simulation step
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert_frame_equal(df, df_expected)



def test__single_compartment_react_main_loop_2_b():
    # VARIABLE steps, with diagnostics

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})
    assert ind.index_to_species == ["A", "B"]

    ### Repeat the simulation from previous test, this time with diagnostics enabled (excessive large time step;
    # so large as to cause a backtracking)
    initial_system = np.array([10., 50.])
    sim = ReactionSimulator(system=initial_system, species_index_map=ind, enable_diagnostics=True,
                            reaction_registry=rxns, method="forward_euler", preset="fast")  # Note the preset
    sim.system_time = 666

    new_count, step_recommended = sim._single_compartment_react_main_loop(time_step=0.8, variable_steps=True,
                                                                          step_count=13, n_steps=1000) #, explain_variable_steps=(-1, 1000)

    assert new_count == 14
    assert math.isclose(step_recommended, 0.0186624)   # "Go much smaller" in next round
    assert np.allclose(sim.system.system_state, [11.306368, 48.693632])
    assert np.allclose(sim.previous_system, [10., 50.])
    assert math.isclose(sim.system_time, 666.0186624)    # A step of just 0.0186624 was actually taken


    # Verify the diagnostic data: part 1 - the "diagnostic_rxn_data"
    assert type(sim.diagnostics.rxn_data.diagnostic_rxn_data) is dict
    assert len(sim.diagnostics.rxn_data.diagnostic_rxn_data) == 1
    coll_tab = sim.diagnostics.rxn_data.diagnostic_rxn_data[0]   # For reaction 0 (the only reaction we have)
    assert type(coll_tab) is CollectionTabular
    df = coll_tab.get_dataframe()
    with pd.option_context('display.max_columns', None):    # The configuration changes automatically revert after exiting the 'with' block
        print(df)

    rows_expected = [
                    {"START_TIME": 666, "time_step": 0.8, "aborted": True, "Delta A": np.nan, "Delta B": np.nan, "rate": -70., "caption": "aborted: neg. conc. in `B`"},
                    {"START_TIME": 666, "time_step": 0.4, "aborted": True, "Delta A": 28., "Delta B": -28., "rate": -70., "caption": "aborted: excessive norm value(s)"},
                    {"START_TIME": 666, "time_step": 0.24, "aborted": True, "Delta A": 16.8, "Delta B": -16.8, "rate": -70., "caption": "aborted: excessive norm value(s)"},
                    {"START_TIME": 666, "time_step": 0.144, "aborted": True, "Delta A": 10.08, "Delta B": -10.08, "rate": -70., "caption": "aborted: excessive norm value(s)"},
                    {"START_TIME": 666, "time_step": 0.0864, "aborted": True, "Delta A": 6.048, "Delta B": -6.048, "rate": -70., "caption": "aborted: excessive norm value(s)"},
                    {"START_TIME": 666, "time_step": 0.05184, "aborted": True, "Delta A": 3.6288, "Delta B": -3.6288, "rate": -70., "caption": "aborted: excessive norm value(s)"},
                    {"START_TIME": 666, "time_step": 0.031104, "aborted": True, "Delta A": 2.17728, "Delta B": -2.17728, "rate": -70., "caption": "aborted: excessive norm value(s)"},
                    {"START_TIME": 666, "time_step": 0.0186624, "aborted": False, "Delta A": 1.306368, "Delta B": -1.306368, "rate": -70., "caption": ""}
                    ]
    df_expected = pd.DataFrame(rows_expected)
    assert_frame_equal(df, df_expected)



def test_reaction_step_common_fixed_step_1():
    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})
    assert ind.index_to_species == ["A", "B"]

    system = np.array([10, 50])      # 10. goes to "A" and 50. to "B"
    sim = ReactionSimulator(system=system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler")
    assert sim.system.get_conc_as_dict() == {"A": 10, "B": 50}

    result = sim.reaction_step_common_fixed_step(delta_time=0.1)
    assert np.allclose(result, [7, -7])     # The increment vector
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.1)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})

    sim.system.set_system_state([10, 50])     # Reset the system state, and repeat identically
    result = sim.reaction_step_common_fixed_step(delta_time=0.1)
    assert np.allclose(result, [7, -7])     # The increment vector
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.1)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})

    sim.system.set_system_state([10, 50])     # Reset the system state, and repeat with double time interval
    result = sim.reaction_step_common_fixed_step(delta_time=0.2)
    assert np.allclose(result, [14, -14])     # The increment vector has doubled
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.2)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})

    sim.system.set_system_state([10, 50])     # Reset the system state,
                                        # and repeat with large time interval just short of triggering error
    result = sim.reaction_step_common_fixed_step(delta_time=0.7142857)
    assert np.allclose(result, [49.999999, -49.999999])     # The increment vector is a hair away from making [B] negative
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.7142857)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    sim.system.set_system_state([10, 50])     # Reset the system state
    with pytest.raises(ExcessiveTimeStepHard) as ex:      # Excessive time step that would make [B] negative
        sim.reaction_step_common_fixed_step(delta_time=0.8)
    msg = ex.value.details
    #details = ex.value.details
    #assert details["function"] == "reaction_step_common_fixed_step"
    #assert details["delta_time"] == 0.8
    #assert details["system_time"] == 0
    #assert details["rate"] == -70.0
    #assert details["rxn_index"] == 0
    #assert details["previous_function"] == "_validate_increment"
    assert msg == """reaction_step_common_fixed_step(): unable to complete the reaction step.  System time is now: 0 
Try REDUCING the time step, or switching to variable time steps. 
DETAILS received from the low-level function `_validate_increment()`: 
      The attempted time step (0.8) would lead to a NEGATIVE concentration in the species `B` from the reaction `A <-> B` (rxn # 0)
      Baseline concentration value of `B` : 50 at system time 0; requested change (NOT carried out): -56, due to a reaction rate of -70"""



    sim.system.set_system_state([10, 50])     # Reset the system state
    sim.method = "heun"
    with pytest.raises(ExcessiveTimeStepHard) as ex:      # Excessive time step that would make [B] negative
        sim.reaction_step_common_fixed_step(delta_time=0.8)
    msg = ex.value.details
    #details = ex.value.details
    #assert details["function"] == "reaction_step_common_fixed_step"
    #assert details["delta_time"] == 0.8
    #assert details["rate"] == -70.0
    #assert details["previous_function"] == "heun_single_rxn"
    assert msg == """reaction_step_common_fixed_step(): unable to complete the reaction step.  System time is now: 0 
Try REDUCING the time step, or switching to variable time steps. 
DETAILS received from the low-level function `heun_single_rxn()`: 
heun_single_rxn(): excessive time step (0.8), leading to negative concentrations"""


def test_reaction_step_common_fixed_step_2():
    species_registry = SpeciesRegistry()

    # Reaction A + B <-> C , with mass-action kinetics
    # Based on experiment "1D/reactions/reaction4"
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants=["A" , "B"], products="C",
                      reaction_model="mass action",
                      kinetic_parameters={"kF": 5., "kR": 2.})

    #ind = SpeciesIndexMap({"A", "B", "C"})     # Now automatically created

    system = np.array([10, 50, 20])
    sim = ReactionSimulator(system=system,
                            reaction_registry=rxns, method="forward_euler")
    assert sim.system.get_conc_as_dict() == {"A": 10, "B": 50, "C": 20}

    # The species index was automatically created
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C']
    assert sim.species_index_map.species_to_index == {"A": 0, "B": 1, "C": 2}


    result = sim.reaction_step_common_fixed_step(delta_time=0.002)
    assert np.allclose(result, [-4.92, -4.92, 4.92])     # The increment vector
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.002)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: 2460}) # 5. * 10 * 50 - 2. * 20

    result = sim.reaction_step_common_fixed_step(delta_time=0.00406504)
    assert np.allclose(result, [-9.9999984, -9.9999984,  9.9999984])   # A hair from making [A] negative
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.00406504)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: 2460})

    with pytest.raises(ExcessiveTimeStepHard) as ex:
        sim.reaction_step_common_fixed_step(delta_time=0.0040651)    # A step so large that it would make [A] negative
    msg = ex.value.details
    #details = ex.value.details
    #assert details["function"] == "reaction_step_common_fixed_step"
    #assert math.isclose(details["delta_time"], 0.0040651)
    #assert details["system_time"] == 0
    #assert details["rate"] == 2460
    #assert details["rxn_index"] == 0
    #assert details["previous_function"] == "_validate_increment"
    assert msg == """reaction_step_common_fixed_step(): unable to complete the reaction step.  System time is now: 0 
Try REDUCING the time step, or switching to variable time steps. 
DETAILS received from the low-level function `_validate_increment()`: 
      The attempted time step (0.0040651) would lead to a NEGATIVE concentration in the species `A` from the reaction `A + B <-> C` (rxn # 0)
      Baseline concentration value of `A` : 10 at system time 0; requested change (NOT carried out): -10.0001, due to a reaction rate of 2460"""


    rxns.clear_reactions_data()       # Re-start with a blank slate of reactions
    # No change to index
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C']
    assert sim.species_index_map.species_to_index == {"A": 0, "B": 1, "C": 2}

    # A <-> B and X <-> Y, both with mass-action kinetics
    rxns.add_reaction(reactants="A", products="B",
                    reaction_model="mass action", kinetic_parameters={"kF": 300., "kR": 2.})
    rxns.add_reaction(reactants="X", products="Y",
                    reaction_model="mass action")       # Extraneous reaction that doesn't participate,
                                                        # because the initial concentrations of X and Y are 0

    with pytest.raises(Exception):
        sim.system.set_system_state([10., 50., 0, 0])       # Trying to set 4 values, but the specie index is 3 values

    sim.system.set_conc({"A": 10, "B": 50})             # This will trigger a synchronization
    # Check the expanded index
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C', 'X', 'Y']
    assert sim.species_index_map.species_to_index == {"A": 0, "B": 1, "C": 2, "X": 3, "Y": 4}

    sim.system.set_system_state([10., 50., 0, 0, 0])

    # No change
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C', 'X', 'Y']
    assert sim.species_index_map.species_to_index == {"A": 0, "B": 1, "C": 2, "X": 3, "Y": 4}

    assert sim.system.get_conc_as_dict() == {"A": 10, "B": 50, "C": 0, "X": 0, "Y": 0}


    result = sim.reaction_step_common_fixed_step(delta_time=0.002)
    # 10 * 300 * .002 - 50 * 2 * .002 = 5.8
    assert np.allclose(result, [-5.8,  5.8, 0,  0, 0])   # X and Y aren't affected by this reaction; hence, 0 change


    # Add the reaction we saw earlier, A + B <-> C, and use different initial concentrations
    rxns.add_reaction(reactants=["A" , "B"], products="C",
                    reaction_model="mass action",
                    kinetic_parameters={"kF": 5., "kR": 2.})
    assert rxns.number_of_reactions() == 3
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C', 'X', 'Y']
    assert sim.species_index_map.species_to_index == {"A": 0, "B": 1, "C": 2, "X": 3, "Y": 4}


    sim.system.set_system_state([10, 50, 20, 0, 0])
    assert sim.system.get_conc_as_dict() == {"A": 10, "B": 50, "C": 20, "X": 0, "Y": 0}    # The reaction X <-> Y will continue to remain irrelevant

    # We saw in earlier runs, with our initial concentrations, that
    # over a delta_time=0.02, one reaction causes a change in [A] of -4.92,
    # and the other a change of -5.8
    # Individually, neither is problematic, given that [A] is initially 10,
    # but combined (-10.72) they would make [A] negative!
    with pytest.raises(ExcessiveTimeStepHard) as ex:
        sim.reaction_step_common_fixed_step(delta_time=0.002)
    msg = ex.value.details
    #details = ex.value.details
    #assert details["function"] == "reaction_step_common_fixed_step"
    #assert details["previous_function"] == "attempt_reaction_step"
    #assert math.isclose(details["delta_time"], 0.002)
    assert msg == """reaction_step_common_fixed_step(): unable to complete the reaction step.  System time is now: 0 
Try REDUCING the time step, or switching to variable time steps. 
DETAILS received from the low-level function `attempt_reaction_step()`: 
      The attempted time step (0.002) would lead to a NEGATIVE concentration 
      in one or more of the species (for instance `A`, of index 0), from the COMBINED reactions.
      Baseline concentration values: [10 50 20  0  0] at system time 0; requested changes (NOT carried out): [-10.72   0.88   4.92   0.     0.  ]"""




def test_reaction_step_common_variable_step_1():
    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thru the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})
    assert ind.index_to_species == ["A", "B"]

    initial_system = np.array([10, 50])
    sim = ReactionSimulator(system=initial_system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler", preset="fast")  # Note the preset

    # Start with a step neither too small nor too large (in the context of the "fast" preset)
    incr, step_taken, step_recommended = sim.reaction_step_common_variable_step(delta_time=0.02)
    assert math.isclose(step_taken, 0.02)           # Done as suggested
    assert math.isclose(step_recommended, 0.02)     # "Stay on course"
    assert np.allclose(incr, [1.4, -1.4])

    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.02)
    assert sim.reaction_step_diagnostics.decision_data == {'action': 'stay', 'operation': 'stay',
                                                           'step_factor': 1,
                                                           'applicable_norms': 'ALL'}
    assert compare_dicts(sim.reaction_step_diagnostics.norms, {'norm_A': 0.98, 'norm_B': 0.14})
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # New simulation, with somewhat large time step (but not so large as to cause hard errors)
    initial_system = np.array([10, 50])
    sim = ReactionSimulator(system=initial_system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler", preset="fast")  # Note the preset

    incr, step_taken, step_recommended = sim.reaction_step_common_variable_step(delta_time=0.04)
    #print(incr, step_taken, step_recommended)
    assert math.isclose(step_taken, 0.024)          # Smaller than suggested [but 20% longer than in our previous test run, above]
    assert math.isclose(step_recommended, 0.0192)   # "Go even smaller" in next round (80% of current value)
    assert np.allclose(incr, [1.68, -1.68])

    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.024)
    assert sim.reaction_step_diagnostics.decision_data == {'action': 'high', 'operation': 'downshift',
                                                           'step_factor': 0.8,
                                                           'applicable_norms': ['norm_A']}
    assert compare_dicts(sim.reaction_step_diagnostics.norms, {'norm_A': 1.4112, 'norm_B': 0.168})
    assert math.isclose(step_recommended, step_taken * sim.reaction_step_diagnostics.decision_data['step_factor'])
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})



    # New simulation, with excessively small time step
    initial_system = np.array([10, 50])
    sim = ReactionSimulator(system=initial_system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler", preset="fast")    # Note the preset

    incr, step_taken, step_recommended = sim.reaction_step_common_variable_step(delta_time=0.01)    # , explain_variable_steps=[-1,1]
    #print(incr, step_taken, step_recommended)
    assert math.isclose(step_taken, 0.01)           # Done as suggested
    assert math.isclose(step_recommended, .015)     # "Go larger smaller" in next round (50% more of current value)
    assert np.allclose(incr, [0.7, -0.7])

    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.01)
    assert sim.reaction_step_diagnostics.decision_data == {'action': 'low', 'operation': 'upshift',
                                                           'step_factor': 1.5,
                                                           'applicable_norms': 'ALL'}
    assert compare_dicts(sim.reaction_step_diagnostics.norms, {'norm_A': 0.245, 'norm_B':0.07})
    assert math.isclose(step_recommended, step_taken * sim.reaction_step_diagnostics.decision_data['step_factor'])
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})



def test_reaction_step_common_variable_step_2(capsys):
    # Capture informational output in the case of a variable step that gets "rewound" and re-done
    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})
    assert ind.index_to_species == ["A", "B"]

    # New simulation, with slightly larger time step
    initial_system = np.array([10, 50])
    sim = ReactionSimulator(system=initial_system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler", preset="fast")

    sim.system.set_system_state(initial_system)
    incr, step_taken, step_recommended = sim.reaction_step_common_variable_step(delta_time=0.04, explain_variable_steps=[-1,1])
    captured = capsys.readouterr()  # Capture the standard output and standard error, from the previous function call

    assert math.isclose(step_taken, 0.024)          # Smaller than suggested [but 20% longer than in our previous test run, above]
    assert math.isclose(step_recommended, 0.0192)   # "Go even smaller" in nest round (80% of current value)
    assert np.allclose(incr, [1.68, -1.68])

    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.024)
    assert sim.reaction_step_diagnostics.decision_data == {'action': 'high', 'operation': 'downshift',
                                                           'step_factor': 0.8,
                                                           'applicable_norms': ['norm_A']}
    assert compare_dicts(sim.reaction_step_diagnostics.norms, {'norm_A': 1.4112, 'norm_B': 0.168})
    assert math.isclose(step_recommended, step_taken * sim.reaction_step_diagnostics.decision_data['step_factor'])
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})

    output_expected = """
(STEP 1 aborted) SYSTEM TIME 0 : Examining Conc. changes due to tentative Δt=0.04 ...
    Previous:  None
    Baseline:  [10 50]
    Deltas:    [ 2.8 -2.8]
    Norms:     { 'norm_A': 3.92 }
    Thresholds:    
                   norm_A : low 0.8 | high 1.2 | abort 1.7 | (VALUE 3.92)
                   norm_B :  (skipped; not needed)
    Step Factors:     {'upshift': 1.5, 'downshift': 0.8, 'abort': 0.6, 'error': 0.5}
    => Action: 'ABORT'  ('abort' with step size factor of 0.6)
       * INFO: the tentative time step (0.04) leads to a value of ['norm_A'] > its ABORT threshold:
       -> will backtrack, and re-do step with a SMALLER Δt, x0.6 (now set to 0.024) [Step started at t=0, and will rewind there]

(STEP 1 completed) SYSTEM TIME 0 : Examining Conc. changes due to tentative Δt=0.024 ...
    Previous:  None
    Baseline:  [10 50]
    Deltas:    [ 1.68 -1.68]
    Norms:     { 'norm_A': 1.4112, 'norm_B': 0.168 }
    Thresholds:    
                   norm_A : low 0.8 | high 1.2 | (VALUE 1.4112) | abort 1.7
                   norm_B : low 0.15 | (VALUE 0.168) | high 0.8 | abort 1.8
    Step Factors:     {'upshift': 1.5, 'downshift': 0.8, 'abort': 0.6, 'error': 0.5}
    => Action: 'HIGH'  ('downshift' with step size factor of 0.8)
       INFO: COMPLETED STEP NORMALLY and MADE INTERVAL SMALLER, multiplied by 0.8 (set to 0.0192) at the next round, because ['norm_A'] is high
    [The current step started at System Time: 0 , and will continue to 0.024]
"""
    assert output_expected == captured.out



def test_attempt_reaction_step(capsys):
    # NO diagnostics.  Instead, capture informational output

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})

    system = np.array([10, 50])
    sim = ReactionSimulator(system=system, species_index_map=ind, reaction_registry=rxns, method="forward_euler")

    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.1, variable_steps=False) # FIXED steps
    assert np.allclose(delta_conc, [7, -7])
    assert math.isclose(rec_next_step, 0.1)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Reset and re-run (this time with "heun")
    sim.system.set_system_state([10, 50])
    sim.method = "heun"
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.1, variable_steps=False)
    assert np.allclose(delta_conc, [5.25, -5.25])
    assert math.isclose(rec_next_step, 0.1)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Reset and re-run (back to "forward_euler", but this time with variable step, and a smaller step)
    sim.system.set_system_state([10, 50])
    sim.method = "forward_euler"
    sim.adaptive_steps.use_adaptive_preset(preset="fast")
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.02, variable_steps=True)
    assert np.allclose(delta_conc, [1.4, -1.4])
    assert math.isclose(rec_next_step, 0.02)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Reset and re-run (this time with a printed explanation of variable time steps)
    sim.system.set_system_state([10, 50])
    sim.method = "forward_euler"
    sim.adaptive_steps.use_adaptive_preset(preset="fast")
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.02, variable_steps=True,
                                                          explain_variable_steps=(-1, 1), step_counter=1)

    captured = capsys.readouterr()  # Capture the standard output and standard error, from the previous function call
    assert "(STEP 1 completed) SYSTEM TIME 0 : Examining Conc. changes due to tentative Δt=0.02 ..." in captured.out
    assert "    Previous:  None" in captured.out
    assert "    Baseline:  [10 50]" in captured.out
    assert "    Deltas:    [ 1.4 -1.4]" in captured.out
    assert "    Norms:     { 'norm_A': 0.98, 'norm_B': 0.14 }" in captured.out
    assert "    Thresholds:" in captured.out
    assert "                   norm_A : low 0.8 | (VALUE 0.98) | high 1.2 | abort 1.7" in captured.out
    assert "                   norm_B : (VALUE 0.14) | low 0.15 | high 0.8 | abort 1.8" in captured.out
    assert "    => Action: 'STAY'  (with step size factor of 1)" in captured.out
    assert "       INFO: COMPLETED STEP NORMALLY - we're inside the target range of all norms.  No change to step size." in captured.out
    assert "    [The current step started at System Time: 0 , and will continue to 0.02]" in captured.out

    assert np.allclose(delta_conc, [1.4, -1.4])
    assert math.isclose(rec_next_step, 0.02)    # Stayed on course (no change)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Reset and re-run (again with a printed explanation of variable time steps)
    sim.system.set_system_state([10, 50])
    sim.method = "forward_euler"
    sim.adaptive_steps.use_adaptive_preset(preset="fast")
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.025, variable_steps=True,
                                                          explain_variable_steps=(-1, 1), step_counter=1)

    captured = capsys.readouterr()  # Capture the standard output and standard error, from the previous function call
    assert "(STEP 1 completed) SYSTEM TIME 0 : Examining Conc. changes due to tentative Δt=0.025 ..." in captured.out
    assert "    Previous:  None" in captured.out
    assert "    Baseline:  [10 50]" in captured.out
    assert "    Deltas:    [ 1.75 -1.75]" in captured.out
    assert "    Norms:     { 'norm_A': 1.5312, 'norm_B': 0.175 }" in captured.out
    assert "    Thresholds:" in captured.out
    assert "                   norm_A : low 0.8 | high 1.2 | (VALUE 1.5312) | abort 1.7" in captured.out
    assert "                   norm_B : low 0.15 | (VALUE 0.175) | high 0.8 | abort 1.8" in captured.out
    assert "    Step Factors:     {'upshift': 1.5, 'downshift': 0.8, 'abort': 0.6, 'error': 0.5}" in captured.out
    assert "    => Action: 'HIGH'  ('downshift' with step size factor of 0.8)" in captured.out
    assert "       INFO: COMPLETED STEP NORMALLY and MADE INTERVAL SMALLER, multiplied by 0.8 (set to 0.02) at the next round, because ['norm_A'] is high" in captured.out
    assert "    [The current step started at System Time: 0 , and will continue to 0.025]" in captured.out

    assert np.allclose(delta_conc, [1.75, -1.75])
    assert math.isclose(rec_next_step, 0.02)    # Slowed down!
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Reset and re-run (again with a printed explanation of variable time steps)
    sim.system.set_system_state([10, 50])
    sim.method = "forward_euler"
    sim.adaptive_steps.use_adaptive_preset(preset="fast")
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.01, variable_steps=True,
                                                          explain_variable_steps=(-1, 1), step_counter=1)

    captured = capsys.readouterr()  # Capture the standard output and standard error, from the previous function call
    assert "(STEP 1 completed) SYSTEM TIME 0 : Examining Conc. changes due to tentative Δt=0.01 ..." in captured.out
    assert "    Previous:  None" in captured.out
    assert "    Baseline:  [10 50]" in captured.out
    assert "    Deltas:    [ 0.7 -0.7]" in captured.out
    assert "    Norms:     { 'norm_A': 0.245, 'norm_B': 0.07 }" in captured.out
    assert "    Thresholds:" in captured.out
    assert "                   norm_A : (VALUE 0.245) | low 0.8 | high 1.2 | abort 1.7" in captured.out
    assert "                   norm_B : (VALUE 0.07) | low 0.15 | high 0.8 | abort 1.8" in captured.out
    assert "    Step Factors:     {'upshift': 1.5, 'downshift': 0.8, 'abort': 0.6, 'error': 0.5}" in captured.out
    assert "    => Action: 'LOW'  ('upshift' with step size factor of 1.5)" in captured.out
    assert "       INFO: COMPLETED STEP NORMALLY and MADE INTERVAL LARGER, multiplied by 1.5 (set to 0.015) at the next round, because all norms are low" in captured.out
    assert "    [The current step started at System Time: 0 , and will continue to 0.01]" in captured.out

    assert np.allclose(delta_conc, [0.7, -0.7])
    assert math.isclose(rec_next_step, 0.015)    # Speed up (larger step)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})



def test_attempt_reaction_step_2_a():
    # With DIAGNOSTICS enabled, for FIXED reaction step

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})

    system = np.array([10, 50])
    sim = ReactionSimulator(system=system, species_index_map=ind, reaction_registry=rxns, method="forward_euler",
                            enable_diagnostics=True)   # Diagnostics are enabled at instantiation time
    sim.system_time = 88
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.1, variable_steps=False)     # FIXED time step
    assert np.allclose(delta_conc, [7, -7])
    assert math.isclose(rec_next_step, 0.1)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Verify the diagnostic data: part 1 - the "diagnostic_rxn_data", created by single_step_single_rxn()
    assert type(sim.diagnostics.rxn_data.diagnostic_rxn_data) is dict
    assert len(sim.diagnostics.rxn_data.diagnostic_rxn_data) == 1
    coll_tab = sim.diagnostics.rxn_data.diagnostic_rxn_data[0]
    assert type(coll_tab) is CollectionTabular
    df = coll_tab.get_dataframe()
    row_expected = {"START_TIME": 88, "time_step": 0.1, "aborted": False, "Delta A": 7.0, "Delta B": -7.0, "rate": -70., "caption": ""}
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert compare_pandas(df_expected, df, disregard_order=True)

    # Verify the diagnostic data: part 2 - the "diagnostic_decisions_data", created by attempt_reaction_step()
    assert type(sim.diagnostics.step_decision_data.diagnostic_decisions_data) is CollectionTabular
    df = sim.diagnostics.step_decision_data.diagnostic_decisions_data.get_dataframe()
    assert type(df) is pd.DataFrame
    assert len(df) == 1
    row_expected = {"START_TIME": 88, "Delta A": 7.0, "Delta B": -7.0, "caption": ""}
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert compare_pandas(df_expected, df, disregard_order=True)



def test_attempt_reaction_step_2_b():
    # With DIAGNOSTICS enabled, for VARIABLE reaction step

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})

    system = np.array([10, 50])
    sim = ReactionSimulator(system=system, species_index_map=ind, reaction_registry=rxns, method="forward_euler",
                            enable_diagnostics=True)   # Diagnostics are enabled at instantiation time
    sim.system_time = 99
    sim.adaptive_steps.use_adaptive_preset(preset="fast")
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.02, variable_steps=True)     # VARIABLE time step
    assert np.allclose(delta_conc, [1.4, -1.4])
    assert math.isclose(rec_next_step, 0.02)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Verify the diagnostic data: part 1 - the "diagnostic_rxn_data", created by single_step_single_rxn()
    assert type(sim.diagnostics.rxn_data.diagnostic_rxn_data) is dict
    assert len(sim.diagnostics.rxn_data.diagnostic_rxn_data) == 1
    coll_tab = sim.diagnostics.rxn_data.diagnostic_rxn_data[0]
    assert type(coll_tab) is CollectionTabular
    df = coll_tab.get_dataframe()
    row = {"START_TIME": 99, "time_step": 0.02, "aborted": False, "Delta A": 1.4, "Delta B": -1.4, "rate": -70., "caption": ""}
    df_expected = pd.DataFrame(row, index=[0])
    assert_frame_equal(df, df_expected)


    # Verify the diagnostic data: part 2 - the "diagnostic_decisions_data", created by attempt_reaction_step(), this time with
    #   extra fields resulting from the VARIABLE step
    assert type(sim.diagnostics.step_decision_data.diagnostic_decisions_data) is CollectionTabular
    df = sim.diagnostics.step_decision_data.diagnostic_decisions_data.get_dataframe()
    assert type(df) is pd.DataFrame
    assert len(df) == 1
    # More fields are present now, because of the variable step
    row = {"START_TIME": 99, "Delta A": 1.4, "Delta B": -1.4,
           "norm_A": 0.98, "norm_B": 0.14, "norm_C": None, "norm_D": None, "action": "OK (stay)", "step_factor": 1, "time_step": 0.02,
           "caption": ""}
    df_expected = pd.DataFrame(row, index=[0])
    assert_frame_equal(df, df_expected)



def test_single_step_all_rxns_1():
    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_registry=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    sim = ReactionSimulator(system=[10, 50], reaction_registry=rxns)
    assert np.allclose(sim.system.system_state , [10., 50.])
    assert sim.species_index_map.index_to_species == ["A", "B"]
    assert sim.species_index_map.species_to_index == {"A": 0, "B": 1}

    result = sim.single_step_all_rxns(delta_time=0.1)
    assert np.allclose(result, [7, -7])
    assert result[0] == - result[1]         # From the stoichiometry
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Reset and re-run
    sim.system.set_system_state([10, 50])
    sim.method = "heun"
    result = sim.single_step_all_rxns(delta_time=0.1)
    assert np.allclose(result, [5.25, -5.25])
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    rxns.clear_reactions_data()       # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , hypothetically with 1st-order kinetics for each species
    # Based on experiment "1D/reactions/reaction2"
    #rxns = ReactionRegistry(species_data=species_registry)
    #sim.reaction_registry = rxns
    rxns.add_reaction(reactants="A", products=[(3,"B")],
                      reaction_model="custom",
                      kinetic_parameters={"kF": 5., "kR": 2., "rate_function": Custom_Model.kinetic_rate_first_order})
    sim.method = "forward_euler"
    result = sim.single_step_all_rxns(delta_time=0.1)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]      # From the stoichiometry


    rxns.clear_reactions_data()       # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B ,hypothetically with 1st-order kinetics for each species
    # Based on experiment "1D/reactions/reaction3"
    rxns.add_reaction(reactants=[(2,"A")], products=[(3,"B")],
                      reaction_model="custom",
                      kinetic_parameters={"kF": 5., "kR": 2., "rate_function": Custom_Model.kinetic_rate_first_order})
    sim.method = "forward_euler"
    result = sim.single_step_all_rxns(delta_time=0.1)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test_single_step_all_rxns_2():
    species_registry = SpeciesRegistry(["A", "B", "C"])
    rxns = ReactionRegistry(species_registry=species_registry)

    # Unimolecular elementary reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    rxns.add_reaction(reactants="A", products="B",
                      reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})

    sim = ReactionSimulator(system=[10, 50], reaction_registry=rxns, method="forward_euler")
    assert sim.species_index_map.index_to_species == ["A", "B"]         # Species "C" not participating in this reaction
    assert sim.species_index_map.species_to_index == {"A": 0, "B": 1}

    result = sim.single_step_all_rxns(delta_time=0.1)
    assert np.allclose(result, [ 7. , -7.])     # Species "C" not participating in this reaction; not present in system state
    assert result[0] == - result[1]             # From the stoichiometry


    rxns.clear_reactions_data()       # Re-start with a blank slate of reactions
    # Synthesis reaction A + B <-> C , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction4"
    rxns.add_reaction(reactants=["A" , "B"], products="C",
                      reaction_model="mass action", kinetic_parameters={"kF": 5., "kR": 2.})

    sim.system.set_conc({"A": 10, "B": 50, "C": 20})
    # The species indexing got expanded, to accommodate `C`
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C']                     # Notice the reaction-wise sorting
    assert sim.species_index_map.species_to_index == {'A': 0, 'B': 1, 'C': 2}

    result = sim.single_step_all_rxns(delta_time=0.002)
    assert np.allclose(result, [-4.92, -4.92, 4.92])
    assert result[0] == result[1]           # From the stoichiometry
    assert result[1] == - result[2]         # From the stoichiometry


def test_single_step_all_rxns_3():
    species_registry = SpeciesRegistry()
    rxns = ReactionRegistry(species_registry=species_registry)

    # Reaction A <-> 2C + D , HYPOTHETICALLY with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction5"
    rxns.add_reaction(reactants=["A"], products=[(2, "C") , ("D")],
                    reaction_model="custom",
                    kinetic_parameters={"kF": 5., "kR": 2., "rate_function": Custom_Model.kinetic_rate_first_order})

    sim = ReactionSimulator(system=[4., 7., 2.], reaction_registry=rxns, method="forward_euler")

    result = sim.single_step_all_rxns(delta_time=0.05)
    assert np.allclose(result, [0.4 , -0.8 , -0.4])
    assert result[0] == - result[1] /2    # From the stoichiometry
    assert result[0] == - result[2]       # From the stoichiometry


def test_single_step_all_rxns_4():
    species_registry = SpeciesRegistry()
    rxns = ReactionRegistry(species_registry=species_registry)

    # Reaction 2A + 5B <-> 4C + 3D , HYPOTHETICALLY with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction6"
    rxns.add_reaction(reactants=[(2,"A") , (5,"B")], products=[(4,"C") , (3,"D")],
                    reaction_model="custom",
                    kinetic_parameters={"kF": 5., "kR": 2., "rate_function": Custom_Model.kinetic_rate_first_order})

    sim = ReactionSimulator(system=[4., 7., 5., 2.], reaction_registry=rxns, method="forward_euler")

    result = sim.single_step_all_rxns(delta_time=0.001)
    assert np.allclose(result, [-0.24 , -0.6 , 0.48, 0.36])
    assert  np.allclose(result[0] /2 , result[1] /5)    # From the stoichiometry
    assert  np.allclose(result[1] /5 , -result[2] /4)   # From the stoichiometry
    assert  np.allclose(result[2] /4 , result[3] /3)    # From the stoichiometry


def test_single_step_all_rxns_5():
    species_registry = SpeciesRegistry()
    rxns = ReactionRegistry(species_registry=species_registry)

    # Reaction  2A <-> B , with mass-action kinetics
    # Based on experiment "1D/reactions/reaction7"
    rxns.add_reaction(reactants=[(2, "A")], products="B",
                    reaction_model="mass action", kinetic_parameters={"kF": 5., "kR": 2.})

    sim = ReactionSimulator(system=[3., 5.], reaction_registry=rxns, method="forward_euler")

    result = sim.single_step_all_rxns(delta_time=0.02)
    assert np.allclose(result, [-1.4 , 0.7])
    assert np.allclose(result[0] /2 , -result[1])      # From the stoichiometry


def test_single_step_all_rxns_6():
    species_registry = SpeciesRegistry()
    rxns = ReactionRegistry(species_registry=species_registry)

    # Coupled reactions A + B <-> C  and  C + D <-> E , each with mass-action kinetics
    # Based on experiment "1D/reactions/reaction8"
    rxns.add_reaction(reactants=["A", "B"], products="C",
                      reaction_model="mass action", kinetic_parameters={"kF": 5., "kR": 2.})

    rxns.add_reaction(reactants=["C", "D"], products="E",
                      reaction_model="mass action", kinetic_parameters={"kF": 8., "kR": 4.})

    assert rxns.number_of_reactions() == 2

    sim = ReactionSimulator(system=[3., 5., 1., 0.4, 0.1], reaction_registry=rxns, method="forward_euler")

    result = sim.single_step_all_rxns(delta_time=0.02)
    assert np.allclose(result, [-1.46 , -1.46  , 1.404 , -0.056 ,  0.056])
    assert np.allclose(result[0] , result[1])                   # From the stoichiometry
    assert np.allclose(result[3] , -result[4])                  # From the stoichiometry
    assert np.allclose(result[0] + result[4], -result[2])       # From the stoichiometry
                                                                # The increase in [A] and [E] combined
                                                                # must match the decrease in [C]



def test_single_step_single_rxn_1():

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created as an independent object)
    rxn_defn = ReactionDefinition(reactants="A", products="B",
                                  species_registry=species_registry, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})
    rxn_sim = rxn_defn.sim_reactions[0]
    assert rxn_sim.analytic_solution_family == "ONE_TO_ONE"

    ind = SpeciesIndexMap()
    species_id_set =  rxn_sim.stoichiometry.get_all_species_ids()
    ind.add_species(species_id_set)

    system = np.array([10, 50])
    sim = ReactionSimulator(method="forward_euler", system=system, species_index_map=ind)
    assert sim.species_index_map.index_to_species == ['A', 'B']
    assert sim.species_index_map.species_to_index == {'A' :0, 'B' :1}
    assert np.allclose(sim.system.system_state, [10, 50])

    increment_vector = np.zeros(2, dtype='d')
    result = sim.single_step_single_rxn(increment_vector=increment_vector, delta_time=0.1,
                                        rxn=rxn_sim, rxn_index=0)
    assert np.allclose(increment_vector, [7, -7])
    assert math.isclose(result, -70)


    system = np.array([10, 50])
    sim = ReactionSimulator(method="heun", system=system, species_index_map=ind)

    increment_vector = np.zeros(2, dtype='d')
    result = sim.single_step_single_rxn(increment_vector=increment_vector, delta_time=0.1,
                                        rxn=rxn_sim, rxn_index=0)
    assert np.allclose(increment_vector, [5.25, -5.25])
    assert math.isclose(result, -70)


    system = np.array([10, 50])
    sim = ReactionSimulator(method="analytic",
                            system=system, species_index_map=ind)

    increment_vector = np.zeros(2, dtype='d')
    result = sim.single_step_single_rxn(increment_vector=increment_vector, delta_time=0.1,
                                        rxn=rxn_sim, rxn_index=0)
    assert np.allclose(increment_vector, [5.508570764023133, -5.508570764023133])
    assert math.isclose(result, -70)


def test_single_step_single_rxn_2():
    # Here we also test the diagnostics

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B
    rxn_defn = ReactionDefinition(reactants="A", products="B",
                                  species_registry=species_registry, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})
    rxn_sim = rxn_defn.sim_reactions[0]

    reaction_registry = ReactionRegistry(species_registry=species_registry)
    reaction_registry.register_reaction(rxn_defn)

    ind = SpeciesIndexMap()
    species_id_set =  rxn_sim.stoichiometry.get_all_species_ids()
    ind.add_species(species_id_set)

    diagnostics = Diagnostics(reactions=reaction_registry, species_index_map=ind)

    system = np.array([10, 50])
    sim = ReactionSimulator(method="forward_euler", system=system, species_index_map=ind,
                            enable_diagnostics=True, diagnostics=diagnostics)

    increment_vector = np.zeros(2, dtype='d')
    sim.system_time = 123
    result = sim.single_step_single_rxn(increment_vector=increment_vector, delta_time=0.1,
                                        rxn=rxn_sim, rxn_index=0)
    assert np.allclose(increment_vector, [7, -7])
    assert math.isclose(result, -70)

    assert len(sim.diagnostics.rxn_data.diagnostic_rxn_data) == 1

    collection = sim.diagnostics.rxn_data.diagnostic_rxn_data[0]
    assert type(collection) is CollectionTabular
    df = collection.get_dataframe()
    assert type(df) is pd.DataFrame

    row = {"START_TIME": 123, "time_step": 0.1, "aborted": False, "Delta A": 7.0, "Delta B": -7.0, "rate": -70.0, "caption": ""}
    df_expected = pd.DataFrame(row, index=[0])

    assert compare_pandas(df_expected, df)



    # Re-start with a fresh Diagnostics object; this time we'll take an EXCESSIVE step size
    diagnostics = Diagnostics(reactions=reaction_registry, species_index_map=ind)

    system = np.array([10, 50])
    sim = ReactionSimulator(method="forward_euler", system=system, species_index_map=ind,
                            enable_diagnostics=True, diagnostics=diagnostics)

    increment_vector = np.zeros(2, dtype='d')
    sim.system_time = 666
    with pytest.raises(ExcessiveTimeStepHard) as ex:      # Excessive time step that would make [B] negative
        sim.single_step_single_rxn(increment_vector=increment_vector, delta_time=0.8,
                                  rxn=rxn_sim, rxn_index=0)
    details = ex.value.details
    assert details.get("function") == "_validate_increment"
    assert details.get("delta_time") == 0.8
    assert details.get("caption") == "aborted: neg. conc. in `B` from rxn # 0"
    assert details.get("system_time") == 666
    assert details.get("rate") == -70
    assert details.get("rxn_index") == 0


    # TODO: move the part commented out below to the higher layer that catches the Exception
    """
    #print(sim.diagnostics)
    # part 1: diagnostic_rxn_data
    assert type(sim.diagnostics.rxn_data.diagnostic_rxn_data) is dict
    assert len(sim.diagnostics.rxn_data.diagnostic_rxn_data) == 1
    collection = sim.diagnostics.rxn_data.diagnostic_rxn_data[0]
    assert type(collection) is CollectionTabular
    df = collection.get_dataframe()
    assert type(df) is pandas.DataFrame
    row = {"START_TIME": 666, "time_step": 0.8, "aborted": True, "Delta A": np.nan, "Delta B": np.nan,
           "rate": -70.0, "caption": "aborted: neg. conc. in `B`"}
    df_expected = pd.DataFrame(row, index=[0])
    assert compare_pandas(df_expected, df)

    # part 2: diagnostic_decisions_data
    assert len(sim.diagnostics.step_decision_data.diagnostic_decisions_data) == 1
    collection = sim.diagnostics.step_decision_data.diagnostic_decisions_data
    assert type(collection) is CollectionTabular
    df = collection.get_dataframe()
    assert type(df) is pandas.DataFrame
    row = {"START_TIME": 666, "action": "ABORT", "caption": "", "time_step": 0.8}
    df_expected = pd.DataFrame(row, index=[0])
    assert compare_pandas(df_expected, df)
    """



def test_validate_increment():
    species_registry = SpeciesRegistry()
    rxns = ReactionRegistry(species_registry=species_registry)

    # A <-> B
    rxns.add_reaction(reactants="A", products="B",
                    reaction_model="mass action", kinetic_parameters={"kF": 1., "kR": 0.2})

    sim = ReactionSimulator(system=np.array([19, 3]), reaction_registry=rxns)
    sim.reaction_step_diagnostics = ReactionStep(delta_time=0.97)

    r = rxns.get_reaction(0)

    sim._validate_increment(delta_conc=50., baseline_conc=10.,
                            rxn=r, rxn_index=0, rxn_rate=123, species_id="B")
    sim._validate_increment(delta_conc=-9.99, baseline_conc=10.,
                            rxn=r, rxn_index=0, rxn_rate=123, species_id="B")

    with pytest.raises(ExcessiveTimeStepHard) as ex:      # Would lead to a negative concentration
        sim._validate_increment(delta_conc=-10.001, baseline_conc=10.,
                               rxn=r, rxn_index=0, rxn_rate=123, species_id="B")
    details = ex.value.details
    assert details["function"] == "_validate_increment"
    assert math.isclose(details["delta_time"], 0.97)
    assert details["caption"] == 'aborted: neg. conc. in `B` from rxn # 0'
    assert details["system_time"] == 0
    assert details["rate"] == 123
    assert details["rxn_index"] == 0
    assert details["message"] == """      The attempted time step (0.97) would lead to a NEGATIVE concentration in the species `B` from the reaction `A <-> B` (rxn # 0)
      Baseline concentration value of `B` : 10 at system time 0; requested change (NOT carried out): -10.001, due to a reaction rate of 123"""

    # TODO: test that diagnostic data, if enabled, gets saved as needed



def test__fetch_concs_for_rxn():
    species_registry = SpeciesRegistry()
    rxns = ReactionRegistry(species_registry=species_registry)

    # D -> F
    rxns.add_reaction(reactants="D", products="F",
                      reaction_model="mass action")

    sim = ReactionSimulator(system=np.array([19, 3]), reaction_registry=rxns)
    assert sim.species_index_map.index_to_species == ['D', 'F']
    assert sim.species_index_map.species_to_index == {'D': 0, 'F': 1}

    r = rxns.get_reaction(0)         # D -> F
    result = sim._fetch_concs_for_rxn(rxn=r)
    assert result == {"D": 19, "F": 3}
    assert np.allclose(sim.system.system_state, np.array([19, 3]))


    # A + F -> C
    rxns.add_reaction(reactants=["A", "F"], products="C",
                      reaction_model="mass action")
    sim.system.set_conc({"A": 12, "C": 31})
    assert sim.species_index_map.index_to_species == ['D', 'F', 'A', 'C']   # Notice how 'A' and 'C' got added, in sorted order
    assert sim.species_index_map.species_to_index == {'D': 0, 'F': 1, 'A': 2, 'C': 3}
    assert np.allclose(sim.system.system_state, np.array([19, 3, 12, 31]))   # `A` got over-written; system array got expanded

    r = rxns.get_reaction(1)         # A + F -> C
    result = sim._fetch_concs_for_rxn(rxn=r)
    assert result == {"A": 12, "C": 31, "F": 3}


    # C -> A
    r_uni = ReactionDefinition(reactants="C", products="A",
                               species_registry=species_registry,
                               reaction_model="mass action")
    rxns.register_reaction(r_uni)
    assert sim.species_index_map.index_to_species == ['D', 'F', 'A', 'C']
    assert sim.species_index_map.species_to_index == {'D': 0, 'F': 1, 'A': 2, 'C': 3}
    assert np.allclose(sim.system.system_state, np.array([19, 3, 12, 31]))

    r = rxns.get_reaction(2)        # C -> A
    result = sim._fetch_concs_for_rxn(rxn=r)
    assert result == {"A": 12, "C": 31}


    # C + D -> B
    r_syn = ReactionDefinition(reactants=["C", "D"], products="B",
                               species_registry=species_registry,
                               reaction_model="mass action", autoregister_species=True)
    rxns.register_reaction(r_syn)
    assert np.allclose(sim.system.system_state, [19,  3, 12, 31])
    sim.system.set_conc({"B": 1})
    assert sim.species_index_map.index_to_species == ['D', 'F', 'A', 'C', 'B']
    assert sim.species_index_map.species_to_index == {'D': 0, 'F': 1, 'A': 2, 'C': 3, 'B': 4}
    assert np.allclose(sim.system.system_state, [19,  3, 12, 31, 1])  # The set_system_conc_by_dict() op expanded the system array

    r = rxns.get_reaction(3)        # C + D -> B
    result = sim._fetch_concs_for_rxn(rxn=r, )
    assert result == {"B": 1, "C": 31, "D": 19}



def test_dispatcher_single_rxn():

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B
    rxn_defn = ReactionDefinition(reactants="A", products="B",
                                  species_registry=species_registry, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})
    rxn_sim = rxn_defn.sim_reactions[0]
    assert rxn_sim.analytic_solution_family == "ONE_TO_ONE"

    sim = ReactionSimulator(method="forward_euler")
    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 10, "B": 50}, delta_time=0.1)
    assert compare_dicts(delta_conc, {'A': 7, 'B': -7})
    assert math.isclose(rate, -70)

    sim = ReactionSimulator(method="heun")
    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 10, "B": 50}, delta_time=0.1)
    assert compare_dicts(delta_conc, {'A': 5.25, 'B': -5.25})
    assert math.isclose(rate, -70)

    sim = ReactionSimulator(method="analytic")
    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 10, "B": 50}, delta_time=0.1)
    assert compare_dicts(delta_conc, {'A': 5.508570764023133, 'B': -5.508570764023133}) # A lot closer to Heun than to forward Euler!
    assert math.isclose(rate, -70)



def test_analytic_solver_single_rxn_1():

    sim = ReactionSimulator(method="analytic")

    # Reaction : A <-> B
    sim.add_reaction(reactants="A", products="B",
                     reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})
    rxn_sim = sim.reaction_registry.get_reaction(0)
    assert rxn_sim.model.reversible
    assert rxn_sim.analytic_solution_family == "ONE_TO_ONE"

    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 10, "B": 50}, delta_time=0.1)
    assert compare_dicts(delta_conc, {'A': 5.508570764023133, 'B': -5.508570764023133})
    assert math.isclose(rate, -70)      # 3. * 10 - 2. * 50.


    # Fresh start with a new reaction A -> B (irreversible)
    sim.reaction_registry.clear_reactions_data()
    sim.add_reaction(reactants="A", products="B",
                     reaction_model="mass action", kinetic_parameters={"kF": 3.})
    rxn_sim = sim.reaction_registry.get_reaction(0)
    assert not rxn_sim.model.reversible
    assert rxn_sim.analytic_solution_family == "ONE_TO_ONE"

    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 80, "B": 20}, delta_time=0.31739)
    assert compare_dicts(delta_conc, {'A': -49.127783572982025, 'B': 49.127783572982025})
    assert math.isclose(rate, 240)     # 3. * 80


def test_analytic_solver_single_rxn_2():

    sim = ReactionSimulator(method="analytic")

    # Reaction A + B <-> P
    sim.add_reaction(reactants=["A", "B"], products="P",
                     reaction_model="mass action", kinetic_parameters={"kF": 5., "kR": 2.})
    rxn_sim = sim.reaction_registry.get_reaction(0)
    assert rxn_sim.model.reversible
    assert rxn_sim.analytic_solution_family == "TWO_TO_ONE"

    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 10, "B": 50, "P":20}, delta_time=0.0001)
    assert compare_dicts(delta_conc, {'A': -0.24233229989273042, 'B': -0.24233229989273042, 'P': 0.24233229989273042})
    assert math.isclose(rate, 2460)     # 5. * 10 * 50 - 2. * 20


    # Fresh start with a new reaction A + B -> P  (irreversible)
    sim.reaction_registry.clear_reactions_data()
    sim.add_reaction(reactants=["A", "B"], products="P",
                     reaction_model="mass action", kinetic_parameters={"kF": 5.})
    rxn_sim = sim.reaction_registry.get_reaction(0)
    assert not rxn_sim.model.reversible
    assert rxn_sim.analytic_solution_family == "TWO_TO_ONE"

    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 10, "B": 50, "P":20}, delta_time=0.01)
    assert compare_dicts(delta_conc, {'A': -8.8871974442945, 'B': -8.8871974442945, 'P': 8.8871974442945})
    assert math.isclose(rate, 2500)     # 5. * 10 * 50 - 2. * 20


    # Fresh start with a new reaction P <-> A + B  (simply the same as the previous reversible reaction A + B <-> P , as long as kF and kR are reversed)
    sim.reaction_registry.clear_reactions_data()
    sim.add_reaction(reactants="P", products=["A", "B"],
                     reaction_model="mass action", kinetic_parameters={"kF": 2., "kR": 5.}) # Notice the reversal of kF and kR
    rxn_sim = sim.reaction_registry.get_reaction(0)
    assert rxn_sim.model.reversible
    assert rxn_sim.analytic_solution_family == "ONE_TO_TWO"

    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 10, "B": 50, "P":20}, delta_time=0.0001)
    assert compare_dicts(delta_conc, {'A': -0.24233229989273042, 'B': -0.24233229989273042, 'P': 0.24233229989273042})  # Same as before
    assert math.isclose(rate, -2460)    # The opposite of what was in the reversed reaction


    # Fresh start with a new reaction R -> P + Q  (irreversible)
    sim.reaction_registry.clear_reactions_data()
    sim.add_reaction(reactants="R", products=["P", "Q"],
                     reaction_model="mass action", kinetic_parameters={"kF": 3.})
    rxn_sim = sim.reaction_registry.get_reaction(0)
    assert not rxn_sim.model.reversible
    assert rxn_sim.analytic_solution_family == "ONE_TO_TWO"

    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"R": 80, "P": 20, "Q":100}, delta_time=0.31739)
    assert compare_dicts(delta_conc, {'R': -49.127783572982025, 'P': 49.127783572982025, 'Q': 49.127783572982025})
    assert math.isclose(rate, 240)        # 3. * 80


def test_analytic_solver_single_rxn_3():
    sim = ReactionSimulator(method="analytic")

    # Reaction 2 A -> P  (irreversible)
    sim.add_reaction(reactants=["A", "A"], products="P",
                     reaction_model="mass action", kinetic_parameters={"kF": 5.})
    rxn_sim = sim.reaction_registry.get_reaction(0)
    assert not rxn_sim.model.reversible
    assert rxn_sim.analytic_solution_family == "TWO_TO_ONE"

    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 30, "P":20}, delta_time=0.01)
    assert compare_dicts(delta_conc, {'A': -18., 'P': 18.})
    assert math.isclose(rate, 4500)     # 5. * 30^2



def test_forward_euler_single_rxn_1():
    sr = SpeciesRegistry(ids=["A", "B"])

    # Reaction : A <-> B
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                                  reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})
    assert len(rxn_defn.sim_reactions) == 1
    sim_rxn = rxn_defn.sim_reactions[0]

    result = ReactionSimulator.forward_euler_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.1)
    assert result[0] == {'A': 7, 'B': -7}       # Reactant [A] is increasing, and product [B] is decreasing
    assert result[1] == -70     # Rate = 3. * 10. - 2. * 50 .  Reaction is progressing in the reverse direction

    result = ReactionSimulator.forward_euler_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.8)
    assert result[0] == {'A': 56, 'B': -56}              # Note: these increments would make [B] negative!
    assert result[1] == -70


    # Reaction : A -> B  (no reverse reaction)
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                                  reaction_model="mass action", kinetic_parameters={"kF": 3.})
    sim_rxn = rxn_defn.sim_reactions[0]

    result = ReactionSimulator.forward_euler_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.1)
    assert result[0] == {'A': -3, 'B': 3}
    assert result[1] == 30    # Rate = 3. * 10.     Reaction is now forward

    result = ReactionSimulator.forward_euler_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.4)
    assert result[0] == {'A': -12, 'B': 12}         # Note: these increments would make [A] negative!
    assert result[1] == 30


    # Reaction : A + B <-> C
    rxn_defn = ReactionDefinition(reactants=["A" , "B"], products="C", species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 5., "kR": 2})
    sim_rxn = rxn_defn.sim_reactions[0]

    C_0 = {"A": 10, "B": 50, "C": 20}
    delta_time = 0.002
    delta, rate = ReactionSimulator.forward_euler_single_rxn(rxn=sim_rxn, conc_init=C_0, delta_time=delta_time)
    assert delta == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert rate == 5 * 10 * 50 - 2 * 20        # 2460 , i.e.  kF [A] [B] - kR [C]
    for k in delta:     # Loop over the keys
        assert delta[k] == rate * 0.002 * np.sign(sim_rxn.stoichiometry.vector[k])


    # Reaction : C <-> A + B
    rxn_defn = ReactionDefinition(reactants="C", products=["A" , "B"], species_registry=sr,
                                  reaction_model="mass action", kinetic_parameters={"kF": 2., "kR": 5.})
    sim_rxn = rxn_defn.sim_reactions[0]

    result = ReactionSimulator.forward_euler_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50, "C": 20}, delta_time=0.002)
    assert result[0] == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert result[1] == -2460


def test_forward_euler_single_rxn_2():

    # Reaction: # E + S <-> ES* -> E + P, with SingleSubstrateMechanism model
    sr = SpeciesRegistry(ids=["S", "P", "E"])
    rxn_defn = ReactionDefinition(id=85, reactants=["S", "E"], products=["P", "E"], species_registry=sr,
                             reaction_model="single substrate mechanism",
                             kinetic_parameters={"k1_F": 18, "k1_R": 100, "k2_F": 49})

    initial_conc = {"E": 1, "S": 20, "P": 0, "ES*": 0}
    dt = 0.002

    # Look at each derived reaction individually
    rxn_sim_tuple = rxn_defn.sim_reactions
    assert len(rxn_sim_tuple) == 2
    rxn_sim_1, rxn_sim_2 = rxn_sim_tuple

    incr_dict_1, rate_1 = ReactionSimulator.forward_euler_single_rxn(rxn=rxn_sim_1, conc_init=initial_conc, delta_time=dt)
    assert rate_1 == 360       # 18 * 1 * 20 - 100 * 0
    assert incr_dict_1 == {'E': -0.72, 'S': -0.72, 'ES*': 0.72}


    incr_dict_2, rate_2 = ReactionSimulator.forward_euler_single_rxn(rxn=rxn_sim_2, conc_init=initial_conc, delta_time=dt)
    assert rate_2 == 0          # 49 * 0
    assert incr_dict_2 == {'E': 0, 'ES*': 0, 'P': 0}


    # Simulating individually each of the 2 sub-reactions (i.e. manually expanding the original enzymatic reaction) will give the same results
    upstream_rxn_defn = ReactionDefinition(reactants=["E", "S"], products="ES*", species_registry=sr,
                             reaction_model="mass action",
                             kinetic_parameters={"kF": 18, "kR": 100})
    upstream_rxn_sim = upstream_rxn_defn.sim_reactions[0]
    assert ReactionSimulator.forward_euler_single_rxn(rxn=upstream_rxn_sim, conc_init=initial_conc, delta_time=dt) \
                    == ( {'E': -0.72, 'S': -0.72, 'ES*': 0.72} , 360 )

    downstream_rxn_defn = ReactionDefinition(reactants="ES*", products=["E", "P"], species_registry=sr,
                             reaction_model="mass action",
                             kinetic_parameters={"kF": 49})
    downstream_rxn_sim = downstream_rxn_defn.sim_reactions[0]
    assert ReactionSimulator.forward_euler_single_rxn(rxn=downstream_rxn_sim, conc_init=initial_conc,delta_time=dt) \
                    == ( {'E': 0, 'ES*': 0, 'P': 0} , 0 )



    # Manually advance the system state by the previous time step, plus one more
    conc = initial_conc

    # Update the system concentrations (thus advancing the simulation)
    update_concentrations(conc, incr_dict_1)
    update_concentrations(conc, incr_dict_2)

    assert conc == {'E': 0.28, 'S': 19.28, 'P': 0.0, 'ES*': 0.72}

    incr_dict_1, rate_1 = ReactionSimulator.forward_euler_single_rxn(rxn=upstream_rxn_sim, conc_init=conc, delta_time=dt)
    assert math.isclose(rate_1, 25.1712)    # 18 * 0.28 * 19.28 - 100 * 0.72
    expected_incr = {'E': -0.0503424, 'S': -0.0503424, 'ES*': 0.0503424}    # Delta_conc = 25.1712 * 0.002 = 0.0503424
    compare_dicts(incr_dict_1, expected_incr)

    incr_dict_2, rate_2 = ReactionSimulator.forward_euler_single_rxn(rxn=downstream_rxn_sim, conc_init=initial_conc, delta_time=dt)
    assert math.isclose(rate_2, 35.28)      # 49 * 0.72
    expected_incr = {'ES*': -0.07056, 'E': 0.07056, 'P': 0.07056}           # Delta_conc = 35.28 * 0.002 = 0.07056
    compare_dicts(incr_dict_2, expected_incr)

    # Update the system concentrations (thus advancing the simulation)
    update_concentrations(conc, incr_dict_1)
    update_concentrations(conc, incr_dict_2)

    expected = {'E': 0.3002176, 'S': 19.2296576, 'P': 0.07056, 'ES*': 0.6997824}
    compare_dicts(conc, expected)


def test_forward_euler_single_rxn_3():
    # Reaction: # A + B -> C + D , with custom reaction model
    sr = SpeciesRegistry(ids=["A", "B", "C", "D"])
    rxn_defn = ReactionDefinition(id=49, reactants=["A", "B"], products=["C", "D"], species_registry=sr,
                             reaction_model="custom",
                             kinetic_parameters={"kF": 10, "rate_function": Custom_Model.kinetic_rate_first_order})
    #print(rxn_defn.describe(concise=False))

    sim_rxn = rxn_defn.sim_reactions[0]

    initial_conc = {"A": 2, "B": 4, "C": 5, "D": 3}
    dt = 0.1

    incr_dict, rate = ReactionSimulator.forward_euler_single_rxn(rxn=sim_rxn, conc_init=initial_conc, delta_time=dt)
    assert rate == 80       #  10. * 2 * 4   (no reverse reaction)
    assert incr_dict == {'A': -8, 'B': -8, 'C': 8, 'D': 8}      # 80 * 0.1 = 8

    # Make reversible
    sim_rxn.set_parameters({"kR": 2})
    incr_dict, rate = ReactionSimulator.forward_euler_single_rxn(rxn=sim_rxn, conc_init=initial_conc, delta_time=dt)
    assert rate == 50      # 80 - 2. * 5 * 3
    assert incr_dict == {'A': -5, 'B': -5, 'C': 5, 'D': 5}      # 50 * 0.1 = 5



def test_heun_single_rxn_1():
    sr = SpeciesRegistry(ids=["A", "B"])

    # Reaction : A <-> B
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                                  reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})
    assert len(rxn_defn.sim_reactions) == 1
    sim_rxn = rxn_defn.sim_reactions[0]

    result = ReactionSimulator.heun_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.1)
    """
    As seen in test_forward_euler_single_rxn_1():
        Start rate = -70
        Delta conc = {'A': 7, 'B': -7}
        So, Euler-based final conc:  {'A': 17, 'B': 43}
        Then: final rate = 3. * 17. - 2. * 43 = -35
        and: average rate = (-70 - 35)/2 = -52.5
    """
    assert result[0] == {'A': 52.5*0.1, 'B': -52.5*0.1}
    assert result[1] == -70     # Start Rate


    with pytest.raises(ExcessiveTimeStepHard) as ex:      # Excessive time step that would make [B], as computed by forward Euler, negative
        ReactionSimulator.heun_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.8)

    details = ex.value.details
    assert details.get("function") == "heun_single_rxn"
    assert details.get("delta_time") == 0.8
    assert details.get("rate") == -70
    assert details.get("message") == "heun_single_rxn(): excessive time step (0.8), " \
                                     "leading to negative concentrations"

    with pytest.raises(ExcessiveTimeStepHard) as ex:
        # Excessive time step (0.7), leading to final rate of 175 which, when averaged with the initial rate of -70 would flip the reaction's direction
        # Euler final_conc:  {'A': 59, 'B': 1}   [A] = 10 + (-1) * -70 * 0.7   ;  [A] = 50 + 1 * -70 * 0.7
        # Final rate =  3. * 59 - 2. * 1 = 175
        ReactionSimulator.heun_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.7)

    details = ex.value.details
    assert details.get("function") == "heun_single_rxn"
    assert details.get("delta_time") == 0.7
    assert details.get("rate") == -70
    assert math.isclose(details.get("rate_final"), 175)
    assert details.get("message") == "heun_single_rxn(): excessive time step (0.7), " \
                                     "leading to final rate of 175 which, when averaged with the " \
                                     "initial rate of -70.0 would flip the reaction's direction"

    with pytest.raises(ExcessiveTimeStepHard) as ex:
        # excessive time step (0.4001), leading to final rate of 70.0353 which, when averaged with the initial rate of -70 would flip the reaction's direction
        ReactionSimulator.heun_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.4001)

    details = ex.value.details
    assert details.get("function") == "heun_single_rxn"
    assert details.get("delta_time") == 0.4001
    assert details.get("rate") == -70
    assert math.isclose(details.get("rate_final"), 70.035)
    assert details.get("message") == "heun_single_rxn(): excessive time step (0.4001), " \
                                     "leading to final rate of 70.035 which, when averaged with the " \
                                     "initial rate of -70.0 would flip the reaction's direction"


    result = ReactionSimulator.heun_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.39)
    # Euler final_conc:  {'A': 37.3, 'B': 22.7}   [A] = 10 + (-1) * (-70) * 0.39   ;  [A] = 50 + 1 * -70 * 0.39
    # Final rate =  3. * 37.3 - 2. * 22.7 = 66.5
    # Heun's rate = (-70 + 66.5)/2 = -1.75
    assert result[0] == {'A': 0.6825, 'B': -0.6825}     #  delta [A] = (-1) * (-1.75) * 0.39
    assert result[1] == -70


    return       # TODO: continue
    # Reaction : A -> B  (no reverse reaction)
    rxn_defn = ReactionDefinition(reactants="A", products="B", species_registry=sr,
                                  reaction_model="mass action", kinetic_parameters={"kF": 3.})
    sim_rxn = rxn_defn.sim_reactions[0]

    result = ReactionSimulator.heun_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.1)
    assert result[0] == {'A': -3, 'B': 3}
    assert result[1] == 30    # Rate = 3. * 10.     Reaction is now forward

    result = ReactionSimulator.heun_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50}, delta_time=0.4)
    assert result[0] == {'A': -12, 'B': 12}         # Note: these increments would make [A] negative!
    assert result[1] == 30


    # Reaction : A + B <-> C
    rxn_defn = ReactionDefinition(reactants=["A" , "B"], products="C", species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action", kinetic_parameters={"kF": 5., "kR": 2})
    sim_rxn = rxn_defn.sim_reactions[0]

    C_0 = {"A": 10, "B": 50, "C": 20}
    delta_time = 0.002
    delta, rate = ReactionSimulator.heun_single_rxn(rxn=sim_rxn, conc_init=C_0, delta_time=delta_time)
    assert delta == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert rate == 5 * 10 * 50 - 2 * 20        # 2460 , i.e.  kF [A] [B] - kR [C]
    for k in delta:     # Loop over the keys
        assert delta[k] == rate * 0.002 * np.sign(sim_rxn.stoichiometry.vector[k])


    # Reaction : C <-> A + B
    rxn_defn = ReactionDefinition(reactants="C", products=["A" , "B"], species_registry=sr,
                                  reaction_model="mass action", kinetic_parameters={"kF": 2., "kR": 5.})
    sim_rxn = rxn_defn.sim_reactions[0]

    result = ReactionSimulator.heun_single_rxn(rxn=sim_rxn, conc_init={"A": 10, "B": 50, "C": 20}, delta_time=0.002)
    assert result[0] == {'A': -4.92, 'B': -4.92, 'C': 4.92}
    assert result[1] == -2460




############    Reactions    ############


def test_add_reaction():
    sim = ReactionSimulator()

    sim.add_reaction(reactants="X", products="Y", reaction_model="mass action")
    sim.add_reaction(reactants="C", products=["B", "A"], reaction_model="mass action")

    assert sim.species_index_map.index_to_species == ['X', 'Y', 'A', 'B', 'C']  # Notice the reaction-wise sorting
    assert sim.species_index_map.species_to_index == {'X': 0, 'Y': 1, 'A': 2, 'B': 3, 'C': 4}

    assert sim.species_index_map.index_of('X') == 0
    assert sim.species_index_map.index_of('Y') == 1
    assert sim.species_index_map.index_of('A') == 2
    assert sim.species_index_map.index_of('B') == 3
    assert sim.species_index_map.index_of('C') == 4

    with pytest.raises(Exception):
        sim.species_index_map.index_of('UNKNOWN')


    assert sim.species_index_map.species_at(0) == 'X'
    assert sim.species_index_map.species_at(1) == 'Y'
    assert sim.species_index_map.species_at(2) == 'A'
    assert sim.species_index_map.species_at(3) == 'B'
    assert sim.species_index_map.species_at(4) == 'C'

    with pytest.raises(Exception):
        sim.species_index_map.locate_species_id(5)



def test_indexes_of_active_species():
    sim = ReactionSimulator()

    assert sim.indexes_of_active_species() == []                 # No reactions yet

    sim.add_reaction(reactants="A", products="B", reaction_model="mass action")

    assert sim.species_index_map.index_to_species == ['A', 'B']
    assert sim.species_index_map.species_to_index == {'A': 0, 'B': 1}
    assert sim.indexes_of_active_species() == [0, 1]

    sim.add_reaction(reactants=["B", "X"], products=["C", "X"], reaction_model="mass action")
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C', 'X']              # Notice the newly-added names got sorted alphabetically
    assert sim.species_index_map.species_to_index == {'A': 0, 'B': 1, 'C': 2, 'X': 3}
    assert sim.indexes_of_active_species() == [0, 1, 2, 3]

    sim.add_reaction(reactants="X", products="Y", reaction_model="mass action")
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C', 'X', 'Y']
    assert sim.species_index_map.species_to_index == {'A': 0, 'B': 1, 'C': 2, 'X': 3, 'Y': 4}
    assert sim.indexes_of_active_species() == [0, 1, 2, 3, 4]

    sim.add_reaction(reactants=["A", "B", "Z"], products=["C", "Z"], reaction_model="mass action")
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C', 'X', 'Y', 'Z']
    assert sim.species_index_map.species_to_index == {'A': 0, 'B': 1, 'C': 2, 'X': 3, 'Y': 4, 'Z': 5}
    assert sim.indexes_of_active_species() == [0, 1, 2, 3, 4, 5]   # All



def test_set_temp():
    sim = ReactionSimulator()

    sim.set_temp(123)
    assert sim.temp == 123

    sim.set_temp(0, "C")
    assert np.allclose(sim.temp, 273.15)





##################   RESULT ANALYSIS   ##################


def test_is_in_equilibrium():
    sim = ReactionSimulator()

    with pytest.raises(Exception):
        assert sim.is_in_equilibrium(explain=False, tolerance=2.3)   # We're failing to provide concentrations, either as argument
                                                                     #   or as system concentrations

    # Reaction 0 : A <-> B
    sim.add_reaction(reactants="A", products="B",
                    reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})

    c = {'A': 23.9931640625, 'B': 36.0068359375}
    assert sim.is_in_equilibrium(rxn_index=0, conc=c, explain=False, tolerance=1)
    assert sim.is_in_equilibrium(conc=c, explain=False, tolerance=1)      # Testing ALL reactions

    # Reaction 1 : A <-> F
    sim.add_reaction(reactants=["A"], products=["F"],
                      reaction_model="mass action", kinetic_parameters={"kF": 20., "kR": 2.})

    c = {'A': 3, 'F': 32.999}
    assert sim.is_in_equilibrium(rxn_index=1, conc=c, explain=False, tolerance=10)   # The deviation is just below the 10% tolerance

    with pytest.raises(Exception):
        assert sim.is_in_equilibrium(conc=c, explain=False, tolerance=10) # Testing ALL reactions, but neglected to provide [B]

    c = {'A': 3, 'B': 4.6, 'F': 33.001}   # We're now including all the concentrations we need for both reactions
    assert sim.is_in_equilibrium(conc=c, explain=False, tolerance=10) \
           == {False: [1]}    # The deviation for reaction 1 is just above the 10% tolerance

    assert sim.is_in_equilibrium(conc=c, explain=False, tolerance=2.2) \
           == {False: [0, 1]}      # With a tolerance so low, not reaction meets the criteria

    assert sim.is_in_equilibrium(conc=c, explain=False, tolerance=2.3) \
           == {False: [1]}                # Reaction 0 barely passes with this tolerance (it deviates by 2.222 %)


    sim.system.set_system_state(conc=[3, 4.6, 33.001])    # "A", "B", "F"
    assert sim.get_single_conc('A') == 3
    assert sim.get_single_conc('B') == 4.6
    assert sim.get_single_conc('F') == 33.001

    assert sim.is_in_equilibrium(conc=c, explain=False, tolerance=2.3) \
           == {False: [1]}        # Now using the System concentrations



def test_reaction_in_equilibrium():
    sim = ReactionSimulator()

    # Reaction 0 : A <-> B
    sim.add_reaction(reactants=["A"], products=["B"],
                    reaction_model="mass action", kinetic_parameters={"kF": 3., "kR": 2.})
    c = {'A': 23.9931640625, 'B': 36.0068359375}
    assert sim.reaction_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)

    # Reaction 1 : A <-> F
    sim.add_reaction(reactants=["A"], products=["F"],
                    reaction_model="mass action", kinetic_parameters={"kF": 20., "kR": 2.})
    c = {'A': 3, 'F': 32.999}
    assert sim.reaction_in_equilibrium(rxn_index = 1, conc=c, explain=False, tolerance=10)   # Just below the 10% tolerance
    c = {'A': 3, 'F': 33.001}
    assert not sim.reaction_in_equilibrium(rxn_index = 1, conc=c, explain=False, tolerance=10)  # Just above the 10% tolerance

    # Reaction 2:  A + B <-> C , with 1st-order kinetics for each species
    sim.add_reaction(reactants=["A" , "B"], products="C",
                    reaction_model="mass action", kinetic_parameters={"kF": 5., "kR": 2.})
    c = {'A': 0.29487741, 'B': 40.29487741, 'C': 29.70512259}
    assert sim.reaction_in_equilibrium(rxn_index = 2, conc=c, explain=False, tolerance=1)


    sim.clear_reactions()   # This will reset the reaction count to 0

    # Reaction 0:  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction (i.e. regular mass-action kinetics)
    sim.add_reaction(reactants=[(2, "A")], products=["B"],
                    reaction_model="mass action", kinetic_parameters={"kF": 5., "kR": 2.})
    c = {'A': 1.51554944, 'B': 5.74222528}
    assert sim.reaction_in_equilibrium(rxn_index = 0, conc=c, tolerance=1, explain=False)






#######################    Macromolecules    #######################


def test_update_occupancy():
    sr = SpeciesRegistry(ids=["A", "B", "C", "M1", "M2"])
    mm_data = MacroMolecules(sr)

    mm_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=10)
    mm_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=20)
    mm_data.set_binding_site_affinity(macromolecule="M1", site_number=3, ligand="C", Kd=30)

    mm_data.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="C", Kd=3)
    mm_data.set_binding_site_affinity(macromolecule="M2", site_number=2, ligand="C", Kd=30)
    mm_data.set_binding_site_affinity(macromolecule="M2", site_number=3, ligand="C", Kd=300)

    assert mm_data.get_ligands() == {"A", "B", "C"}


    sim = ReactionSimulator(macromolecules=mm_data)

    # Verify that the ligands got added to the system state
    assert sim.species_index_map.index_to_species == ['A', 'B', 'C']
    assert sim.species_index_map.species_to_index == {"A": 0, "B": 1, "C": 2}

    sim.system.set_macromolecules()
    assert sim.system.macro_system == {"M1": 1, "M2": 1}

    sim.set_conc({"A": 10, "B": 20, "C": 30})

    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=1) , 0.)
    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=2) , 0.)
    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=3) , 0.)

    sim.update_occupancy()

    # All fractional occupancies will be 1/2, because the concentrations of the ligands below
    # exactly match their binding affinities
    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=1) , 0.5)
    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=2) , 0.5)
    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=3) , 0.5)

    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=1) , 0.9)  # Ligand conc is 10x binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=2) , 0.5)  # Ligand conc = binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=3) , 0.1)  # Ligand conc is 1/10 binding affinity


    # Vary the concentration of ligand C, starting with zero
    sim.system.set_single_conc(conc=0, species_id="C")     # No ligand C
    sim.update_occupancy()

    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=1) , 0.5)  # Unaffected (different ligand)
    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=2) , 0.5)  # Unaffected (different ligand)
    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=3) , 0)    # No occupancy in absence of ligand

    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=1) , 0)    # No occupancy in absence of ligand
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=2) , 0)    # No occupancy in absence of ligand
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=3) , 0)    # No occupancy in absence of ligand


    # Very low concentration of ligand C
    sim.system.set_single_conc(conc=0.3, species_id="C")
    sim.update_occupancy()

    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=1) , 0.5)      # Unaffected (different ligand)
    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=2) , 0.5)      # Unaffected (different ligand)
    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=3) , 0.012195) # Ligand conc is 1/100 binding affinity

    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=1) , 0.1)      # Ligand conc is 1/10 binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=2) , 0.012195) # Ligand conc is 1/100 binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=3) , 0.00136986)  # Ligand conc is 1/1000 binding affinity


    # Low concentration of ligand C
    sim.system.set_single_conc(conc=3, species_id="C")
    sim.update_occupancy()

    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=3) , 0.1)      # Ligand conc is 1/10 binding affinity

    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=1) , 0.5)      # Ligand conc = binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=2) , 0.1)      # Ligand conc is 1/10 binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=3) , 0.012195) # Ligand conc is 1/100 binding affinity


    # Mid concentration of ligand C
    sim.system.set_single_conc(conc=30, species_id="C")
    sim.update_occupancy()

    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=3) , 0.5)  # Ligand conc = binding affinity

    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=1) , 0.9)  # Ligand conc is 10x binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=2) , 0.5)  # Ligand conc = binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=3) , 0.1)  # Ligand conc is 1/10 binding affinity


    # High concentration of ligand C
    sim.system.set_single_conc(conc=300, species_id="C")
    sim.update_occupancy()

    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=3) , 0.9)       # Ligand conc is 10x binding affinity

    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=1) , 0.9878049)    # Ligand conc is 100x binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=2) , 0.9)          # Ligand conc is 10x binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=3) , 0.5)          # Ligand conc = binding affinity


    # Very High concentration of ligand C
    sim.system.set_single_conc(conc=3000, species_id="C")
    sim.update_occupancy()

    assert np.allclose(sim.system.get_occupancy(macromolecule="M1", site_number=3) , 0.9878049)    # Ligand conc is 100x binding affinity

    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=1) , 0.99863)      # Ligand conc is 1000x binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=2) , 0.9878049)    # Ligand conc is 100x binding affinity
    assert np.allclose(sim.system.get_occupancy(macromolecule="M2", site_number=3) , 0.9)          # Ligand conc is 10x binding affinity

    #print(rxn.get_occupancy(macromolecule="M2", site_number=1))
    #print(rxn.get_occupancy(macromolecule="M2", site_number=2))
    #print(rxn.get_occupancy(macromolecule="M2", site_number=3))



def test_sigmoid():
    sim = ReactionSimulator()

    # When the concentration matches Kd, the occupancy is 1/2 by definition
    assert np.allclose(sim.sigmoid(conc=1., Kd=1.), 0.5)
    assert np.allclose(sim.sigmoid(conc=0.1, Kd=0.1), 0.5)
    assert np.allclose(sim.sigmoid(conc=10., Kd=10.), 0.5)

    # When the concentration is 10% of Kd, the occupancy is about 0.1 ;
    # when the concentration is 10 times Kd, the occupancy is about 0.9
    assert np.allclose(sim.sigmoid(conc=0.1, Kd=1.), 0.1)
    assert np.allclose(sim.sigmoid(conc=10., Kd=1.), 0.9)

    assert np.allclose(sim.sigmoid(conc=10., Kd=100.), 0.1)
    assert np.allclose(sim.sigmoid(conc=1000., Kd=100.), 0.9)

    assert np.allclose(sim.sigmoid(conc=0.001, Kd=0.01), 0.1)
    assert np.allclose(sim.sigmoid(conc=0.1, Kd=0.01), 0.9)

    #print(rxn.sigmoid(conc=0.1, Kd=1.))



def test_logistic():
    sim = ReactionSimulator()

    #print(rxn.logistic(-10))
    assert np.allclose(sim.logistic(-100), 0)
    assert np.allclose(sim.logistic(-10), 0.00004539786)
    assert np.allclose(sim.logistic(0), 0.5)
    assert np.allclose(sim.logistic(10), 0.99995460213)
    assert np.allclose(sim.logistic(100), 1)






########    class VariableTimeSteps    ###########################################################################

def test_set_thresholds__VARIABLE_TIMESTEPS():
    rd = VariableTimeSteps()
    assert rd.thresholds == []

    with pytest.raises(Exception):
        rd.set_thresholds(norm=123)     # Bad norm name

    with pytest.raises(Exception):
        rd.set_thresholds(norm="")      # Bad norm name

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", low=5, high=1)     # Can't have low > high

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", low=5, abort=1)    # Can't have low > abort

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", high=8, abort=1)   # Can't have high > abort

    rd.set_thresholds(norm="norm_A", low=5)                 # Create a new rule, for norm_A
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 5}]

    rd.set_thresholds(norm="norm_A", low=6)                 # Update an existing value
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 6}]

    rd.set_thresholds(norm="norm_A", high=8)                # Add a new value to an existing rule
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 6, 'high': 8}]

    rd.set_thresholds(norm="norm_A", abort=10)               # Add a new value to an existing rule
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 6, 'high': 8, 'abort': 10}]

    # Bad values that violate low < high < abort
    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", low=8)

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", high=6)    # Too small

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", high=10)   # Too big

    with pytest.raises(Exception):
        rd.set_thresholds(norm="norm_A", abort=8)

    assert rd.thresholds == [{'norm': 'norm_A', 'low': 6, 'high': 8, 'abort': 10}]  # Nothing got changed by the failed calls



def test_delete_thresholds():
    rd = VariableTimeSteps()
    assert rd.thresholds == []

    with pytest.raises(Exception):
        rd.delete_thresholds(norm="not found")  # No rule with that norm exists

    rd.set_thresholds(norm="norm_A", low=1, high=2, abort=3)
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 1, 'high': 2, 'abort': 3}]

    rd.delete_thresholds(norm="norm_A", high=True)
    assert rd.thresholds == [{'norm': 'norm_A', 'low': 1, 'abort': 3}]

    with pytest.raises(Exception):
        rd.delete_thresholds(norm="norm_A", high=True)  # Trying to delete a non-existing threshold

    rd.delete_thresholds(norm="norm_A", low=True)
    assert rd.thresholds == [{'norm': 'norm_A', 'abort': 3}]

    rd.delete_thresholds(norm="norm_A", abort=True)
    assert rd.thresholds == []      # Nothing is left of that rule



def test_adjust_timestep():
    var_ts = VariableTimeSteps()

    prev =     np.array([1,   8, 8, 10, 10])
    baseline = np.array([2,   5, 5, 14, 14])
    delta =    np.array([0.5, 1, 4, -2, -5])

    n_chems = len(baseline)     # 5 chemicals (with indexes 0 thru 4)


    normA = var_ts.norm_A(delta_conc=delta)
    assert np.allclose(normA, 1.85)

    normB = var_ts.norm_B(baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normB, 0.8)


    var_ts.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.84)
    var_ts.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.79)
    var_ts.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4, error=0.25)

    indexes_of_active_chemicals = [0]     # To indicate that just the 0-th chemical is to be considered in the norms
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'operation': 'stay', 'step_factor': 1, 'norms': {'norm_A': 0.25, 'norm_B': 0.25}, 'applicable_norms': 'ALL'}

    indexes_of_active_chemicals = [0, 1]     # To indicate that just chemicals with indices 0 and 1 are to be considered in the norms
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'operation': 'stay', 'step_factor': 1, 'norms': {'norm_A': 0.3125, 'norm_B': 0.25}, 'applicable_norms': 'ALL'}

    indexes_of_active_chemicals = [0, 1, 2, 3, 4]       # All the chemicals are to be considered in the norms, from now on

    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,  delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'operation': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85}, 'applicable_norms': ['norm_A']}

    var_ts.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.86)     # normA (1.85) no longer triggers abort, but normB (0.8) still does
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'operation': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_B']}

    var_ts.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.81)    # normB (0.8) no longer triggers abort, but triggers a high.
                                                                        # normA (1.85) triggers a high, too
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'high', 'operation': 'downshift', 'step_factor': 0.5, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_A', 'norm_B']}

    var_ts.set_thresholds(norm="norm_A", low=0.5, high=1.86, abort=1.87)    # normA (1.85) no longer triggers high nor abort
    var_ts.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.79)
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'operation': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_B']}

    var_ts.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.81)    # normB (0.8) no longer triggers abort, but still triggers a high
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'high', 'operation': 'downshift', 'step_factor': 0.5, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_B']}

    var_ts.set_thresholds(norm="norm_B", low=0.08, high=0.81, abort=0.82)    # normB (0.8) no longer triggers high nor abort
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'operation': 'stay', 'step_factor': 1, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'ALL'}

    var_ts.set_thresholds(norm="norm_A", low=1.86, high=1.87, abort=1.88)   # normA (1.85) will now trigger a "low"
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'operation': 'stay', 'step_factor': 1, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'ALL'}
    # We're still on the 'stay' action because we aren't below ALL the thresholds

    var_ts.set_thresholds(norm="norm_B", low=0.81, high=0.82, abort=0.83)   # normB (0.8) will now trigger a "low", too
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'low', 'operation': 'upshift', 'step_factor': 1.2, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'ALL'}

    var_ts.set_thresholds(norm="norm_C", low=1.34, high=2.60, abort=2.61)   # normC (1.333) will still continue to trigger a "low"
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'low', 'operation': 'upshift', 'step_factor': 1.2, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 4/3}, 'applicable_norms': 'ALL'}

    var_ts.set_thresholds(norm="norm_C", low=1.32, high=2.60, abort=2.61)   # normC (1.333) will no longer trigger a "low" - but not a "high" nor an "abort"
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'operation': 'stay', 'step_factor': 1, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 4/3}, 'applicable_norms': 'ALL'}

    var_ts.set_thresholds(norm="norm_C", low=1.3, high=1.32, abort=2.61)   # normC (1.333) will now trigger a "high"
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,  delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'high', 'operation': 'downshift', 'step_factor': 0.5, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 4/3}, 'applicable_norms': ['norm_C']}

    var_ts.set_thresholds(norm="norm_C", low=1, high=1.31, abort=1.32)   # normC (1.333) will now trigger an "abort"
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'operation': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 4/3}, 'applicable_norms': ['norm_C']}

    var_ts.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.79)    # normB (0.8) will now trigger an "abort" before we can even get to normC
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'operation': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': ['norm_B']}

    var_ts.set_thresholds(norm="norm_A", low=0.5, high=1.0, abort=1.84)    # normA (1.85) will now trigger an "abort" before we can even get to normB
    result = var_ts.adjust_timestep(n_chems=n_chems, indexes_of_active_chemicals=indexes_of_active_chemicals,
                                delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'operation': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85}, 'applicable_norms': ['norm_A']}

    # TODO: also test with the other norms



def test_relative_significance():
    rd = VariableTimeSteps()

    # Assess the relative significance of various quantities
    #       relative to a baseline value of 10
    assert rd.relative_significance(1, 10) == "S"
    assert rd.relative_significance(4.9, 10) == "S"
    assert rd.relative_significance(5.1, 10) == "C"
    assert rd.relative_significance(10, 10) == "C"
    assert rd.relative_significance(19.9, 10) == "C"
    assert rd.relative_significance(20.1, 10) == "L"
    assert rd.relative_significance(137423, 10) == "L"



def test_compute_reaction_quotient():
    # Reaction : A <-> B
    q = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"A")], product_data=[(1,"B")],
                                                 conc={'A': 24., 'B': 36.}, explain=False)
    assert np.allclose(q, 1.5)
    q, formula = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"A")], product_data=[(1,"B")],
                                                          conc={'A': 24., 'B': 36.}, explain=True)
    assert np.allclose(q, 1.5)
    assert formula == '[B] / [A]'


    # Reaction : 2 A + B <-> C
    q, formula = ReactionKinetics.compute_reaction_quotient(reactant_data=[(2, "A"), (1, "B")], product_data=["C"],
                                                          conc={'A': 24., 'B': 36., 'C': 10}, explain=True)

    assert np.allclose(q, 10 / (24.**2 * 36))
    assert formula == '[C] / ( [A]^2 [B])'


    # Reaction: R + S <-> P + Q
    q, formula = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"R"), "S"],
                                                          product_data=["P", (1,"Q")],
                                                          conc={'R': 10., 'S': 4., 'P': 5., 'Q': 20.}, explain=True)
    assert np.allclose(q, 2.5)    #  (5 * 20) / (10 * 4)
    assert formula == '([P][Q]) / ([R][S])'

    # Reaction: R + 3 S <-> 2 P + Q
    q, formula = ReactionKinetics.compute_reaction_quotient(reactant_data=[(1,"R"), (3,"S")],
                                                          product_data=[(2,"P"), (1,"Q")],
                                                          conc={'R': 10., 'S': 4., 'P': 5., 'Q': 20.}, explain=True)
    assert np.allclose(q, 0.78125)    #  (5**2 * 20) / (10 * 4**3)
    assert formula == '( [P]^2 [Q]) / ([R] [S]^3 )'





#################   INDIVIDUAL NORMS   #################

def test_norms():
    rd = VariableTimeSteps()

    prev =     np.array([1,   8, 8, 10, 10])
    baseline = np.array([2,   5, 5, 14, 14])
    delta =    np.array([0.5, 1, 4, -2, -5])

    normA = rd.norm_A(delta_conc=delta)
    assert np.allclose(normA, 1.85)

    normB = rd.norm_B(baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normB, 0.8)

    normC = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normC, 4/3)

    normD = rd.norm_D(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normD, 0.7833333333333333)



def test_norm_A():
    rd = VariableTimeSteps()

    delta_conc = np.array([1, 4])
    result = rd.norm_A(delta_conc)
    assert np.allclose(result, 4.25)

    delta_conc = np.array([.5, 2])
    result = rd.norm_A(delta_conc)
    assert np.allclose(result, 1.0625)

    delta_conc = np.array([.5, 2, 1])
    result = rd.norm_A(delta_conc)
    assert np.allclose(result, 0.5833333333)

    delta_conc = np.array([0.5, 1, 4, -2, -5])
    result = rd.norm_A(delta_conc)
    assert np.allclose(result, 1.85)



def test_norm_B():
    rd = VariableTimeSteps()

    delta = np.array([4, -1])
    base = np.array([10, 2])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.5)

    base = np.array([10, 0])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.4)     # The zero baseline concentration is disregarded

    delta = np.array([300,        6, 0, -1])
    base = np.array([0.00000001, 10, 0,  2])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.6)

    delta = np.array([300,   6, 0, -1])
    base = np.array([0.001, 10, 0,  2])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 300000)

    base = np.array([2,   5, 5, 14, 14])
    delta = np.array([0.5, 1, 4, -2, -5])
    result = rd.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.8)

    with pytest.raises(Exception):
        rd.norm_B(baseline_conc=np.array([1, 2, 3]), delta_conc=delta) # Too many entries in array

    with pytest.raises(Exception):
        rd.norm_B(baseline_conc=base, delta_conc=np.array([1]))        # Too few entries in array



def test_norm_C():
    rd = VariableTimeSteps()

    result = rd.norm_C(prev_conc=np.array([1]), baseline_conc=np.array([2]), delta_conc=np.array([0.5]))
    assert result == 0

    result = rd.norm_C(prev_conc=np.array([8]), baseline_conc=np.array([5]), delta_conc=np.array([1]))
    assert result == 0

    result = rd.norm_C(prev_conc=np.array([8]), baseline_conc=np.array([5]), delta_conc=np.array([4]))
    assert np.allclose(result, 4/3)

    result = rd.norm_C(prev_conc=np.array([10]), baseline_conc=np.array([14]), delta_conc=np.array([-2]))
    assert result == 0

    result = rd.norm_C(prev_conc=np.array([10]), baseline_conc=np.array([14]), delta_conc=np.array([-5]))
    assert np.allclose(result, 5/4)

    prev =     np.array([1,   8, 8, 10, 10])
    baseline = np.array([2,   5, 5, 14, 14])
    delta =    np.array([0.5, 1, 4, -2, -5])
    result = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(result, 4/3)

    # A scenario where the 'prev' and 'baseline' values are almost identical
    prev = np.append(prev, 3)
    baseline = np.append(baseline, 2.999999999)
    delta = np.append(delta, 8)
    result = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(result, 4/3)

    # A scenario where the 'delta' dwarfs the change between 'prev' and 'baseline'
    prev = np.append(prev, 10)
    baseline = np.append(baseline, 10.05)
    delta = np.append(delta, -9)
    result = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(result, 4/3)

    prev = np.append(prev, 10)
    baseline = np.append(baseline, 10.2)
    delta = np.append(delta, -9)
    result = rd.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(result, 45)



def test_norm_D():
    rd = VariableTimeSteps()

    prev =     np.array([ 12.96672432,  31.10067726,  55.93259842,  44.72389482, 955.27610518])
    baseline = np.array([ 12.99244738,  31.04428765,  55.96326497,  43.91117372, 956.08882628])
    delta =    np.array([-2.56160549,   -0.03542113,   2.59702662,   1.36089356,  -1.36089356])
    result = rd.norm_D(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(result, 37.64942285399873)
