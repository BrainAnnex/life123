import math
import pandas as pd
from pandas.testing import assert_frame_equal
import pytest
import numpy as np
from life123.uniform_compartment import UniformCompartment
from life123.reaction_kinetics import ReactionKinetics
from life123.reaction_simulator import ReactionSimulator, AnalyticReactionSolver, VariableTimeSteps
from life123.species_registry import SpeciesRegistry
from life123.reactions import ReactionDefinition
from life123.reaction_registry import ReactionRegistry
from life123.kinetics import Custom_Model
from life123.reaction_simulator import ExcessiveTimeStepHard, ExcessiveTimeStepSoft
from life123.species_index_map import SpeciesIndexMap
from life123.diagnostics import Diagnostics
from life123.history import HistoryUniformConcentration, HistoryReactionRate
from life123.collections import CollectionTabular
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





########    class ReactionSimulator    ###########################################################################



def test_single_compartment_react():

    # Test based on experiment "cycles_1"
    species_registry = SpeciesRegistry(ids=["A", "B", "C", "E_high", "E_low"])
    rxns = ReactionRegistry(species_data=species_registry)

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


    # Assign an array index to all the species we're dealing with
    all_species = rxns.get_species_in_any_reaction(sort=True)

    ind = SpeciesIndexMap(all_species)
    assert ind.index_to_species == ['A', 'B', 'C', 'E_high', 'E_low']

    uc = UniformCompartment(species_data=species_registry, index_species=ind)
    initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}
    uc.set_conc(conc=initial_conc, snapshot=True)
    #print(uc.system)
    assert np.allclose(uc.system, [ 100.,    0.,    0., 1000. ,   0.])

    sim = ReactionSimulator(system=uc.system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler",
                            diagnostics_enabled=True)
    sim.uniform_compartment = uc

    uc.reaction_simulator = sim
    uc.diagnostics = sim.diagnostics
    #print(sim.system)
    assert np.allclose(sim.system, [ 100.,    0.,    0., 1000. ,   0.])


    sim.single_compartment_react(initial_step=0.0005, target_end_time=0.0035, variable_steps=False)
    assert np.allclose(sim.system, uc.system)

    run1 = uc.get_system_conc()

    assert np.allclose(sim.system_time, 0.0035)
    assert np.allclose(run1, [9.69252541e+01, 3.05696280e+00, 1.77831454e-02, 9.99980686e+02, 1.93144884e-02])
    #assert sim.diagnostics.explain_time_advance(return_times=True, silent=True) == \
    #           ([0.0, 0.0035], [0.0005])       TODO: fix



def test__single_compartment_react_main_loop_1():
    # FIXED steps

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_data=species_registry)
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
    assert np.allclose(sim.system, [17., 43.])
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
    sim.system = np.array([10., 50.])
    sim.system_time = 0
    sim.rate_history = HistoryReactionRate(active=True)         # We're not instantiating a new "ReactionSimulator"; so, we must reset some values
    sim.conc_history = HistoryUniformConcentration(active=True)
    # Pick a time step just below 50/70 (where 70 is the initial reaction rate),
    # to bring the second species concentration to almost zero
    new_count, recommended_next_step = sim._single_compartment_react_main_loop(time_step=0.7142857142, variable_steps=False,
                                                                               step_count=5, n_steps=1000)
    assert new_count == 6
    assert math.isclose(recommended_next_step, 0.7142857142)
    assert np.allclose(sim.system, [60., 0.])       # This time step was so large that it converted all B to A !
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
    sim.system = np.array([10., 50.])
    sim.system_time = 20
    # Pick a time step just above the previous time, to tip [B] into negative values
    with pytest.raises(ExcessiveTimeStepHard) as ex:
        sim._single_compartment_react_main_loop(time_step=0.714285715, variable_steps=False,
                                                step_count=1, n_steps=1000)
    details = ex.value.details
    assert details["function"] == "reaction_step_common_fixed_step"
    assert math.isclose(details["delta_time"], 0.714285715)
    assert details["caption"] == 'aborted: neg. conc. in `B` from rxn # 0'
    assert details["system_time"] == 20
    assert details["rate"] == -70.0
    assert details["rxn_index"] == 0
    assert details["previous_function"] == "_validate_increment"
    assert details["message"] == """reaction_step_common_fixed_step(): unable to complete the reaction step.  Try REDUCING the time step, or switching to variable time steps. 
DETAILS: 
      The tentative time step (0.714286) would lead to a NEGATIVE concentration in the species `B` from the reaction `A <-> B` (rxn # 0)
      Baseline concentration value of `B` : 50 at system time 20; requested change (NOT carried out): -50"""

    # The following values are all unchanged
    assert sim.system_time == 20
    assert np.allclose(sim.system, [10., 50.])
    assert np.allclose(sim.previous_system, [10., 50.])



def test__single_compartment_react_main_loop_2_a():
    # VARIABLE steps

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_data=species_registry)
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
    assert np.allclose(sim.system, [11.4, 48.6])
    assert np.allclose(sim.previous_system, [10., 50.])
    assert math.isclose(sim.system_time, 8.02)

    # Inspect the `sim.reaction_step_diagnostics` data object
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.02)
    assert sim.reaction_step_diagnostics.number_neg_concs == 0
    assert sim.reaction_step_diagnostics.number_soft_aborts == 0
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
    assert np.allclose(sim.system, [11.68, 48.32])
    assert np.allclose(sim.previous_system, [10., 50.])
    assert math.isclose(sim.system_time, 17.024)    # A step of just 0.024 was actually taken

    # Inspect the `sim.reaction_step_diagnostics` data object
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.024)
    assert sim.reaction_step_diagnostics.number_neg_concs == 0
    assert sim.reaction_step_diagnostics.number_soft_aborts == 0
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
    assert np.allclose(sim.system, [11.306368, 48.693632])
    assert np.allclose(sim.previous_system, [10., 50.])
    assert math.isclose(sim.system_time, 666.0186624)    # A step of just 0.0186624 was actually taken

    # Inspect the `sim.reaction_step_diagnostics` data object
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.0186624)    # 0.8 * 0.5 * (0.6)^6, from multiple re-tries
    assert sim.reaction_step_diagnostics.number_neg_concs == 0
    assert sim.reaction_step_diagnostics.number_soft_aborts == 0
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
    rxns = ReactionRegistry(species_data=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})
    assert ind.index_to_species == ["A", "B"]

    ### Repeat the simulation from previous test, this time with diagnostics enabled (excessive large time step; so large as to cause a backtracking)
    initial_system = np.array([10., 50.])
    sim = ReactionSimulator(system=initial_system, species_index_map=ind, diagnostics_enabled=True,
                            reaction_registry=rxns, method="forward_euler", preset="fast")  # Note the preset
    sim.system_time = 666

    new_count, step_recommended = sim._single_compartment_react_main_loop(time_step=0.8, variable_steps=True,
                                                                          step_count=13, n_steps=1000, explain_variable_steps=(-1, 1000)) #

    assert new_count == 14
    assert math.isclose(step_recommended, 0.0186624)   # "Go much smaller" in next round
    assert np.allclose(sim.system, [11.306368, 48.693632])
    assert np.allclose(sim.previous_system, [10., 50.])
    assert math.isclose(sim.system_time, 666.0186624)    # A step of just 0.0186624 was actually taken


    # Verify the diagnostic data: part 1 - the "diagnostic_rxn_data"
    assert type(sim.diagnostics.diagnostic_rxn_data) is dict
    assert len(sim.diagnostics.diagnostic_rxn_data) == 1
    coll_tab = sim.diagnostics.diagnostic_rxn_data[0]   # For reaction 0 (the only reaction we have)
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



def test_reaction_step_common_fixed_step():
    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_data=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})
    assert ind.index_to_species == ["A", "B"]

    system = np.array([10, 50])
    sim = ReactionSimulator(system=system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler")

    result = sim.reaction_step_common_fixed_step(delta_time=0.1)
    assert np.allclose(result, [7, -7])
    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.1)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    sim.system = np.array([10, 50])     # Reset the system state
    with pytest.raises(ExcessiveTimeStepHard) as ex:      # Excessive time step that would make [B] negative
        sim.reaction_step_common_fixed_step(delta_time=0.8)
    details = ex.value.details
    assert details["function"] == "reaction_step_common_fixed_step"
    assert details["delta_time"] == 0.8
    assert details["system_time"] == 0
    assert details["rate"] == -70.0
    assert details["rxn_index"] == 0
    assert details["previous_function"] == "_validate_increment"
    assert details["message"] == """reaction_step_common_fixed_step(): unable to complete the reaction step.  Try REDUCING the time step, or switching to variable time steps. 
DETAILS: 
      The tentative time step (0.8) would lead to a NEGATIVE concentration in the species `B` from the reaction `A <-> B` (rxn # 0)
      Baseline concentration value of `B` : 50 at system time 0; requested change (NOT carried out): -56"""



    sim.system = np.array([10, 50])     # Reset the system state
    sim.method = "heun"
    with pytest.raises(ExcessiveTimeStepHard) as ex:      # Excessive time step that would make [B] negative
        sim.reaction_step_common_fixed_step(delta_time=0.8)
    details = ex.value.details
    assert details["function"] == "reaction_step_common_fixed_step"
    assert details["delta_time"] == 0.8
    assert details["rate"] == -70.0
    assert details["previous_function"] == "heun_single_rxn"
    assert details["message"] == """reaction_step_common_fixed_step(): unable to complete the reaction step.  Try REDUCING the time step, or switching to variable time steps. 
DETAILS: 
heun_single_rxn(): excessive time step (0.8), leading to negative concentrations"""

    #for k, v in details.items():
    #    print(k, " : ", v)



def test_reaction_step_common_variable_step_1():
    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_data=species_registry)
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
    assert sim.reaction_step_diagnostics.number_neg_concs == 0
    assert sim.reaction_step_diagnostics.number_soft_aborts == 0
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
    assert sim.reaction_step_diagnostics.number_neg_concs == 0
    assert sim.reaction_step_diagnostics.number_soft_aborts == 0
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
    assert sim.reaction_step_diagnostics.number_neg_concs == 0
    assert sim.reaction_step_diagnostics.number_soft_aborts == 0
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
    rxns = ReactionRegistry(species_data=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})
    assert ind.index_to_species == ["A", "B"]

    # New simulation, with slightly larger time step
    initial_system = np.array([10, 50])
    sim = ReactionSimulator(system=initial_system, species_index_map=ind,
                            reaction_registry=rxns, method="forward_euler", preset="fast")

    sim.system = initial_system
    incr, step_taken, step_recommended = sim.reaction_step_common_variable_step(delta_time=0.04, explain_variable_steps=[-1,1])
    captured = capsys.readouterr()  # Capture the standard output and standard error, from the previous function call

    assert math.isclose(step_taken, 0.024)          # Smaller than suggested [but 20% longer than in our previous test run, above]
    assert math.isclose(step_recommended, 0.0192)   # "Go even smaller" in nest round (80% of current value)
    assert np.allclose(incr, [1.68, -1.68])

    assert math.isclose(sim.reaction_step_diagnostics.delta_time, 0.024)
    assert sim.reaction_step_diagnostics.number_neg_concs == 0
    assert sim.reaction_step_diagnostics.number_soft_aborts == 0
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
    rxns = ReactionRegistry(species_data=species_registry)
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
    sim.system = np.array([10, 50])
    sim.method = "heun"
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.1, variable_steps=False)
    assert np.allclose(delta_conc, [5.25, -5.25])
    assert math.isclose(rec_next_step, 0.1)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Reset and re-run (back to "forward_euler", but this time with variable step, and a smaller step)
    sim.system = np.array([10, 50])
    sim.method = "forward_euler"
    sim.adaptive_steps.use_adaptive_preset(preset="fast")
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.02, variable_steps=True)
    assert np.allclose(delta_conc, [1.4, -1.4])
    assert math.isclose(rec_next_step, 0.02)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Reset and re-run (this time with a printed explanation of variable time steps)
    sim.system = np.array([10, 50])
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
    sim.system = np.array([10, 50])
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
    sim.system = np.array([10, 50])
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
    rxns = ReactionRegistry(species_data=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})

    system = np.array([10, 50])
    sim = ReactionSimulator(system=system, species_index_map=ind, reaction_registry=rxns, method="forward_euler",
                            diagnostics_enabled=True)   # Diagnostics are enabled at instantiation time
    sim.system_time = 88
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.1, variable_steps=False)     # FIXED time step
    assert np.allclose(delta_conc, [7, -7])
    assert math.isclose(rec_next_step, 0.1)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Verify the diagnostic data: part 1 - the "diagnostic_rxn_data", created by single_step_single_rxn()
    assert type(sim.diagnostics.diagnostic_rxn_data) is dict
    assert len(sim.diagnostics.diagnostic_rxn_data) == 1
    coll_tab = sim.diagnostics.diagnostic_rxn_data[0]
    assert type(coll_tab) is CollectionTabular
    df = coll_tab.get_dataframe()
    row_expected = {"START_TIME": 88, "time_step": 0.1, "aborted": False, "Delta A": 7.0, "Delta B": -7.0, "rate": -70., "caption": ""}
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert compare_pandas(df_expected, df, disregard_order=True)

    # Verify the diagnostic data: part 2 - the "diagnostic_decisions_data", created by attempt_reaction_step()
    assert type(sim.diagnostics.diagnostic_decisions_data) is CollectionTabular
    df = sim.diagnostics.diagnostic_decisions_data.get_dataframe()
    assert type(df) is pd.DataFrame
    assert len(df) == 1
    row_expected = {"START_TIME": 88, "Delta A": 7.0, "Delta B": -7.0, "caption": ""}
    df_expected = pd.DataFrame(row_expected, index=[0])
    assert compare_pandas(df_expected, df, disregard_order=True)



def test_attempt_reaction_step_2_b():
    # With DIAGNOSTICS enabled, for VARIABLE reaction step

    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_data=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})

    system = np.array([10, 50])
    sim = ReactionSimulator(system=system, species_index_map=ind, reaction_registry=rxns, method="forward_euler",
                            diagnostics_enabled=True)   # Diagnostics are enabled at instantiation time
    sim.system_time = 99
    sim.adaptive_steps.use_adaptive_preset(preset="fast")
    delta_conc, rec_next_step = sim.attempt_reaction_step(delta_time=0.02, variable_steps=True)     # VARIABLE time step
    assert np.allclose(delta_conc, [1.4, -1.4])
    assert math.isclose(rec_next_step, 0.02)
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Verify the diagnostic data: part 1 - the "diagnostic_rxn_data", created by single_step_single_rxn()
    assert type(sim.diagnostics.diagnostic_rxn_data) is dict
    assert len(sim.diagnostics.diagnostic_rxn_data) == 1
    coll_tab = sim.diagnostics.diagnostic_rxn_data[0]
    assert type(coll_tab) is CollectionTabular
    df = coll_tab.get_dataframe()
    row = {"START_TIME": 99, "time_step": 0.02, "aborted": False, "Delta A": 1.4, "Delta B": -1.4, "rate": -70., "caption": ""}
    df_expected = pd.DataFrame(row, index=[0])
    assert_frame_equal(df, df_expected)


    # Verify the diagnostic data: part 2 - the "diagnostic_decisions_data", created by attempt_reaction_step(), this time with
    #   extra fields resulting from the VARIABLE step
    assert type(sim.diagnostics.diagnostic_decisions_data) is CollectionTabular
    df = sim.diagnostics.diagnostic_decisions_data.get_dataframe()
    assert type(df) is pd.DataFrame
    assert len(df) == 1
    # More fields are present now, because of the variable step
    row = {"START_TIME": 99, "Delta A": 1.4, "Delta B": -1.4,
           "norm_A": 0.98, "norm_B": 0.14, "norm_C": None, "norm_D": None, "action": "OK (stay)", "step_factor": 1, "time_step": 0.02,
           "caption": ""}
    df_expected = pd.DataFrame(row, index=[0])
    assert_frame_equal(df, df_expected)



def test_single_step_all_rxns():
    species_registry = SpeciesRegistry()

    # Reaction : A <-> B  (created thu the ReactionRegistry object)
    rxns = ReactionRegistry(species_data=species_registry)
    rxns.add_reaction(reactants="A", products="B", reaction_model="mass action",
                      kinetic_parameters={"kF": 3., "kR": 2.})

    ind = SpeciesIndexMap({"A", "B"})

    system = np.array([10, 50])
    sim = ReactionSimulator(system=system, species_index_map=ind, reaction_registry=rxns)
    result = sim.single_step_all_rxns(delta_time=0.1)
    assert np.allclose(result, [7, -7])
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    # Reset and re-run
    sim.system = np.array([10, 50])
    sim.method = "heun"
    result = sim.single_step_all_rxns(delta_time=0.1)
    assert np.allclose(result, [5.25, -5.25])
    assert compare_dicts(sim.reaction_step_diagnostics.system_rxn_rates, {0: -70.0})


    #print(result)   #TODO: in progress



def test_single_step_single_rxn():

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
    assert np.allclose(sim.system, [10, 50])

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
    sim = ReactionSimulator(method="analytic", exact=True,
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

    reaction_registry = ReactionRegistry(species_data=species_registry)
    reaction_registry.register_reaction(rxn_defn)

    ind = SpeciesIndexMap()
    species_id_set =  rxn_sim.stoichiometry.get_all_species_ids()
    ind.add_species(species_id_set)

    diagnostics = Diagnostics(reactions=reaction_registry, species_to_index=ind.species_to_index)

    system = np.array([10, 50])
    sim = ReactionSimulator(method="forward_euler", system=system, species_index_map=ind,
                            diagnostics_enabled=True, diagnostics=diagnostics)

    increment_vector = np.zeros(2, dtype='d')
    sim.system_time = 123
    result = sim.single_step_single_rxn(increment_vector=increment_vector, delta_time=0.1,
                                        rxn=rxn_sim, rxn_index=0)
    assert np.allclose(increment_vector, [7, -7])
    assert math.isclose(result, -70)

    assert len(sim.diagnostics.diagnostic_rxn_data) == 1

    collection = sim.diagnostics.diagnostic_rxn_data[0]
    assert type(collection) is CollectionTabular
    df = collection.get_dataframe()
    assert type(df) is pd.DataFrame

    row = {"START_TIME": 123, "time_step": 0.1, "aborted": False, "Delta A": 7.0, "Delta B": -7.0, "rate": -70.0, "caption": ""}
    df_expected = pd.DataFrame(row, index=[0])

    assert compare_pandas(df_expected, df)



    # Re-start with a fresh Diagnostics object; this time we'll take an EXCESSIVE step size
    diagnostics = Diagnostics(reactions=reaction_registry, species_to_index=ind.species_to_index)

    system = np.array([10, 50])
    sim = ReactionSimulator(method="forward_euler", system=system, species_index_map=ind,
                            diagnostics_enabled=True, diagnostics=diagnostics)

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
    assert type(sim.diagnostics.diagnostic_rxn_data) is dict
    assert len(sim.diagnostics.diagnostic_rxn_data) == 1
    collection = sim.diagnostics.diagnostic_rxn_data[0]
    assert type(collection) is CollectionTabular
    df = collection.get_dataframe()
    assert type(df) is pandas.DataFrame
    row = {"START_TIME": 666, "time_step": 0.8, "aborted": True, "Delta A": np.nan, "Delta B": np.nan,
           "rate": -70.0, "caption": "aborted: neg. conc. in `B`"}
    df_expected = pd.DataFrame(row, index=[0])
    assert compare_pandas(df_expected, df)

    # part 2: diagnostic_decisions_data
    assert len(sim.diagnostics.diagnostic_decisions_data) == 1
    collection = sim.diagnostics.diagnostic_decisions_data
    assert type(collection) is CollectionTabular
    df = collection.get_dataframe()
    assert type(df) is pandas.DataFrame
    row = {"START_TIME": 666, "action": "ABORT", "caption": "", "time_step": 0.8}
    df_expected = pd.DataFrame(row, index=[0])
    assert compare_pandas(df_expected, df)
    """



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

    sim = ReactionSimulator(method="analytic", exact=True)
    delta_conc, rate = sim.dispatcher_single_rxn(rxn=rxn_sim, conc_init={"A": 10, "B": 50}, delta_time=0.1)
    assert compare_dicts(delta_conc, {'A': 5.508570764023133, 'B': -5.508570764023133}) # A lot closer to Heun than to forward Euler!
    assert math.isclose(rate, -70)



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






########    class AnalyticalReactionSolver    ###########################################################################


def test_exact_advance_unimolecular_reversible():
    # Reaction A <-> P
    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0)
    assert np.allclose(p, 10.)      # No change

    incr = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0, incremental=True)
    assert np.allclose(incr, 0)     # No change


    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.005)
    assert np.allclose(p, 11.08636387)

    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.2)
    assert np.allclose(p, 37.81330458845654)

    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.31739)
    assert np.allclose(p, 45.)

    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=1.12)
    assert np.allclose(p, 53.837294)

    incr = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=1.12, incremental=True)
    assert np.allclose(incr, 53.837294-10)      # P(t) - P0


    p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=100.)
    assert np.allclose(p, 54.)

    equil = ReactionKinetics._compute_equilibrium_conc_first_order(kF=3., kR=2., a=1, A0=80., p=1, P0=10.)
    assert np.allclose(p, equil["P"])


def test_exact_advance_unimolecular_reversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law,
    # at the middle point of 3 sampled points
    h = 0.0005
    t_start = 0.2
    p1 = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=t_start)
    p2 = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=t_start + h)
    a2 = 90 - p2    # From mass conservation
    p3 = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=t_start + 2 * h)
    #print(p2)
    gradient = np.gradient([p1, p2, p3], h)
    derivative_p2 = gradient[1]
    rate_at_mid_point = 3 * a2 - 2 * p2     # kF * a2 - kR * p2
    assert np.allclose(derivative_p2, rate_at_mid_point)


def test_exact_advance_unimolecular_reversible_3():
    # Compare the exact solution against a fine-grained forward-Euler approximation
    sr = SpeciesRegistry()
    rxn_defn = ReactionDefinition(reactants="A", products="P",
                                  species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF": 3., "kR": 2.})
    rxn = rxn_defn.sim_reactions[0]

    t_final = 0.2
    n_steps = 50000
    t_step = t_final / n_steps
    a = 80
    p = 10
    for i in range(n_steps):
        increment_dict, _ = \
            rxn.step_simulation(delta_time=t_step, conc_dict={"A": a, "P": p})
        delta_p = increment_dict["P"]
        p += delta_p
        a -= delta_p

    exact_p = AnalyticReactionSolver.exact_advance_unimolecular_reversible(kF=3., kR=2., A0=80., P0=10., t=0.2)   # 37.81330458845654
    assert np.allclose(p, exact_p)



def test_exact_advance_unimolecular_irreversible():
    # Reaction A -> P
    # Compare against the REVERSIBLE reaction solver with zero reverse rate constant
    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0)
    assert np.allclose(p, 10.)      # No change

    incr = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0, incremental=True)
    assert np.allclose(incr, 0)     # No change

    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0.005)
    assert np.allclose(p, 11.191044831754994)

    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0.31739)
    assert np.allclose(p, 59.127783572982025)

    incr = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=0.31739, incremental=True)
    assert np.allclose(incr, p-10)      # P(t) - P0

    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=1.12)
    assert np.allclose(p, 87.22117928442091)

    p = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=100.)
    assert np.allclose(p, 90)

    incr = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=100., incremental=True)
    assert np.allclose(incr, 80)


def test_exact_advance_unimolecular_irreversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law, at the middle point
    h = 0.0005
    t_start = 0.3
    p1 = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=t_start)
    p2 = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=t_start + h)
    a2 = 90 - p2    # From mass conservation
    p3 = AnalyticReactionSolver.exact_advance_unimolecular_irreversible(kF=3., A0=80., P0=10., t=t_start + 2 * h)

    gradient = np.gradient([p1, p2, p3], h)
    derivative_p2 = gradient[1]
    rate_at_mid_point = 3 * a2   # kF * a2
    assert np.allclose(derivative_p2, rate_at_mid_point)



def test_exact_advance_synthesis_reversible():
    # Reaction A + B <-> P

    # General case A0 != B0

    Dp = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0, incremental=True)
    assert np.allclose(Dp, 0.)     # No change

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0, incremental=False)
    assert np.allclose(P_t, 20.)    # No change

    # The comparison values below are from Octave, v. 10.3.0  ;  for example:
    '''
    lsode_options("absolute tolerance", 1e-12);
    lsode_options("relative tolerance", 1e-10);
    
    kf = 5.0;
    kr = 2.0;
    a0 = 10.0;
    b0 = 50.0;
    p0 = 20.0;
    
    p_init = p0;
    
    t = [0, 0.0001, 0.0003, 0.0004, 0.0005, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.008, 0.01, 1.];
    
    function dpdt = bimol_ode(p, t, kf, kr, a0, b0, p0)
      dpdt = kf * (a0 - p + p0) .* (b0 - p + p0) - kr * p;
    end
    
    p = lsode(@(p,t) bimol_ode(p,t,kf,kr,a0,b0,p0), p_init, t);
    
    format long
    p
    '''
    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0001, incremental=False)
    assert np.allclose(P_t, 20.24233230049967)

    Dp = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0001, incremental=True)
    assert np.allclose(Dp, 0.24233230049967)    # This is a DELTA value

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0003, incremental=False)
    assert np.allclose(P_t, 20.70580470925962)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0004, incremental=False)
    assert np.allclose(P_t, 20.92746192779109)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0005, incremental=False)
    assert np.allclose(P_t, 21.14272449451510)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0007, incremental=False)
    assert np.allclose(P_t, 21.55497320836003)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.001, incremental=False)
    assert np.allclose(P_t, 22.13079644707523)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.0015, incremental=False)
    assert np.allclose(P_t, 22.98939523899655)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.002, incremental=False)
    assert np.allclose(P_t, 23.73870896747882)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.003, incremental=False)
    assert np.allclose(P_t, 24.97201963453988)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.005, incremental=False)
    assert np.allclose(P_t, 26.68110034281166)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.008, incremental=False)
    assert np.allclose(P_t, 28.12354983857080)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=0.01, incremental=False)
    assert np.allclose(P_t, 28.66884983321557)

    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=1., incremental=False)
    assert np.allclose(P_t, 29.70512259124693)


    # Verify reaching equilibrium at a large t
    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=10., incremental=False)
    eq = ReactionKinetics._compute_equilibrium_conc_first_order(kF=5., kR=2., a=1, A0=10., b=1, B0=50., p=1, P0=20.)
    assert np.allclose(P_t, eq["P"])    # 29.705122591242464


    # Special case A0 = B0
    P_t = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=35, B0=35, P0=20., t=0.002, incremental=False)
    assert np.allclose(P_t, 28.99659131142920)



def test_exact_advance_synthesis_reversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law,
    # at the middle point of 3

    # Case where `A` is the limiting reagent
    A0=10
    B0=50
    P0=20
    kF=5
    kR=2

    #print(ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0))
    #Equilibrium at {'A': 0.2948774087575341, 'B': 40.294877408757536, 'P': 29.705122591242464}

    h = 0.00001
    t_start = 0.002
    times = t_start + h * np.arange(3)

    # Sample at 3 closely-spaced points
    p = [AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
         for t in times
         ]
    #print(p)    # [23.738708979042038, 23.752698461282, 23.766650993990154]

    p_middle = p[1]
    a_middle = A0 - (p_middle - P0)    # From mass conservation
    b_middle = B0 - (p_middle - P0)    # From mass conservation

    #print(a_middle, b_middle)        # 6.247301538717998 46.247301538718
    gradient = np.gradient(p, h)
    derivative_middle = gradient[1]
    rate_at_mid_point = kF * a_middle * b_middle - kR * p_middle    # kF [A] [B] - kR [P]
    assert np.allclose(derivative_middle, rate_at_mid_point)


    # Case where `B` is the limiting reagent
    A0=60
    B0=30
    P0=10
    kF=7
    kR=3

    #print(ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0))
    #Equilibrium at {'A': 30.553318635734474, 'B': 0.5533186357344739, 'P': 39.44668136426553}

    h = 0.00001
    t_start = 0.002
    times = t_start + h * np.arange(3)

    # Sample at 3 closely-spaced points
    p = [AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
         for t in times
         ]
    #print(p)    # [25.243842573299553, 25.289220064086344, 25.334407842935416]
    p_middle = p[1]
    a_middle = A0 - (p_middle - P0)    # From mass conservation
    b_middle = B0 - (p_middle - P0)    # From mass conservation

    #print(a_middle, b_middle)   44.710779935913656 14.710779935913656
    gradient = np.gradient(p, h)
    derivative_middle = gradient[1]
    rate_at_mid_point = kF * a_middle * b_middle - kR * p_middle  # kF [A] [B] - kR [P]
    assert np.allclose(derivative_middle, rate_at_mid_point)


    # Case where A0 = B0 (and reaction mostly in reverse)
    A0=8
    B0=8
    P0=30
    kF=5
    kR=8

    #print(ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0))
    #Equilibrium at {'A': 7.038367176906169, 'B': 7.038367176906169, 'P': 30.96163282309383}

    h = 0.00005
    t_start = 0.01
    times = t_start + h * np.arange(3)

    # Sample at 3 closely-spaced points
    p = [AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
         for t in times
         ]
    #print(p)    # [30.536666653790505, 30.538373794434957, 30.54007389774443]
    p_middle = p[1]
    a_middle = A0 - (p_middle - P0)    # From mass conservation
    b_middle = B0 - (p_middle - P0)    # From mass conservation

    #print(a_middle, b_middle, p)    # 7.461626205565043 7.461626205565043
    gradient = np.gradient(p, h)
    derivative_middle = gradient[1]
    rate_at_mid_point = kF * a_middle * b_middle - kR * p_middle  # kF [A] [B] - kR [P]
    assert np.allclose(derivative_middle, rate_at_mid_point)



def test_exact_advance_synthesis_reversible_3():
    # Compare the exact solution against a fine-grained forward Euler approximation
    # Reaction A + B <-> P
    sr = SpeciesRegistry()
    rxn_defn = ReactionDefinition(reactants=["A", "B"], products="P",
                                  species_registry=sr, autoregister_species=True,
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF": 5., "kR": 2.})
    rxn = rxn_defn.sim_reactions[0]

    t_final = 0.001
    n_steps = 1800
    t_step = t_final / n_steps
    a = 10
    b = 50
    p = 20

    #print(ReactionKinetics._compute_equilibrium_conc_first_order(kF=5., kR=2., a=1, A0=a, b=1, B0=b, p=1, P0=p))
    #Equilibrium at {'A': 0.2948774087575341, 'B': 40.294877408757536, 'P': 29.705122591242464}

    for i in range(n_steps):
        increment_dict, _ = \
            rxn.step_simulation(delta_time=t_step, conc_dict={"A": a, "B": b, "P": p})
        delta_p = increment_dict["P"]
        p += delta_p
        a -= delta_p
        b -= delta_p

    #print(p)    # 22.130930189693114   Value at t_final, from the fine-grained forward Euler approximation
    exact_p = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=t_final, incremental=False)   # 22.130796453845
    assert np.allclose(p, exact_p)


    # Let's do a second identical round, to reach t_final * 2
    for i in range(n_steps):
        increment_dict, _ = \
            rxn.step_simulation(delta_time=t_step, conc_dict={"A": a, "B": b, "P": p})
        delta_p = increment_dict["P"]
        p += delta_p
        a -= delta_p
        b -= delta_p

    #print(p)    # 23.738906198053474  Value at t_final * 2, from the fine-grained forward Euler approximation
    exact_p = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=2., A0=10., B0=50., P0=20., t=t_final * 2, incremental=False)
    assert np.allclose(p, exact_p)



def test_exact_advance_synthesis_irreversible():
    # Reaction A + B -> P

    # We'll start with `A` as the limiting reagent
    # No change at time 0
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0, incremental=False)
    assert np.allclose(P_t, 20.)

    Delta_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0, incremental=True)
    assert np.allclose(Delta_t, 0)

    # Compare against the REVERSIBLE reaction solver with a near-zero reverse rate constant : not exactly zero, because its
    # implementation reverts to the irreversible solver when the reverse rate constant is very close to zero
    eps = 0.00001
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0.01, incremental=False)
    c = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=10., B0=50., P0=20., t=0.01, incremental=False)
    assert np.allclose(P_t, c)  # 28.8871974442945

    Delta_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0.01, incremental=True)
    D_c = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=10., B0=50., P0=20., t=0.01, incremental=True)
    assert np.allclose(Delta_t, D_c)

    for t  in [0.0003, 0.0004, 0.5, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.008, 0.01, 1.]:
        P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=t, incremental=False)
        c = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=10., B0=50., P0=20., t=t, incremental=False)
        assert np.allclose(P_t, c)


    # A time so large that the reaction has gone to completion (`A` being the limiting reagent)
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=3., incremental=False)
    assert np.allclose(P_t, 30.)    # All `A`(the limiting reagent) converted to `P`

    # A value of t so ridiculously large that an OverflowError is caused (but caught) in the internal math.exp usage
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=200., incremental=False)
    assert np.allclose(P_t, 30.)    # All `A`(the limiting reagent) converted to `P`


    # Now, let `B` be the limiting reagent at large times
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=5., P0=20., t=3., incremental=False)
    assert np.allclose(P_t, 25.)    # All `B`(the limiting reagent) converted to `P`

    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=5., P0=20., t=200., incremental=False)
    assert np.allclose(P_t, 25.)    # All `B`(the limiting reagent) converted to `P`


    # If either A0 or B0 is zero, the reaction doesn't proceed
    Delta_P = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=0, B0=50., P0=20., t=1., incremental=True)
    assert np.allclose(Delta_P, 0)

    Delta_P = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=0, P0=20., t=1., incremental=True)
    assert np.allclose(Delta_P, 0)


    # Case A0 = B0  (equivalently, 2 A -> P)

    # No change at time 0
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30, P0=20, t=0, incremental=False)
    assert np.allclose(P_t, 20)     # No change

    for t  in [0.0003, 0.0004, 0.5, 0.0007, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.008, 0.01, 1.]:
        P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30, P0=20, t=t, incremental=False)
        P_t_approx = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30.000001, P0=20, t=t, incremental=False)
        c = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=5., kR=eps, A0=30, B0=30, P0=20, t=t, incremental=False)
        assert np.allclose(P_t, c)
        assert np.allclose(P_t, P_t_approx)

    # A time so large that the reaction has gone to completion
    P_t = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=30, B0=30, P0=20, t=1000., incremental=False)
    assert np.allclose(P_t, 50.)    # The reagents fully converted to `P`


def test_exact_advance_synthesis_irreversible_2():
    # Verify that the actual (numerically computed) reaction rate matches what's expected from the rate law,
    # at the middle point of 3
    h = 0.00001
    t_start = 0.003
    p1 = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50, P0=20., t=t_start)
    p2 = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50, P0=20., t=t_start + h)
    a2 = 10 - (p2 - 20)    # From mass conservation
    b2 = 50 - (p2 - 20)    # From mass conservation
    p3 = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50, P0=20., t=t_start + 2 * h)
    #print(a2, b2, p2)
    gradient = np.gradient([p1, p2, p3], h)
    derivative_p2 = gradient[1]
    rate_at_mid_point = 5 * a2 * b2   # kF * a2 * b2
    assert np.allclose(derivative_p2, rate_at_mid_point)


def test_exact_advance_synthesis_irreversible_3():
    # Compare the exact solution against a fine-grained forward Euler approximation
    # Reaction A + B -> P
    sr = SpeciesRegistry()
    rxn_defn = ReactionDefinition(reactants=["A", "B"], products="P",
                                  species_registry= sr, autoregister_species=True,
                                  reaction_model="mass action",
                                  kinetic_parameters={"kF": 5., "kR": 0})
    rxn = rxn_defn.sim_reactions[0]

    t_final = 0.01
    n_steps = 10000
    t_step = t_final / n_steps
    a = 10
    b = 50
    p = 20
    for i in range(n_steps):
        increment_dict, _ = \
            rxn.step_simulation(delta_time=t_step, conc_dict={"A": a, "B": b, "P": p})
        delta_p = increment_dict["P"]
        p += delta_p
        a -= delta_p
        b -= delta_p

    exact_p = AnalyticReactionSolver.exact_advance_synthesis_irreversible(kF=5., A0=10., B0=50., P0=20., t=0.01, incremental=False)   # 28.8871974442945
    assert np.allclose(p, exact_p)



def test_approx_solution_synthesis_rxn():
    # Reaction A + B <-> P
    kF=5.
    kR=2.
    A0=10.
    B0=50.
    P0=20.

    with pytest.raises(Exception):
        AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=0, kR=kR, A0=A0, B0=B0, P0=P0, t=0)   # kF cannot be zero

    p = AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=0)
    assert np.allclose(p, P0)      # No change from P0

    #equil = ReactionKinetics._compute_equilibrium_conc_first_order(kF=kF, kR=kR, a=1, A0=A0, p=1, P0=P0, b=1, B0=B0)
    #print(equil)       # {'A': 0.2948774087575341, 'B': 40.294877408757536, 'P': 29.705122591242464}

    p = AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=0.1)
    assert np.allclose(p, 29.705122591242464)       # Reaching the equilibrium value after a sufficiently long time


    # Process a whole array of times in one call
    t = np.array([0, 0.000864, 0.001555, 0.009850, 0.067400])
    result = AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=5., kR=2., A0=10., B0=50., P0=20., t=t)
    assert np.allclose(result, [20., 22.11257633, 23.46603241, 29.11416425, 29.70512254])

    # Compare against exact solution
    for t in [0, 0.000864, 0.001555, 0.009850, 0.067400]:
        p_approx = AnalyticReactionSolver.approx_solution_synthesis_rxn(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
        p_exact = AnalyticReactionSolver.exact_advance_synthesis_reversible(kF=kF, kR=kR, A0=A0, B0=B0, P0=P0, t=t)
        assert abs(p_approx - p_exact) / p_exact < 0.02     # Less than 2% discrepancy



def test_approx_solution_synthesis_rxn_ALT():
    pass    # TODO





########    class VariableTimeSteps    ###########################################################################

def test_set_thresholds():
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
