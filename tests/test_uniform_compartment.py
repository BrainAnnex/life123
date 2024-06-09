import pytest
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from src.modules.chemicals.chem_data import ChemData
from src.modules.reactions.uniform_compartment import UniformCompartment
from src.modules.movies.movies import MovieTabular



def test_set_conc():
    #TODO: test the snapshot argument
    chem_data = ChemData(names=["A", "B", "C"])
    rxn = UniformCompartment(chem_data)

    with pytest.raises(Exception):
        rxn.set_conc(conc=[1, 2, 3, 4])         # Wrong number of entries
    with pytest.raises(Exception):
        rxn.set_conc(conc=[1., -2., -3.])       # Negative values
    with pytest.raises(Exception):
        rxn.set_conc(conc=(10., 20., -0.01))    # Negative values

    rxn.set_conc(conc=[1., 2., 3.])
    assert np.allclose(rxn.system, [1., 2., 3.])

    rxn.set_conc(conc=(10., 20., 30.))
    assert np.allclose(rxn.system, [10., 20., 30.])

    rxn.set_conc(conc={"B": 100})
    assert np.allclose(rxn.system, [10., 100., 30.])


    rxn = UniformCompartment(chem_data)
    rxn.set_conc(conc={"C": 3})
    assert np.allclose(rxn.system, [0., 0., 3.])

    rxn.set_conc(conc={"C": 8, "A": 1, "B": 5})
    assert np.allclose(rxn.system, [1., 5., 8.])

    with pytest.raises(Exception):
        rxn.set_conc(conc={"A": -0.01})



def test_get_system_conc():
    chem_data = ChemData(names=["A", "B", "C"])
    rxn = UniformCompartment(chem_data)
    rxn.set_conc(conc=(10., 20., 30.))

    result = rxn.get_system_conc()
    assert np.allclose(result, [10., 20., 30.])



def test_get_chem_conc():
    chem_data = ChemData(names=["A", "B"])
    dyn = UniformCompartment(chem_data)
    dyn.set_conc(conc=(10., 20.))

    assert np.allclose(dyn.get_chem_conc("A"), 10.)
    assert np.allclose(dyn.get_chem_conc("B"), 20.)
    with pytest.raises(Exception):
        dyn.get_chem_conc("Unknown")    # Non-existent chemical



def test_get_conc_dict():
    chem_data = ChemData(names=["A", "B", "C", "D"])
    rxn = UniformCompartment(chem_data)
    rxn.set_conc(conc=(100, 200, 300, 400))

    result = rxn.get_conc_dict()
    assert result == {"A": 100, "B": 200, "C": 300, "D": 400}

    result = rxn.get_conc_dict(species=["D", "A"])
    assert result == {"A": 100, "D": 400}

    result = rxn.get_conc_dict(species=("C",))      # Tuple with 1 element
    assert result == {"C": 300}

    with pytest.raises(Exception):
        rxn.get_conc_dict(species="C")                      # Wrong data type
        rxn.get_conc_dict(system_data=np.array([1, 2]))     # Wrong number of entries

    result = rxn.get_conc_dict(system_data=np.array([1, 2, 3, 4]))
    assert result == {"A": 1, "B": 2, "C": 3, "D": 4}

    result = rxn.get_conc_dict(species=["B"], system_data=np.array([1, 2, 3, 4]))
    assert result == {"B": 2}



def test_specify_steps():
    rxn = UniformCompartment(None)

    with pytest.raises(Exception):
        # Too few arguments
        rxn.specify_steps()
        rxn.specify_steps(total_duration=15.1)
        rxn.specify_steps(time_step=0.2)
        rxn.specify_steps(n_steps=30)
        # Too many arguments
        rxn.specify_steps(total_duration=15.1, time_step=0.2, n_steps=30)

    assert rxn.specify_steps(time_step=0.5, n_steps=24) == (0.5, 24)
    assert rxn.specify_steps(total_duration=12.0, time_step=0.5) == (0.5, 24)
    assert rxn.specify_steps(total_duration=12.0, n_steps=24) == (0.5, 24)



def test_single_compartment_react():

    # Test based on experiment "cycles_1"
    chem_data = ChemData(names=["A", "B", "C", "E_high", "E_low"])

    # Reaction A <-> B, mostly in forward direction (favored energetically)
    chem_data.add_reaction(reactants="A", products="B",
                           forward_rate=9., reverse_rate=3.)

    # Reaction B <-> C, also favored energetically
    chem_data.add_reaction(reactants="B", products="C",
                           forward_rate=8., reverse_rate=4.)

    # Reaction C + E_High <-> A + E_Low, also favored energetically, but kinetically slow
    chem_data.add_reaction(reactants=["C" , "E_high"], products=["A", "E_low"],
                           forward_rate=1., reverse_rate=0.2)

    initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}

    dynamics = UniformCompartment(chem_data=chem_data)
    dynamics.set_conc(conc=initial_conc, snapshot=False)

    dynamics.set_diagnostics()

    dynamics.set_step_factors(error=0.5)    # Will be used by an excessive first step
                                            # leading to a hard abort

    dynamics.single_compartment_react(initial_step=0.0010, target_end_time=0.0035, variable_steps=False)

    run1 = dynamics.get_system_conc()

    assert np.allclose(run1, [9.69124339e+01, 3.06982519e+00, 1.77408783e-02, 9.99985757e+02, 1.42431633e-02])
    assert np.allclose(dynamics.system_time, 0.0035)
    assert dynamics.explain_time_advance(return_times=True, silent=True) == \
               ([0.0, 0.002, 0.0025, 0.0035],
                [0.001, 0.0005, 0.001])

    # The above computation automatically took 2 normal steps, 1 half-size step and 1 normal step;
    # now repeat the process manually

    dynamics2 = UniformCompartment(chem_data=chem_data)
    dynamics2.set_conc(conc=initial_conc, snapshot=False)

    dynamics2.set_diagnostics()
    dynamics2.single_compartment_react(initial_step=0.0010, n_steps=2, variable_steps=False)
    dynamics2.single_compartment_react(initial_step=0.0005, n_steps=1, variable_steps=False)
    dynamics2.single_compartment_react(initial_step=0.0010, n_steps=1, variable_steps=False)
    run2 = dynamics.get_system_conc()
    assert np.allclose(run2, run1)      # Same result as before
    assert np.allclose(dynamics2.system_time, 0.0035)
    assert dynamics2.explain_time_advance(return_times=True) == \
           ([0.0, 0.002, 0.0025, 0.0035],
            [0.001, 0.0005, 0.001])
    # The time advance is now different, because multiple calls to single_compartment_react() break up the counting
    # (just a convention)



def test_single_reaction_fixed_step():
    chem_data = ChemData(names=["A", "B"])
    rxn = UniformCompartment(chem_data)

    rxn.set_conc(conc=[10., 50.], snapshot=False)

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [ 7. , -7.])
    assert result[0] == - result[1]         # From the stoichiometry


    rxn.clear_reactions()       # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction2"
    chem_data.add_reaction(reactants="A", products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]      # From the stoichiometry


    rxn.clear_reactions()       # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction3"
    chem_data.add_reaction(reactants=[(2,"A",1)], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test_single_reaction_variable_step_1():
    chem_data = ChemData(names=["A", "B"])
    rxn = UniformCompartment(chem_data)

    rxn.set_conc(conc=[10., 50.], snapshot=False)

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [ 7. , -7.])
    assert result[0] == - result[1]         # From the stoichiometry


    rxn.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction2"
    chem_data.add_reaction(reactants=["A"], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]      # From the stoichiometry


    rxn.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction3"
    chem_data.add_reaction(reactants=[(2,"A",1)], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test_single_reaction_step_2():
    chem_data = ChemData(names=["A", "B", "C"])
    rxn = UniformCompartment(chem_data)

    rxn.set_conc(conc=[10., 50., 20.], snapshot=False)

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # # Based on experiment "reactions_single_compartment/react_1"
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [ 7. , -7. , 0.])    # Chemical "C" not participating in this reaction; its delta conc. is 0
    assert result[0] == - result[1]         # From the stoichiometry


    chem_data.clear_reactions_data()   # Re-start with a blank slate of reactions
    # Reaction A + B <-> C , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction4"
    chem_data.add_reaction(reactants=[("A") , ("B")], products=[("C")],
                     forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.002)
    assert np.allclose(result, [-4.92, -4.92, 4.92])
    assert result[0] == result[1]           # From the stoichiometry
    assert result[1] == - result[2]         # From the stoichiometry



def test_single_reaction_step_3():
    chem_data = ChemData(names=["A", "C", "D"])
    rxn = UniformCompartment(chem_data)

    rxn.set_conc(conc=[4., 7., 2.], snapshot=False)

    # Reaction A <-> 2C + D , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction5"
    chem_data.add_reaction(reactants=[("A")], products=[(2, "C", 1) , ("D")],
                     forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.05)
    assert np.allclose(result, [0.4 , -0.8 , -0.4])
    assert result[0] == - result[1] /2    # From the stoichiometry
    assert result[0] == - result[2]       # From the stoichiometry



def test_single_reaction_step_4():
    chem_data = ChemData(names=["A", "B", "C", "D"])
    rxn = UniformCompartment(chem_data)

    rxn.set_conc(conc=[4., 7., 5., 2.], snapshot=False)

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction6"
    chem_data.add_reaction(reactants=[(2,"A",1) , (5,"B",1)], products=[(4,"C",1) , (3,"D",1)],
                     forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.001)
    assert np.allclose(result, [-0.24 , -0.6 , 0.48, 0.36])
    assert  np.allclose(result[0] /2 , result[1] /5)    # From the stoichiometry
    assert  np.allclose(result[1] /5 , -result[2] /4)   # From the stoichiometry
    assert  np.allclose(result[2] /4 , result[3] /3)    # From the stoichiometry



def test_single_reaction_step_5():
    chem_data = ChemData(names=["A", "B"])
    rxn = UniformCompartment(chem_data)

    rxn.set_conc(conc=[3., 5.], snapshot=False)

    # Reaction  2A <-> B , with 2nd-order kinetics in forward reaction, and 1st-order in reverse.
    # Based on experiment "1D/reactions/reaction7"
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)

    result = rxn._reaction_elemental_step(delta_time=0.02)
    assert np.allclose(result, [-1.4 , 0.7])
    assert np.allclose(result[0] /2 , -result[1])      # From the stoichiometry



def test_single_reaction_step_6():
    chem_data = ChemData(names=["A", "B", "C", "D", "E"])
    rxn = UniformCompartment(chem_data)

    rxn.set_conc(conc=[3., 5., 1., 0.4, 0.1], snapshot=False)

    # Coupled reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction8"
    chem_data.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)
    assert chem_data.number_of_reactions() == 2

    result = rxn._reaction_elemental_step(delta_time=0.02)
    assert np.allclose(result, [-1.46 , -1.46  , 1.404 , -0.056 ,  0.056])
    assert np.allclose(result[0] , result[1])                   # From the stoichiometry
    assert np.allclose(result[3] , -result[4])                   # From the stoichiometry
    assert np.allclose(result[0] + result[4], -result[2])       # From the stoichiometry
                                                                # The increase in [A] and [E] combined
                                                                # must match the decrease in [C]


def test_single_compartment_react_variable_steps_1():
    # Based on experiment "variable_steps_1"

    # Initialize the system
    chem_data = ChemData(names=["U", "X", "S"])

    # Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
    chem_data.add_reaction(reactants=[(2, "S", 1)], products="U",
                           forward_rate=8., reverse_rate=2.)

    # Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
    chem_data.add_reaction(reactants="S", products="X",
                           forward_rate=6., reverse_rate=3.)

    dynamics = UniformCompartment(chem_data=chem_data)
    dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})

    dynamics.set_thresholds(norm="norm_A", low=0.25, high=0.64, abort=1.44)
    dynamics.set_thresholds(norm="norm_B", low=None, high=None, abort=None)
    dynamics.set_step_factors(abort=0.5, downshift=0.5, upshift=2.0)

    dynamics.single_compartment_react(initial_step=0.01, target_end_time=0.2,
                                      variable_steps=True)

    df = dynamics.get_history()
    #print(df)
    assert len(df) == 23

    assert np.allclose(df.iloc[0][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.0000,  50.000000,  100.000000,  0.000000])

    assert np.allclose(df.iloc[1][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.0050,  49.500000,  98.500000,  2.500000], rtol=1e-03)     # Notice the halved step size

    assert np.allclose(df.iloc[2][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.0075,  49.302500,  97.798750,  3.596250], rtol=1e-03)

    assert np.allclose(df.iloc[22][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.2050,  55.600598,   69.127620,  19.671183], rtol=1e-03)

    #print(df.iloc[22][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'))



def test_single_compartment_correct_neg_conc_1():
    # Based on "Run 3" of experiment "negative_concentrations_1

    # Initialize the system
    chem_data = ChemData(names=["U", "X", "S"])

    # Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
    chem_data.add_reaction(reactants=[(2, "S", 1)], products="U",
                           forward_rate=8., reverse_rate=2.)

    # Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
    chem_data.add_reaction(reactants="S", products="X",
                           forward_rate=6., reverse_rate=3.)

    dynamics = UniformCompartment(chem_data=chem_data)
    dynamics.set_conc(conc={"U": 50., "X": 100., "S": 0.})

    dynamics.set_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

    dynamics.thresholds = [{"norm": "norm_A", "low": 0.5, "high": 0.8, "abort": 1.44},
                           {"norm": "norm_B", "low": 0.08, "high": 0.5, "abort": 1.5}]
    dynamics.step_factors = {"upshift": 1.5, "downshift": 0.5, "abort": 0.5, "error": 0.5}

    dynamics.single_compartment_react(initial_step=0.1, target_end_time=0.8, variable_steps=False)
    # Note: negative concentrations that would arise from the given step size, get automatically intercepted - and
    #       the step sizes get reduced as needed

    df = dynamics.get_history()
    #print(df)
    assert len(df) == 11

    assert np.allclose(df.iloc[0][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.0000,  50.000000,  100.000000,  0.000000])

    assert np.allclose(df.iloc[1][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.1,  40, 70, 50], rtol=1e-03)

    assert np.allclose(df.iloc[2][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.15,  56, 74.50, 13.5], rtol=1e-03)        # Notice the halved step size

    assert np.allclose(df.iloc[3][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.25,  55.6, 60.25, 28.55], rtol=1e-03)        # Back to the requested step size

    assert np.allclose(df.iloc[5][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.45,  58.7, 45.1465, 37.4535], rtol=1e-03)

    assert np.allclose(df.iloc[6][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.5,  67.8114, 49.610575, 14.766625], rtol=1e-03)        # Notice another instance of halved step size

    assert np.allclose(df.iloc[7][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.6,  66.062420, 43.587378, 24.287783], rtol=1e-03)        # Back to the requested step size

    assert np.allclose(df.iloc[10][['SYSTEM TIME', 'U', 'X', 'S']].to_numpy(dtype='float16'),
                       [0.9,  76.895206, 44.446655,	1.762933], rtol=1e-03)        # Final step




###########################  LOWER-LEVEL METHODS  ###########################


def test_norm_A():
    uc = UniformCompartment()

    delta_conc = np.array([1, 4])
    result = uc.norm_A(delta_conc)
    assert np.allclose(result, 4.25)

    delta_conc = np.array([.5, 2])
    result = uc.norm_A(delta_conc)
    assert np.allclose(result, 1.0625)

    delta_conc = np.array([.5, 2, 1])
    result = uc.norm_A(delta_conc)
    assert np.allclose(result, 0.5833333333)

    delta =    np.array([0.5, 1, 4, -2, -5])
    result = uc.norm_A(delta_conc=delta)
    assert np.allclose(result, 1.85)



def test_norm_B():
    uc = UniformCompartment()

    delta = np.array([4, -1])
    base = np.array([10, 2])
    result = uc.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.5)

    base = np.array([10, 0])
    result = uc.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.4)     # The zero baseline concentration is disregarded

    delta = np.array([300,        6, 0, -1])
    base = np.array([0.00000001, 10, 0,  2])
    result = uc.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.6)

    delta = np.array([300,   6, 0, -1])
    base = np.array([0.001, 10, 0,  2])
    result = uc.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 300000)

    base = np.array([2,   5, 5, 14, 14])
    delta = np.array([0.5, 1, 4, -2, -5])
    result = uc.norm_B(baseline_conc=base, delta_conc=delta)
    assert np.allclose(result, 0.8)

    with pytest.raises(Exception):
        uc.norm_B(baseline_conc=np.array([1, 2, 3]), delta_conc=delta) # Too many entries in array

    with pytest.raises(Exception):
        uc.norm_B(baseline_conc=base, delta_conc=np.array([1]))        # Too few entries in array



def test_norm_C():
    uc = UniformCompartment()

    result = uc.norm_C(prev_conc=np.array([1]), baseline_conc=np.array([2]), delta_conc=np.array([0.5]))
    assert result == 0

    result = uc.norm_C(prev_conc=np.array([8]), baseline_conc=np.array([5]), delta_conc=np.array([1]))
    assert result == 0

    result = uc.norm_C(prev_conc=np.array([8]), baseline_conc=np.array([5]), delta_conc=np.array([4]))
    assert np.isclose(result, 4/3)

    result = uc.norm_C(prev_conc=np.array([10]), baseline_conc=np.array([14]), delta_conc=np.array([-2]))
    assert result == 0

    result = uc.norm_C(prev_conc=np.array([10]), baseline_conc=np.array([14]), delta_conc=np.array([-5]))
    assert np.isclose(result, 5/4)

    prev =     np.array([1,   8, 8, 10, 10])
    baseline = np.array([2,   5, 5, 14, 14])
    delta =    np.array([0.5, 1, 4, -2, -5])
    result = uc.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3 + 5/4)

    # A scenario where the 'prev' and 'baseline' values are almost identical
    prev = np.append(prev, 3)
    baseline = np.append(baseline, 2.999999999)
    delta = np.append(delta, 8)
    result = uc.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3 + 5/4)

    # A scenario where the 'delta' dwarfs the change between 'prev' and 'baseline'
    prev = np.append(prev, 10)
    baseline = np.append(baseline, 10.05)
    delta = np.append(delta, -9)
    result = uc.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3 + 5/4)

    prev = np.append(prev, 10)
    baseline = np.append(baseline, 10.1)
    delta = np.append(delta, -9)
    result = uc.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.isclose(result, 4/3 + 5/4 + 90)



def test_adjust_timestep():
    uc = UniformCompartment(names=["C1", "C2", "C3", "C4", "C5"])
    prev =     np.array([1,   8, 8, 10, 10])
    baseline = np.array([2,   5, 5, 14, 14])
    delta =    np.array([0.5, 1, 4, -2, -5])

    normA = uc.norm_A(delta_conc=delta)
    assert np.allclose(normA, 1.85)

    normB = uc.norm_B(baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normB, 0.8)

    normC = uc.norm_C(prev_conc=prev, baseline_conc=baseline, delta_conc=delta)
    assert np.allclose(normC, 2.583333333333333)        # 4/3 + 5/4

    uc.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.84)
    uc.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.79)
    uc.set_step_factors(upshift=1.2, downshift=0.5, abort=0.4, error=0.25)

    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85}, 'applicable_norms': 'norm_A'}

    uc.set_thresholds(norm="norm_A", low=0.5, high=0.8, abort=1.86)     # normA (1.85) no longer triggers abort
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'high', 'step_factor': 0.5, 'norms': {'norm_A': 1.85}, 'applicable_norms': 'norm_A'}

    uc.set_thresholds(norm="norm_A", low=0.5, high=1.86, abort=1.87)    # normA (1.85) no longer triggers high nor abort
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'norm_B'}

    uc.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.81)    # normB (0.8) no longer triggers abort
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'high', 'step_factor': 0.5, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'norm_B'}

    uc.set_thresholds(norm="norm_B", low=0.08, high=0.81, abort=0.81)    # normB (0.8) no longer triggers high nor abort
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'step_factor': 1, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'ALL'}

    uc.set_thresholds(norm="norm_A", low=1.86, high=1.87, abort=1.87)   # normA (1.85) will now trigger a "low"
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'step_factor': 1, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'ALL'}
                      # We're still on the 'stay' action because we aren't below ALL the thresholds

    uc.set_thresholds(norm="norm_B", low=0.81, high=0.82, abort=0.83)   # normB (0.8) will now trigger a "low", too
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'low', 'step_factor': 1.2, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'ALL'}

    uc.set_thresholds(norm="norm_C", low=2.59, high=2.60, abort=2.61)   # normC (2.583) will still continue to trigger a "low"
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'low', 'step_factor': 1.2, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 2.583333333333333}, 'applicable_norms': 'ALL'}

    uc.set_thresholds(norm="norm_C", low=2.58, high=2.60, abort=2.61)   # normC (2.583) will no longer trigger a "low" - but not a "high" nor an "abort"
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'stay', 'step_factor': 1, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 2.583333333333333}, 'applicable_norms': 'ALL'}

    uc.set_thresholds(norm="norm_C", low=2.57, high=2.58, abort=2.61)   # normC (2.583) will now trigger a "high"
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'high', 'step_factor': 0.5, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 2.583333333333333}, 'applicable_norms': 'norm_C'}

    uc.set_thresholds(norm="norm_C", low=2.56, high=2.57, abort=2.58)   # normC (2.583) will now trigger an "abort"
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8, 'norm_C': 2.583333333333333}, 'applicable_norms': 'norm_C'}

    uc.set_thresholds(norm="norm_B", low=0.08, high=0.5, abort=0.79)    # normB (0.8) will now trigger an "abort" before we can even get to normC
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85, 'norm_B': 0.8}, 'applicable_norms': 'norm_B'}

    uc.set_thresholds(norm="norm_A", low=0.5, high=1.0, abort=1.84)    # normA (1.85) will now trigger an "abort" before we can even get to normB
    result = uc.adjust_timestep(delta_conc=delta, baseline_conc=baseline, prev_conc=prev)
    assert result == {'action': 'abort', 'step_factor': 0.4, 'norms': {'norm_A': 1.85}, 'applicable_norms': 'norm_A'}



def test_compute_all_rate_deltas():
    chem_data = ChemData(names=["A", "B", "C", "D"])
    rxn = UniformCompartment(chem_data)

    # Start with reaction A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=20., reverse_rate=2.)  # Reaction 0
    rxn.set_conc(conc=[5., 8., 0, 0], snapshot=False)

    result = rxn.compute_all_reaction_deltas(delta_time=0.5)    # {0: 42.0}
    assert len(result) == 1
    assert np.allclose(result[0], 42.0)

    # Add reaction 2B <-> 3C , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=[(2, "B", 1)], products=[(3, "C", 1)],
                           forward_rate=10., reverse_rate=25.)   # Rxn 1
    rxn.set_conc(conc=[5., 8., 15., 0], snapshot=False)
    result = rxn.compute_all_reaction_deltas(delta_time=0.5) # {0: 42.0, 1: -147.5}
    assert len(result) == 2
    assert np.allclose(result[0], 42.0)
    assert np.allclose(result[1], -147.5)

    # Add reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[(2,"A",1) , (5,"B",1)], products=[(4,"C",1) , (3,"D",1)],
                     forward_rate=5., reverse_rate=2.)          # Rxn 2
    rxn.set_conc(conc=[5., 8., 15.,  7.], snapshot=False)
    result = rxn.compute_all_reaction_deltas(delta_time=0.5)    # {0: 42.0, 1: -147.5, 2: -5.0}
    assert len(result) == 3
    assert np.allclose(result[0], 42.0)
    assert np.allclose(result[1], -147.5)
    assert np.allclose(result[2], -5)


    # Add reaction  2A <-> B , with 2nd-order kinetics in the forward direction
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=3., reverse_rate=2.)
    result = rxn.compute_all_reaction_deltas(delta_time=0.5)    # {0: 42.0, 1: -147.5, 2: -5.0, 3: 29.5}
    assert len(result) == 4
    assert np.allclose(result[0], 42.0)
    assert np.allclose(result[1], -147.5)
    assert np.allclose(result[2], -5)
    assert np.allclose(result[3], 29.5)

    # This time, only process reactions 0 and 2
    result_0_2 = rxn.compute_all_reaction_deltas(delta_time=0.5,
                                                 rxn_list=[0, 2])   # {0: 42.0, 2: -5.0}
    assert len(result_0_2) == 2
    assert np.allclose(result_0_2[0], 42.0)
    assert np.allclose(result_0_2[2], -5.0)


    # FLUSH OUT ALL REACTIONS (to start over)
    rxn.clear_reactions()
    # Start with reaction A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=20., reverse_rate=2.)   # Rxn 0
    rxn.set_conc(conc=[5., 8., 0, 0], snapshot=False)
    result = rxn.compute_all_reaction_deltas(delta_time=0.25)   #  {0: 21.0}
    assert len(result) == 1
    assert np.allclose(result[0], 21.0)



def test_compute_reaction_delta_rate():
    chem_data = ChemData(names=["A", "B", "C", "D"])
    dynamics = UniformCompartment(chem_data)
    dynamics.set_conc(conc=[5., 8., 0, 0], snapshot=False)

    # Reaction A <-> B , with 1st-order kinetics in both directions
    rxn_index = chem_data.add_reaction(reactants="A", products="B", forward_rate=20., reverse_rate=2.)
    result = dynamics.compute_reaction_delta_rate(rxn=chem_data.get_reaction(rxn_index))
    assert np.allclose(result, 20. * 5. - 2. * 8.)

    # Reaction 5A <-> 2B , with 1st-order kinetics in both directions.
    # Same as before, but different stoichiometry (which does NOT influence the result)
    rxn_index = chem_data.add_reaction(reactants=[(5, "A", 1)], products=[(2, "B", 1)],
                                       forward_rate=20., reverse_rate=2.)
    result = dynamics.compute_reaction_delta_rate(rxn=chem_data.get_reaction(rxn_index))
    assert np.allclose(result, 20. * 5. - 2. * 8.)

    # Reaction C <-> D , with 1st-order kinetics in both directions
    rxn_index = chem_data.add_reaction(reactants="C", products="D", forward_rate=20., reverse_rate=2.)
    dynamics.set_conc(conc=[0., 0., 5., 8.], snapshot=False)
    result = dynamics.compute_reaction_delta_rate(rxn=chem_data.get_reaction(rxn_index))
    assert np.allclose(result, 20. * 5. - 2. * 8.)

    # Reaction 2B <-> 3C , with 1st-order kinetics in both directions
    rxn_index = chem_data.add_reaction(reactants=[(2, "B", 1)], products=[(3, "C", 1)],
                                       forward_rate=10., reverse_rate=25.)
    dynamics.set_conc(conc=[0., 8., 15., 0.], snapshot=False)
    result = dynamics.compute_reaction_delta_rate(rxn=chem_data.get_reaction(rxn_index))
    assert np.allclose(result,  10. * 8. - 25. * 15.)

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    rxn_index = chem_data.add_reaction(reactants=[(2,"A",1) , (5,"B",1)], products=[(4,"C",1) , (3,"D",1)],
                                 forward_rate=5., reverse_rate=2.)
    dynamics.set_conc(conc=[3.5, 9., 11., 7.], snapshot=False)
    result = dynamics.compute_reaction_delta_rate(rxn=chem_data.get_reaction(rxn_index))
    assert np.allclose(result,  5. * 3.5 * 9. - 2. * 11. * 7.)

    # Reaction  2A <-> B , with 2nd-ORDER kinetics in the forward direction
    rxn_index = chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    dynamics.set_conc(conc=[4.5, 6., 0., 0.], snapshot=False)
    result = dynamics.compute_reaction_delta_rate(rxn=chem_data.get_reaction(rxn_index))
    assert np.allclose(result, 5. * 4.5 **2 - 2. * 6.)

    # Reaction  B <-> 2C , with 2nd-ORDER kinetics in the reverse direction
    rxn_index = chem_data.add_reaction(reactants=[("B")], products=[(2, "C", 2)], forward_rate=4., reverse_rate=2.)
    dynamics.set_conc(conc=[0., 5., 4, 0.], snapshot=False)
    result = dynamics.compute_reaction_delta_rate(rxn=chem_data.get_reaction(rxn_index))
    assert np.allclose(result, 4. * 5. - 2. * 4. **2)



def test_is_in_equilibrium():
    chem_data = ChemData(names=["A", "B", "C", "D", "E", "F"])
    rxn = UniformCompartment(chem_data)

    # Reaction 0 : A <-> B
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    c = {'A': 23.9931640625, 'B': 36.0068359375}
    assert rxn.is_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)
    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=1)      # Testing ALL reactions

    # Reaction 1 : A <-> F
    chem_data.add_reaction(reactants=["A"], products=["F"], forward_rate=20, reverse_rate=2.)
    c = {'A': 3, 'F': 32.999}
    assert rxn.is_in_equilibrium(rxn_index=1, conc=c, explain=False, tolerance=10)   # The deviation is just below the 10% tolerance

    with pytest.raises(Exception):
        assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=10) # Testing ALL reactions, but neglected to provide [B]

    c = {'A': 3, 'B': 4.6, 'F': 33.001}   # We're now including all the concentrations we need for both reactions
    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=10) \
           == {False: [1]}    # The deviation for reaction 1 is just above the 10% tolerance

    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=2.2) \
           == {False: [0, 1]}      # With a tolerance so low, not reaction meets the criteria

    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=2.3) \
           == {False: [1]}                # Reaction 0 barely passes with this tolerance (it deviates by 2.222 %)

    with pytest.raises(Exception):
        assert rxn.is_in_equilibrium(explain=False, tolerance=2.3)  # We're failing to provide concentrations

    rxn.set_conc(conc=[3, 4.6, 0, 0, 0, 33.001])    # The concentrations are in the same order as the declared chemicals
    assert rxn.is_in_equilibrium(conc=c, explain=False, tolerance=2.3) \
           == {False: [1]}        # Now using the System concentrations



def test_reaction_in_equilibrium():
    chem_data = ChemData(names=["A", "B", "C", "D", "E", "F"])
    rxn = UniformCompartment(chem_data)

    # Reaction 0 : A <-> B
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    c = {'A': 23.9931640625, 'B': 36.0068359375}
    assert rxn.reaction_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)

    # Reaction 1 : A <-> F
    chem_data.add_reaction(reactants=["A"], products=["F"], forward_rate=20, reverse_rate=2.)
    c = {'A': 3, 'F': 32.999}
    assert rxn.reaction_in_equilibrium(rxn_index = 1, conc=c, explain=False, tolerance=10)   # Just below the 10% tolerance
    c = {'A': 3, 'F': 33.001}
    assert not rxn.reaction_in_equilibrium(rxn_index = 1, conc=c, explain=False, tolerance=10)  # Just above the 10% tolerance

    # Reaction 2 : A <-> 3B
    chem_data.add_reaction(reactants=["A"], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)
    c = {'A': 14.54545455, 'B': 36.36363636}
    assert rxn.reaction_in_equilibrium(rxn_index = 2, conc=c, explain=False, tolerance=1)

    # Reaction 3:  2A <-> 3B
    chem_data.add_reaction(reactants=[(2,"A",1)], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)
    c = {'A': 16.25, 'B': 40.625}
    assert rxn.reaction_in_equilibrium(rxn_index = 3, conc=c, explain=False, tolerance=1)

    # Reaction 4:  A + B <-> C , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=["A" , "B"], products=[("C")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 0.29487741, 'B': 40.29487741, 'C': 29.70512259}
    assert rxn.reaction_in_equilibrium(rxn_index = 4, conc=c, explain=False, tolerance=1)

    # Reaction 5:  A <-> 2C + D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=["A"], products=[(2, "C", 1) , ("D")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 4.31058733, 'C': 6.37882534, 'D': 1.68941267}
    assert rxn.reaction_in_equilibrium(rxn_index = 5, conc=c, explain=False, tolerance=1)

    # Reaction 6:  2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    chem_data.add_reaction(reactants=[(2,"A",1) , (5,"B",1)], products=[(4,"C",1) , (3,"D",1)],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 2.80284552, 'B': 4.00711381, 'C': 7.39430896, 'D': 3.79573172}
    assert rxn.reaction_in_equilibrium(rxn_index = 6, conc=c, explain=False, tolerance=1)


    rxn.clear_reactions()   # This will reset the reaction count to 0

    # Reaction 0:  2A <-> B , with 1st-order kinetics in both directions
    chem_data.add_reaction(reactants=[(2, "A", 1)], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 2.16928427, 'B': 5.41535786}
    assert rxn.reaction_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)

    # Reaction 1:  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction
    chem_data.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 1.51554944, 'B': 5.74222528}
    assert rxn.reaction_in_equilibrium(rxn_index = 1, conc=c, tolerance=1, explain=False)

    # Reaction 2:  A + B <-> C + D
    chem_data.add_reaction(reactants=["A", "B"], products=["C", "D"], forward_rate=5., reverse_rate=2.)

    # with zero concentrations of a reactant, it won't be an equilibrium...
    c = {'A': 15.2, 'B': 0, 'C': 21.3, 'D': 4.1}
    assert not rxn.reaction_in_equilibrium(rxn_index = 2, conc=c, tolerance=10, explain=False)

    # likewise with zero concentrations of a reaction product, it won't be an equilibrium...
    c = {'A': 15.2, 'B': 21.3, 'C': 0 , 'D': 4.1}
    assert not rxn.reaction_in_equilibrium(rxn_index = 2, conc=c, tolerance=10, explain=False)

    # but with zero concentrations of both a reactant and a reaction product, the reaction is stuck in place -
    # so, in equilibrium
    c = {'A': 15.2, 'B': 0, 'C': 0 , 'D': 4.1}
    assert rxn.reaction_in_equilibrium(rxn_index = 2, conc=c, tolerance=1, explain=False)



def test_validate_increment():
    chem_data = ChemData(names=["A", "B", "C"])
    rxn = UniformCompartment(chem_data)

    chem_data.add_reaction(reactants=["A"], products=["B"])

    rxn.validate_increment(delta_conc=50., baseline_conc=10., rxn_index=0, species_index=2, delta_time=0.02)
    rxn.validate_increment(delta_conc=-9.99, baseline_conc=10., rxn_index=0, species_index=2, delta_time=0.02)

    with pytest.raises(Exception):      # Would lead to a negative concentration
        rxn.validate_increment(delta_conc=-10.01, baseline_conc=10., rxn_index=0, species_index=2, delta_time=0.02)



def test_stoichiometry_checker():
    chem = ChemData(names=["A", "B", "C", "D"])
    rxn = UniformCompartment(chem)

    chem.add_reaction(reactants=["A"], products=["B"])          # Reaction 0:   A <--> B

    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 50, 0, 0]), conc_arr_after=np.array([10, 40, 0, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 50, 0, 0]), conc_arr_after=np.array([10, 39.9, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 0, 0, 0]), conc_arr_after=np.array([90, 10, 0, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=[(2, "A")], products=["B"])         # Reaction 0:   2A <--> B

    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 0, 0, 0]), conc_arr_after=np.array([80, 10, 0, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 0, 0, 0]), conc_arr_after=np.array([80, 10.1, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 50, 0, 0]), conc_arr_after=np.array([10, 45, 0, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=["A", "B"], products=["C"])         # Reaction 0:   A + B <--> C

    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 40, 10, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 40, 9.9, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=["A", (3, "B")], products=["C"])     # Reaction 0:   A + 3B <--> C

    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 10, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 9.9, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=["A", (3, "B")], products=[(4, "C")])                   # Reaction 0:   A + 3B <--> 4C

    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 40, 0]))
    assert not rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 39.9, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=[(2, "A"), (3, "B")], products=[(4, "C"), (5, "D")])     # Reaction 0:   2A + 3B <--> 4C + 5D
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 100, 100, 100]), conc_arr_after=np.array([120, 130, 60, 50]))
    assert not rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 100, 100, 100]), conc_arr_after=np.array([120.1, 130, 60, 50]))
    assert rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 100, 100, 100]), conc_arr_after=np.array([80, 70, 140, 150]))
    assert not rxn.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 100, 100, 100.1]), conc_arr_after=np.array([80, 70, 140, 150]))



def test_sigmoid():
    rxn = UniformCompartment(None)

    # When the concentration matches Kd, the occupancy is 1/2 by definition
    assert np.allclose(rxn.sigmoid(conc=1., Kd=1.), 0.5)
    assert np.allclose(rxn.sigmoid(conc=0.1, Kd=0.1), 0.5)
    assert np.allclose(rxn.sigmoid(conc=10., Kd=10.), 0.5)

    # When the concentration is 10% of Kd, the occupancy is about 0.1 ;
    # when the concentration is 10 times Kd, the occupancy is about 0.9
    assert np.allclose(rxn.sigmoid(conc=0.1, Kd=1.), 0.1)
    assert np.allclose(rxn.sigmoid(conc=10., Kd=1.), 0.9)

    assert np.allclose(rxn.sigmoid(conc=10., Kd=100.), 0.1)
    assert np.allclose(rxn.sigmoid(conc=1000., Kd=100.), 0.9)

    assert np.allclose(rxn.sigmoid(conc=0.001, Kd=0.01), 0.1)
    assert np.allclose(rxn.sigmoid(conc=0.1, Kd=0.01), 0.9)

    #print(rxn.sigmoid(conc=0.1, Kd=1.))



def test_logistic():
    rxn = UniformCompartment(None)

    #print(rxn.logistic(-10))
    assert np.allclose(rxn.logistic(-100), 0)
    assert np.allclose(rxn.logistic(-10), 0.00004539786)
    assert np.allclose(rxn.logistic(0), 0.5)
    assert np.allclose(rxn.logistic(10), 0.99995460213)
    assert np.allclose(rxn.logistic(100), 1)



def test_set_macromolecules():
    chem = ChemData(names=["A", "B"])
    chem.add_macromolecules(["M1", "M2"])
    chem.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=3)
    chem.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=5)
    chem.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="B", Kd=11)
    chem.set_binding_site_affinity(macromolecule="M2", site_number=3, ligand="B", Kd=102)

    rxn = UniformCompartment(chem)
    rxn.set_macromolecules()    # By default, set counts to 1 for all the registered macromolecules

    assert rxn.macro_system == {"M1": 1, "M2": 1}
    assert len(rxn.macro_system_state) == 2
    assert rxn.macro_system_state["M1"] == {1: ("A", 0.), 2: ("B", 0.)}
    assert rxn.macro_system_state["M2"] == {1: ("B", 0.), 3: ("B", 0.)}


    # Over-write the previous settings
    rxn.set_macromolecules({"M2": 4})

    assert rxn.macro_system == {"M2": 4}
    assert len(rxn.macro_system_state) == 1
    assert rxn.macro_system_state["M2"] == {1: ("B", 0.), 3: ("B", 0.)}


    with pytest.raises(Exception):
        rxn.set_macromolecules({"M999": 2})        # Unknown macromolecule

    chem.add_macromolecules("M999")

    # Over-write the previous settings
    rxn.set_macromolecules({"M999": 2})
    assert rxn.macro_system == {"M999": 2}
    assert len(rxn.macro_system_state) == 1
    assert rxn.macro_system_state["M999"] == {}     # No binding sites were registered

    chem.set_binding_site_affinity(macromolecule="M999", site_number=12, ligand="A", Kd=4.5)
    assert rxn.macro_system_state["M999"] == {}     # The system state has not been updated yet
    # Over-write the previous settings
    rxn.set_macromolecules({"M999": 2})
    assert rxn.macro_system_state["M999"] == {12: ("A", 0.)}



def test_set_occupancy():
    chem_data = ChemData(names=["A", "B"])
    chem_data.add_macromolecules(["M1", "M2"])

    rxn = UniformCompartment(chem_data)

    with pytest.raises(Exception):
        # Occupancy out of range
        rxn.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=1.1)

    with pytest.raises(Exception):
        # Occupancy out of range
        rxn.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=-0.3)

    with pytest.raises(Exception):
        # Unknown macromolecule
        rxn.set_occupancy(macromolecule="Unknown", site_number=1, fractional_occupancy=0.5)

    with pytest.raises(Exception):
        # No binding sites are defined on macromolecule M1
        rxn.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=0.5)

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="B", Kd=6.)
    rxn.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=0.5)

    assert rxn.macro_system == {"M1": 1}
    assert rxn.macro_system_state == {"M1": {1: ("B", 0.5)}}

    assert rxn.get_occupancy(macromolecule="M1", site_number=1) == 0.5



def test_get_occupancy():
    chem_data = ChemData(names=["A", "B", "C"])
    chem_data.add_macromolecules(["M1", "M2"])

    rxn = UniformCompartment(chem_data)

    with pytest.raises(Exception):
        # The system state for macromolecules has not been set yet
        rxn.get_occupancy(macromolecule="M1", site_number=1)

    rxn.set_macromolecules()

    with pytest.raises(Exception):
        # No occupancy data yet set for macromolecule `Unknown`
        rxn.get_occupancy(macromolecule="Unknown", site_number=1)

    with pytest.raises(Exception):
        # No occupancy data yet set for site number 1 of macromolecule `M1`
        rxn.get_occupancy(macromolecule="M1", site_number=1)

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="C", Kd=23.5)
    rxn.set_occupancy(macromolecule="M1", site_number=1, fractional_occupancy=0.5)

    assert rxn.macro_system == {"M1": 1, "M2": 1}
    assert rxn.macro_system_state == {"M1": {1: ("C", 0.5)}, "M2": {} }

    assert rxn.get_occupancy(macromolecule="M1", site_number=1) == 0.5



def test_update_occupancy():
    chem_data = ChemData(names=["A", "B", "C"])
    chem_data.add_macromolecules(["M1", "M2"])

    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=1, ligand="A", Kd=10)
    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=2, ligand="B", Kd=20)
    chem_data.set_binding_site_affinity(macromolecule="M1", site_number=3, ligand="C", Kd=30)

    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=1, ligand="C", Kd=3)
    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=2, ligand="C", Kd=30)
    chem_data.set_binding_site_affinity(macromolecule="M2", site_number=3, ligand="C", Kd=300)


    rxn = UniformCompartment(chem_data)

    rxn.set_macromolecules()
    assert rxn.macro_system == {"M1": 1, "M2": 1}

    rxn.set_conc([10, 20, 30])      # Set the concentrations, respectively, of A, B and C

    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=1) , 0.)
    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=2) , 0.)
    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=3) , 0.)

    rxn.update_occupancy()

    # All fractional occupancies will be 1/2, because the concentrations of the ligands below
    # exactly match their binding affinities
    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=1) , 0.5)
    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=2) , 0.5)
    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=3) , 0.5)

    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=1) , 0.9)  # Ligand conc is 10x binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=2) , 0.5)  # Ligand conc = binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=3) , 0.1)  # Ligand conc is 1/10 binding affinity


    # Vary the concentration of ligand C, starting with zero
    rxn.set_single_conc(conc=0, species_name="C")     # No ligand C
    rxn.update_occupancy()

    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=1) , 0.5)  # Unaffected (different ligand)
    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=2) , 0.5)  # Unaffected (different ligand)
    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=3) , 0)    # No occupancy in absence of ligand

    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=1) , 0)    # No occupancy in absence of ligand
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=2) , 0)    # No occupancy in absence of ligand
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=3) , 0)    # No occupancy in absence of ligand


    # Very low concentration of ligand C
    rxn.set_single_conc(conc=0.3, species_name="C")
    rxn.update_occupancy()

    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=1) , 0.5)      # Unaffected (different ligand)
    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=2) , 0.5)      # Unaffected (different ligand)
    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=3) , 0.012195) # Ligand conc is 1/100 binding affinity

    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=1) , 0.1)      # Ligand conc is 1/10 binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=2) , 0.012195) # Ligand conc is 1/100 binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=3) , 0.00136986)  # Ligand conc is 1/1000 binding affinity


    # Low concentration of ligand C
    rxn.set_single_conc(conc=3, species_name="C")
    rxn.update_occupancy()

    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=3) , 0.1)      # Ligand conc is 1/10 binding affinity

    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=1) , 0.5)      # Ligand conc = binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=2) , 0.1)      # Ligand conc is 1/10 binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=3) , 0.012195) # Ligand conc is 1/100 binding affinity


    # Mid concentration of ligand C
    rxn.set_single_conc(conc=30, species_name="C")
    rxn.update_occupancy()

    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=3) , 0.5)  # Ligand conc = binding affinity

    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=1) , 0.9)  # Ligand conc is 10x binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=2) , 0.5)  # Ligand conc = binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=3) , 0.1)  # Ligand conc is 1/10 binding affinity


    # High concentration of ligand C
    rxn.set_single_conc(conc=300, species_name="C")
    rxn.update_occupancy()

    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=3) , 0.9)       # Ligand conc is 10x binding affinity

    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=1) , 0.9878049)    # Ligand conc is 100x binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=2) , 0.9)          # Ligand conc is 10x binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=3) , 0.5)          # Ligand conc = binding affinity


    # Very High concentration of ligand C
    rxn.set_single_conc(conc=3000, species_name="C")
    rxn.update_occupancy()

    assert np.allclose(rxn.get_occupancy(macromolecule="M1", site_number=3) , 0.9878049)    # Ligand conc is 100x binding affinity

    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=1) , 0.99863)      # Ligand conc is 1000x binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=2) , 0.9878049)    # Ligand conc is 100x binding affinity
    assert np.allclose(rxn.get_occupancy(macromolecule="M2", site_number=3) , 0.9)          # Ligand conc is 10x binding affinity

    #print(rxn.get_occupancy(macromolecule="M2", site_number=1))
    #print(rxn.get_occupancy(macromolecule="M2", site_number=2))
    #print(rxn.get_occupancy(macromolecule="M2", site_number=3))



def test_explain_time_advance():
    rxn = UniformCompartment(None)

    with pytest.raises(Exception):
        rxn.explain_time_advance(return_times=True)      # Diagnostics weren't first enabled

    rxn.set_diagnostics()

    with pytest.raises(Exception):
        rxn.explain_time_advance(return_times=True)      # No diagnostic data yet present

    # Start out with uniform steps
    rxn.diagnostic_conc_data.store(par=20.,
                                   data_snapshot={"time_step": 100.})

    assert rxn.explain_time_advance(return_times=True, silent=True) is None


    rxn.diagnostic_conc_data.store(par=30.,
                                   data_snapshot={"time_step": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)    # TODO: also test the returned step sizes
    assert np.allclose(result, [20., 30.])

    rxn.diagnostic_conc_data.store(par=40.,
                                   data_snapshot={"time_step": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40.])

    # Switching to smaller step
    rxn.diagnostic_conc_data.store(par=45.,
                                   data_snapshot={"time_step": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 45.])

    rxn.diagnostic_conc_data.store(par=50.,
                                   data_snapshot={"time_step": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50.])

    # Switching to larger step
    rxn.diagnostic_conc_data.store(par=70.,
                                   data_snapshot={"time_step": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70.])

    # Yet larger
    rxn.diagnostic_conc_data.store(par=95.,
                                   data_snapshot={"time_step": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95.])

    # Smaller again
    rxn.diagnostic_conc_data.store(par=96.,
                                   data_snapshot={"time_step": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95., 96.])

    rxn.diagnostic_conc_data.store(par=97.,
                                   data_snapshot={"time_step": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95., 97.])

    rxn.diagnostic_conc_data.store(par=98.,
                                   data_snapshot={"time_step": 100.})
    result, _ = rxn.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95., 98.])

    #print(rxn.diagnostic_data_baselines.get())
    #print(result)



def test__delta_names():
    chem_data = ChemData(names=["A", "B", "X"])
    dyn = UniformCompartment(chem_data)

    assert dyn._delta_names() == ["Delta A", "Delta B", "Delta X"]



def test__delta_conc_dict():
    chem_data = ChemData(names=["A", "B", "X"])
    dyn = UniformCompartment(chem_data)

    assert dyn._delta_conc_dict(np.array([10, 20, 30])) == \
           {"Delta A": 10, "Delta B": 20, "Delta X": 30}

    with pytest.raises(Exception):
        dyn._delta_conc_dict(np.array([10, 20, 30, 40]))    # One element too many



def test_save_diagnostic_rxn_data():
    chem_data = ChemData(names=["A", "B", "C", "X"])
    # Add 3 reactions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A"], products=["X"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A", "B"], products=["X"], forward_rate=5., reverse_rate=2.)

    dyn = UniformCompartment(chem_data)
    assert len(dyn.diagnostic_rxn_data) == 0

    dyn.system_time = 100
    with pytest.raises(Exception):
        dyn.save_diagnostic_rxn_data(rxn_index=0, time_step=4,
                                     increment_vector_single_rxn=np.array([2, -2])) # Wrong size of Numpy array

    # Add data for reaction 0
    dyn.save_diagnostic_rxn_data(rxn_index=0, time_step=4,
                                 increment_vector_single_rxn=np.array([2, -2, 0, 0]))

    assert len(dyn.diagnostic_rxn_data) == 3    # Same as the number of reactions

    diagnostic_data_rxn_0 = dyn.diagnostic_rxn_data[0]
    diagnostic_data_rxn_1 = dyn.diagnostic_rxn_data[1]
    diagnostic_data_rxn_2 = dyn.diagnostic_rxn_data[2]

    assert (type(diagnostic_data_rxn_0)) == MovieTabular
    assert (type(diagnostic_data_rxn_1)) == MovieTabular
    assert (type(diagnostic_data_rxn_2)) == MovieTabular

    df_0 = diagnostic_data_rxn_0.get_dataframe()
    df_1 = diagnostic_data_rxn_1.get_dataframe()
    df_2 = diagnostic_data_rxn_2.get_dataframe()

    expected_df_0 = pd.DataFrame([[100, 2, -2, 0, 0, 4, ""]],
                                 columns = ["START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    assert_frame_equal(df_0, expected_df_0, check_dtype=False)

    assert df_1.empty
    assert df_2.empty

    # Add data for reaction 1
    dyn.save_diagnostic_rxn_data(rxn_index=1, time_step=4,
                                 increment_vector_single_rxn=np.array([7, 0, 0, -7]))

    df_0 = diagnostic_data_rxn_0.get_dataframe()
    df_1 = diagnostic_data_rxn_1.get_dataframe()
    df_2 = diagnostic_data_rxn_2.get_dataframe()

    assert_frame_equal(df_0, expected_df_0, check_dtype=False)  # Nothing was done to df_0 by processing reaction index 1

    expected_df_1 = pd.DataFrame([[100, 7, 0, 0, -7, 4, ""]],
                                 columns = ["START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    assert_frame_equal(df_1, expected_df_1, check_dtype=False)

    assert df_2.empty

    # Add data for reaction 2
    dyn.save_diagnostic_rxn_data(rxn_index=2, time_step=4,
                                 increment_vector_single_rxn=np.array([-8, -8, 0, 8]),
                                 caption="I'm a caption")

    df_0 = diagnostic_data_rxn_0.get_dataframe()
    df_1 = diagnostic_data_rxn_1.get_dataframe()
    df_2 = diagnostic_data_rxn_2.get_dataframe()

    assert_frame_equal(df_0, expected_df_0, check_dtype=False)  # Nothing was done to df_0 by processing reaction index 1
    assert_frame_equal(df_1, expected_df_1, check_dtype=False)  # Nothing was done to df_1 by processing reaction index 1

    expected_df_2 = pd.DataFrame([[100, -8, -8, 0, 8, 4, "I'm a caption"]],
                                 columns = ["START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    #TODO: test adding multiple entries for any reaction



def test_get_diagnostic_rxn_data():
    chem_data = ChemData(names=["A", "B", "C", "X"])
    # Add 3 reactions
    chem_data.add_reaction(reactants=["A"], products=["B"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A"], products=["X"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A", "B"], products=["X"], forward_rate=5., reverse_rate=2.)

    dyn = UniformCompartment(chem_data)

    dyn.system_time = 100

    # Add data for reaction 0
    dyn.save_diagnostic_rxn_data(rxn_index=0, time_step=4,
                                 increment_vector_single_rxn=np.array([2, -2, 0, 0]))

    df_0 = dyn.get_diagnostic_rxn_data(rxn_index=0, print_reaction=False)

    expected_df_0 = pd.DataFrame([[100, 2, -2, 0, 0, 4, ""]],
                                 columns = ["START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    assert_frame_equal(df_0, expected_df_0, check_dtype=False)

    df_1 = dyn.get_diagnostic_rxn_data(rxn_index=1, print_reaction=False)
    assert df_1.empty
    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False)
    assert df_2.empty

    # Add data for reaction 1
    dyn.save_diagnostic_rxn_data(rxn_index=1, time_step=4,
                                 increment_vector_single_rxn=np.array([7, 0, 0, -7]))

    df_1 = dyn.get_diagnostic_rxn_data(rxn_index=1, print_reaction=False)

    expected_df_1 = pd.DataFrame([[100, 7, 0, 0, -7, 4, ""]],
                                 columns = ["START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    assert_frame_equal(df_1, expected_df_1, check_dtype=False)

    df_0 = dyn.get_diagnostic_rxn_data(rxn_index=0, print_reaction=False)
    assert_frame_equal(df_0, expected_df_0, check_dtype=False)      # No change made to df_0 from the last step

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False)
    assert df_2.empty

    # Add data for reaction 2
    dyn.save_diagnostic_rxn_data(rxn_index=2, time_step=4,
                                 increment_vector_single_rxn=np.array([-8, -8, 0, 8]),
                                 caption="I'm a caption")

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False, tail=1) # With just one row, tail=1 won't make a difference

    expected_df_2 = pd.DataFrame([[100, -8, -8, 0, 8, 4, "I'm a caption"]],
                                 columns = ["START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False, t=50) # With just one row, the time selector won't matter

    expected_df_2 = pd.DataFrame([[50, 100, -8, -8, 0, 8, 4, "I'm a caption"]],
                                 columns = ["search_value", "START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])

    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    # Add a 2nd data row for reaction 2
    dyn.system_time = 104
    dyn.save_diagnostic_rxn_data(rxn_index=2, time_step=4,
                                 increment_vector_single_rxn=np.array([-11, -11, 0, 11]),
                                 caption="2nd row")

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False)
    expected_df_2 = pd.DataFrame([[100, -8, -8, 0, 8, 4, "I'm a caption"] , [104, -11, -11, 0, 11, 4, "2nd row"]],
                                 columns = ["START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False, tail=2)   # The full dataset, again
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False, head=2)   # The full dataset, again
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False, head=1)   # Just the first row
    expected_df_2 = pd.DataFrame([[100, -8, -8, 0, 8, 4, "I'm a caption"]],
                                 columns = ["START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False, tail=1)   # Just the last row
    expected_df_2 = pd.DataFrame([[104, -11, -11, 0, 11, 4, "2nd row"]],
                                 columns = ["START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    expected_df_2.index = [1]   # To conform to the original index, which the tail parameter doesn't alter
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False, t=150)   # The row closest in time will be the last row
    expected_df_2 = pd.DataFrame([[150, 104, -11, -11, 0, 11, 4, "2nd row"]],
                                 columns = ["search_value", "START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    expected_df_2.index = [1]   # To conform to the original index, which the t parameter doesn't alter
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = dyn.get_diagnostic_rxn_data(rxn_index=2, print_reaction=False, t=30)   # The row closest in time will be the first row
    expected_df_2 = pd.DataFrame([[30, 100, -8, -8, 0, 8, 4, "I'm a caption"]],
                                 columns = ["search_value", "START_TIME", "Delta A", "Delta B", "Delta C", "Delta X", "time_step", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)
