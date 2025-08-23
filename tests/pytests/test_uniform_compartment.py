import pytest
import numpy as np
from life123 import ChemData, UniformCompartment, Reactions



def test_constructor():
    names = ["A", "B", "C"]
    chem_data = ChemData(names=names)

    with pytest.raises(Exception):
        UniformCompartment(chem_data=chem_data, names=names)

    rxns = Reactions(chem_data=chem_data)

    with pytest.raises(Exception):
        UniformCompartment(reactions=rxns, names=names)

    alt_chem_data = ChemData()
    with pytest.raises(Exception):
        UniformCompartment(reactions=rxns, chem_data=alt_chem_data) # inconsistent


    uc = UniformCompartment(reactions=rxns)
    assert uc.reactions == rxns
    assert uc.chem_data == chem_data
    assert uc.chem_data.get_all_labels() == names

    uc = UniformCompartment(chem_data=chem_data)
    assert uc.chem_data == chem_data
    assert uc.chem_data.get_all_labels() == names
    assert uc.reactions.number_of_reactions() == 0

    uc = UniformCompartment(names=names)
    assert uc.chem_data.get_all_labels() == names
    assert uc.reactions.number_of_reactions() == 0

    uc = UniformCompartment(reactions=rxns, chem_data=chem_data)
    assert uc.reactions == rxns
    assert uc.chem_data == chem_data
    assert uc.chem_data.get_all_labels() == names



def test_set_conc():
    #TODO: test the snapshot argument
    chem_data = ChemData(names=["A", "B", "C"])
    uc = UniformCompartment(chem_data=chem_data)

    with pytest.raises(Exception):
        uc.set_conc(conc=[1, 2, 3, 4])         # Wrong number of entries
    with pytest.raises(Exception):
        uc.set_conc(conc=[1., -2., -3.])       # Negative values
    with pytest.raises(Exception):
        uc.set_conc(conc=(10., 20., -0.01))    # Negative values

    uc.set_conc(conc=[1., 2., 3.])
    assert np.allclose(uc.system, [1., 2., 3.])

    uc.set_conc(conc=(10., 20., 30.))
    assert np.allclose(uc.system, [10., 20., 30.])

    uc.set_conc(conc={"B": 100})
    assert np.allclose(uc.system, [10., 100., 30.])


    uc = UniformCompartment(chem_data=chem_data)
    uc.set_conc(conc={"C": 3})
    assert np.allclose(uc.system, [0., 0., 3.])

    uc.set_conc(conc={"C": 8, "A": 1, "B": 5})
    assert np.allclose(uc.system, [1., 5., 8.])

    with pytest.raises(Exception):
        uc.set_conc(conc={"A": -0.01})



def test_get_system_conc():
    chem_data = ChemData(names=["A", "B", "C"])
    uc = UniformCompartment(chem_data=chem_data)
    uc.set_conc(conc=(10., 20., 30.))

    result = uc.get_system_conc()
    assert np.allclose(result, [10., 20., 30.])



def test_get_chem_conc():
    chem_data = ChemData(names=["A", "B"])
    uc = UniformCompartment(chem_data=chem_data)
    uc.set_conc(conc=(10., 20.))

    assert np.allclose(uc.get_chem_conc("A"), 10.)
    assert np.allclose(uc.get_chem_conc("B"), 20.)
    with pytest.raises(Exception):
        uc.get_chem_conc("Unknown")    # Non-existent chemical



def test_get_conc_dict():
    chem_data = ChemData(names=["A", "B", "C", "D"])
    rxn = UniformCompartment(chem_data=chem_data)
    rxn.set_conc(conc=(100, 200, 300, 400))

    result = rxn.get_conc_dict()
    assert result == {"A": 100, "B": 200, "C": 300, "D": 400}

    result = rxn.get_conc_dict(species=["D", "A"])
    assert result == {"A": 100, "D": 400}

    result = rxn.get_conc_dict(species=("C",))      # Tuple with 1 element
    assert result == {"C": 300}

    with pytest.raises(Exception):
        rxn.get_conc_dict(species="C")                      # Wrong data type

    with pytest.raises(Exception):
        rxn.get_conc_dict(system_data=np.array([1, 2]))     # Wrong number of entries

    result = rxn.get_conc_dict(system_data=np.array([1, 2, 3, 4]))
    assert result == {"A": 1, "B": 2, "C": 3, "D": 4}

    result = rxn.get_conc_dict(species=["B"], system_data=np.array([1, 2, 3, 4]))
    assert result == {"B": 2}





##########################################################################################################

def test_specify_steps():
    uc = UniformCompartment()

    with pytest.raises(Exception):
        # Too few arguments
        uc.specify_steps()
        uc.specify_steps(total_duration=15.1)
        uc.specify_steps(time_step=0.2)
        uc.specify_steps(n_steps=30)
        # Too many arguments
        uc.specify_steps(total_duration=15.1, time_step=0.2, n_steps=30)

    assert uc.specify_steps(time_step=0.5, n_steps=24) == (0.5, 24)
    assert uc.specify_steps(total_duration=12.0, time_step=0.5) == (0.5, 24)
    assert uc.specify_steps(total_duration=12.0, n_steps=24) == (0.5, 24)



def test_single_compartment_react():

    # Test based on experiment "cycles_1"
    chem_data = ChemData(names=["A", "B", "C", "E_high", "E_low"])
    rxns = Reactions(chem_data=chem_data)

    # Reaction A <-> B, mostly in forward direction (favored energetically)
    rxns.add_reaction(reactants="A", products="B",
                      forward_rate=9., reverse_rate=3.)

    # Reaction B <-> C, also favored energetically
    rxns.add_reaction(reactants="B", products="C",
                      forward_rate=8., reverse_rate=4.)

    # Reaction C + E_High <-> A + E_Low, also favored energetically, but kinetically slow
    rxns.add_reaction(reactants=["C" , "E_high"], products=["A", "E_low"],
                      forward_rate=1., reverse_rate=0.2)

    initial_conc = {"A": 100., "B": 0., "C": 0., "E_high": 1000., "E_low": 0.}

    uc = UniformCompartment(reactions=rxns, enable_diagnostics=True)

    uc.set_conc(conc=initial_conc, snapshot=True)


    uc.single_compartment_react(initial_step=0.0005, target_end_time=0.0035, variable_steps=False)

    run1 = uc.get_system_conc()

    assert np.allclose(uc.system_time, 0.0035)
    assert np.allclose(run1, [9.69252541e+01, 3.05696280e+00, 1.77831454e-02, 9.99980686e+02, 1.93144884e-02])
    assert uc.diagnostics.explain_time_advance(return_times=True, silent=True) == \
               ([0.0, 0.0035], [0.0005])


    # Now repeat the process, step-by-step

    uc2 = UniformCompartment(reactions=rxns, enable_diagnostics=True)

    uc2.set_conc(conc=initial_conc, snapshot=True)

    for _ in range(7):
        uc2.single_compartment_react(initial_step=0.0005, n_steps=1, variable_steps=False)

    run2 = uc.get_system_conc()
    assert np.allclose(run2, run1)      # Same result as before
    assert np.allclose(uc2.system_time, 0.0035)

    print(uc2.diagnostics.explain_time_advance(return_times=True))
    assert uc.diagnostics.explain_time_advance(return_times=True, silent=True) == \
               ([0.0, 0.0035], [0.0005])



def test_reaction_step_common_fixed_step_1():
    uc = UniformCompartment(names=["A", "B"])

    uc.set_conc(conc=[10., 50.], snapshot=False)

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    uc.add_reaction(reactants="A", products="B", forward_rate=3., reverse_rate=2.)

    result = uc.reaction_step_common_fixed_step(delta_time=0.1)
    assert np.allclose(result, [ 7. , -7.])     # The increment vector

    # Reset the system concentration to their original values, and repeat
    uc.set_conc(conc=[10., 50.], snapshot=False)
    result = uc.reaction_step_common_fixed_step(delta_time=0.1)
    assert np.allclose(result, [ 7. , -7.])

    uc.set_conc(conc=[10., 50.], snapshot=False)
    result = uc.reaction_step_common_fixed_step(delta_time=0.2)
    assert np.allclose(result, [ 14. , -14.])     # The increment vector has doubled

    uc.set_conc(conc=[10., 50.], snapshot=False)
    result = uc.reaction_step_common_fixed_step(delta_time=0.7142857)
    assert np.allclose(result, [49.999999, -49.999999])     # The increment vector is a hair from making [B] negative

    uc.set_conc(conc=[10., 50.], snapshot=False)
    with pytest.raises(Exception):
        uc.reaction_step_common_fixed_step(delta_time=0.71429)  # A step so large that it would make [B] negative



def test_reaction_step_common_fixed_step_2():

    uc = UniformCompartment(names=["A", "B", "C"])

    uc.set_conc(conc=[10., 50., 20.], snapshot=False)

    # Reaction A + B <-> C , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction4"
    uc.add_reaction(reactants=["A" , "B"], products="C",
                    forward_rate=5., reverse_rate=2.)

    result = uc.reaction_step_common_fixed_step(delta_time=0.002)
    assert np.allclose(result, [-4.92, -4.92, 4.92])

    result = uc.reaction_step_common_fixed_step(delta_time=0.00406504)
    assert np.allclose(result, [-9.9999984, -9.9999984,  9.9999984])   # A hair from making [A] negative

    with pytest.raises(Exception):
        uc.reaction_step_common_fixed_step(delta_time=0.0040651)    # A step so large that it would make [A] negative


    # Now let's consider a different system, with a reaction A <-> B , with 1st-order kinetics in both directions
    uc = UniformCompartment(names=["A", "B", "C"])
    uc.add_reaction(reactants="A", products="B", forward_rate=300., reverse_rate=2.)

    uc.set_conc(conc=[10., 50., 20.], snapshot=False)

    result = uc.reaction_step_common_fixed_step(delta_time=0.002)
    # 10 * 300 * .002 - 50 * 2 * .002 = 5.8
    assert np.allclose(result, [-5.8,  5.8,  0.])   # C isn't affected by this reaction; hence, 0 change


    # Add the reaction we saw earlier, A + B <-> C, and reset the concentrations
    uc.add_reaction(reactants=["A" , "B"], products="C",
                    forward_rate=5., reverse_rate=2.)
    uc.set_conc(conc=[10., 50., 20.], snapshot=False)

    # We now have 2 reaction
    assert uc.number_of_reactions() == 2

    # We saw in earlier runs, with our initial concentrations, that
    # over a delta_time=0.02, one reaction causes a change in [A] of -4.92,
    # and the other a change of -5.8
    # Individually, neither is problematic, given that [A] is initially 10,
    # but combined (-10.72) they would make [A] negative!
    with pytest.raises(Exception):
        uc.reaction_step_common_fixed_step(delta_time=0.002)



def test__reaction_elemental_step():
    uc = UniformCompartment(names=["A", "B"])

    uc.set_conc(conc=[10., 50.], snapshot=False)

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    uc.add_reaction(reactants="A", products="B", forward_rate=3., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [ 7. , -7.])
    assert result[0] == - result[1]         # From the stoichiometry


    uc.clear_reactions()       # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction2"
    uc.add_reaction(reactants="A", products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]      # From the stoichiometry


    uc.clear_reactions()       # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction3"
    uc.add_reaction(reactants=[(2,"A",1)], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test__reaction_elemental_step_1():
    uc = UniformCompartment(names=["A", "B"])

    uc.set_conc(conc=[10., 50.], snapshot=False)

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # Based on experiment "reactions_single_compartment/react_1"
    uc.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [ 7. , -7.])
    assert result[0] == - result[1]         # From the stoichiometry


    uc.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction2"
    uc.add_reaction(reactants=["A"], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [5. , -15.])
    assert -3 * result[0] == result[1]      # From the stoichiometry


    uc.clear_reactions()   # Re-start with a blank slate of reactions
    # Reaction 2A <-> 3B , with 1st-order kinetics in both directions.
    # Based on experiment "1D/reactions/reaction3"
    uc.add_reaction(reactants=[(2,"A",1)], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [10., -15.])
    assert result[0]/2 == - result[1] /3   # From the stoichiometry



def test__reaction_elemental_step_2():
    uc = UniformCompartment(names=["A", "B", "C"])

    uc.set_conc(conc=[10., 50., 20.], snapshot=False)

    # Reaction A <-> B , with 1st-order kinetics in both directions.
    # # Based on experiment "reactions_single_compartment/react_1"
    uc.add_reaction(reactants="A", products="B", forward_rate=3., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.1)
    assert np.allclose(result, [ 7. , -7. , 0.])    # Chemical "C" not participating in this reaction; its delta conc. is 0
    assert result[0] == - result[1]         # From the stoichiometry


    uc.reactions.clear_reactions_data()   # Re-start with a blank slate of reactions
    # Reaction A + B <-> C , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction4"
    uc.add_reaction(reactants=["A" , "B"], products="C",
                     forward_rate=5., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.002)
    assert np.allclose(result, [-4.92, -4.92, 4.92])
    assert result[0] == result[1]           # From the stoichiometry
    assert result[1] == - result[2]         # From the stoichiometry



def test__reaction_elemental_step_3():
    uc = UniformCompartment(names=["A", "C", "D"])

    uc.set_conc(conc=[4., 7., 2.], snapshot=False)

    # Reaction A <-> 2C + D , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction5"
    uc.add_reaction(reactants=[("A")], products=[(2, "C", 1) , ("D")],
                    forward_rate=5., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.05)
    assert np.allclose(result, [0.4 , -0.8 , -0.4])
    assert result[0] == - result[1] /2    # From the stoichiometry
    assert result[0] == - result[2]       # From the stoichiometry



def test__reaction_elemental_step_4():
    uc = UniformCompartment(names=["A", "B", "C", "D"])

    uc.set_conc(conc=[4., 7., 5., 2.], snapshot=False)

    # Reaction 2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction6"
    uc.add_reaction(reactants=[(2,"A",1) , (5,"B",1)], products=[(4,"C",1) , (3,"D",1)],
                     forward_rate=5., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.001)
    assert np.allclose(result, [-0.24 , -0.6 , 0.48, 0.36])
    assert  np.allclose(result[0] /2 , result[1] /5)    # From the stoichiometry
    assert  np.allclose(result[1] /5 , -result[2] /4)   # From the stoichiometry
    assert  np.allclose(result[2] /4 , result[3] /3)    # From the stoichiometry



def test__reaction_elemental_step_5():
    uc = UniformCompartment(names=["A", "B"])

    uc.set_conc(conc=[3., 5.], snapshot=False)

    # Reaction  2A <-> B , with 2nd-order kinetics in forward reaction, and 1st-order in reverse.
    # Based on experiment "1D/reactions/reaction7"
    uc.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)

    result = uc._reaction_elemental_step(delta_time=0.02)
    assert np.allclose(result, [-1.4 , 0.7])
    assert np.allclose(result[0] /2 , -result[1])      # From the stoichiometry



def test__reaction_elemental_step_6():
    uc = UniformCompartment(names=["A", "B", "C", "D", "E"])

    uc.set_conc(conc=[3., 5., 1., 0.4, 0.1], snapshot=False)

    # Coupled reactions A + B <-> C  and  C + D <-> E , with 1st-order kinetics for each species.
    # Based on experiment "1D/reactions/reaction8"
    uc.add_reaction(reactants=["A", "B"], products=["C"], forward_rate=5., reverse_rate=2.)
    uc.add_reaction(reactants=["C", "D"], products=["E"], forward_rate=8., reverse_rate=4.)
    assert uc.number_of_reactions() == 2

    result = uc._reaction_elemental_step(delta_time=0.02)
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

    uc = UniformCompartment(chem_data=chem_data, preset=None)

    # Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
    uc.add_reaction(reactants=[(2, "S", 1)], products="U",
                           forward_rate=8., reverse_rate=2.)

    # Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
    uc.add_reaction(reactants="S", products="X",
                           forward_rate=6., reverse_rate=3.)
     
    uc.set_conc(conc={"U": 50., "X": 100., "S": 0.})

    uc.adaptive_steps.set_thresholds(norm="norm_A", low=0.25, high=0.64, abort=1.44)
    uc.adaptive_steps.set_step_factors(abort=0.5, downshift=0.5, upshift=2.0)

    uc.single_compartment_react(initial_step=0.01, target_end_time=0.2,
                                      variable_steps=True)

    df = uc.get_history()
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



def test_single_compartment_correct_neg_conc():
    # Based on "Run 3" of experiment "negative_concentrations_1

    # Initialize the system
    chem_data = ChemData(names=["U", "X", "S"])

    uc = UniformCompartment(chem_data=chem_data)

    # Reaction 2 S <-> U , with 1st-order kinetics for all species (mostly forward)
    uc.add_reaction(reactants=[(2, "S", 1)], products="U",
                    forward_rate=8., reverse_rate=2.)

    # Reaction S <-> X , with 1st-order kinetics for all species (mostly forward)
    uc.add_reaction(reactants="S", products="X",
                    forward_rate=6., reverse_rate=3.)
    
    uc.set_conc(conc={"U": 50., "X": 100., "S": 0.})

    uc.enable_diagnostics()       # To save diagnostic information about the call to single_compartment_react()

    uc.single_compartment_react(initial_step=0.25, n_steps=1, variable_steps=False)

    assert np.allclose(uc.system_time, 0.25)
    assert np.allclose(uc.system, [ 25.,  25., 125.])

    with pytest.raises(Exception):
        # This step would make [S] negative
        uc.single_compartment_react(initial_step=0.25, n_steps=1, variable_steps=False)

    # Nothing has changed, since that last step wasn't actually taken
    assert np.allclose(uc.system_time, 0.25)
    assert np.allclose(uc.system, [ 25.,  25., 125.])

    # A smaller step saves the day!
    uc.single_compartment_react(initial_step=0.03, n_steps=1, variable_steps=False)
    assert np.allclose(uc.system_time, 0.28)    # 0.25 + 0.03
    assert np.allclose(uc.system, [53.5,  45.25, 47.75])





###########################  LOWER-LEVEL METHODS  ###########################


def test_is_in_equilibrium():
    chem_data = ChemData(names=["A", "B", "C", "D", "E", "F"])
    uc = UniformCompartment(chem_data=chem_data)

    # Reaction 0 : A <-> B
    uc.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    c = {'A': 23.9931640625, 'B': 36.0068359375}
    assert uc.is_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)
    assert uc.is_in_equilibrium(conc=c, explain=False, tolerance=1)      # Testing ALL reactions

    # Reaction 1 : A <-> F
    uc.add_reaction(reactants=["A"], products=["F"], forward_rate=20, reverse_rate=2.)
    c = {'A': 3, 'F': 32.999}
    assert uc.is_in_equilibrium(rxn_index=1, conc=c, explain=False, tolerance=10)   # The deviation is just below the 10% tolerance

    with pytest.raises(Exception):
        assert uc.is_in_equilibrium(conc=c, explain=False, tolerance=10) # Testing ALL reactions, but neglected to provide [B]

    c = {'A': 3, 'B': 4.6, 'F': 33.001}   # We're now including all the concentrations we need for both reactions
    assert uc.is_in_equilibrium(conc=c, explain=False, tolerance=10) \
           == {False: [1]}    # The deviation for reaction 1 is just above the 10% tolerance

    assert uc.is_in_equilibrium(conc=c, explain=False, tolerance=2.2) \
           == {False: [0, 1]}      # With a tolerance so low, not reaction meets the criteria

    assert uc.is_in_equilibrium(conc=c, explain=False, tolerance=2.3) \
           == {False: [1]}                # Reaction 0 barely passes with this tolerance (it deviates by 2.222 %)

    with pytest.raises(Exception):
        assert uc.is_in_equilibrium(explain=False, tolerance=2.3)  # We're failing to provide concentrations

    uc.set_conc(conc=[3, 4.6, 0, 0, 0, 33.001])    # The concentrations are in the same order as the declared chemicals
    assert uc.is_in_equilibrium(conc=c, explain=False, tolerance=2.3) \
           == {False: [1]}        # Now using the System concentrations



def test_reaction_in_equilibrium():
    chem_data = ChemData(names=["A", "B", "C", "D", "E", "F"])
    uc = UniformCompartment(chem_data=chem_data)

    # Reaction 0 : A <-> B
    uc.add_reaction(reactants=["A"], products=["B"], forward_rate=3., reverse_rate=2.)
    c = {'A': 23.9931640625, 'B': 36.0068359375}
    assert uc.reaction_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)

    # Reaction 1 : A <-> F
    uc.add_reaction(reactants=["A"], products=["F"], forward_rate=20, reverse_rate=2.)
    c = {'A': 3, 'F': 32.999}
    assert uc.reaction_in_equilibrium(rxn_index = 1, conc=c, explain=False, tolerance=10)   # Just below the 10% tolerance
    c = {'A': 3, 'F': 33.001}
    assert not uc.reaction_in_equilibrium(rxn_index = 1, conc=c, explain=False, tolerance=10)  # Just above the 10% tolerance

    # Reaction 2 : A <-> 3B
    uc.add_reaction(reactants=["A"], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)
    c = {'A': 14.54545455, 'B': 36.36363636}
    assert uc.reaction_in_equilibrium(rxn_index = 2, conc=c, explain=False, tolerance=1)

    # Reaction 3:  2A <-> 3B
    uc.add_reaction(reactants=[(2,"A",1)], products=[(3,"B",1)], forward_rate=5., reverse_rate=2.)
    c = {'A': 16.25, 'B': 40.625}
    assert uc.reaction_in_equilibrium(rxn_index = 3, conc=c, explain=False, tolerance=1)

    # Reaction 4:  A + B <-> C , with 1st-order kinetics for each species
    uc.add_reaction(reactants=["A" , "B"], products=[("C")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 0.29487741, 'B': 40.29487741, 'C': 29.70512259}
    assert uc.reaction_in_equilibrium(rxn_index = 4, conc=c, explain=False, tolerance=1)

    # Reaction 5:  A <-> 2C + D , with 1st-order kinetics for each species
    uc.add_reaction(reactants=["A"], products=[(2, "C", 1) , ("D")],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 4.31058733, 'C': 6.37882534, 'D': 1.68941267}
    assert uc.reaction_in_equilibrium(rxn_index = 5, conc=c, explain=False, tolerance=1)

    # Reaction 6:  2A + 5B <-> 4C + 3D , with 1st-order kinetics for each species
    uc.add_reaction(reactants=[(2,"A",1) , (5,"B",1)], products=[(4,"C",1) , (3,"D",1)],
                     forward_rate=5., reverse_rate=2.)
    c = {'A': 2.80284552, 'B': 4.00711381, 'C': 7.39430896, 'D': 3.79573172}
    assert uc.reaction_in_equilibrium(rxn_index = 6, conc=c, explain=False, tolerance=1)


    uc.clear_reactions()   # This will reset the reaction count to 0

    # Reaction 0:  2A <-> B , with 1st-order kinetics in both directions
    uc.add_reaction(reactants=[(2, "A", 1)], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 2.16928427, 'B': 5.41535786}
    assert uc.reaction_in_equilibrium(rxn_index = 0, conc=c, explain=False, tolerance=1)

    # Reaction 1:  2A <-> B , NOW WITH 2nd-order kinetics in the forward direction
    uc.add_reaction(reactants=[(2, "A", 2)], products=["B"], forward_rate=5., reverse_rate=2.)
    c = {'A': 1.51554944, 'B': 5.74222528}
    assert uc.reaction_in_equilibrium(rxn_index = 1, conc=c, tolerance=1, explain=False)

    # Reaction 2:  A + B <-> C + D
    uc.add_reaction(reactants=["A", "B"], products=["C", "D"], forward_rate=5., reverse_rate=2.)

    # with zero concentrations of a reactant, it won't be an equilibrium...
    c = {'A': 15.2, 'B': 0, 'C': 21.3, 'D': 4.1}
    assert not uc.reaction_in_equilibrium(rxn_index = 2, conc=c, tolerance=10, explain=False)

    # likewise with zero concentrations of a reaction product, it won't be an equilibrium...
    c = {'A': 15.2, 'B': 21.3, 'C': 0 , 'D': 4.1}
    assert not uc.reaction_in_equilibrium(rxn_index = 2, conc=c, tolerance=10, explain=False)

    # but with zero concentrations of both a reactant and a reaction product, the reaction is stuck in place -
    # so, in equilibrium
    c = {'A': 15.2, 'B': 0, 'C': 0 , 'D': 4.1}
    assert uc.reaction_in_equilibrium(rxn_index = 2, conc=c, tolerance=1, explain=False)



def test_validate_increment():
    uc = UniformCompartment(names=["A", "B", "C"])

    uc.add_reaction(reactants="A", products="B")

    uc.validate_increment(delta_conc=50., baseline_conc=10., rxn_index=0, species_index=2, delta_time=0.02)
    uc.validate_increment(delta_conc=-9.99, baseline_conc=10., rxn_index=0, species_index=2, delta_time=0.02)

    with pytest.raises(Exception):      # Would lead to a negative concentration
        uc.validate_increment(delta_conc=-10.001, baseline_conc=10., rxn_index=0, species_index=2, delta_time=0.97)

    # TODO: test that diagnostic data, if enabled, gets saved as needed



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

    rxn = UniformCompartment(chem_data=chem)
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

    rxn = UniformCompartment(chem_data=chem_data)

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

    rxn = UniformCompartment(chem_data=chem_data)

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


    rxn = UniformCompartment(chem_data=chem_data)

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
