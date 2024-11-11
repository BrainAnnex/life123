import pytest
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from life123 import ChemData, MovieTabular
from life123.diagnostics import Diagnostics



def test_explain_time_advance():
    diag = Diagnostics(ChemData())   # Argument isn't actually used, but it's required

    # Start out with uniform steps
    diag.diagnostic_conc_data.store(par=20.,
                                   data_snapshot={"time_step": 100.})

    assert diag.explain_time_advance(return_times=True, silent=True) is None


    diag.diagnostic_conc_data.store(par=30.,
                                    data_snapshot={"time_step": 100.})
    result, _ = diag.explain_time_advance(return_times=True, silent=True)    # TODO: also test the returned step sizes
    assert np.allclose(result, [20., 30.])

    diag.diagnostic_conc_data.store(par=40.,
                                    data_snapshot={"time_step": 100.})
    result, _ = diag.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40.])

    # Switching to smaller step
    diag.diagnostic_conc_data.store(par=45.,
                                    data_snapshot={"time_step": 100.})
    result, _ = diag.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 45.])

    diag.diagnostic_conc_data.store(par=50.,
                                    data_snapshot={"time_step": 100.})
    result, _ = diag.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50.])

    # Switching to larger step
    diag.diagnostic_conc_data.store(par=70.,
                                    data_snapshot={"time_step": 100.})
    result, _ = diag.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70.])

    # Yet larger
    diag.diagnostic_conc_data.store(par=95.,
                                    data_snapshot={"time_step": 100.})
    result, _ = diag.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95.])

    # Smaller again
    diag.diagnostic_conc_data.store(par=96.,
                                    data_snapshot={"time_step": 100.})
    result, _ = diag.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95., 96.])

    diag.diagnostic_conc_data.store(par=97.,
                                    data_snapshot={"time_step": 100.})
    result, _ = diag.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95., 97.])

    diag.diagnostic_conc_data.store(par=98.,
                                    data_snapshot={"time_step": 100.})
    result, _ = diag.explain_time_advance(return_times=True, silent=True)
    assert np.allclose(result, [20., 40., 50., 70., 95., 98.])

    #print(rxn.diagnostic_data_baselines.get())
    #print(result)



def test__delta_names():
    chem_data = ChemData(names=["A", "B", "X"])
    diag = Diagnostics(chem_data)

    assert diag._delta_names() == ["Delta A", "Delta B", "Delta X"]



def test__delta_conc_dict():
    chem_data = ChemData(names=["A", "B", "X"])
    diag = Diagnostics(chem_data)

    assert diag._delta_conc_dict(np.array([10, 20, 30])) == \
           {"Delta A": 10, "Delta B": 20, "Delta X": 30}

    with pytest.raises(Exception):
        diag._delta_conc_dict(np.array([10, 20, 30, 40]))    # One element too many



def test_save_diagnostic_rxn_data():
    chem_data = ChemData(names=["A", "B", "C", "X"])
    # Add 3 reactions
    chem_data.add_reaction(reactants="A", products="B", forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants="A", products="X", forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A", "B"], products="X", forward_rate=5., reverse_rate=2.)

    diag = Diagnostics(chem_data)
    assert len(diag.diagnostic_rxn_data) == 0

    with pytest.raises(Exception):
        diag.save_rxn_data(rxn_index=0, system_time=100, time_step=4,
                           increment_vector_single_rxn=np.array([2, -2]))     # Wrong size of Numpy array

    # Add data for reaction index 0
    diag.save_rxn_data(rxn_index=0, system_time=100, time_step=4,
                       increment_vector_single_rxn=np.array([2, -2, 0, 0]))

    assert len(diag.diagnostic_rxn_data) == 1

    diagnostic_data_rxn_0 = diag.diagnostic_rxn_data[0]

    assert (type(diagnostic_data_rxn_0)) == MovieTabular

    df_0 = diagnostic_data_rxn_0.get_dataframe()

    expected_df_0 = pd.DataFrame([[100, 4, False, 2, -2, ""]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B",  "caption"])
    assert_frame_equal(df_0, expected_df_0, check_dtype=False)


    # Add data for reaction index 1
    diag.save_rxn_data(system_time=100, rxn_index=1, time_step=4,
                       increment_vector_single_rxn=np.array([7, 0, 0, -7]))

    assert len(diag.diagnostic_rxn_data) == 2       # 2 reactions added so far
    df_0 = diagnostic_data_rxn_0.get_dataframe()
    diagnostic_data_rxn_1 = diag.diagnostic_rxn_data[1]
    df_1 = diagnostic_data_rxn_1.get_dataframe()

    assert_frame_equal(df_0, expected_df_0, check_dtype=False)  # Nothing was done to df_0 by processing reaction index 1

    expected_df_1 = pd.DataFrame([[100, 4, False, 7, -7, ""]],
                                 columns = ["START_TIME",  "time_step", "aborted", "Delta A", "Delta X","caption"])
    assert_frame_equal(df_1, expected_df_1, check_dtype=False)


    # Add data for reaction index 2
    diag.save_rxn_data(system_time=100, rxn_index=2, time_step=4,
                       increment_vector_single_rxn=np.array([-8, -8, 0, 8]),
                       caption="I'm a caption")

    assert len(diag.diagnostic_rxn_data) == 3       # 3 reactions added so far
    diagnostic_data_rxn_2 = diag.diagnostic_rxn_data[2]
    df_2 = diagnostic_data_rxn_2.get_dataframe()

    assert_frame_equal(df_0, expected_df_0, check_dtype=False)  # Nothing was done to df_0 by processing reaction index 2
    assert_frame_equal(df_1, expected_df_1, check_dtype=False)  # Nothing was done to df_1 by processing reaction index 2

    expected_df_2 = pd.DataFrame([[100, 4, False, -8, -8, 8, "I'm a caption"]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B", "Delta X", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)


    # Add more data for reaction index 0
    diag.save_rxn_data(rxn_index=0, system_time=104, time_step=6,
                       increment_vector_single_rxn=np.array([-1, 1, 0, 0]), caption="my comment")

    assert len(diag.diagnostic_rxn_data) == 3       # Still 3 reactions

    assert_frame_equal(df_1, expected_df_1, check_dtype=False)  # Nothing was done to df_1 by processing reaction index 0
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)  # Nothing was done to df_2 by processing reaction index 0


    df_0 = diagnostic_data_rxn_0.get_dataframe()
    expected_df_0.loc[len(expected_df_0)] = [104, 6, False, -1, 1, "my comment"]    # To add a row to the existing df
    assert_frame_equal(df_0, expected_df_0, check_dtype=False)


    # Add a 3rd entry for reaction index 0, this time with rate information
    diag.save_rxn_data(rxn_index=0, system_time=110, time_step=12,
                       increment_vector_single_rxn=np.array([-4, 4, 0, 0]),
                       caption="start recording rate", rate=3)


    df_0 = diagnostic_data_rxn_0.get_dataframe()
    expected_df_0 = pd.DataFrame([[100,  4,  False, 2, -2, "", np.nan],
                                  [104,  6, False, -1,  1, "my comment", np.nan],
                                  [110, 12, False, -4,  4, "start recording rate", 3.0]
                                  ],
                                  columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B", "caption", "rate"])
                                  # Notice "retroactively" addind NaN's to the earlier rows that didn't save a rate value
    assert_frame_equal(df_0, expected_df_0, check_dtype=False)

    #TODO: more test adding multiple entries for any reaction



def test_save_diagnostic_aborted_rxns():
    chem_data = ChemData(names=["A", "B", "C", "X"])
    # Add 3 reactions
    chem_data.add_reaction(reactants="A", products="B", forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants="A", products="X", forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A", "B"], products="X", forward_rate=5., reverse_rate=2.)

    diag = Diagnostics(chem_data)

    diag.save_diagnostic_aborted_rxns(system_time=100, time_step=40, caption="aborted step")

    assert len(diag.diagnostic_rxn_data) == 3       # All 3 reactions got set

    df_0 = diag.get_rxn_data(rxn_index=0, print_reaction=False)
    expected_df_0 = pd.DataFrame([[100, 40, True, np.nan, np.nan, "aborted step"]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B",  "caption"])
    assert_frame_equal(df_0, expected_df_0, check_dtype=False)

    df_1 = diag.get_rxn_data(rxn_index=1, print_reaction=False)
    expected_df_1 = pd.DataFrame([[100, 40, True, np.nan, np.nan, "aborted step"]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta X",  "caption"])
    assert_frame_equal(df_1, expected_df_1, check_dtype=False)

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False)
    expected_df_2 = pd.DataFrame([[100, 40, True, np.nan, np.nan, np.nan, "aborted step"]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B",  "Delta X",  "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)



def test_get_diagnostic_rxn_data():
    chem_data = ChemData(names=["A", "B", "C", "X"])
    # Add 3 reactions
    chem_data.add_reaction(reactants="A", products="B", forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A"], products=["X"], forward_rate=5., reverse_rate=2.)
    chem_data.add_reaction(reactants=["A", "B"], products=["X"], forward_rate=5., reverse_rate=2.)

    diag = Diagnostics(chem_data)


    # Add data for reaction index 0
    diag.save_rxn_data(system_time=100, rxn_index=0, time_step=4,
                       increment_vector_single_rxn=np.array([2, -2, 0, 0]))

    df_0 = diag.get_rxn_data(rxn_index=0, print_reaction=False)

    expected_df_0 = pd.DataFrame([[100, 4, False, 2, -2, ""]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B", "caption"])
    assert_frame_equal(df_0, expected_df_0, check_dtype=False)

    df_1 = diag.get_rxn_data(rxn_index=1, print_reaction=False)
    assert df_1 is None
    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False)
    assert df_2 is None


    # Add data for reaction index 1
    diag.save_rxn_data(system_time=100, rxn_index=1, time_step=4,
                       increment_vector_single_rxn=np.array([7, 0, 0, -7]))

    df_1 = diag.get_rxn_data(rxn_index=1, print_reaction=False)

    expected_df_1 = pd.DataFrame([[100, 4, False, 7, -7, ""]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta X", "caption"])
    assert_frame_equal(df_1, expected_df_1, check_dtype=False)

    df_0 = diag.get_rxn_data(rxn_index=0, print_reaction=False)
    assert_frame_equal(df_0, expected_df_0, check_dtype=False)      # No change made to df_0 from the last step

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False)
    assert df_2 is None


    # Add data for reaction index 2
    diag.save_rxn_data(system_time=100, rxn_index=2, time_step=4,
                       increment_vector_single_rxn=np.array([-8, -8, 0, 8]),
                       caption="I'm a caption")

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False, tail=1) # With just one row, tail=1 won't make a difference

    expected_df_2 = pd.DataFrame([[100, 4, False, -8, -8,  8,  "I'm a caption"]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B", "Delta X", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False, t=50) # With just one row, the time selector won't matter

    expected_df_2 = pd.DataFrame([[50, 100, 4, False, -8, -8, 8, "I'm a caption"]],
                                 columns = ["search_value", "START_TIME", "time_step", "aborted", "Delta A", "Delta B", "Delta X", "caption"])

    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    # Add a 2nd data row for reaction 2
    diag.save_rxn_data(system_time=104, rxn_index=2, time_step=4,
                       increment_vector_single_rxn=np.array([-11, -11, 0, 11]),
                       caption="2nd row")

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False)
    expected_df_2 = pd.DataFrame([[100, 4, False, -8, -8, 8, "I'm a caption"] , [104, 4, False, -11, -11, 11, "2nd row"]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B", "Delta X", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False, tail=2)   # The full dataset, again
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False, head=2)   # The full dataset, again
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False, head=1)   # Just the first row
    expected_df_2 = pd.DataFrame([[100, 4, False, -8, -8, 8, "I'm a caption"]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B", "Delta X", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False, tail=1)   # Just the last row
    expected_df_2 = pd.DataFrame([[104, 4, False, -11, -11, 11, "2nd row"]],
                                 columns = ["START_TIME", "time_step", "aborted", "Delta A", "Delta B", "Delta X", "caption"])
    expected_df_2.index = [1]   # To conform to the original index, which the tail parameter doesn't alter
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False, t=150)   # The row closest in time will be the last row
    expected_df_2 = pd.DataFrame([[150, 104, 4, False, -11, -11, 11, "2nd row"]],
                                 columns = ["search_value", "START_TIME", "time_step", "aborted", "Delta A", "Delta B", "Delta X", "caption"])
    expected_df_2.index = [1]   # To conform to the original index, which the t parameter doesn't alter
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)

    df_2 = diag.get_rxn_data(rxn_index=2, print_reaction=False, t=30)   # The row closest in time will be the first row
    expected_df_2 = pd.DataFrame([[30, 100, 4, False, -8, -8, 8, "I'm a caption"]],
                                 columns = ["search_value", "START_TIME", "time_step", "aborted", "Delta A", "Delta B", "Delta X", "caption"])
    assert_frame_equal(df_2, expected_df_2, check_dtype=False)



def test_stoichiometry_checker():
    chem = ChemData(names=["A", "B", "C", "D"])
    diag = Diagnostics(chem)

    chem.add_reaction(reactants="A", products="B")          # Reaction 0:   A <--> B

    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 50, 0, 0]), conc_arr_after=np.array([10, 40, 0, 0]))
    assert not diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 50, 0, 0]), conc_arr_after=np.array([10, 39.9, 0, 0]))
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 0, 0, 0]), conc_arr_after=np.array([90, 10, 0, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=[(2, "A")], products=["B"])         # Reaction 0:   2A <--> B

    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 0, 0, 0]), conc_arr_after=np.array([80, 10, 0, 0]))
    assert not diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 0, 0, 0]), conc_arr_after=np.array([80, 10.1, 0, 0]))
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 50, 0, 0]), conc_arr_after=np.array([10, 45, 0, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=["A", "B"], products=["C"])         # Reaction 0:   A + B <--> C

    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 40, 10, 0]))
    assert not diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 40, 9.9, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=["A", (3, "B")], products=["C"])     # Reaction 0:   A + 3B <--> C

    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 10, 0]))
    assert not diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 9.9, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=["A", (3, "B")], products=[(4, "C")])                   # Reaction 0:   A + 3B <--> 4C

    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 40, 0]))
    assert not diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 50, 0, 0]), conc_arr_after=np.array([90, 20, 39.9, 0]))


    chem.clear_reactions_data()
    chem.add_reaction(reactants=[(2, "A"), (3, "B")], products=[(4, "C"), (5, "D")])     # Reaction 0:   2A + 3B <--> 4C + 5D
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([0, 0, 0, 0]), conc_arr_after=np.array([0, 0, 0, 0]))
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 100, 100, 100]), conc_arr_after=np.array([120, 130, 60, 50]))
    assert not diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 100, 100, 100]), conc_arr_after=np.array([120.1, 130, 60, 50]))
    assert diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 100, 100, 100]), conc_arr_after=np.array([80, 70, 140, 150]))
    assert not diag.stoichiometry_checker(rxn_index=0, conc_arr_before=np.array([100, 100, 100, 100.1]), conc_arr_after=np.array([80, 70, 140, 150]))
