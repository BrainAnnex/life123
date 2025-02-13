import pytest
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from life123 import CollectionTabular, CollectionArray, Collection



###############  For class CollectionTabular  ###############

def test_store():
    m = CollectionTabular()

    d = {"A": 1, "B": 2, "C": 3}
    d_original = d
    m.store(par=10, data_snapshot=d, caption="first entry")  # Add a snapshot
    assert d == d_original
    assert len(m) == 1
    assert str(m) == "`CollectionTabular` object with 1 snapshot(s) parametrized by `SYSTEM TIME`"
    row = list(m.collection_df.iloc[0])                 # By row index
    assert row == [10, 1, 2, 3, 'first entry']

    m.store(par=20, data_snapshot={"A": 10, "B": 20, "C": 30}, caption="second entry")  # Add a snapshot
    assert len(m) == 2
    assert str(m) == "`CollectionTabular` object with 2 snapshot(s) parametrized by `SYSTEM TIME`"
    row = list(m.collection_df.iloc[0])                 # By row index
    assert row == [10, 1, 2, 3, 'first entry']
    row = list(m.collection_df.iloc[1])                 # By row index
    assert row == [20, 10, 20, 30, 'second entry']

    m.store(par=30, data_snapshot={"A": -1, "B": -2, "C": -3})      # Add a snapshot
    assert len(m) == 3
    assert str(m) == "`CollectionTabular` object with 3 snapshot(s) parametrized by `SYSTEM TIME`"
    row = list(m.collection_df.iloc[0])                 # By row index
    assert row == [10, 1, 2, 3, 'first entry']
    row = list(m.collection_df.iloc[1])                 # By row index
    assert row == [20, 10, 20, 30, 'second entry']
    row = list(m.collection_df.iloc[2])                 # By row index
    assert row == [30, -1, -2, -3, '']

    m.store(par=40, data_snapshot={"A": 111, "B": 222}, caption="notice that C is missing")  # Add a snapshot
    assert len(m) == 4
    assert str(m) == "`CollectionTabular` object with 4 snapshot(s) parametrized by `SYSTEM TIME`"
    df = m.collection_df
    data_values = [{"SYSTEM TIME": 10, "A": 1,   "B": 2,  "C": 3,  "caption": "first entry"},
                   {"SYSTEM TIME": 20, "A": 10,  "B": 20, "C": 30, "caption": "second entry"},
                   {"SYSTEM TIME": 30, "A": -1,  "B": -2, "C": -3, "caption": ""},
                   {"SYSTEM TIME": 40, "A": 111, "B": 222,         "caption": "notice that C is missing"}
                  ]
    expected = pd.DataFrame(data_values)
    assert df.equals(expected)

    m.store(par=50, data_snapshot={"A": 8, "B": 88, "C": 888, "D": 1}, caption="notice the newly-appeared D")  # Add a snapshot
    assert len(m) == 5
    assert str(m) == "`CollectionTabular` object with 5 snapshot(s) parametrized by `SYSTEM TIME`"
    df = m.collection_df

    data_values.append({"SYSTEM TIME": 50, "A": 8, "B": 88, "C": 888, "D": 1, "caption": "notice the newly-appeared D"})
    expected = pd.DataFrame(data_values)
    assert_frame_equal(df, expected, check_dtype=False)    # To allow for slight discrepancies in floating-point
                                                           # (since int's get converted to floats in columns with Nan's)
    """
       SYSTEM TIME    A    B      C                      caption    D
    0           10    1    2    3.0                  first entry  NaN
    1           20   10   20   30.0                 second entry  NaN
    2           30   -1   -2   -3.0                               NaN
    3           40  111  222    NaN     notice that C is missing  NaN
    4           50    8   88  888.0  notice the newly-appeared D  1.0
    """



def test_get_dataframe():
    m = CollectionTabular()

    # Same data as used in test_CollectionTabular(), except that the `SYSTEM TIME` parameter now has floats values
    m.store(par=10,   data_snapshot={"A": 1, "B": 2, "C": 3}, caption="first entry")
    m.store(par=12.4, data_snapshot={"A": 10, "B": 20, "C": 30}, caption="second entry")
    m.store(par=33.1, data_snapshot={"A": -1, "B": -2, "C": -3})
    m.store(par=40,   data_snapshot={"A": 111, "B": 222}, caption="notice that C is missing")
    m.store(par=50.5, data_snapshot={"A": 8, "B": 88, "C": 888, "D": 1}, caption="notice the newly-appeared D")

    """
       SYSTEM TIME      A    B      C                      caption    D
    0           10      1    2    3.0                  first entry  NaN
    1           12.4   10   20   30.0                 second entry  NaN
    2           33.1   -1   -2   -3.0                               NaN
    3           40    111  222    NaN     notice that C is missing  NaN
    4           50.5    8   88  888.0  notice the newly-appeared D  1.0
    """
    # Check the extraction of the last row
    df_last_row = m.get_dataframe(tail=1)
    #print("\n", df_last_row)

    data_values = [{"SYSTEM TIME": 50.5, "A": 8, "B": 88, "C": 888, "caption": "notice the newly-appeared D", "D": 1}]
    expected_df = pd.DataFrame(data_values, index=[4])
    #print("\n", expected_df)

    assert_frame_equal(df_last_row, expected_df, check_dtype=False) # To allow for slight discrepancies in floating-point
                                                                    # (since int's get converted to floats in columns with Nan's)

    # Check the extraction of the last 2 rows
    df_last_2_rows = m.get_dataframe(tail=2)

    data_values = [ {"SYSTEM TIME": 40, "A": 111, "B": 222, "C": np.nan, "caption": "notice that C is missing"},
                    {"SYSTEM TIME": 50.5, "A": 8, "B": 88,  "C": 888,    "caption": "notice the newly-appeared D", "D": 1}]
    expected_df = pd.DataFrame(data_values, index=[3, 4])

    assert_frame_equal(df_last_2_rows, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point


    # Check the extraction of a row by value SEARCH (using a value for "SYSTEM TIME" a tad smaller than what is in the dataframe)
    df_extracted_row = m.get_dataframe(search_col="SYSTEM TIME", search_val=33.099)

    data_values = [{"search_value": 33.099, "SYSTEM TIME": 33.1, "A": -1, "B": -2, "C": -3.0, "caption": "", "D": np.nan}]
    expected_df = pd.DataFrame(data_values, index=[2])
    #print("\n", expected_df)
    assert_frame_equal(df_extracted_row, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point


    # Check the extraction of a row by value (this time using a slightly larger value for "SYSTEM TIME" than what is in the dataframe)
    df_extracted_row = m.get_dataframe(search_col="SYSTEM TIME", search_val=33.1234)
    #print("\n", df_extracted_row)

    data_values = [{"search_value": 33.1234, "SYSTEM TIME": 33.1, "A": -1, "B": -2, "C": -3.0, "caption": "", "D": np.nan}]
    expected_df = pd.DataFrame(data_values, index=[2])
    assert_frame_equal(df_extracted_row, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point


    # Check the extraction of a group of row by value-range filtering
    df_filtered = m.get_dataframe(search_col="SYSTEM TIME", val_start=35)        # This corresponds to the last 2 rows

    data_values = [ {"SYSTEM TIME": 40, "A": 111, "B": 222, "C": np.nan, "caption": "notice that C is missing"},
                    {"SYSTEM TIME": 50.5, "A": 8, "B": 88,  "C": 888,    "caption": "notice the newly-appeared D", "D": 1}]
    expected_df = pd.DataFrame(data_values, index=[3, 4])

    assert_frame_equal(df_filtered, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point



def test_get_2():
    m = CollectionTabular()

    m.store(par=10, data_snapshot={"A": 1, "B": 2, }, caption="first entry")  # Add a snapshot
    row = list(m.collection_df.iloc[0])                 # By row index, reaching into the internal data structure
    assert row == [10, 1, 2, 'first entry']

    df_returned_copy = m.get_dataframe(return_copy=True)
    expected_df = pd.DataFrame([{"SYSTEM TIME": 10, "A": 1, "B": 2, "caption": "first entry"}])
    assert_frame_equal(df_returned_copy, expected_df)

    # Now we mess around with the returned value
    df_returned_copy.loc[0, "A"] = 999

    row = list(m.collection_df.iloc[0])
    assert row == [10, 1, 2, 'first entry']     # The internal data structure is untouched


    df_returned_view = m.get_dataframe(return_copy=False)
    expected_df = pd.DataFrame([{"SYSTEM TIME": 10, "A": 1, "B": 2, "caption": "first entry"}])
    assert_frame_equal(df_returned_view, expected_df)

    # Now we mess around with the returned value
    df_returned_view.loc[0, "A"] = 999

    row = list(m.collection_df.iloc[0])
    assert row == [10, 999, 2, 'first entry']     # The internal data structure got messed up as well




def test_set_caption_last_snapshot():
    m = CollectionTabular()
    m.store(par=100, data_snapshot={"A": 1, "B": 2, "C": 3}, caption="first entry")
    m.store(par=200, data_snapshot={"A": 10, "B": 20, "C": 30})
    m.set_caption_last_snapshot("End of experiment")
    #print(m.collection_df)
    last_row = list(m.collection_df.loc[1])
    assert last_row == [200, 10, 20, 30, 'End of experiment']



def test_set_field_last_snapshot():
    m = CollectionTabular()
    m.store(par=100, data_snapshot={"A": 1, "B": 2}, caption="first entry")     # Add a 1st row

    m.set_field_last_snapshot("B", 22)
    last_row = list(m.collection_df.loc[0])
    assert last_row == [100, 1, 22, 'first entry']

    m.set_field_last_snapshot("X", 99)
    last_row = list(m.collection_df.loc[0])
    assert list(m.collection_df.columns) == ["SYSTEM TIME", "A", "B", "caption", "X"]   # New column present
    assert last_row == [100, 1, 22, 'first entry', 99]

    m.store(par=200, data_snapshot={"A": -1, "B": -2})      # Add a 2nd row

    m.set_field_last_snapshot("Y", -123)

    assert list(m.collection_df.columns) == ["SYSTEM TIME", "A", "B", "caption", "X", "Y"]   # New column present

    data_values = [ {"SYSTEM TIME": 100, "A": 1,  "B": 22, "caption": "first entry", "X": 99,     "Y": np.nan},
                    {"SYSTEM TIME": 200, "A": -1, "B": -2, "caption": "",            "X": np.nan, "Y": -123}]
    expected_df = pd.DataFrame(data_values)

    assert_frame_equal(m.collection_df, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point



def test_update_last_snapshot():
    m = CollectionTabular()
    m.store(par=100, data_snapshot={"A": 1, "B": 2}, caption="first entry")     # Add a 1st row

    m.update_last_snapshot({"A": 11, "B": 22})
    last_row = list(m.collection_df.loc[0])
    assert last_row == [100, 11, 22, 'first entry']

    m.update_last_snapshot({"A": 9, "X": 99})
    last_row = list(m.collection_df.loc[0])
    assert list(m.collection_df.columns) == ["SYSTEM TIME", "A", "B", "caption", "X"]   # New column present
    assert last_row == [100, 9, 22, 'first entry', 99]

    m.store(par=200, data_snapshot={"A": -1, "B": -2})      # Add a 2nd row

    m.update_last_snapshot({"Y": -123})

    assert list(m.collection_df.columns) == ["SYSTEM TIME", "A", "B", "caption", "X", "Y"]   # New column present

    data_values = [ {"SYSTEM TIME": 100, "A": 9,  "B": 22, "caption": "first entry", "X": 99,     "Y": np.nan},
                    {"SYSTEM TIME": 200, "A": -1, "B": -2, "caption": "",            "X": np.nan, "Y": -123}]
    expected_df = pd.DataFrame(data_values)

    assert_frame_equal(m.collection_df, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point





###############  For class CollectionArray  ###############

def test_CollectionArray():
    m = CollectionArray()

    m.store(par=10, data_snapshot=np.array([1., 2., 3.]), caption="first entry")
    assert m.parameters == [10]
    assert m.captions == ["first entry"]
    assert m.snapshot_shape == (3,)
    assert np.allclose(m.data_arr, [[1., 2., 3.]])
    assert len(m) == 1
    assert str(m) == "CollectionArray object with 1 snapshot(s) parametrized by `SYSTEM TIME`"

    m.store(par=20, data_snapshot=np.array([10., 11., 12.]), caption="second entry")
    assert m.parameters == [10, 20]
    assert m.captions == ["first entry", "second entry"]
    assert m.snapshot_shape == (3,)
    assert np.allclose(m.data_arr, [[1., 2., 3.],
                                    [10., 11., 12.]])
    assert len(m) == 2
    assert str(m) == "CollectionArray object with 2 snapshot(s) parametrized by `SYSTEM TIME`"

    m.store(par=30, data_snapshot=np.array([-10., -11., -12.]))
    assert m.parameters == [10, 20, 30]
    assert m.captions == ["first entry", "second entry", ""]
    assert m.snapshot_shape == (3,)
    assert np.allclose(m.data_arr, [[1., 2., 3.],
                                    [10., 11., 12.],
                                    [-10., -11., -12.]])
    assert len(m) == 3
    assert str(m) == "CollectionArray object with 3 snapshot(s) parametrized by `SYSTEM TIME`"

    with pytest.raises(Exception):
        m.store(par=666, data_snapshot=np.array([1., 2., 3., 4., 5.]),
                caption="doesn't conform to shape of earlier entries!")


    # Again, in higher dimensions (and using methods to fetch the attributes)
    m_2D = CollectionArray(parameter_name="a,b values")

    m_2D.store(par={"a": 4., "b": 12.3},
               data_snapshot=np.array([[1., 2., 3.],
                                       [10., 11., 12.]]))
    assert m_2D.get_parameters() == [{"a": 4., "b": 12.3}]
    assert m_2D.get_captions() == [""]
    assert m_2D.get_shape() == (2, 3)
    assert np.allclose(m_2D.data_arr, [[1., 2., 3.],
                                       [10., 11., 12.]])
    assert len(m_2D) == 1
    assert str(m_2D) == "CollectionArray object with 1 snapshot(s) parametrized by `a,b values`"

    m_2D.store(par={"a": 400., "b": 123},
               data_snapshot=np.array([[-1., -2., -3.],
                                       [-10., -11., -12.]]
                                      ),
               caption="2nd matrix")
    assert m_2D.get_parameters() == [{"a": 4., "b": 12.3}, {"a": 400., "b": 123}]
    assert m_2D.get_captions() == ["", "2nd matrix"]
    assert m_2D.get_shape() == (2, 3)
    expected = np.array([
                            [[1., 2., 3.],
                             [10., 11., 12.]]
                            ,
                            [[-1., -2., -3.],
                             [-10., -11., -12.]]
                        ])
    assert np.allclose(m_2D.data_arr, expected)
    assert len(m_2D) == 2
    assert str(m_2D) == "CollectionArray object with 2 snapshot(s) parametrized by `a,b values`"

    with pytest.raises(Exception):
        m_2D.store(par=666, data_snapshot=np.array([10., 11., 12.]), caption="wrong shape for the data!")




###############  For class CollectionGeneral  ###############

def test_Collection():
    m = Collection()

    m.store(par=10, data_snapshot={"c1": 1, "c2": 2}, caption="first entry")
    assert len(m) == 1
    assert str(m) == "Collection object with 1 snapshot(s) parametrized by `SYSTEM TIME`"
    assert m.get_collection() == [(10, {"c1": 1, "c2": 2}, "first entry")]

    m.store(par=20, data_snapshot=[999, 111], caption="data snapshots can be anything")
    assert len(m) == 2
    assert str(m) == "Collection object with 2 snapshot(s) parametrized by `SYSTEM TIME`"
    data = m.get_collection()
    assert len(data) == 2
    assert data[0] == (10, {"c1": 1, "c2": 2}, "first entry")
    assert data[1] == (20, [999, 111], "data snapshots can be anything")

    m.store(par=(1,2), data_snapshot="001001101")
    assert len(m) == 3
    assert str(m) == "Collection object with 3 snapshot(s) parametrized by `SYSTEM TIME`"
    data = m.get_collection()
    assert len(data) == 3
    assert data[0] == (10, {"c1": 1, "c2": 2}, "first entry")
    assert data[1] == (20, [999, 111], "data snapshots can be anything")
    assert data[2] == ((1,2), "001001101", "")

    assert m.get_captions() == ["first entry", "data snapshots can be anything", ""]
    assert m.get_parameters() == [10, 20, (1,2)]
    assert m.get_data() == [{"c1": 1, "c2": 2} , [999, 111] , "001001101"]
