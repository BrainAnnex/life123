import pytest
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from life123 import CollectionTabular, CollectionArray, Collection
from tests.utilities.comparisons import *



###############  For class CollectionTabular  ###############

def test_store():
    ct = CollectionTabular()

    d = {"A": 1, "B": 2, "C": 3}
    ct.store(par=10, data_snapshot=d, caption="first entry")  # Add a snapshot
    assert d == {"A": 1, "B": 2, "C": 3}        # Unchanged
    assert len(ct) == 1
    assert str(ct) == "`CollectionTabular` object with 1 snapshot(s) parametrized by `SYSTEM TIME`.  To access, use its get_dataframe() method"
    data = ct.get_raw_data()
    assert data == [{"A": 1, "B": 2, "C": 3, "SYSTEM TIME" :10, "caption": "first entry"}]

    ct.store(par=20, data_snapshot={"A": 10, "B": 20, "C": 30}, caption=None)  # Add a snapshot (no caption)
    assert len(ct) == 2
    assert str(ct) == "`CollectionTabular` object with 2 snapshot(s) parametrized by `SYSTEM TIME`.  To access, use its get_dataframe() method"
    data = ct.get_raw_data()
    assert data == [    {"A": 1, "B": 2, "C": 3,    "SYSTEM TIME" :10, "caption": "first entry"},
                        {"A": 10, "B": 20, "C": 30, "SYSTEM TIME" :20}
                   ]

    ct.store(par=30, data_snapshot={"A": -1, "B": -2, "C": -3})      # Add a snapshot (blank caption)
    assert len(ct) == 3
    assert str(ct) == "`CollectionTabular` object with 3 snapshot(s) parametrized by `SYSTEM TIME`.  To access, use its get_dataframe() method"
    data = ct.get_raw_data()
    assert data == [    {"A": 1, "B": 2, "C": 3,    "SYSTEM TIME" :10, "caption": "first entry"},
                        {"A": 10, "B": 20, "C": 30, "SYSTEM TIME" :20},
                        {"A": -1, "B": -2, "C": -3, "SYSTEM TIME" :30, "caption": ""}
                   ]



def test_get_dataframe_1():
    ct = CollectionTabular()

    ct.store(par=10, data_snapshot={"A": 1, "B": 2, "C": 3}, caption="first entry")
    df = ct.get_dataframe()
    assert list(df.columns) == ["SYSTEM TIME", "A", "B", "C", "caption"]
    row = list(df.iloc[0])                      # By row index
    assert row == [10, 1, 2, 3, 'first entry']

    ct.store(par=20, data_snapshot={"A": 10, "B": 20, "C": 30}, caption=None)
    df = ct.get_dataframe()
    assert list(df.columns) == ["SYSTEM TIME", "A", "B", "C", "caption"]

    row = list(df.iloc[0])                      # The previous row
    assert row == [10, 1, 2, 3, 'first entry']
    row = list(df.iloc[1])                      # The new row
    np.testing.assert_equal(row, [20, 10, 20, 30, np.nan])

    ct.store(par=30, data_snapshot={"A": -1, "B": -2, "C": -3})      # Add a snapshot (blank caption)
    df = ct.get_dataframe()
    assert list(df.columns) == ["SYSTEM TIME", "A", "B", "C", "caption"]
    row = list(df.iloc[2])
    assert row == [30, -1, -2, -3, ""]

    ct.store(par=40, data_snapshot={"A": 111, "B": 222}, caption="notice that C is missing")  # Add a snapshot
    df = ct.get_dataframe()
    assert list(df.columns) == ["SYSTEM TIME", "A", "B", "C", "caption"]
    row = list(df.iloc[3])
    np.testing.assert_equal(row, [40, 111, 222, np.nan, "notice that C is missing"])

    # Add a snapshot with an extra field, a boolean
    ct.store(par=50, data_snapshot={"A": 8, "B": 88, "C": 888, "aborted": True}, caption="notice the newly-appeared field")
    df = ct.get_dataframe()
    assert list(df.columns) == ["SYSTEM TIME", "A", "B", "C", "aborted", "caption"]
    row = list(df.iloc[4])
    np.testing.assert_equal(row, [50, 8, 88, 888, True, "notice the newly-appeared field"])

    data =  [
                {"SYSTEM TIME": 10, "A": 1, "B": 2, "C": 3, "caption": "first entry"},
                {"SYSTEM TIME": 20, "A": 10, "B": 20, "C": 30},
                {"SYSTEM TIME": 30, "A": -1, "B": -2, "C": -3, "caption": ""},
                {"SYSTEM TIME": 40, "A": 111, "B": 222, "caption": "notice that C is missing"},
                {"SYSTEM TIME": 50, "A": 8, "B": 88, "C": 888, "aborted": True, "caption": "notice the newly-appeared field"}
            ]
    df_expected = pd.DataFrame(data)
    #print(df)
    #print(df_expected)

    assert compare_pandas(df, df_expected, disregard_order=True)



def test_get_dataframe_2():
    ct = CollectionTabular()

    # Now the `SYSTEM TIME` parameter now has floats values
    ct.store(par=10,   data_snapshot={"A": 1, "B": 2, "C": 3}, caption="first entry")
    ct.store(par=12.4, data_snapshot={"A": 10, "B": 20, "C": 30}, caption="second entry")
    ct.store(par=33.1, data_snapshot={"A": -1, "B": -2, "C": -3})
    ct.store(par=40,   data_snapshot={"A": 111, "B": 222}, caption="notice that C is missing")
    ct.store(par=50.5, data_snapshot={"A": 8, "B": 88, "C": 888, "D": 1}, caption="notice the newly-appeared D")

    """
       SYSTEM TIME      A    B      C     D                    caption
    0           10      1    2    3.0   NaN                first entry  
    1           12.4   10   20   30.0   NaN               second entry  
    2           33.1   -1   -2   -3.0   NaN                            
    3           40    111  222    NaN   NaN   notice that C is missing  
    4           50.5    8   88  888.0   1.0  notice the newly-appeared  
    """
    # Check the extraction of the last row
    df_last_row = ct.get_dataframe(tail=1)
    #print("\n", df_last_row)

    data_values = [{"SYSTEM TIME": 50.5, "A": 8, "B": 88, "C": 888, "D": 1, "caption": "notice the newly-appeared D"}]
    expected_df = pd.DataFrame(data_values, index=[4])
    #print("\n", expected_df)

    assert_frame_equal(df_last_row, expected_df, check_dtype=False) # To allow for slight discrepancies in floating-point
                                                                    # (since int's get converted to floats in columns with Nan's)

    # Check the extraction of the last 2 rows
    df_last_2_rows = ct.get_dataframe(tail=2)

    data_values = [ {"SYSTEM TIME": 40, "A": 111, "B": 222, "C": np.nan, "D": np.nan,"caption": "notice that C is missing"},
                    {"SYSTEM TIME": 50.5, "A": 8, "B": 88,  "C": 888,    "D": 1,     "caption": "notice the newly-appeared D"}]
    expected_df = pd.DataFrame(data_values, index=[3, 4])

    assert_frame_equal(df_last_2_rows, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point


    # Check the extraction of a row by value SEARCH (using a value for "SYSTEM TIME" a tad smaller than what is in the dataframe)
    df_extracted_row = ct.get_dataframe(search_col="SYSTEM TIME", search_val=33.099)

    data_values = [{"search_value": 33.099, "SYSTEM TIME": 33.1, "A": -1, "B": -2, "C": -3.0, "D": np.nan, "caption": ""}]
    expected_df = pd.DataFrame(data_values, index=[2])
    #print("\n", expected_df)
    assert_frame_equal(df_extracted_row, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point


    # Check the extraction of a row by value (this time using a slightly larger value for "SYSTEM TIME" than what is in the dataframe)
    df_extracted_row = ct.get_dataframe(search_col="SYSTEM TIME", search_val=33.1234)
    #print("\n", df_extracted_row)

    data_values = [{"search_value": 33.1234, "SYSTEM TIME": 33.1, "A": -1, "B": -2, "C": -3.0, "D": np.nan, "caption": ""}]
    expected_df = pd.DataFrame(data_values, index=[2])
    assert_frame_equal(df_extracted_row, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point


    # Check the extraction of a group of row by value-range filtering
    df_filtered = ct.get_dataframe(search_col="SYSTEM TIME", val_start=35)        # This corresponds to the last 2 rows

    data_values = [ {"SYSTEM TIME": 40, "A": 111, "B": 222, "C": np.nan, "D": np.nan, "caption": "notice that C is missing"},
                    {"SYSTEM TIME": 50.5, "A": 8, "B": 88,  "C": 888,    "D": 1,      "caption": "notice the newly-appeared D"}]
    expected_df = pd.DataFrame(data_values, index=[3, 4])

    assert_frame_equal(df_filtered, expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point



def test_set_caption_last_snapshot():
    m = CollectionTabular()
    m.store(par=100, data_snapshot={"A": 1, "B": 2, "C": 3}, caption="first entry")
    m.store(par=200, data_snapshot={"A": 10, "B": 20, "C": 30})
    m.set_caption_last_snapshot("End of experiment")
    last_row = list(m.get_dataframe().loc[1])
    assert last_row == [200, 10, 20, 30, 'End of experiment']



def test_set_field_last_snapshot():
    m = CollectionTabular()
    m.store(par=100, data_snapshot={"A": 1, "B": 2}, caption="first entry")     # Add a 1st row

    m.set_field_last_snapshot("B", 22)
    last_row = list(m.get_dataframe().loc[0])
    assert last_row == [100, 1, 22, 'first entry']

    m.set_field_last_snapshot("X", 99)
    last_row = list(m.get_dataframe().loc[0])
    assert list(m.get_dataframe().columns) == ["SYSTEM TIME", "A", "B", "X", "caption"]   # New column present
    assert last_row == [100, 1, 22, 99, 'first entry']

    m.store(par=200, data_snapshot={"A": -1, "B": -2})      # Add a 2nd row

    m.set_field_last_snapshot("Y", -123)

    assert list(m.get_dataframe().columns) == ["SYSTEM TIME", "A", "B", "X", "Y", "caption"]   # New column present

    data_values = [ {"SYSTEM TIME": 100, "A": 1,  "B": 22, "X": 99,     "Y": np.nan, "caption": "first entry"},
                    {"SYSTEM TIME": 200, "A": -1, "B": -2, "X": np.nan, "Y": -123,   "caption": ""}]
    expected_df = pd.DataFrame(data_values)

    assert_frame_equal(m.get_dataframe(), expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point



def test_update_last_snapshot():
    m = CollectionTabular()
    m.store(par=100, data_snapshot={"A": 1, "B": 2}, caption="first entry")     # Add a 1st row

    m.update_last_snapshot({"A": 11, "B": 22})
    last_row = list(m.get_dataframe().loc[0])
    assert last_row == [100, 11, 22, 'first entry']

    m.update_last_snapshot({"A": 9, "X": 99})
    last_row = list(m.get_dataframe().loc[0])
    assert list(m.get_dataframe().columns) == ["SYSTEM TIME", "A", "B", "X", "caption"]   # New column present
    assert last_row == [100, 9, 22, 99, 'first entry']

    m.store(par=200, data_snapshot={"A": -1, "B": -2})      # Add a 2nd row

    m.update_last_snapshot({"Y": -123})

    assert list(m.get_dataframe().columns) == ["SYSTEM TIME", "A", "B", "X", "Y", "caption"]   # New column present

    data_values = [ {"SYSTEM TIME": 100, "A": 9,  "B": 22, "X": 99,     "Y": np.nan, "caption": "first entry"},
                    {"SYSTEM TIME": 200, "A": -1, "B": -2, "X": np.nan, "Y": -123,   "caption": ""}]
    expected_df = pd.DataFrame(data_values)

    assert_frame_equal(m.get_dataframe(), expected_df, check_dtype=False)  # To allow for slight discrepancies in floating-point






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
