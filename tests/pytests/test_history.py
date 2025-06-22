import pytest
import numpy as np
import pandas as pd
from pandas.testing import assert_frame_equal
from life123 import CollectionTabular, CollectionArray, Collection
from life123 import HistoryBinConcentration, HistoryReactionRate, HistoryUniformConcentration
from life123.history import History



###############  For class History  ###############



###############  For class HistoryBinConcentration  ###############

def test_bin_history():
    hist = HistoryBinConcentration(active=False)
    hist.enable_history()

    data_snapshot = \
            { 5: {"A": 1, "B": 3},
              8: {"A": 4, "B": 2}
            }
    hist.save_snapshot(system_time=0, data_snapshot=data_snapshot)

    data_snapshot = \
            { 5: {"A": 2, "B": 7},
              8: {"A": 5, "B": 0}
            }
    hist.save_snapshot(system_time=100, data_snapshot=data_snapshot, caption="2nd entry")

    #l = hist.history.get_collection()
    #print(l)

    with pytest.raises(Exception):
        hist.bin_history(bin_address="bad address")

    assert hist.bin_history(bin_address=666) == "" \
        "No historical concentration data available for bin 666. History collecting IS enabled for ALL bins - but was it enabled PRIOR to running the simulation?"


    df = hist.bin_history(bin_address=5)
    expected = pd.DataFrame([[0, 1, 3], [100, 2, 7]], columns = ["SYSTEM TIME", "A", "B"])  # FULL MATRIX, specified by row
    assert df.equals(expected)

    df = hist.bin_history(bin_address=8)
    expected = pd.DataFrame([[0, 4, 2], [100, 5, 0]], columns = ["SYSTEM TIME", "A", "B"])  # FULL MATRIX, specified by row
    assert df.equals(expected)

    df = hist.bin_history(bin_address=8, include_captions=True)
    expected = pd.DataFrame([[0, 4, 2, ""], [100, 5, 0, "2nd entry"]], columns = ["SYSTEM TIME", "A", "B", "caption"])
    assert df.equals(expected)


    data_snapshot = \
            { 5: {"A": 10, "B": 14},
              8: {"A": 18, "B": 20}
            }
    hist.save_snapshot(system_time=200, data_snapshot=data_snapshot, caption="one more entry")

    df = hist.bin_history(bin_address=8)
    expected = pd.DataFrame([[0, 4, 2], [100, 5, 0], [200, 18, 20]],
                            columns = ["SYSTEM TIME", "A", "B"])
    assert df.equals(expected)

    df = hist.bin_history(bin_address=8, downsize=2)    # Include every other row, starting with the initial one
    expected = pd.DataFrame([[0, 4, 2], [200, 18, 20]],
                            columns = ["SYSTEM TIME", "A", "B"])
    assert df.equals(expected)

    df = hist.bin_history(bin_address=8, downsize=3)    # Every 3rd row, but always include the last one
    expected = pd.DataFrame([[0, 4, 2], [200, 18, 20]],
                            columns = ["SYSTEM TIME", "A", "B"])
    assert df.equals(expected)


    data_snapshot = \
            { 5: {"A": 60, "B": 42},
              8: {"A": 38, "B": 2}
            }
    hist.save_snapshot(system_time=500, data_snapshot=data_snapshot, caption="last entry")

    df = hist.bin_history(bin_address=8, downsize=3)    # Every 3rd row
    expected = pd.DataFrame([[0, 4, 2], [500, 38, 2]],
                            columns = ["SYSTEM TIME", "A", "B"])
    assert df.equals(expected)

    df = hist.bin_history(bin_address=8, include_captions=True, downsize=3)    # Every 3rd row
    expected = pd.DataFrame([[0, 4, 2, ""], [500, 38, 2, "last entry"]],
                            columns = ["SYSTEM TIME", "A", "B", "caption"])
    assert df.equals(expected)

    df = hist.bin_history(bin_address=8, include_captions=True, downsize=1000)    # We again end up with the first (0-th) and last rows
    assert df.equals(expected)



###############  For class HistoryReactionRate  ###############



###############  For class HistoryUniformConcentration  ###############