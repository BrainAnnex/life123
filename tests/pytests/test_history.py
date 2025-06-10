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
            { 5: {"A": 1.3, "B": 3.9},
              8: {"A": 4.6, "B": 2.7}
            }

    hist.save_snapshot(system_time=0, data_snapshot=data_snapshot)

    data_snapshot = \
            { 5: {"A": 2.3, "B": 3.1},
              8: {"A": 0.6, "B": 0.7}
            }

    hist.save_snapshot(system_time=0, data_snapshot=data_snapshot)

    l = hist.history.get_collection()

    print(l)

    df = hist.bin_history(bin_address=5)
    print(df)

    df = hist.bin_history(bin_address=8)
    print(df)

    '''
    Try this Pandas df structure;
    SYSTEM TIME | bin_5_A | bin_5_B | bin_8_A | bin_8_B
    
    OR:
    
    SYSTEM TIME | A_bin_5 | A_bin_8 | B_bin_5 | B_bin_8
    '''





###############  For class HistoryReactionRate  ###############



###############  For class HistoryUniformConcentration  ###############