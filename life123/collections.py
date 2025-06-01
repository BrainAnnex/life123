# 3 CLASSES: CollectionTabular, CollectionArray, Collection

import pandas as pd
import numpy as np


class CollectionTabular:
    """
    A "tabular collection" is a Pandas dataframe
    built up from a sequence of "snapshots" of data that's in  the form of a python dictionary
    (representing a list of values and their corresponding names),
    such as the state of the system or of parts thereof.

    Each data "snapshots" is taken at different times,
    or results from varying some parameter.

    Each snapshot - incl. its parameter values and optional captions -
    will constitute a "row" in a tabular format

    MAIN DATA STRUCTURE for "tabular" collections:
        A Pandas dataframe
    """

    def __init__(self, parameter_name="SYSTEM TIME"):
        """
        :param parameter_name:  A label explaining what the snapshot parameter is.
                                Typically it's "SYSTEM TIME" (default), but could be anything
                                Used as the Pandas column name for the
                                parameter value attached to the various snapshot captures
        """
        self.parameter_name = parameter_name

        self.collection_df = pd.DataFrame()      # Empty Pandas dataframe



    def __len__(self) -> int:
        """
        Return the number of snapshots in the collection

        :return:    An integer
        """
        return len(self.collection_df)



    def __str__(self):
        return f"`CollectionTabular` object with {len(self.collection_df)} snapshot(s) parametrized by `{self.parameter_name}`"



    def store(self, par, data_snapshot: dict, caption="") -> None:
        """
        Save up the given data snapshot, alongside the specified parameter value and optional caption

        EXAMPLE :
                store(par=8., data_snapshot={"A": 1., "B": 2.}, caption="State immediately before injection of 2nd reactant")

        :param par:             Typically, the System Time - but could be any value that parametrizes the snapshots
        :param data_snapshot:   A dict of data to preserve for later use;
                                    it's acceptable to contain new fields not used in previous calls
                                    (in that case, the dataframe will add new columns automatically - and NaN values
                                     will appear in earlier rows)
        :param caption:         [OPTIONAL] String to describe the snapshot.
                                    Use None to avoid including that column (if it already exists in the dataframe, it'll appear as NaN)
        :return:                None (the object variable "self.collection_df" will get updated)
        """
        assert type(data_snapshot) == dict, \
            "CollectionTabular.store(): The argument `data_snapshot` must be a python dictionary"

        if self.collection_df.empty:            # No Pandas dataframe was yet started
            self.collection_df = pd.DataFrame(data_snapshot, index=[0])     # Form the initial Pandas dataframe (zero refers to the initial row)
                                                                            #   from the given data dictionary
            self.collection_df.insert(0, self.parameter_name, par)          # Add a column at the beginning, for the snapshot parameter
            if caption is not None:
                self.collection_df["caption"] = caption                     # Add a column at the end, for the caption
        else:                                   # The Pandas dataframe was already started
            d = data_snapshot.copy()                    # Make a copy, to avoid altering the passed dict
            d[self.parameter_name] = par                # Expand the snapshot dict
            if caption is not None:
                d["caption"] = caption                  # Expand the snapshot dict

            self.collection_df = pd.concat([self.collection_df, pd.DataFrame([d])], ignore_index=True)    # Append new row to dataframe
            # Note: we cannot do an in-place addition of a new row to an existing dataframe,
            #       because this new row might contain fields not
            #       yet present in the dataframe
            #       TODO: look into in-place addition when the keys of the dict are all among existing df columns



    def get_dataframe(self, head=None, tail=None,
                      val_start=None, val_end=None,
                      search_col=None, search_val=None, return_copy=True) -> pd.DataFrame:
        """
        Return the main data structure (a Pandas dataframe) 
        - or a part thereof (in which case a column named "search_value" is inserted to the left.)

        Optionally, limit the dataframe to a specified numbers of rows at the end,
        or just return row(s) corresponding to a specific search value(s) in the specified column
        - i.e. the row(s) with the CLOSEST value to the requested one(s).

        IMPORTANT:  if multiple options to restrict the dataset are present, only one is carried out;
                    the priority is:  1) head,  2) tail,  3) filtering,  3) search

        :param head:        [OPTIONAL] Integer.  If provided, only show the first several rows;
                                as many as specified by that number.
        :param tail:        [OPTIONAL] Integer.  If provided, only show the last several rows;
                                as many as specified by that number.
                                If the "head" argument is passed, this argument will get ignored

        :param val_start:  [OPTIONAL] Perform a FILTERING using the start value in the the specified column
                                - ASSUMING the dataframe is ordered by that value (e.g. a system time)
        :param val_end:    [OPTIONAL] FILTER by end value.
                                Either one or both of start/end values may be provided

        :param search_col:  [OPTIONAL] String with the name of a column in the dataframe,
                                against which to match the value below
        :param search_val:  [OPTIONAL] Number, or list/tuple of numbers, with value(s)
                                to search in the above column

        :param return_copy: [OPTIONAL] If True (default), the returned dataframe is guaranteed to be a (deep) copy -
                                so that modifying it won't affect the internal dataframe

        :return:            A Pandas dataframe, with all or some of the rows
                                that were stored in the main data structure.
                                If a search was requested, insert a column named "search_value" to the left
        """
        # The main data structure (a Pandas dataframe), with the "saved snapshots", is available as self.collection_df
        if return_copy:
            df = self.collection_df.copy()  # Note: some of the operations below also make a copy
        else:
            df = self.collection_df

        if head is not None:
            return df.head(head)    # This request is given top priority

        if tail is not None:
            return df.tail(tail)    # This request is given next-highest priority

        if search_col is None:
            return df               # Without a search column, neither filtering nor search are possible;
                                    #   so, return everything


        # If we get thus far, we're doing a SEARCH or FILTERING, and were given a column name;
        # we'll first look into whether we have a FILTERING request

        if (val_start is not None) and (val_end is not None):
            return df[df[search_col].between(val_start, val_end)]

        if val_start is not None:
            return df[df[search_col] >= val_start]

        if val_end is not None:
            return df[df[search_col] <= val_end]


        # If we get thus far, we're doing a SEARCH, and were given a column name

        if search_val is not None:
            # Perform a lookup of a value, or set of values,
            # in the specified column
            if (type(search_val) == tuple) or (type(search_val) == list):
                # Looking up a group of values
                lookup = pd.merge_asof(pd.DataFrame({'search_value':search_val}), df,
                                       right_on=search_col, left_on='search_value',
                                       direction='nearest')
                # Note: the above operation will insert a column named "search_value" to the left
                return lookup
            else:
                # Looking up a single value lookup
                index = df[search_col].sub(search_val).abs().idxmin()   # Locate index of row with smallest
                                                                        # absolute value of differences from the search value
                lookup = df.iloc[index : index+1]               # Select 1 row
                lookup.insert(0, 'search_value', [search_val])  # Insert a column named "search_value" to the left
                return lookup


        # In the absence of a passed search_val, return the full dataset
        return df



    def clear_dataframe(self) -> None:
        """
        Do a clean start

        :return:    None
        """
        self.collection_df = pd.DataFrame()      # Empty Pandas dataframe



    def set_caption_last_snapshot(self, caption: str) -> None:
        """
        Set the caption field of the last (most recent) snapshot to the given value.
        Any previous value gets over-written

        :param caption: String containing a caption to write into the last (most recent) snapshot
        :return:        None
        """
        index = len(self.collection_df) - 1
        self.collection_df.loc[index, "caption"] = caption



    def set_field_last_snapshot(self, field_name: str, field_value) -> None:
        """
        Set the specified field of the last (most recent) snapshot to the given value.
        Any previous value gets over-written.
        If the specified field name is not already one of the columns in the underlying
        data frame, a new column by that name gets added; any previous rows will have the value NaN
        assigned to that column

        :param field_name:  Name of field of interest
        :param field_value: Value to write into the above field for the last (most recent) snapshot
        :return:            None
        """
        last_row_index = len(self.collection_df) - 1
        self.collection_df.loc[last_row_index, field_name] = field_value



    def update_last_snapshot(self, update_values: dict) -> None:
        """
        Set some fields of the last (most recent) snapshot to the given values.
        Any previous value gets over-written.
        If any field name is not already among the columns in the underlying
        data frame, a new column by that name gets added; any previous rows will have the value NaN
        assigned to that column

        :param update_values:   Dict whose keys are the names of the columns to update
        :return:                None
        """
        last_row_index = len(self.collection_df) - 1
        self.collection_df.loc[last_row_index, update_values.keys()] = update_values.values()





###############################################################################################################

class CollectionArray:
    """
    Use this structure if your "snapshots" (data to add to the cumulative collection) are Numpy arrays,
    of any dimension - but always retaining that same dimension.

    Usually, the snapshots will be dumps of the entire system state, or parts thereof, but could be anything.
    Typically, each snapshot is taken at a different time (for example, to create a history), but could also
    be the result of varying some parameter(s)

    DATA STRUCTURE:
        A Numpy array 1 dimension larger than that of the snapshots

        EXAMPLE: if the snapshots are the 1-d numpy arrays [1., 2., 3.] and [10., 20., 30.]
                        then the internal structure will be the matrix
                        [[1., 2., 3.],
                         [10., 20., 30.]]
    """

    def __init__(self, parameter_name="SYSTEM TIME"):
        """
        :param parameter_name:  A label explaining what the snapshot parameter is.
                                Typically it's "SYSTEM TIME" (default), but could be anything
        """
        self.parameter_name = parameter_name

        self.data_arr = None           # A Numpy Array
        self.snapshot_shape = None
        self.parameters = []
        self.captions = []



    def __len__(self) -> int:
        """
        Return the number of snapshots in the collection

        :return:    An integer
        """
        return self.data_arr.shape[0]


    def __str__(self):
        return f"CollectionArray object with {self.__len__()} snapshot(s) parametrized by `{self.parameter_name}`"



    def store(self, par, data_snapshot: np.array, caption = "") -> None:
        """
        Save up the given data snapshot, and its associated parameters and optional caption

        EXAMPLES:
                store(par = 8., data_snapshot = np.array([1., 2., 3.]), caption = "State after injection of 2nd reactant")
                store(par = {"a": 4., "b": 12.3}, data_snapshot = np.array([1., 2., 3.]))

        :param par:             Typically, the System Time - but could be anything that parametrizes the snapshots
                                    (e.g., a dictionary, or any desired data structure.)
                                    It doesn't have to remain consistent, but it's probably good practice to keep it so
        :param data_snapshot:   A Numpy array, of any shape - but must keep that same shape across snapshots
        :param caption:         OPTIONAL string to describe the snapshot
        :return:                None
        """
        assert type(data_snapshot) == np.ndarray, \
            "CollectionArray.store(): The argument `data_snapshot` must be a dictionary whenever a 'tabular' collection is created"


        if self.data_arr is None:   # If this is the first call to this function
            self.data_arr = np.array([data_snapshot])
            self.snapshot_shape = data_snapshot.shape
        else:                       # this is a later call, to expand existing stored data
            assert data_snapshot.shape == self.snapshot_shape, \
                f"CollectionArray.store(): The argument `data_snapshot` must have the same shape across calls, namely {self.snapshot_shape}"
            new_arr = [data_snapshot]
            self.data_arr = np.concatenate((self.data_arr, new_arr), axis=0)    # "Stack up" along the first axis

        self.parameters.append(par)
        self.captions.append(caption)



    def get_array(self) -> np.array:
        """
        Return the main data structure - the Numpy Array

        :return:    A Numpy Array with the main data structure
        """
        return self.data_arr


    def get_parameters(self) -> list:
        """
        Return all the parameter values

        :return:    A list with the parameter values
        """
        return self.parameters


    def get_captions(self) -> [str]:
        """
        Return all the captions

        :return:    A list with the captions
        """
        return self.captions


    def get_shape(self) -> tuple:
        """

        :return:    A tuple with the shape of the snapshots
        """
        return self.snapshot_shape




###############################################################################################################

class Collection:
    """
    A "Collection" is a list of snapshots that the user wants to preserve,
    such as the state of the entire system, or of parts thereof,
    either taken at different times,
    or resulting from varying some parameter(s)

    This class accept data in ARBITRARY formats.
    If your data is Numpy arrays, you may use the more specialized class "CollectionArray"

    MAIN DATA STRUCTURE:
        A list of triplets.
        Each triplet is of the form (parameter value, caption, snapshot_data)
            1) The "parameter" is typically time, but could be anything.
               (a descriptive meaning of this parameter is stored in the object attribute "parameter_name")
            2) "snapshot_data" can be anything of interest, typically a clone of some data element.
            3) "caption" is just a string with an optional label.

        If the "parameter" is time, it's assumed to be in increasing order

        EXAMPLE:
            [
                (0., DATA_STRUCTURE_1, "Initial state"),
                (8., DATA_STRUCTURE_2, "State immediately after injection of 2nd reactant")
            ]
    """

    def __init__(self, parameter_name="SYSTEM TIME"):
        """
        :param parameter_name:  A label explaining what the snapshot parameter is.
                                Typically it's "SYSTEM TIME" (default), but could be anything
        """
        self.parameter_name = parameter_name

        self.data = []     # List of triples



    def __len__(self):
        """
        Return the number of snapshots in the collection

        :return:    An integer
        """
        return len(self.data)



    def __str__(self):
        return f"Collection object with {len(self.data)} snapshot(s) parametrized by `{self.parameter_name}`"



    def store(self, par, data_snapshot, caption = "") -> None:
        """
        Save up the given data snapshot

        EXAMPLE:
                store(par = 8.,
                      data_snapshot = {"c1": 10., "c2": 20.},
                      caption = "State immediately before injection of 2nd reactant")

                store(par = {"a": 4., "b": 12.3},
                     data_snapshot = [999., 111.],
                     caption = "Parameter is a dict; data is a list")

        IMPORTANT:  if passing a variable pointing to an existing mutable structure (such as a list, dict, object)
                    make sure to first *clone* it, to preserve it as it!

        :param par:             Typically, the System Time - but could be anything that parametrizes the snapshots
                                    (e.g., a dictionary, or any desired data structure.)
                                    It doesn't have to remain consistent, but it's probably good practice to keep it so
        :param data_snapshot:   Data in any format (such as a Numpy array, or an object)
        :param caption:         OPTIONAL string to describe the snapshot
        :return:                None
        """
        self.data.append((par, data_snapshot, caption))



    def get_collection(self) -> list:
        """
        Return the main data structure - the list of snapshots, with their attributes

        :return:
        """
        return self.data



    def get_data(self) -> list:
        """
        Return a list of all the data snapshots

        :return:    A list of all the snapshots
        """
        #TODO: maybe offer a way to only extract a portion
        return [triplet[1] for triplet in self.data]


    def get_parameters(self) -> list:
        """
        Return a list of all the parameter values

        :return:    A list with all the parameter values
        """
        return [triplet[0] for triplet in self.data]


    def get_captions(self) -> [str]:
        """
        Return a list of all the captions

        :return:    A list with all the captions
        """
        return [triplet[2] for triplet in self.data]
