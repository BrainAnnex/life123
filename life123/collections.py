# 3 CLASSES: CollectionTabular, CollectionArray, Collection

import pandas as pd
import numpy as np


class CollectionTabular:
    """
    A "tabular collection" is a Pandas dataframe
    built up from a sequence of "snapshots" of data that's in the form of a python dictionary
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

        self.snapshot_list = []                 # The data being accumulated;
                                                # each element represents a "data snapshot" (which is a dict)



    def __len__(self) -> int:
        """
        Return the number of snapshots in the collection

        :return:    An integer
        """
        return len(self.snapshot_list)



    def __str__(self) -> str:
        return f"`CollectionTabular` object with {self.__len__()} snapshot(s)" \
        f" parametrized by `{self.parameter_name}`.  To access, use its get_dataframe() method"



    def get_raw_data(self):
        return self.snapshot_list



    def store(self, par, data_snapshot :dict, caption="") -> None:
        """
        Save up the given data snapshot, alongside the specified parameter value and optional caption.
        NOTE: the caption - or any other field - may also be changed
              by a later call to `set_caption_last_snapshot()`

        EXAMPLE :
                store(par=8., data_snapshot={"A": 1., "B": 2.}, caption="State immediately before injection of 2nd reactant")

        :param par:             Typically, the System Time - but could be any value that parametrizes the snapshots
        :param data_snapshot:   A dict of data to preserve for later use; it does NOT get modified.
                                    It's acceptable for it to contain new fields not used in previous calls
                                    (in that case, the dataframe will add new columns automatically - and NaN values
                                     will appear in earlier rows)
        :param caption:         [OPTIONAL] String to describe the snapshot.
                                    Use None to avoid including that column (if it already exists in the dataframe, it'll appear as NaN)

        :return:                None (the object variable "self.collection" will get updated)
        """
        #print(f"CollectionTabular.store(): par={par} | data_snapshot={data_snapshot} | caption={caption}")
        assert type(data_snapshot) == dict, \
            "CollectionTabular.store(): The argument `data_snapshot` must be a python dictionary"

        d = data_snapshot.copy()                    # Make a copy, to avoid altering the passed dict

        d[self.parameter_name] = par                # Expand the snapshot dict
        if caption is not None:
            d["caption"] = caption                  # Expand the snapshot dict

        self.snapshot_list.append(d)



    def get_dataframe(self, head=None, tail=None,
                      val_start=None, val_end=None,
                      search_col=None, search_val=None, return_copy=True) -> pd.DataFrame:
        """
        Return the main data structure (a Pandas dataframe) 
        - or a part thereof (in which case a column named "search_value"
                             is inserted to the left into the result.)

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

        :param return_copy: TODO: OBSOLETE  - [OPTIONAL] If True (default), the returned dataframe is guaranteed to be a (deep) copy -
                                so that modifying it won't affect the internal dataframe

        :return:            A Pandas dataframe, with all or some of the rows
                                that were stored in the main data structure.
                                If a search was requested, insert a column named "search_value" to the left
        """
        # Turn the data into a Pandas dataframe
        df = pd.DataFrame(self.snapshot_list)

        # Make the column `self.parameter_name` to become the first (leftmost)
        col_name = self.parameter_name
        df.insert(0, col_name, df.pop(col_name))    # Extracts the column and inserts it at position 0

        # Make the column "caption", if present, to become the last (rightmost)
        col_name = "caption"
        if col_name in df:
            df[col_name] = df.pop(col_name)


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
        self.snapshot_list = []



    def set_caption_last_snapshot(self, caption :str) -> None:
        """
        Set the caption field of the last (most recent) snapshot to the given value.
        Any previous value gets over-written

        :param caption: String containing a caption to write into the last (most recent) snapshot
        :return:        None
        """
        # TODO: it'd be ideal not to have to depend on this...
        #print(f"*** SETTING CAPTION TO: '{caption}'")
        last_entry_index = len(self.snapshot_list) - 1
        snapshot_dict = self.snapshot_list[last_entry_index]
        snapshot_dict["caption"] = caption



    def set_field_last_snapshot(self, field_name :str, field_value) -> None:
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
        # TODO: it'd be ideal not to have to depend on this...
        #print(f"*** SETTING field `{field_name}` TO: '{field_value}'")
        last_entry_index = len(self.snapshot_list) - 1
        snapshot_dict = self.snapshot_list[last_entry_index]
        snapshot_dict[field_name] = field_value



    def update_last_snapshot(self, update_values :dict) -> None:
        """
        Set some fields of the last (most recent) snapshot to the given values.
        Any previous value gets over-written.
        If any field name is not already among the columns in the underlying
        data frame, a new column by that name gets added; any previous rows will have the value NaN
        assigned to that column

        :param update_values:   Dict whose keys are the names of the columns to update
        :return:                None
        """
        # TODO: it'd be ideal not to have to depend on this...
        last_entry_index = len(self.snapshot_list) - 1
        snapshot_dict = self.snapshot_list[last_entry_index]
        snapshot_dict |= update_values      # update a dictionary in place






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
