# 3 CLASSES: MovieTabular, MovieArray, MovieGeneral

import pandas as pd
import numpy as np
from typing import Union


class MovieTabular:
    """
    A "tabular movie" is a list of snapshots of data that comes as a list or one-dimensional array,
    such as something based on the state of the system or of parts thereof,
    either taken at different times,
    or resulting from varying some parameter(s)

    Each snapshot - incl. its parameter values and optional captions -
    will constitute a "row" in a tabular format

    MAIN DATA STRUCTURE for "tabular" mode:
        A Pandas dataframe
    """

    def __init__(self, parameter_name="SYSTEM TIME"):
        """
        :param parameter_name:  A label explaining what the snapshot parameter is.
                                Typically it's "SYSTEM TIME" (default), but could be anything
                                Used as the Pandas column name for the
                                    parameter value attached to the snapshot captures
        """
        self.parameter_name = parameter_name

        self.movie = None      # Pandas dataframe



    def __len__(self) -> int:
        """
        Return the number of snapshots comprising the movie
        :return:    An integer
        """
        return len(self.movie)



    def __str__(self):
        return f"MovieTabular object with {len(self.movie)} snapshot(s) parametrized by `{self.parameter_name}`"



    def store(self, par, data_snapshot: dict, caption="") -> None:
        """
        Save up the given data snapshot, alongside the specified parameter value and optional caption

        EXAMPLE :
                store(8., {"A": 1., "B": 2.}, "State immediately before injection of 2nd reagent")

        :param par:             Typically, the System Time - but could be anything that parametrizes the snapshots
                                    (e.g., a dictionary, or any desired data structure.)
                                    It doesn't have to remain consistent, but it's probably good practice to keep it so
        :param data_snapshot:   A dict of data to preserve for later use
        :param caption:         OPTIONAL string to describe the snapshot
        :return:                None (the object variable "self.movie" will get updated)
        """
        if self.movie is None:     # No Pandas dataframe was yet started
            assert type(data_snapshot) == dict, \
                "MovieTabular.store(): The argument `data_snapshot` must be a dictionary"

            self.movie = pd.DataFrame(data_snapshot, index=[0])     # Form the initial Pandas dataframe (zero refers to the initial row)
            self.movie.insert(0, self.parameter_name, par)          # Add a column at the beginning
            self.movie["caption"] = caption                         # Add a column at the end
        else:                       # The Pandas dataframe was already started
            data_snapshot[self.parameter_name] = par                # Expand the snapshot dict
            data_snapshot["caption"] = caption                      # Expand the snapshot dict
            self.movie = pd.concat([self.movie, pd.DataFrame([data_snapshot])], ignore_index=True)    # Append new row to dataframe



    def get_movie(self, search_col=None, search_val=None, val_start=None, val_end=None, tail=None) -> pd.DataFrame:
        """
        Return the main data structure (a Pandas dataframe) 
        - or a part thereof (in which case insert a column named "search_value" to the left.)

        Optionally, limit the dataframe to a specified numbers of rows at the end,
        or just return row(s) corresponding to a specific search value(s) in the specified column
        - i.e. the row(s) with the CLOSEST value to the requested one(s).

        :param search_col:  (OPTIONAL) String with the name of a column in the dataframe,
                                against which to match the value below
        :param search_val:  (OPTIONAL) Number, or list/tuple of numbers, with value(s)
                                to search in the above column

        :param val_start:  (OPTIONAL) Perform a FILTERING using the start value in the the specified column
                                - ASSUMING the dataframe is ordered by that value (e.g. a system time)
        :param val_end:    (OPTIONAL) FILTER by end value.
                                Either one or both of start/end values may be provided

        :param tail:        (OPTIONAL) Integer.  If provided, only show the last several rows;
                                as many as specified by that number.

        NOTE: if multiple options to restrict the dataset are present, only one is carried out;
              the priority is:  1) tail,  2) filtering,  3) search

        :return:            A Pandas dataframe, with all or some of the rows
                                that were stored in the main data structure.
                                If a search was requested, insert a column named "search_value" to the left
        """
        df = self.movie     # The main data structure (a Pandas dataframe), with the "saved snapshots"

        if tail is not None:
            return df.tail(tail)    # This request is given top priority

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



    def set_caption_last_snapshot(self, caption: str) -> None:
        """
        Set the caption field of the last (most recent) snapshot to the given value.
        Any previous value gets over-written

        :param caption:
        :return:        None
        """
        index = len(self.movie) - 1
        self.movie.loc[index, "caption"] = caption




###############################################################################################################

class MovieArray:
    """
    Use this structure if your "snapshots" (data to add to the cumulative collection) are Numpy arrays,
    of any dimension - but always retaining that same dimension.

    Usually, the snapshots will be dump of the entire system state, or parts thereof, but could be anything.
    Typically, each snapshot is taken at a different time (for example, to create a history), but could also
    be the result of varying some parameter(s)

    DATA STRUCTURE:
        A Numpy array of 1 dimension larger than that of the snapshots

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

        self.movie = None           # A Numpy Array
        self.snapshot_shape = None
        self.parameters = []
        self.captions = []



    def __len__(self) -> int:
        """
        Return the number of snapshots comprising the movie
        :return:    An integer
        """
        return self.movie.shape[0]


    def __str__(self):
        return f"MovieArray object with {self.__len__()} snapshot(s) parametrized by `{self.parameter_name}`"



    def store(self, par, data_snapshot: np.array, caption = "") -> None:
        """
        Save up the given data snapshot, and its associated parameters and optional caption

        EXAMPLES:
                store(par = 8., data_snapshot = np.array([1., 2., 3.]), caption = "State after injection of 2nd reagent")
                store(par = {"a": 4., "b": 12.3}, data_snapshot = np.array([1., 2., 3.]))

        :param par:             Typically, the System Time - but could be anything that parametrizes the snapshots
                                    (e.g., a dictionary, or any desired data structure.)
                                    It doesn't have to remain consistent, but it's probably good practice to keep it so
        :param data_snapshot:   A Numpy array, of any shape - but must keep that same shape across snapshots
        :param caption:         OPTIONAL string to describe the snapshot
        :return:                None
        """
        assert type(data_snapshot) == np.ndarray, \
            "MovieArray.store(): The argument `data_snapshot` must be a dictionary whenever a 'tabular' movie is created"


        if self.movie is None:    # If this is the first call to this function
            self.movie = np.array([data_snapshot])
            self.snapshot_shape = data_snapshot.shape
        else:                   # this is a later call, to expand existing stored data
            assert data_snapshot.shape == self.snapshot_shape, \
                f"MovieArray.store(): The argument `data_snapshot` must have the same shape across calls, namely {self.snapshot_shape}"
            new_arr = [data_snapshot]
            self.movie = np.concatenate((self.movie, new_arr), axis=0)    # "Stack up" along the first axis

        self.parameters.append(par)
        self.captions.append(caption)



    def get_array(self) -> np.array:
        """
        Return the main data structure - the Numpy Array

        :return:    A Numpy Array with the main data structure
        """
        return self.movie


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

class MovieGeneral:
    """
    A "general movie" is a list of snapshots of any values that the user wants to preserve,
    such as the state of the entire system, or of parts thereof,
    either taken at different times,
    or resulting from varying some parameter(s)

    This class accept data in arbitrary formats

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
                (8., DATA_STRUCTURE_2, "State immediately after injection of 2nd reagent")
            ]
    """

    def __init__(self, parameter_name="SYSTEM TIME"):
        """
        :param parameter_name:  A label explaining what the snapshot parameter is.
                                Typically it's "SYSTEM TIME" (default), but could be anything
        """
        self.parameter_name = parameter_name

        self.movie = []     # List of triples



    def __len__(self):
        return len(self.movie)



    def __str__(self):
        return f"MovieGeneral object with {len(self.movie)} snapshot(s) parametrized by `{self.parameter_name}`"



    def store(self, par, data_snapshot, caption = "") -> None:
        """
        Save up the given data snapshot

        EXAMPLE:
                store(par = 8.,
                      data_snapshot = {"c1": 10., "c2": 20.},
                      caption = "State immediately before injection of 2nd reagent")

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
        self.movie.append((par, data_snapshot, caption))



    def get_movie(self) -> list:
        """
        Return the main data structure - the list of snapshots, with their attributes

        :return:
        """
        return self.movie



    def get_data(self) -> list:
        """
        Return a list of all the snapshots

        :return:    A list of all the snapshots
        """
        return [triplet[1] for triplet in self.movie]


    def get_parameters(self) -> list:
        """
        Return all the parameter values

        :return:    A list with all the parameter values
        """
        return [triplet[0] for triplet in self.movie]


    def get_captions(self) -> [str]:
        """
        Return all the captions

        :return:    A list with all the captions
        """
        return [triplet[2] for triplet in self.movie]
