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



    def __len__(self):
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
        :return:                None
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



    def get(self, tail=None) -> pd.DataFrame:
        """
        Return the main data structure - a Pandas dataframe

        :param tail:    If an integer value is provided, only show the last several rows,
                            as many as specified by that number
        :return:        A Pandas dataframe
        """
        if tail is None:
            return self.movie
        else:
            return self.movie.tail(tail)



    def set_caption_last_snapshot(self, caption: str):
        """
        Set the caption field of the last (most recent) snapshot to the given value.
        Any previous value gets over-written

        :param caption:
        :return:
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



    def __len__(self):
        # Return the number of snapshots comprising the movie
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



    def get(self) -> list:
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
