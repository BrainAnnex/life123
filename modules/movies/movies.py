import pandas as pd
import numpy as np
from typing import Union


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
            2) "caption" is just a string with an optional label.
            3) "snapshot_data" can be anything of interest, typically a clone of some data element.

        If the "parameter" is time, it's assumed to be in increasing order

        EXAMPLE:
            [
                (0., "Initial state", DATA_STRUCTURE_1),
                (8., "State immediately after injection of 2nd reagent", DATA_STRUCTURE_2)
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
        return f"General-type Movie object with {len(self.movie)} snapshots parametrized by `{self.parameter_name}`"



    def store(self, pars, data_snapshot, caption = "") -> None:
        """
        Save up the given data snapshot

        EXAMPLE:
                store(8., NUMPY_ARRAY_2, "State immediately before injection of 2nd reagent")

        IMPORTANT:  if passing a variable pointing to an existing mutable structure (such as a list, dict, object)
                    make sure to first *clone* it, to preserve it as it!

        :param pars:            Typical, the System Time... but could be anything that parametrizes the snapshots
        :param data_snapshot:   Data in any format (such as a Numpy array, or an object)
        :param caption:         OPTIONAL string to describe the snapshot
        :return:                None
        """
        self.movie.append( (pars, caption, data_snapshot) )



    def get(self) -> list:
        """
        Return the main data structure - the list of snapshots, with their attributes

        :return:
        """
        return self.movie




###############################################################################

class MovieTabular:
    """
    A "tabular movie" is a list of snapshots of data that comes as a list or one-dimensional array,
    such as something based on the state of the system or of parts thereof,
    either taken at different times,
    or resulting from varying some parameter(s)

    Each snapshot will constitute a "row" in a tabular format

    MAIN DATA STRUCTURE for "tabular" mode:
        A Pandas dataframe
    """

    def __init__(self, parameter_name="SYSTEM TIME", tabular=True):
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
        return f"Tabular Movie object with {len(self.movie)} snapshots parametrized by `{self.parameter_name}`"



    def store(self, pars, data_snapshot: dict, caption = "") -> None:
        """
        Save up the given data snapshot

        EXAMPLE :
                store(8., {"A": 1., "B": 2.}, "State immediately before injection of 2nd reagent")

        :param pars:            Typical, the System Time... but could be anything that parametrizes the snapshots
        :param data_snapshot:   A dict
        :param caption:         OPTIONAL string to describe the snapshot
        :return:                None
        """
        if self.movie is None:     # No Pandas dataframe was yet started
            assert type(data_snapshot) == dict, \
                "Movie.store(): The argument `data_snapshot` must be a dictionary whenever a 'tabular' movie is created"
            self.movie = pd.DataFrame(data_snapshot, index=[0])            # Form the initial Pandas dataframe
            self.movie[self.parameter_name] = pars                         # Add a column
            if caption:
                self.movie["caption"] = caption                            # Add a column
        else:
            data_snapshot[self.parameter_name] = pars                   # Expand the snapshot dict
            if caption:
                data_snapshot["caption"] = caption                      # Expand the snapshot dict
            self.movie = pd.concat([self.movie, pd.DataFrame([data_snapshot])], ignore_index=True)    # Append new row to dataframe



    def get(self) -> pd.DataFrame:
        """
        Return the main data structure - a Pandas dataframe

        :return:
        """
        return self.movie



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



    def store(self, pars, data_snapshot: np.array, caption = "") -> None:
        """
        Save up the given data snapshot, and its associated parameters and optional caption

        EXAMPLES:
                store(8., np.array([1., 2., 3.]), caption = "State immediately after injection of 2nd reagent")
                store({"a": 4., "b": 12.3}, np.array([1., 2., 3.]))

        :param pars:            Typical, the System Time - but could be anything that parametrizes the snapshots
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

        self.parameters.append(pars)
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





#############    *** BEING PHASED OUT ***    ##############################

class Movie:
    """
    A "movie" is a list of snapshots
    of the state of the entire system, or of parts thereof, or any values that the user wants to preserve,
    either taken at different times,
    or resulting from varying some parameter(s)

    2 modalities are allowed (to be picked at object instantiation):
        A - a "tabular" mode.  Straightforward and convenient; it lends itself to handy Pandas dataframes
        B - a "non-tabular" mode.  More complex: to preserve data in arbitrary formats

    TODO: maybe split into 2 separate classes?  (IN-PROGRESS!)


    MAIN DATA STRUCTURE for "tabular" mode:
        A Pandas dataframe


    MAIN DATA STRUCTURE for "non-tabular" mode:
        A list of triplets.
        Each triplet is of the form (parameter value, caption, snapshot_data)
            1) The "parameter" is typically time, but could be anything.
               (a descriptive meaning of this parameter is stored in the object attribute "parameter_name")
            2) "caption" is just a string with an optional label.
            3) "snapshot_data" can be anything of interest, typically a clone of some data element.

        If the "parameter" is time, it's assumed to be in increasing order

        EXAMPLE:
            [
                (0., "Initial state", NUMPY_ARRAY_1),
                (8., "State immediately after injection of 2nd reagent", NUMPY_ARRAY_2)
            ]
    """


    def __init__(self, parameter_name="SYSTEM TIME", tabular=True):
        """

        :param parameter_name:  Only used in the "tabular" mode, as the Pandas column name for the
                                    parameter value attached to any particular snapshot capture
        :param tabular:         A flag indicating whether the "tabular" format will be used in this object.
                                    (Explained in the notes above)
        """
        self.parameter_name = parameter_name
        self.tabular = tabular

        self.movie = []
        self.df = None      # It will only be applicable if tabular is True



    def __len__(self):
        return len(self.movie)



    def __str__(self):
        return f"Movie with {len(self.movie)} snapshots of type `{self.parameter_name}`"



    def store(self, pars, data_snapshot, caption = "") -> None:
        """
        Save up the given data snapshot

        EXAMPLE in tabular mode:
                store(8., {"A": 1., "B": 2.}, "State immediately after injection of 2nd reagent")
        EXAMPLE not using tabular mode:
                store(8., NUMPY_ARRAY_2, "State immediately after injection of 2nd reagent")

        IMPORTANT:  if passing a variable pointing to an existing mutable structure (such as a list, dict, object)
                    make sure to first *clone* it, to preserve it as it!

        :param pars:            Typical, the System Time... but could be anything that parametrizes the snapshots
        :param data_snapshot:   If using a tabular system, this must be a dict;
                                    otherwise, it can be anything
        :param caption:         OPTIONAL string to describe the snapshot
        :return:                None
        """
        if not self.tabular:
            self.movie.append( (pars, caption, data_snapshot) )
        else:
            if self.df is None:     # No Pandas dataframe was yet started
                assert type(data_snapshot) == dict, \
                    "Movie.store(): The argument `data_snapshot` must be a dictionary whenever a 'tabular' movie is created"
                self.df = pd.DataFrame(data_snapshot, index=[0])            # Form the initial Pandas dataframe
                self.df[self.parameter_name] = pars                         # Add a column
                if caption:
                    self.df["caption"] = caption                            # Add a column
            else:
                data_snapshot[self.parameter_name] = pars                   # Expand the snapshot dict
                if caption:
                    data_snapshot["caption"] = caption                      # Expand the snapshot dict
                self.df = pd.concat([self.df, pd.DataFrame([data_snapshot])], ignore_index=True)    # Append new row to dataframe



    def get(self) -> Union[list, pd.DataFrame]:
        """
        Return the main data structure - the list of snapshots, with their attributes, OR a Pandas dataframe

        :return:
        """
        if self.tabular:
            return self.df
        else:
            return self.movie
