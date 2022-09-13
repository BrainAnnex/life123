import pandas as pd


class Movie:
    """
    A "movie" is a list of snapshots
    of the state of the entire system, or of parts thereof, or any values that the user wants to preserve,
    either taken at different times,
    or resulting from varying some parameter(s)

    2 modalities are allowed (to be picked at object instantiation):
        A - a "tabular" mode.  Straightforward and convenient; it lends itself to handy Pandas dataframes
        B - a "non-tabular" mode.  More complex: to preserve data in arbitrary formats


    MAIN DATA STRUCTURE for "tabular" mode:
        A Pandas frame


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

        :param parameter_name:
        :param tabular:
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

        EXAMPLE with tabular mode:
                append(8., {"A": 1., "B": 2.}, "State immediately after injection of 2nd reagent")
        EXAMPLE without tabular mode:
                append(8., NUMPY_ARRAY_2, "State immediately after injection of 2nd reagent")

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
                    "Movie.append(): The argument `data_snapshot` must be a dictionary whenever a 'tabular' movie is created"
                self.df = pd.DataFrame(data_snapshot, index=[0])            # Form the initial Pandas dataframe
                self.df[self.parameter_name] = pars                         # Add a column
                if caption:
                    self.df["caption"] = caption                            # Add a column
            else:
                data_snapshot[self.parameter_name] = pars                   # Expand the snapshot dict
                if caption:
                    data_snapshot["caption"] = caption                      # Expand the snapshot dict
                self.df = pd.concat([self.df, pd.DataFrame([data_snapshot])], ignore_index=True)    # Append new row to dataframe



    def get(self) -> list:
        """
        Return the main data structure - the list of snapshots, with their attributes

        :return:
        """
        if self.tabular:
            return self.df
        else:
            return self.movie
