class Movies:
    """
    A "movie" is a list of snapshots of the state of the entire system,
    or of parts thereof,
    either taken at different times,
    or as a results of varying some parameter(s)

    MAIN DATA STRUCTURE:
        A list of triplets.
        Each triplet is of the form (parameter, caption, snapshot_data)
            The "parameter" is typically time, but could be anything.
            (a descriptive meaning of this parameter is stored in the object attribute "parameter_name"
            "caption" is just a string with a label.
            "snapshot_data" can be anything of interest, typically a clone of some data element.

        If the "parameter" is time, it's assumed to be in increasing order

    EXAMPLE:
        [
            (0., "State at t=0.", NUMPY_ARRAY_1),
            (8., "State at t=8.", NUMPY_ARRAY_2)
        ]
    """


    def __init__(self, parameter_name="TIME"):
        self.movie = []
        self.parameter_name = parameter_name



    def __len__(self):
        return len(self.movie)



    def __str__(self):
        return f"Movie with {len(self.movie)} snapshots of type {self.parameter_name}"



    def append(self, pars, caption:str, data_snapshot) -> None:
        """
        EXAMPLE:  append(8., "State at t=8.", NUMPY_ARRAY_2)

        :param pars:
        :param caption:
        :param data_snapshot:
        :return:                None
        """
        self.movie.append( (pars, caption, data_snapshot) )



    def get(self) -> list:
        """
        Return the main data structure - the list of snapshots, with their attributes

        :return:
        """
        return self.movie
