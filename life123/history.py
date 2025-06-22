import os
import csv
import pandas as pd
from life123.collections import CollectionTabular, Collection



class History:

    def __init__(self, active=False, chem_labels=None, frequency=1):
        """
        For the management of historical data of various types,
        each of which relies on a child class.
        This base (parent) class is NOT meant for the end user
        - and should NOT be directly instantiated;
        it provides common properties and methods.

        :param active:      Flag indicating whether the history is enabled
        :param chem_labels: List of chemicals for which history is to be kept; use None to mean all
        :param frequency:   How many simulation cycles to wait until taking another data snapshot
        """
        self.history = None                 # Object set by the children classes.
                                            # The type, depending on the child class, will be one of:
                                            # CollectionTabular, CollectionArray, Collection

        self.active = active                # Flag indicating whether the history is enabled
        self.capture_frequency = frequency

        self.restrict_chemicals = chem_labels   # None means "no restriction" (i.e. ALL chemicals)

        self.step_of_last_snapshot = -1     # An out-of-range value

        self.initial_caption = None

        self.log_file = None




    def is_enabled(self):
        """
        Return True if the history-keeping is enabled

        :return:
        """
        return self.active



    def get_history(self):
        """

        :return:    Object of one of the following types:
                        CollectionTabular, CollectionArray, Collection
        """
        return self.history



    def enable_history_common(self, frequency=1, chem_labels=None) -> None:
        """
        Activate the keeping of historical data (the specifics
        will depend on the child class)

        :param frequency:   How many simulation cycles to wait until taking another data snapshot
        :param chem_labels:
        :return:            None
        """
        self.active = True
        self.capture_frequency = frequency
        self.restrict_chemicals = chem_labels



    def disable_history(self) -> None:
        """
        Inactivate any further captures related to this history.
        Nothing is deleted: neither the history nor its parameters;
        this history-keeping may later be resumed (possibly with different parameters)
        by a call to enable_history()

        :return:    None
        """
        self.active = False



    def to_capture(self, step_count=None, extra=False) -> bool:
        """
        Return True if various criteria indicate that a snapshot should be taken at the give step_count;
        False otherwise

        Note that when step_count is 0 there's no capture, unless the frequency is 1

        :param step_count:  None means always capture
        :param extra:       [OPTIONAL] If True, it means that this is a special extra capture;
                                the capture frequency will NOT considered in the
                                decision about saving it, but an extra check will be performed
                                to make sure it's not a duplicate of the earlier capture
        :return:
        """
        if step_count is None:
            return True

        if not self.active:
            return False

        if extra and (step_count == self.step_of_last_snapshot):
            return False     # In the case of special extra snapshots, make sure we didn't already take it

        if (self.capture_frequency == 1) or extra:
            return True

        if (step_count+1) % self.capture_frequency == 0:
            return True

        return False



    def start_csv_log(self, log_file :str) -> None:
        """
        Register the specified filename for all future CSV logs.
        Any existing file by that name will get over-written.

        :param log_file:Name of a file to use for the CSV logs
        :return:        None
        """
        # TODO: validate log_file
        self.log_file = log_file

        if os.path.exists(log_file):
            os.remove(log_file)
            print(f"-> CSV-format output will be LOGGED into the file '{log_file}' . An existing file by that name was over-written")
        else:
            print(f"-> CSV-format output will be LOGGED into the file '{log_file}'")



    def save_snapshot_common(self, system_time, data_snapshot, step_count=None, caption="") -> str:
        """

        :param system_time:
        :param data_snapshot:   Data format will vary in different child classes
        :param step_count:
        :param caption:         [OPTIONAL] String to save alongside this snapshot; if not provided,
                                    the object property self.initial_caption is used, if that was set
        :return:                The caption that was actually used
        """
        if (not caption) and self.initial_caption:
            caption = self.initial_caption
            self.initial_caption = None             # Reset value

        # The following is to avoid unsightly NaN's and floating-point numbers
        if step_count is None:
            data_snapshot["step"] = ""
        else:
            data_snapshot["step"] = str(step_count)

        self.history.store(par=system_time,
                           data_snapshot=data_snapshot, caption=caption)

        if step_count is not None:
            self.step_of_last_snapshot = step_count

        return caption



    def set_caption_last_snapshot(self, caption) -> None:
        """
        Add a caption to the very last entry in the system history.
        No action taken is caption is None

        :param caption:
        :return:        None
        """
        if caption is not None:
            self.history.set_caption_last_snapshot(caption)




##########################################################################################################################

class HistoryBinConcentration(History):
    """
    For the management of historical data for bin concentrations
    """

    def __init__(self, bins=None, *args, **kwargs):
        """

        :param bins:    List of bin addresses for which history is to be kept; use None to mean all
        :param args:
        :param kwargs:
        """
        super().__init__(*args, **kwargs)          # Invoke the constructor of its parent class, passing all the available args

        self.restrict_bins = bins                               # Bin address, or list of them,
                                                                #   or None (meaning ALL bins, with no restriction)
        self.history = Collection(parameter_name="SYSTEM TIME") # To store user-selected snapshots of quantities of interest
        # Note: an alternative (untried) to using general Collections for the history
        #       might be to use CollectionTabular tabular instead,
        #       with columns such as:  SYSTEM TIME | bin_5_A | bin_5_B | bin_8_A | bin_8_B | caption



    def enable_history(self, frequency=1, chem_labels=None, bins=None) -> None:
        """
        Request history capture, with the specified parameters.
        If history was already enabled, this function can be used to alter its capture parameters.

        :param frequency:   [OPTIONAL] How many simulation cycles to wait until taking another data snapshot
        :param chem_labels: [OPTIONAL] List of chemicals to include in the history;
                                if None (default), include them all.
        :param bins:        [OPTIONAL] Bin address, or list of them;
                                if None (default), include them all.
        :return:            None
        """
        self.enable_history_common(frequency=frequency, chem_labels=chem_labels)
        self.restrict_bins = bins



    def save_snapshot(self, system_time, data_snapshot, step_count=None, caption="") -> None:
        """

           EXAMPLE of data_snapshot (for 1D) systems:
                { 5: {"A": 1.3, "B": 3.9},
                  8: {"A": 4.6, "B": 2.7}
                }

           EXAMPLE of data_snapshot (for 2D) systems:
                { (0,0): {"A": 1.3, "B": 3.9},
                  (3,5): {"A": 4.6, "B": 2.7}
                }
        :param system_time:     The time of this data snapshot
        :param data_snapshot:
        :param step_count:
        :param caption:         [OPTIONAL] String to save alongside this snapshot
        :return:                None
        """
        self.save_snapshot_common(system_time=system_time, data_snapshot=data_snapshot, step_count=step_count,
                                  caption=caption)



    def _diagnose_history_problem(self, bin_address) -> str:
        """
        Return an error message elaborating on the problem, and possible causes,
        of lacking historical concentration data for the given bin

        :param bin_address: A single bin address.  EXAMPLES, in 1D: 8   In 2D : (3,3)
        :return:            A string elaborating on the problem, and possible causes,
                                of lacking historical concentration data for the given bin
        """
        if not self.is_enabled():
            return "History collecting is NOT enabled - you need to enable it PRIOR to running the simulation. Use enable_history()"

        msg = f"No historical concentration data available for bin {bin_address}. "

        if self.restrict_bins is None:
             msg += f"History collecting IS enabled for ALL bins - but was it enabled PRIOR to running the simulation?"
        elif bin_address in self.restrict_bins:
            msg += f"History collecting IS indeed enabled for this bin - but was it enabled PRIOR to running the simulation?"
        else:
            msg += f"History recording is NOT currently enabled for this bin; it is only enabled for bins {self.restrict_bins}"

        return  msg



    def bin_history(self, bin_address, include_captions=False, downsize=1):
        """
        Return the history at the given bin, as a Pandas dataframe.
        The first column is "SYSTEM TIME", and the other ones are the various chemicals for which
        history had been enabled.
        If no historical data is located, an informational string is returned instead

        :param bin_address:         A single bin address.  EXAMPLES, in 1D: 8   In 2D : (3,3)
        :param include_captions:    [OPTIONAL] If True, the captions are returned as an extra "caption" column at the end
        :param downsize:            [OPTIONAL] Pare down the returned history.
                                        If different from 1, include only every n-th entry, where n = downsize,
                                        starting with the initial (0-th) one.  The final row is also always included, regardless.
                                        Note: handy when the user unwisely collected an excessive amount of history!
        :return:                    A Pandas data frame,
                                        or a string with diagnostic suggestions if no historical data is present
        """
        assert type(bin_address) == int or type(bin_address) == tuple, \
            "bin_history(): Argument `bin_address` must be an integer (for 1D) of a tuple (for 2+ D)"

        snapshot_list = self.history.get_data()
        if snapshot_list == []:
            return self._diagnose_history_problem(bin_address)

        first_snapshot = snapshot_list[0]          # EXAMPLE: {(3, 3): {'A': 0.0, 'B': 0.0, 'C': 0.0}}
        first_conc_data = first_snapshot.get(bin_address)   # EXAMPLE: {'A': 0.0, 'B': 0.0, 'C': 0.0}
        if first_conc_data is None:
            return self._diagnose_history_problem(bin_address)

        chem_labels = list(first_conc_data)                 # Turn the dict's keys into list.  EXAMPLE: ['A', 'B', 'C']

        if include_captions:
            cols = ['SYSTEM TIME'] + chem_labels + ['caption']
        else:
            cols = ['SYSTEM TIME'] + chem_labels

        # EXAMPLE of cols:  ['SYSTEM TIME', 'A', 'B', 'C']
        #       or:         ['SYSTEM TIME', 'A', 'B', 'C', 'caption']

        history_list = self.history.get_collection()
        last_index = len(history_list) - 1

        data_matrix = []        # A list of lists
        for i, row in enumerate(history_list):   # Traverse the returned list
            if (downsize != 1) and (i != last_index) and (i % downsize) != 0:
                continue

            t , data, caption = row                 # Unpack the system time, concentration data and caption
            conc_dict = data[bin_address]               # EXAMPLE: {'A': 3.4, 'B': 0.8}
            conc_list = list(conc_dict.values())        # EXAMPLE: [3.4, 0.8]
            if include_captions:
                data_row = [t] + conc_list + [caption]
            else:
                data_row = [t] + conc_list

            data_matrix.append(data_row)

        df = pd.DataFrame(data_matrix, columns = cols)  # FULL MATRIX, specified by row

        return df




##########################################################################################################################

class HistoryReactionRate(History):
    """
    For the management of historical reaction-rate data for Uniform Compartments
    """

    def __init__(self, rxns=None, *args, **kwargs):
        """

        :param rxns:    List of reactions (identified by their index) for which history is to be kept;
                            use None to mean all
        :param args:
        :param kwargs:
        """
        super().__init__(*args, **kwargs)          # Invoke the constructor of its parent class, passing all the available args
        self.restrict_rxns = rxns
        self.history = CollectionTabular(parameter_name="SYSTEM TIME")  # To store user-selected snapshots of quantities of interest



    def enable_history(self, frequency=1, chem_labels=None, rxns=None) -> None:
        """
        Request history capture, with the specified parameters.
        If history was already enabled, this function can be used to alter its capture parameters.

        :param frequency:   [OPTIONAL] How many simulation cycles to wait until taking another data snapshot
        :param chem_labels: [OPTIONAL] List of chemicals to include in the history;
                                if None (default), include them all.
        :param rxns:
        :return:            None
        """
        self.enable_history_common(frequency=frequency, chem_labels=chem_labels)
        self.restrict_rxns = rxns



    def save_snapshot(self, system_time, data_snapshot, step_count=None, caption="") -> None:
        """

           EXAMPLE of data_snapshot:
                {"rxn1_rate": 6.3, "rxn2_rate": 14.3}

        :param system_time:
        :param data_snapshot:
        :param step_count:
        :param caption:         [OPTIONAL] String to save alongside this snapshot
        :return:                None
        """
        self.save_snapshot_common(system_time=system_time, data_snapshot=data_snapshot, step_count=step_count,
                                  caption=caption)




##########################################################################################################################

class HistoryUniformConcentration(History):
    """
    For the management of historical concentration data for Uniform Compartments
    """

    def __init__(self, *args, **kwargs):
        """

        :param args:
        :param kwargs:
        """
        super().__init__(*args, **kwargs)          # Invoke the constructor of its parent class, passing all the available args
        self.history = CollectionTabular(parameter_name="SYSTEM TIME")  # To store user-selected snapshots of quantities of interest



    def enable_history(self, frequency=1, chem_labels=None) -> None:
        """
        Request history capture, with the specified parameters.
        If history was already enabled, this function can be used to alter its capture parameters.

        :param frequency:   [OPTIONAL] How many simulation cycles to wait until taking another data snapshot
        :param chem_labels: [OPTIONAL] List of chemicals to include in the history;
                                if None (default), include them all.
        :return:            None
        """
        self.enable_history_common(frequency=frequency, chem_labels=chem_labels)



    def save_snapshot(self, system_time, data_snapshot, step_count=None, caption="") -> None:
        """

        :param system_time:
        :param data_snapshot:   EXAMPLE: {"A": 1.3, "B": 4.9}
        :param step_count:
        :param caption:         [OPTIONAL] String to save alongside this snapshot
        :return:                None
        """
        caption = self.save_snapshot_common(system_time=system_time, data_snapshot=data_snapshot, step_count=step_count,
                                            caption=caption)

        if self.log_file is None:
            return

        # If we get thus far, we have a log file to manage
        if os.path.exists(self.log_file):
            new_file = False
        else:
            new_file = True

        # Open a CSV file in write mode
        with open(self.log_file, mode="a", newline="") as fh:
            fieldnames = ["SYSTEM TIME"] + list(data_snapshot.keys()) + ["caption"]
            writer = csv.DictWriter(fh, fieldnames = fieldnames)

            data_snapshot["SYSTEM TIME"] = system_time
            data_snapshot["caption"] = caption

            # Write the header (only of a newly-created file)
            if new_file:
                writer.writeheader()

            # Write the dictionary as a row
            writer.writerow(data_snapshot)
