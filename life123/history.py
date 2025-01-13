import os
import csv
import pandas as pd
from life123.collections import CollectionTabular, Collection



class History:

    def __init__(self, active=False, chem_labels=None, frequency=1):
        """
        For the management of historical data of various types,
        each of which relies on a child class.
        This parent class - which should NOT be directly instantiated -
        provides common properties and methods.

        :param active:
        :param chem_labels:
        :param frequency:
        """
        self.history = None                 # Object set by the children classes
                                            # The type, depending on the child class, will be one of:
                                            # CollectionTabular, CollectionArray, Collection
        self.active = active
        self.capture_frequency = frequency

        self.restrict_chemicals = chem_labels

        self.step_of_last_snapshot = -1     # An out-of-range value

        self.initial_caption = None

        self.log_file = None



    def get_history(self):
        """

        :return:    Object of one of the following types:
                        CollectionTabular, CollectionArray, Collection
        """
        return self.history



    def enable_history_common(self, frequency=1, chem_labels=None) -> None:
        """

        :param frequency:
        :param chem_labels:
        :return:            None
        """
        self.active = True
        self.capture_frequency = frequency
        self.restrict_chemicals = chem_labels



    def disable_history(self):
        self.active = False



    def to_capture(self, step_count=None, extra=False) -> bool:
        """
        Return True if various criteria indicate that a snapshot should be taken at the give step_count;
        False otherwise

        Note that when step_count is 0 there's no capture, unless the frequency is 1

        :param step_count:  None means always recapture
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
                                    the object property self.initial_caption is used, if set
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

class HistoryUniformConcentration(History):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)          # Invoke the constructor of its parent class, passing all the available args
        self.history = CollectionTabular(parameter_name="SYSTEM TIME")  # To store user-selected snapshots of quantities of interest



    def enable_history(self, frequency=1, chem_labels=None) -> None:
        """
        Request history capture, with the specified parameters.
        If history was already enabled, this function can be used to alter its capture parameters.

        :param frequency:
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





##########################################################################################################################

class HistoryBinConcentration(History):

    def __init__(self, bins=None, *args, **kwargs):
        super().__init__(*args, **kwargs)          # Invoke the constructor of its parent class, passing all the available args
        self.restrict_bins = bins                               # Bin address, or list of them, or None
        self.history = Collection(parameter_name="SYSTEM TIME") # To store user-selected snapshots of quantities of interest



    def enable_history(self, frequency=1, chem_labels=None, bins=None) -> None:
        """
        Request history capture, with the specified parameters.
        If history was already enabled, this function can be used to alter its capture parameters.

        :param frequency:
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
                { 6: {"A": 1.3, "B": 3.9},
                  8: {"A": 4.6, "B": 2.7}
                }

           EXAMPLE of data_snapshot (for 2D) systems:
                { (0,0): {"A": 1.3, "B": 3.9},
                  (3,5): {"A": 4.6, "B": 2.7}
                }
        :param system_time:
        :param data_snapshot:
        :param step_count:
        :param caption:         [OPTIONAL] String to save alongside this snapshot
        :return:                None
        """
        self.save_snapshot_common(system_time=system_time, data_snapshot=data_snapshot, step_count=step_count,
                                  caption=caption)



    def bin_history(self, bin_address, include_captions=False):
        """
        Return the history at the given bin, as a Pandas dataframe.
        The first column is "SYSTEM TIME", and the other ones are the various chemicals for which
        history had been enabled.
        If no historical data is located, an informational string is returned

        :param bin_address:         A single bin address.  EXAMPLE, in 2D : (3,3)
        :param include_captions:    If True, the captions are returned as an extra "caption" column at the end
        :return:                    A Pandas data frame, or a string if no historical data is present
        """
        assert type(bin_address) == int or type(bin_address) == tuple, \
            "bin_history(): Argument `bin_address` must be an integer (for 1D) of a tuple (2+ D)"

        snapshot_list = self.history.get_data()
        if snapshot_list == []:
            return f"No concentration historical data available for bin {bin_address}"

        first_snapshot = snapshot_list[0]         # EXAMPLE: {(3, 3): {'A': 0.0, 'B': 0.0, 'C': 0.0}}
        first_conc_data = first_snapshot.get(bin_address)   # EXAMPLE: {'A': 0.0, 'B': 0.0, 'C': 0.0}
        if first_conc_data is None: \
            return f"No concentration historical data available for bin {bin_address}"

        chem_labels = list(first_conc_data)                 # Turn the dict's keys into list.  EXAMPLE: ['A', 'B', 'C']

        if include_captions:
            cols = ['SYSTEM TIME'] + chem_labels + ['caption']
        else:
            cols = ['SYSTEM TIME'] + chem_labels

        df = pd.DataFrame(columns = cols) # EXAMPLE: ['SYSTEM TIME', 'A', 'B', 'C']

        for row in self.history.get_collection():
            t , data, caption = row           # Unpack
            conc_dict = data[bin_address]
            l = len(df)
            df.loc[l] = conc_dict           # Add a new row, with all columns except "SYSTEM TIME"
            df.loc[l, "SYSTEM TIME"] = t    # Set the "SYSTEM TIME" cell for the row just added
            if include_captions:
                df.loc[l, "caption"] = caption

        return df




##########################################################################################################################

class HistoryReactionRate(History):

    def __init__(self, rxns=None, *args, **kwargs):
        super().__init__(*args, **kwargs)          # Invoke the constructor of its parent class, passing all the available args
        self.restrict_rxns = rxns
        self.history = CollectionTabular(parameter_name="SYSTEM TIME")  # To store user-selected snapshots of quantities of interest



    def enable_history(self, frequency=1, chem_labels=None, rxns=None) -> None:
        """
        Request history capture, with the specified parameters.
        If history was already enabled, this function can be used to alter its capture parameters.

        :param frequency:
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
